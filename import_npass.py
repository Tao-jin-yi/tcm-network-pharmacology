"""导入 NPASS 天然产物数据集，扩充成分数据库。

NPASS（https://bidd.group/NPASS/）：Natural Product Activity & Species Source，
约 2 万天然产物 + 实验活性数据（IC50/Ki/EC50 等，来自文献与 PubChem BioAssay）。

流程：
1. 下载 4 个原始 txt 到 npass_data/（已存在则跳过，--force 重新下载）；
2. compound_info + compound_properties 合并 → npass_compounds.json（成分库）；
3. activities + target_info → npass_targets.json（人源蛋白实验活性，每对
   化合物-靶点保留最优活性值）。

用法：
    python import_npass.py                 # 全量
    python import_npass.py --no-activity   # 只建成分库，不处理活性
    python import_npass.py --force         # 重新下载原始文件

生成的 npass_compounds.json / npass_targets.json 由 app.py 自动加载；
原始 txt 不入 git（已加入 .gitignore），删掉后重跑本脚本即可重建。
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import requests

HERE = Path(__file__).parent
DATA_DIR = HERE / "npass_data"
OUT_COMPOUNDS = HERE / "npass_compounds.json"
OUT_TARGETS = HERE / "npass_targets.json"
BASE = "https://bidd.group/NPASS/download/"

FILES = [
    "NPASS_compound_info.txt",
    "NPASS_compound_properties.txt",
    "NPASS_target_info.txt",
    "NPASS_activities.txt",
]

# 只保留这些"可比较的抑制活性"类型（单位 nM），用于展示实验作用靶点
KEPT_ACTIVITY_TYPES = {"IC50", "IC90", "IC95", "IC99", "Ki", "EC50", "ED50",
                       "Potency", "A2", "Kb", "Kd", "MIC"}
HUMAN_TAX_ID = "9606"


def download(name: str, force: bool = False, max_retries: int = 12) -> Path:
    """下载原始文件；支持断点续传（Range），网络中断自动重试。"""
    target = DATA_DIR / name
    if target.exists() and not force:
        print(f"已存在，跳过下载：{name}")
        return target
    DATA_DIR.mkdir(exist_ok=True)
    start = time.time()
    existing = target.stat().st_size if target.exists() else 0
    retries = 0
    total = existing
    while True:
        headers = {"Range": f"bytes={total}-"} if total else {}
        try:
            with requests.get(BASE + name, stream=True,
                              timeout=(20, 120), headers=headers) as resp:
                resp.raise_for_status()
                mode = "ab" if total else "wb"
                with open(target, mode) as fh:
                    for chunk in resp.iter_content(chunk_size=1 << 20):
                        fh.write(chunk)
                        total += len(chunk)
                        if total % (10 << 20) < (1 << 20):
                            print(f"  已下载 {total/1e6:.0f} MB ...", flush=True)
            print(f"  完成：{total/1e6:.1f} MB，用时 {time.time()-start:.0f} 秒",
                  flush=True)
            return target
        except requests.RequestException as exc:
            retries += 1
            if retries > max_retries:
                raise
            print(f"  网络中断（{type(exc).__name__}），10 秒后从 "
                  f"{total/1e6:.1f} MB 处续传（第 {retries}/{max_retries} 次）",
                  flush=True)
            time.sleep(10)


def rows_of(path: Path):
    """逐行解析制表符分隔文件（跳过空行）。"""
    with open(path, encoding="utf-8", errors="ignore") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            yield dict(zip(header, line.split("\t")))


def build_compounds() -> list[dict]:
    """合并成分信息与理化属性表。"""
    info = {r["np_id"]: r for r in rows_of(DATA_DIR / "NPASS_compound_info.txt")}
    props = {r["np_id"]: r for r in rows_of(DATA_DIR / "NPASS_compound_properties.txt")}
    compounds = []
    for np_id, row in info.items():
        p = props.get(np_id, {})
        if p.get("canonical_smiles", "n.a.") in ("", "n.a."):
            continue  # 无结构的条目无法用于成分库
        compounds.append({
            "id": np_id,
            "name": row.get("pref_name", "n.a."),
            "iupac": row.get("iupac_name", ""),
            "chembl_id": row.get("chembl_id", "") if row.get("chembl_id") != "n.a." else "",
            "cid": row.get("pubchem_cid", "") if row.get("pubchem_cid") != "n.a." else "",
            "zinc_id": row.get("zinc_id", "") if row.get("zinc_id") != "n.a." else "",
            "formula": p.get("formula", ""),
            "mw": p.get("mw", ""),
            "xlogp": p.get("xlogp", ""),
            "psa": p.get("psa", ""),
            "hbd": p.get("hbd", ""),
            "hba": p.get("hba", ""),
            "rotb": p.get("rotatable_bond", ""),
            "lipinski_fail": p.get("lipinski_failure", ""),
            "inchikey": p.get("standard_inchi_key", ""),
            "smiles": p.get("canonical_smiles", ""),
        })
    return compounds


def build_targets() -> dict[str, list[dict]]:
    """人源蛋白实验活性：np_id → 靶点列表（每对化合物-靶点保留最优值）。"""
    targets = {r["target_id"]: r for r in rows_of(DATA_DIR / "NPASS_target_info.txt")}
    best: dict[tuple[str, str], dict] = {}
    for row in rows_of(DATA_DIR / "NPASS_activities.txt"):
        t = targets.get(row.get("target_id", ""))
        if not t:
            continue
        if t.get("target_organism_tax_id") != HUMAN_TAX_ID:
            continue
        if t.get("target_type") not in ("Individual Protein", "Protein Complex"):
            continue
        act_type = (row.get("activity_type") or "").strip()
        if act_type not in KEPT_ACTIVITY_TYPES:
            continue
        try:
            value = float(row.get("activity_value", ""))
        except ValueError:
            continue
        unit = (row.get("activity_units") or "").strip()
        if unit != "nM":  # 仅比较 nM 标准单位
            continue
        np_id = row["np_id"]
        uniprot = t.get("uniprot_id", "")
        key = (np_id, uniprot or t["target_id"])
        candidate = {"uni": uniprot, "name": t.get("target_name", ""),
                     "type": act_type, "val": value,
                     "rel": (row.get("activity_relation") or "").strip()}
        old = best.get(key)
        if old is None or candidate["val"] < old["val"]:
            best[key] = candidate

    grouped: dict[str, list[dict]] = {}
    for (np_id, _), rec in best.items():
        grouped.setdefault(np_id, []).append(rec)
    return grouped


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-activity", action="store_true",
                        help="只下载成分表并生成 npass_compounds.json")
    parser.add_argument("--force", action="store_true", help="重新下载原始文件")
    args = parser.parse_args()

    needed = FILES[:3] if args.no_activity else FILES
    for name in needed:
        download(name, force=args.force)

    compounds = build_compounds()
    OUT_COMPOUNDS.write_text(json.dumps(compounds, ensure_ascii=False, indent=0),
                             encoding="utf-8")
    print(f"成分库：{len(compounds)} 条 → {OUT_COMPOUNDS.name} "
          f"({OUT_COMPOUNDS.stat().st_size/1e6:.1f} MB)")

    if not args.no_activity:
        grouped = build_targets()
        OUT_TARGETS.write_text(json.dumps(grouped, ensure_ascii=False, indent=0),
                               encoding="utf-8")
        n_pairs = sum(len(v) for v in grouped.values())
        n_human = sum(1 for c in compounds
                      if c["id"] in grouped and "Homo" not in "")
        print(f"实验活性靶点：{n_pairs} 对（{len(grouped)} 个成分有人源蛋白靶点）"
              f" → {OUT_TARGETS.name} ({OUT_TARGETS.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    sys.exit(main())
