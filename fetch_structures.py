"""用 PubChem 官方接口补全成分的化学结构信息。

PubChem PUG REST（免费、无需 key）：
    https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<名称>/property/<字段>/JSON

获取字段：CanonicalSMILES、分子式、分子量、InChIKey、CID，并拼出 2D 结构图地址。
结果写入 compound_structures.json，供 data_prep.py 使用。

用法：
    python fetch_structures.py            # 只补缺失的成分
    python fetch_structures.py --all      # 重新抓取全部
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, Optional

import requests

from data_prep import COMPOUNDS

HERE = Path(__file__).parent
OUTPUT_FILE = HERE / "compound_structures.json"

API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/{props}/JSON"
PROPS = "SMILES,MolecularFormula,MolecularWeight,InChIKey"
IMG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size=large"


def fetch_pubchem(name: str) -> Optional[Dict]:
    """按化合物英文名查询 PubChem，失败返回 None。"""
    url = API.format(name=requests.utils.quote(name), props=PROPS)
    try:
        resp = requests.get(url, timeout=25)
        if resp.status_code != 200:
            return None
        info = resp.json()["PropertyTable"]["Properties"][0]
        cid = info.get("CID")
        return {
            "smiles": info.get("SMILES") or info.get("CanonicalSMILES"),
            "formula": info.get("MolecularFormula"),
            "mw": info.get("MolecularWeight"),
            "inchikey": info.get("InChIKey"),
            "cid": cid,
            "img_url": IMG.format(cid=cid) if cid else None,
            "source": "PubChem",
        }
    except Exception:
        return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true", help="重新抓取全部成分")
    parser.add_argument("--filtered", action="store_true",
                        help="只抓通过 OB≥30 / DL≥0.18 筛选的活性成分"
                             "（TCMSP 全库成分上万，全抓耗时数小时）")
    args = parser.parse_args()

    compounds = COMPOUNDS
    if args.filtered:
        compounds = {mol_id: comp for mol_id, comp in COMPOUNDS.items()
                     if comp.get("ob", 0) >= 30.0 and comp.get("dl", 0) >= 0.18}
        print(f"筛选后待抓成分 {len(compounds)} 个")

    cache: Dict[str, dict] = {}
    if OUTPUT_FILE.exists() and not args.all:
        cache = json.loads(OUTPUT_FILE.read_text(encoding="utf-8"))

    ok, fail = 0, []
    n_total = len(compounds)
    n_done = 0
    for mol_id, compound in compounds.items():
        n_done += 1
        if mol_id in cache and not args.all:
            continue
        name = compound.get("en") or compound["name"]
        data = fetch_pubchem(name)
        if not data and compound.get("en"):  # 英文名查不到时试中文名
            data = fetch_pubchem(compound["name"])
        if data:
            cache[mol_id] = data
            ok += 1
        else:
            fail.append(f"{mol_id} {compound['name']}({name})")
        # 每 50 个增量落盘并打印进度（中断可续传）
        if n_done % 50 == 0:
            OUTPUT_FILE.write_text(json.dumps(cache, ensure_ascii=False, indent=2),
                                   encoding="utf-8")
            print(f"[{n_done}/{n_total}] 累计 {len(cache)}，本次成功 {ok}",
                  flush=True)
        time.sleep(0.34)  # PubChem 限速：不超过 5 次/秒

    OUTPUT_FILE.write_text(json.dumps(cache, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"已写入 {OUTPUT_FILE}：累计 {len(cache)} 个成分，本次成功 {ok} 个")
    if fail:
        print("未获取到结构的成分：")
        for item in fail:
            print("  -", item)


if __name__ == "__main__":
    main()
