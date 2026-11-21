"""把 TCMSP 靶点蛋白名映射到 UniProt 标准基因符号（疾病层的前置步骤）。

TCMSP 的靶点以蛋白全名记录（如 "Prostaglandin G/H synthase 1"），
网络药理学后续分析（PPI、通路富集）需要标准基因符号（PTGS1）。
本脚本从 tcmsp_full.json 收集去重后的蛋白名，逐个查询 UniProt REST
（限人源 + Swiss-Prot 已审阅），结果增量缓存到 tcmsp_target_map.json：

    { 蛋白名: {"gene": "PTGS1", "accession": "P23219", "protein": "...", ...} }

未匹配的蛋白名记录 {"gene": ""}，避免重复查询。

用法：
    python fetch_target_mapping.py          # 全量映射（约 1500 个名字，10-15 分钟）
    python fetch_target_mapping.py --limit 5
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

HERE = Path(__file__).parent
FULL_FILE = HERE / "tcmsp_full.json"
OUTPUT_FILE = HERE / "tcmsp_target_map.json"

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,id,protein_name,gene_names,organism_name,length"
DELAY = 0.4  # UniProt 公共 API 限速


def collect_protein_names() -> list[str]:
    """从 tcmsp_full.json 收集去重靶点蛋白名（成分-靶点表 + 靶点-疾病表）。"""
    store = json.loads(FULL_FILE.read_text(encoding="utf-8"))
    names: set[str] = set()
    for herb in store.values():
        for row in herb.get("targets", []):
            name = (row.get("target_name") or "").strip()
            if name:
                names.add(name)
        for row in herb.get("diseases", []):
            name = (row.get("target_name") or "").strip()
            if name:
                names.add(name)
    return sorted(names)


def query_uniprot(protein_name: str) -> dict | None:
    """按蛋白全名查 UniProt（人源、已审阅），返回首个（推荐）匹配。"""
    params = {
        "query": f'protein_name:"{protein_name}" AND organism_id:9606 AND reviewed:true',
        "fields": FIELDS,
        "format": "json",
        "size": 1,
    }
    try:
        resp = requests.get(UNIPROT, params=params, timeout=30)
    except Exception:
        return None
    if resp.status_code != 200:
        return None
    results = resp.json().get("results", [])
    if not results:
        return {}
    entry = results[0]
    desc = entry.get("proteinDescription", {})
    protein = ""
    if "recommendedName" in desc:
        protein = desc["recommendedName"]["fullName"]["value"]
    elif desc.get("submittedName"):
        protein = desc["submittedName"][0]["fullName"]["value"]
    gene = ""
    if entry.get("genes"):
        gene = entry["genes"][0].get("geneName", {}).get("value", "")
    return {
        "gene": gene,
        "accession": entry.get("primaryAccession", ""),
        "entry_name": entry.get("uniProtkbId", ""),
        "protein": protein,
        "length": entry.get("sequence", {}).get("length"),
        "organism": entry.get("organism", {}).get("scientificName", ""),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0, help="只映射前 N 个（测试用）")
    parser.add_argument("--workers", type=int, default=4,
                        help="并发线程数（UniProt 公共 API 保持低并发）")
    args = parser.parse_args()

    names = collect_protein_names()
    print(f"待映射靶点蛋白名 {len(names)} 个")

    mapping: dict = {}
    if OUTPUT_FILE.exists():
        mapping = json.loads(OUTPUT_FILE.read_text(encoding="utf-8"))
        print(f"已缓存 {len(mapping)} 个")

    todo = [n for n in names if n not in mapping]
    if args.limit:
        todo = todo[:args.limit]
    print(f"本次查询 {len(todo)} 个（{args.workers} 线程并发）")

    t0 = time.time()
    n_ok = n_fail = 0
    lock = threading.Lock()
    done_count = 0

    def work(name: str) -> None:
        nonlocal n_ok, n_fail, done_count
        info = query_uniprot(name)
        with lock:
            if info is None:      # 网络失败：不缓存，下次重试
                n_fail += 1
            else:
                mapping[name] = info   # 空对象 = UniProt 无匹配（已确认）
                n_ok += 1
            done_count += 1
            if done_count % 25 == 0:
                OUTPUT_FILE.write_text(json.dumps(mapping, ensure_ascii=False),
                                       encoding="utf-8")
            if done_count % 50 == 0 or done_count == len(todo):
                gene_ok = sum(1 for v in mapping.values() if v.get("gene"))
                print(f"[{done_count}/{len(todo)}] 已映射 {gene_ok}/{len(mapping)}"
                      f"  ({time.time() - t0:.0f}s)", flush=True)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = []
        for i, name in enumerate(todo):
            futures.append(pool.submit(work, name))
            time.sleep(DELAY / args.workers)   # 整体保持约 1 请求/秒间隔
        for f in as_completed(futures):
            f.result()

    OUTPUT_FILE.write_text(json.dumps(mapping, ensure_ascii=False), encoding="utf-8")
    gene_ok = sum(1 for v in mapping.values() if v.get("gene"))
    print()
    print(f"完成：成功查询 {n_ok}，网络失败 {n_fail}")
    print(f"映射结果：{gene_ok}/{len(mapping)} 个蛋白名找到基因符号")
    print(f"已保存 {OUTPUT_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
