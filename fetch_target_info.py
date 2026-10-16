"""从 UniProt 补全靶点蛋白信息（蛋白名、登录号、亚细胞定位、功能）。

UniProt REST（免费、无需 key），按基因符号批量查询，限定人源且已审阅（Swiss-Prot）：
    https://rest.uniprot.org/uniprotkb/search?query=...&fields=...&format=json

结果写入 target_info.json，键为基因符号（大写）。

用法：
    python fetch_target_info.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List

import requests

HERE = Path(__file__).parent
NETWORK_FILE = HERE / "network.json"
OUTPUT_FILE = HERE / "target_info.json"

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = ("accession,id,protein_name,gene_names,organism_name,length,"
          "cc_subcellular_location,cc_function")
CHUNK = 25  # 每批发查询的基因数，避免 URL 过长


def gene_symbols() -> List[str]:
    """从 network.json 取出全部靶点基因符号。"""
    elements = json.loads(NETWORK_FILE.read_text(encoding="utf-8"))
    return sorted({n["data"]["name"] for n in elements["nodes"]
                   if n["data"]["label"] == "TARGET"})


def query_chunk(genes: List[str]) -> List[dict]:
    query = "(" + " OR ".join(f"gene_exact:{g}" for g in genes) + ")"
    params = {
        "query": f"{query} AND organism_id:9606 AND reviewed:true",
        "fields": FIELDS,
        "format": "json",
        "size": 500,
    }
    resp = requests.get(UNIPROT, params=params, timeout=40)
    if resp.status_code != 200:
        return []
    return resp.json().get("results", [])


def parse_entry(entry: dict) -> dict:
    protein = ""
    desc = entry.get("proteinDescription", {})
    if "recommendedName" in desc:
        protein = desc["recommendedName"]["fullName"]["value"]
    elif desc.get("submittedName"):
        protein = desc["submittedName"][0]["fullName"]["value"]

    gene = ""
    if entry.get("genes"):
        gene = entry["genes"][0].get("geneName", {}).get("value", "")

    locations, function = [], ""
    for comment in entry.get("comments", []):
        ctype = comment.get("commentType")
        if ctype == "SUBCELLULAR LOCATION":
            for loc in comment.get("subcellularLocations", []):
                value = loc.get("location", {}).get("value")
                if value:
                    locations.append(value)
                topo = loc.get("topology", {}).get("value")
                if topo:
                    locations.append(f"{value}({topo})" if value else topo)
        elif ctype == "FUNCTION" and not function:
            texts = comment.get("texts", [])
            if texts:
                function = texts[0].get("value", "")

    return {
        "accession": entry.get("primaryAccession", ""),
        "entry_name": entry.get("uniProtkbId", ""),
        "protein": protein,
        "gene": gene,
        "organism": entry.get("organism", {}).get("scientificName", ""),
        "length": entry.get("sequence", {}).get("length"),
        "locations": sorted(set(locations)),
        "function": function[:500] + ("…" if len(function) > 500 else ""),
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/"
                       f"{entry.get('primaryAccession', '')}/entry",
    }


def main() -> None:
    genes = gene_symbols()
    info: Dict[str, dict] = {}
    missing: List[str] = []

    for i in range(0, len(genes), CHUNK):
        chunk = genes[i:i + CHUNK]
        for entry in query_chunk(chunk):
            item = parse_entry(entry)
            if item["gene"]:
                info[item["gene"].upper()] = item
        time.sleep(0.2)

    for gene in genes:
        if gene.upper() not in info:
            missing.append(gene)

    OUTPUT_FILE.write_text(json.dumps(info, ensure_ascii=False, indent=2),
                           encoding="utf-8")
    print(f"已写入 {OUTPUT_FILE}：靶点 {len(genes)} 个，匹配到 {len(info)} 个")
    if missing:
        print("UniProt 中未匹配（多为自定义/非标准符号）：", "、".join(missing[:30]),
              f"（共 {len(missing)} 个）")


if __name__ == "__main__":
    main()
