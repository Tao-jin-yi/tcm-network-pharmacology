"""为已抓取的 PubMed 文献补全直达链接（DOI / Europe PMC / 开放获取全文）。

Europe PMC REST 可用 PMID（EXT_ID）查询，返回 DOI、PMCID、是否开放获取（isOpenAccess）
与全文链接列表。本脚本把结果写回 literature.json，供应用直接生成「一键打开」按钮，
解决文献只有一个 PubMed 链接、点开后还要自己找全文的问题。

用法：
    python resolve_links.py            # 补全所有缺失链接的文献
    python resolve_links.py --force    # 全部重查
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict, Optional

import requests

HERE = Path(__file__).parent
LITERATURE_FILE = HERE / "literature.json"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def resolve(pmid: str) -> Dict[str, Optional[str]]:
    """按 PMID 查 Europe PMC，返回 doi / pmcid / 开放获取全文等链接。"""
    try:
        resp = requests.get(EPMC, params={"query": f"EXT_ID:{pmid}",
                                          "resultType": "core", "format": "json"},
                            timeout=25)
        if resp.status_code != 200:
            return {}
        results = resp.json().get("resultList", {}).get("result", [])
    except Exception:
        return {}
    if not results:
        return {}

    entry = results[0]
    urls = entry.get("fullTextUrlList", {}).get("fullTextUrl", [])
    free = [u for u in urls
            if u.get("availability") in ("Open access", "Free", "Open Access")]
    pdf = [u for u in free if u.get("documentStyle") == "pdf"]
    doi = entry.get("doi")
    pmcid = entry.get("pmcid")

    return {
        "doi": doi,
        "pmcid": pmcid,
        "is_oa": entry.get("isOpenAccess") == "Y",
        "doi_url": f"https://doi.org/{doi}" if doi else None,
        "epmc_url": f"https://europepmc.org/article/MED/{pmid}",
        "pmc_url": f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/" if pmcid else None,
        "fulltext_url": (pdf[0].get("url") if pdf
                         else (free[0].get("url") if free else None)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="已解析过的也重新查询")
    args = parser.parse_args()

    if not LITERATURE_FILE.exists():
        print("未找到 literature.json，请先运行 python fetch_literature.py")
        return

    data = json.loads(LITERATURE_FILE.read_text(encoding="utf-8"))
    todo = []
    for herb, articles in data.items():
        for article in articles:
            if args.force or "epmc_url" not in article:
                todo.append((herb, article))

    print(f"待解析 {len(todo)} 篇")
    done = oa = doi_hit = 0
    for i, (_, article) in enumerate(todo, 1):
        pmid = str(article.get("pmid") or "")
        if not pmid:
            continue
        info = resolve(pmid)
        article.update(info)
        if info:
            done += 1
            oa += 1 if info.get("is_oa") else 0
            doi_hit += 1 if info.get("doi") else 0
        if i % 20 == 0:
            print(f"  已处理 {i}/{len(todo)}")
        time.sleep(0.15)

    LITERATURE_FILE.write_text(json.dumps(data, ensure_ascii=False, indent=2),
                               encoding="utf-8")
    print(f"完成：成功 {done} 篇，其中有 DOI {doi_hit} 篇、开放获取 {oa} 篇")


if __name__ == "__main__":
    main()
