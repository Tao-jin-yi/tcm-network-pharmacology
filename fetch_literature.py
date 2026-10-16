"""从 PubMed 抓取中药药理机制文献，并从摘要中提取机制关键词与涉及靶点。

数据源：NCBI E-utilities（免费，无需 key；请遵守 3 次/秒 的限速）。
    esearch → 取 PMID；efetch → 取题目、期刊、年份、摘要。

机制关键词采用词典匹配（通路 / 生物学过程），同时提取摘要里出现的本方靶点基因符号。

用法：
    python fetch_literature.py              # 抓取全部中药
    python fetch_literature.py --herb 丹参   # 只抓一味药
    python fetch_literature.py --max 30      # 每味药最多抓 30 篇
"""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from typing import Dict, List, Set

import requests

HERE = Path(__file__).parent
NETWORK_FILE = HERE / "network.json"
OUTPUT_FILE = HERE / "literature.json"

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# 中药名 →（拉丁学名 / 英文药材名）
HERB_QUERY = {
    "黄芪": "Astragalus membranaceus OR Astragali Radix OR Astragalus",
    "熟地黄": "Rehmannia glutinosa OR Rehmanniae Radix",
    "山茱萸": "Cornus officinalis OR Corni Fructus",
    "山药": "Dioscorea opposita OR Dioscoreae Rhizoma",
    "茯苓": "Poria cocos OR Poria",
    "泽泻": "Alisma orientale OR Alismatis Rhizoma",
    "牡丹皮": "Moutan Cortex OR Paeonia suffruticosa",
    "丹参": "Salvia miltiorrhiza OR Salviae Miltiorrhizae Radix OR Danshen",
    "甘草": "Glycyrrhiza uralensis OR Glycyrrhizae Radix OR licorice",
}

# 机制 / 通路词典：摘要中命中即记为该文献的机制
MECHANISMS = {
    "NF-κB 通路": [r"NF-?\s?kappa\s?B", r"NF-?κB", r"NFKB", r"IκB", r"IKK"],
    "PI3K/AKT 通路": [r"PI3K", r"phosphatidylinositol 3-kinase", r"\bAKT\b", r"PKB"],
    "MAPK 通路": [r"MAPK", r"\bERK\b", r"\bJNK\b", r"p38 MAPK", r"\bp38\b"],
    "JAK/STAT 通路": [r"JAK", r"STAT3?", r"\bSTAT\b"],
    "TLR4/NLRP3 炎症小体": [r"TLR4", r"NLRP3", r"inflammasome", r"Toll-like receptor"],
    "Nrf2/HO-1 抗氧化": [r"Nrf2", r"NFE2L2", r"HO-1", r"heme oxygenase", r"antioxidant"],
    "AMPK/mTOR 自噬": [r"AMPK", r"\bmTOR\b", r"autophagy", r"自噬"],
    "TGF-β/Smad 纤维化": [r"TGF-?β", r"TGF-beta", r"\bSmad\b", r"fibrosis"],
    "Wnt/β-catenin": [r"\bWnt\b", r"beta-?catenin", r"β-catenin"],
    "HIF-1/VEGF 血管新生": [r"HIF-?1", r"\bVEGF\b", r"angiogenesis"],
    "凋亡（线粒体途径）": [r"apoptosis", r"\bBax\b", r"Bcl-?2", r"caspase-?3?", r"凋亡"],
    "铁死亡": [r"ferroptosis", r"GPX4", r"铁死亡"],
    "焦亡": [r"pyroptosis", r"GSDMD", r"焦亡"],
    "氧化应激": [r"oxidative stress", r"\bROS\b", r"氧化应激"],
    "炎症反应": [r"inflammation", r"inflammatory", r"炎症"],
    "肠道菌群": [r"gut microbiota", r"intestinal flora", r"肠道菌群"],
    "免疫调节": [r"immunomodulat", r"immune regulation", r"免疫调节"],
    "糖脂代谢": [r"lipid metabolism", r"glucose metabolism", r"insulin resistance",
                 r"糖脂代谢", r"胰岛素抵抗"],
}


def esearch(query: str, retmax: int) -> List[str]:
    params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": retmax,
              "sort": "relevance"}
    resp = requests.get(ESEARCH, params=params, timeout=30)
    if resp.status_code != 200:
        return []
    return resp.json().get("esearchresult", {}).get("idlist", [])


def efetch(pmids: List[str]) -> List[dict]:
    if not pmids:
        return []
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    resp = requests.get(EFETCH, params=params, timeout=40)
    if resp.status_code != 200:
        return []

    articles = []
    for block in re.findall(r"<PubmedArticle>(.*?)</PubmedArticle>", resp.text, re.S):
        pmid = re.search(r"<PMID[^>]*>(\d+)</PMID>", block)
        title = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", block, re.S)
        journal = re.search(r"<ISOAbbreviation>(.*?)</ISOAbbreviation>", block, re.S)
        year = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", block, re.S)
        abstract_parts = re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", block, re.S)

        def clean(text: str) -> str:
            return re.sub(r"<[^>]+>", "", text).replace("&lt;", "<").replace("&gt;", ">").strip()

        articles.append({
            "pmid": pmid.group(1) if pmid else "",
            "title": clean(title.group(1)) if title else "",
            "journal": clean(journal.group(1)) if journal else "",
            "year": int(year.group(1)) if year else None,
            "abstract": " ".join(clean(p) for p in abstract_parts),
        })
    return articles


def extract_mechanisms(abstract: str) -> List[str]:
    found = []
    for name, patterns in MECHANISMS.items():
        for pattern in patterns:
            if re.search(pattern, abstract, re.I):
                found.append(name)
                break
    return found


# 摘要里常见的非标准写法 → 标准基因符号
SYNONYMS = [
    (r"TNF-?\s?α|TNF-?alpha|TNF-?a\b", "TNF"),
    (r"IL-?\s?6\b|interleukin-?6", "IL6"),
    (r"IL-?\s?1\s?β|IL-?\s?1\s?beta|IL-?\s?1B\b", "IL1B"),
    (r"IL-?\s?10\b", "IL10"),
    (r"IFN-?\s?γ|IFN-?\s?gamma", "IFNG"),
    (r"TGF-?\s?β|TGF-?beta", "TGFB1"),
    (r"MMP-?\s?9\b", "MMP9"),
    (r"MMP-?\s?2\b", "MMP2"),
    (r"COX-?\s?2\b", "PTGS2"),
    (r"iNOS|NOS2", "NOS2"),
    (r"caspase-?\s?3\b", "CASP3"),
    (r"caspase-?\s?9\b", "CASP9"),
    (r"Bcl-?\s?2\b", "BCL2"),
    (r"\bBax\b", "BAX"),
    (r"VEGF-?A\b", "VEGFA"),
    (r"CXCL-?8\b", "CXCL8"),
]


def normalize_symbols(abstract: str) -> str:
    """把摘要中的非标准写法替换为标准基因符号，便于后续匹配。"""
    for pattern, symbol in SYNONYMS:
        abstract = re.sub(pattern, symbol, abstract, flags=re.I)
    return abstract


def extract_targets(abstract: str, targets: Set[str]) -> List[str]:
    """提取摘要中出现的本方靶点（全大写且作为独立词出现）。"""
    text = normalize_symbols(abstract)
    return sorted({gene for gene in targets
                   if re.search(rf"\b{re.escape(gene)}\b", text)})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--herb", help="只抓取指定中药")
    parser.add_argument("--max", type=int, default=20, help="每味药最多抓取篇数")
    parser.add_argument("--reparse", action="store_true",
                        help="不重新抓取，仅用已有摘要重新提取机制与靶点")
    args = parser.parse_args()

    elements = json.loads(NETWORK_FILE.read_text(encoding="utf-8"))
    targets = {n["data"]["name"] for n in elements["nodes"]
               if n["data"]["label"] == "TARGET"}

    if args.reparse and OUTPUT_FILE.exists():
        stored = json.loads(OUTPUT_FILE.read_text(encoding="utf-8"))
        for herb, articles in stored.items():
            for article in articles:
                abstract = article.get("abstract", "")
                article["mechanisms"] = extract_mechanisms(abstract)
                article["targets"] = extract_targets(abstract, targets)
        OUTPUT_FILE.write_text(json.dumps(stored, ensure_ascii=False, indent=2),
                               encoding="utf-8")
        hits = sum(len(a.get("targets", [])) for arts in stored.values() for a in arts)
        print(f"已重新解析：{len(stored)} 味药，命中靶点 {hits} 次")
        return
    herbs = [h for h in (n["data"]["name"] for n in elements["nodes"]
                         if n["data"]["label"] == "HERB") if h in HERB_QUERY]
    if args.herb:
        herbs = [args.herb]

    result: Dict[str, list] = {}
    if OUTPUT_FILE.exists():
        result = json.loads(OUTPUT_FILE.read_text(encoding="utf-8"))

    for herb in herbs:
        query = (f'({HERB_QUERY[herb]}) AND '
                 f'(pharmacology[Title/Abstract] OR mechanism[Title/Abstract] OR '
                 f'pharmacological[Title/Abstract])')
        pmids = esearch(query, args.max)
        articles = efetch(pmids)
        for article in articles:
            abstract = article["abstract"]
            article["mechanisms"] = extract_mechanisms(abstract)
            article["targets"] = extract_targets(abstract, targets)
            article["url"] = f"https://pubmed.ncbi.nlm.nih.gov/{article['pmid']}/"
        result[herb] = articles
        print(f"{herb}: {len(articles)} 篇")
        time.sleep(0.4)

    OUTPUT_FILE.write_text(json.dumps(result, ensure_ascii=False, indent=2),
                           encoding="utf-8")
    total = sum(len(v) for v in result.values())
    print(f"已写入 {OUTPUT_FILE}：{len(result)} 味药，共 {total} 篇文献")


if __name__ == "__main__":
    main()
