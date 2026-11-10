"""中药成分数据库数据层：整合 TCMSP 成分数据与 PubChem 理化性质，生成跨库直达链接。

整合来源
- TCMSP：成分（OB/DL/所属中药/作用靶点）、MOL 编号 —— 来自 network.json；
- PubChem PUG REST：IUPAC 名、精确质量、TPSA、XLogP、氢键供受体、可旋转键、
  同义词、2D 结构图 —— 批量抓取后缓存到 pubchem_props.json，离线可用；
- 跨库直达：PubChem / TCMSP / BATMAN-TCM / ETCM / HERB / SymMap / TCMID / NPASS / ChEMBL。

本模块不依赖 streamlit，便于在脚本与测试中复用。
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from urllib.parse import quote

import requests

HERE = Path(__file__).parent
PROPS_CACHE = HERE / "pubchem_props.json"
PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PROP_FIELDS = ("IUPACName,ExactMass,TPSA,XLogP,HBondDonorCount,"
               "HBondAcceptorCount,RotatableBondCount,HeavyAtomCount")

# 主流中药成分数据库入口（在页面「数据来源与跨库检索」中展示）
EXTERNAL_DBS = [
    ("TCMSP", "https://www.tcmsp.cn/",
     "中药系统药理学数据库与分析平台（本库成分与 OB / DL 数据来源，按 MOL 编号检索）"),
    ("BATMAN-TCM 2.0", "http://bionet.ncpsb.org.cn/batman-tcm/",
     "中药成分-靶点预测与功能模块分析"),
    ("ETCM v2.0", "http://www.tcmip.cn/ETCM/",
     "中药百科、成分-靶点-通路整合数据库"),
    ("HERB", "https://herb.ac.cn/",
     "中药高维信息整合数据库（成分 / 靶点 / 疾病）"),
    ("SymMap", "https://www.symmap.org/",
     "中医证候-分子-靶点整合数据库"),
    ("TCMID", "http://www.medicineresources.org.cn/",
     "中药成分-药材-方剂-疾病整合数据库"),
    ("NPASS", "https://bidd.group/NPASS/",
     "天然产物活性来源数据库（含活性值）"),
    ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/",
     "化合物结构、理化性质、同义词（本库实时富化来源）"),
    ("ChEMBL", "https://www.ebi.ac.uk/chembl/",
     "药物化学与生物活性数据库"),
]


def fetch_pubchem_props(cids: list[int],
                        cache_path: Path = PROPS_CACHE) -> dict[str, dict]:
    """批量抓取 PubChem 计算属性，结果缓存到本地 JSON（离线可复用）。

    PubChem PUG REST 支持一次查询多个 CID（上限约 100），故 53 个成分
    只需一次请求；失败不影响应用（缺失字段显示为 —）。
    """
    cache: dict[str, dict] = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text(encoding="utf-8"))
        except Exception:
            cache = {}

    missing = [str(c) for c in cids if str(c) not in cache]
    if not missing:
        return cache

    for start in range(0, len(missing), 90):
        chunk = missing[start:start + 90]
        url = f"{PUG_BASE}/cid/{','.join(chunk)}/property/{PROP_FIELDS}/JSON"
        try:
            resp = requests.get(url, timeout=30)
        except Exception:
            break
        if resp.status_code == 200:
            for row in resp.json().get("PropertyTable", {}).get("Properties", []):
                cache[str(row["CID"])] = row
            cache_path.write_text(json.dumps(cache, ensure_ascii=False, indent=0),
                                  encoding="utf-8")
        time.sleep(0.3)  # 遵守 PubChem 速率限制（≤5 次/秒）
    return cache


def fetch_synonyms(cid: int) -> list[str]:
    """取 PubChem 同义词（贸易名、俗名等），失败返回空列表。"""
    try:
        resp = requests.get(f"{PUG_BASE}/cid/{cid}/synonyms/JSON", timeout=30)
        if resp.status_code != 200:
            return []
        info = resp.json().get("InformationList", {}).get("Information", [])
        if not info:
            return []
        return [s for s in info[0].get("Synonym", [])][:15]
    except Exception:
        return []


def lipinski_flags(mw, xlogp, hbd, hba) -> list[str]:
    """Lipinski 类药性五规则（Ro5）违规项；缺失或非数值时跳过该规则。"""

    def _num(value):
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    flags = []
    mw_n = _num(mw)
    if mw_n is not None and mw_n > 500:
        flags.append("分子量>500")
    xlogp_n = _num(xlogp)
    if xlogp_n is not None and xlogp_n > 5:
        flags.append("XLogP>5")
    hbd_n = _num(hbd)
    if hbd_n is not None and hbd_n > 5:
        flags.append("氢键供体>5")
    hba_n = _num(hba)
    if hba_n is not None and hba_n > 10:
        flags.append("氢键受体>10")
    return flags


def external_links(comp: dict) -> list[tuple[str, str]]:
    """按成分生成可靠的跨库直达链接（只生成确定有效的 URL）。"""
    links: list[tuple[str, str]] = []
    name_en = quote(str(comp.get("en") or comp.get("name") or ""))
    if comp.get("cid"):
        cid = comp["cid"]
        links.append(("PubChem 页面", f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"))
        links.append(("PubChem SDF", f"{PUG_BASE}/cid/{cid}/SDF"))
    if comp.get("inchikey"):
        links.append(("PubChem 按 InChIKey",
                      f"https://pubchem.ncbi.nlm.nih.gov/search/"
                      f"search.cgi?structure={comp['inchikey']}"))
    if name_en:
        links.append(("ChEMBL 检索",
                      f"https://www.ebi.ac.uk/chembl/search_results/?query={name_en}"))
        links.append(("Google Scholar 文献",
                      f"https://scholar.google.com/scholar?q={name_en}"))
        links.append(("PubMed 文献",
                      f"https://pubmed.ncbi.nlm.nih.gov/?term={name_en}"))
    return links


def compound_target_map(net: dict) -> dict[str, list[str]]:
    """从 TARGETS 边构建 成分ID → 靶点基因符号列表 的映射（全库口径）。"""
    mapping: dict[str, list[str]] = {}
    for edge in net.get("edges", []):
        data = edge.get("data", {})
        if data.get("label") == "TARGETS":
            targets = mapping.setdefault(data["source"], [])
            symbol = str(data["target"]).replace("TGT::", "")
            if symbol not in targets:
                targets.append(symbol)
    return mapping
