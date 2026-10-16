"""作用通路库：通路图（配体 → 受体 → 激酶 → 转录因子 → 效应）与靶点归类。

节点类型 kind：
    LIGAND 配体 / RECEPTOR 受体 / KINASE 激酶 / ADAPTOR 接头蛋白 /
    TF 转录因子 / EFFECT 生物学效应
节点带 gene 字段时，可与本方靶点（基因符号）对接。
"""

from __future__ import annotations

from typing import Dict, List, Optional

# ------------------------------------------------------------------ 通路定义
def _n(node_id: str, label: str, kind: str, gene: str = "") -> dict:
    return {"id": node_id, "label": label, "kind": kind, "gene": gene}


PATHWAYS: Dict[str, dict] = {
    "NF-κB 通路": {
        "desc": "促炎因子结合受体后激活 IKK，降解 IκB，释放 NF-κB 入核，诱导 TNF/IL6/PTGS2 等炎症介质。",
        "nodes": [
            _n("nf_01", "TNF-α", "LIGAND", "TNF"),
            _n("nf_02", "IL-1β", "LIGAND", "IL1B"),
            _n("nf_03", "TNFR / IL-1R", "RECEPTOR"),
            _n("nf_04", "MyD88 / TRADD", "ADAPTOR", "MYD88"),
            _n("nf_05", "IKK 复合物", "KINASE", "IKBKB"),
            _n("nf_06", "IκBα 降解", "ADAPTOR", "NFKBIA"),
            _n("nf_07", "NF-κB (p65)", "TF", "RELA"),
            _n("nf_08", "炎症因子释放", "EFFECT", "IL6"),
            _n("nf_09", "PTGS2 / COX-2", "EFFECT", "PTGS2"),
        ],
        "edges": [("nf_01", "nf_03", "结合"), ("nf_02", "nf_03", "结合"),
                  ("nf_03", "nf_04", "招募"), ("nf_04", "nf_05", "激活"),
                  ("nf_05", "nf_06", "磷酸化"), ("nf_06", "nf_07", "释放"),
                  ("nf_07", "nf_08", "转录诱导"), ("nf_07", "nf_09", "转录诱导")],
    },
    "PI3K/AKT 通路": {
        "desc": "生长因子激活受体酪氨酸激酶 → PI3K → AKT → mTOR，调控存活、增殖与代谢；PTEN 为负调控。",
        "nodes": [
            _n("pi_01", "生长因子 (IGF/EGF)", "LIGAND", "VEGFA"),
            _n("pi_02", "受体酪氨酸激酶", "RECEPTOR", "EGFR"),
            _n("pi_03", "PI3K", "KINASE", "PIK3CA"),
            _n("pi_04", "AKT1", "KINASE", "AKT1"),
            _n("pi_05", "PTEN (抑制)", "ADAPTOR", "PTEN"),
            _n("pi_06", "mTOR", "KINASE", "MTOR"),
            _n("pi_07", "细胞存活 / 增殖", "EFFECT", "MYC"),
            _n("pi_08", "凋亡抑制 (BCL2)", "EFFECT", "BCL2"),
        ],
        "edges": [("pi_01", "pi_02", "结合"), ("pi_02", "pi_03", "激活"),
                  ("pi_03", "pi_04", "激活"), ("pi_05", "pi_04", "抑制"),
                  ("pi_04", "pi_06", "激活"), ("pi_06", "pi_07", "促进"),
                  ("pi_04", "pi_08", "促进")],
    },
    "MAPK 通路": {
        "desc": "RAS-RAF-MEK-ERK 级联把胞外刺激转成核内转录反应，调控增殖与分化；p38/JNK 支线对应激与炎症更敏感。",
        "nodes": [
            _n("mk_01", "胞外刺激 / 应激", "LIGAND"),
            _n("mk_02", "受体", "RECEPTOR", "EGFR"),
            _n("mk_03", "RAS", "ADAPTOR", "HRAS"),
            _n("mk_04", "RAF", "KINASE", "RAF1"),
            _n("mk_05", "MEK", "KINASE", "MAP2K1"),
            _n("mk_06", "ERK / MAPK1", "KINASE", "MAPK1"),
            _n("mk_07", "p38 / JNK", "KINASE", "MAPK14"),
            _n("mk_08", "转录因子活化", "TF", "JUN"),
            _n("mk_09", "增殖 / 炎症反应", "EFFECT", "MYC"),
        ],
        "edges": [("mk_01", "mk_02", "刺激"), ("mk_02", "mk_03", "激活"),
                  ("mk_03", "mk_04", "激活"), ("mk_04", "mk_05", "磷酸化"),
                  ("mk_05", "mk_06", "磷酸化"), ("mk_01", "mk_07", "应激激活"),
                  ("mk_06", "mk_08", "入核"), ("mk_07", "mk_08", "入核"),
                  ("mk_08", "mk_09", "转录")],
    },
    "JAK/STAT 通路": {
        "desc": "IL-6 等细胞因子结合受体后激活 JAK，磷酸化 STAT3 入核，诱导急性期蛋白与抗凋亡基因。",
        "nodes": [
            _n("js_01", "IL-6", "LIGAND", "IL6"),
            _n("js_02", "IL-6R / gp130", "RECEPTOR", "IL6R"),
            _n("js_03", "JAK1/2", "KINASE", "JAK2"),
            _n("js_04", "STAT3", "TF", "STAT3"),
            _n("js_05", "SOCS3 (负反馈)", "ADAPTOR", "SOCS3"),
            _n("js_06", "急性期反应 / 增殖", "EFFECT", "MYC"),
            _n("js_07", "抗凋亡 (BCL2)", "EFFECT", "BCL2"),
        ],
        "edges": [("js_01", "js_02", "结合"), ("js_02", "js_03", "激活"),
                  ("js_03", "js_04", "磷酸化"), ("js_04", "js_05", "诱导"),
                  ("js_05", "js_03", "抑制"), ("js_04", "js_06", "转录"),
                  ("js_04", "js_07", "转录")],
    },
    "TLR4/NLRP3 炎症小体": {
        "desc": "LPS 经 TLR4  priming 上调 pro-IL-1β，NLRP3 炎症小体激活 CASP1，成熟 IL-1β 释放并可致焦亡。",
        "nodes": [
            _n("tl_01", "LPS / DAMP", "LIGAND"),
            _n("tl_02", "TLR4", "RECEPTOR", "TLR4"),
            _n("tl_03", "NF-κB 活化", "TF", "RELA"),
            _n("tl_04", "pro-IL-1β", "LIGAND", "IL1B"),
            _n("tl_05", "NLRP3 炎症小体", "ADAPTOR", "NLRP3"),
            _n("tl_06", "CASP1", "KINASE", "CASP1"),
            _n("tl_07", "IL-1β / IL-18 成熟", "LIGAND", "IL18"),
            _n("tl_08", "焦亡 / 炎症", "EFFECT", "GSDMD"),
        ],
        "edges": [("tl_01", "tl_02", "识别"), ("tl_02", "tl_03", "激活"),
                  ("tl_03", "tl_04", "转录"), ("tl_01", "tl_05", "激活"),
                  ("tl_05", "tl_06", "招募"), ("tl_06", "tl_07", "剪切成熟"),
                  ("tl_06", "tl_08", "执行")],
    },
    "Nrf2/HO-1 抗氧化": {
        "desc": "氧化应激使 KEAP1 释放 Nrf2，Nrf2 入核结合 ARE，诱导 HO-1、NQO1 等抗氧化酶。",
        "nodes": [
            _n("nr_01", "ROS / 氧化应激", "LIGAND"),
            _n("nr_02", "KEAP1", "ADAPTOR", "KEAP1"),
            _n("nr_03", "Nrf2 (NFE2L2)", "TF", "NFE2L2"),
            _n("nr_04", "ARE 元件", "ADAPTOR"),
            _n("nr_05", "HO-1 (HMOX1)", "EFFECT", "HMOX1"),
            _n("nr_06", "NQO1 / SOD / CAT", "EFFECT", "NQO1"),
            _n("nr_07", "抗氧化 / 细胞保护", "EFFECT", "CAT"),
        ],
        "edges": [("nr_01", "nr_02", "修饰"), ("nr_02", "nr_03", "释放"),
                  ("nr_03", "nr_04", "结合"), ("nr_04", "nr_05", "转录"),
                  ("nr_04", "nr_06", "转录"), ("nr_05", "nr_07", "效应"),
                  ("nr_06", "nr_07", "效应")],
    },
    "AMPK/mTOR 自噬": {
        "desc": "能量应激激活 AMPK，抑制 mTOR 并激活 ULK1，启动自噬流，参与代谢与炎症调控。",
        "nodes": [
            _n("am_01", "能量应激 / AMP↑", "LIGAND"),
            _n("am_02", "AMPK", "KINASE", "PRKAA1"),
            _n("am_03", "mTORC1", "KINASE", "MTOR"),
            _n("am_04", "ULK1 复合物", "KINASE", "ULK1"),
            _n("am_05", "自噬体形成 (LC3)", "EFFECT", "MAP1LC3B"),
            _n("am_06", "自噬流 / 代谢改善", "EFFECT"),
        ],
        "edges": [("am_01", "am_02", "激活"), ("am_02", "am_03", "抑制"),
                  ("am_03", "am_04", "抑制"), ("am_02", "am_04", "激活"),
                  ("am_04", "am_05", "启动"), ("am_05", "am_06", "完成")],
    },
    "TGF-β/Smad 纤维化": {
        "desc": "TGF-β 结合受体磷酸化 Smad2/3，与 Smad4 入核，驱动胶原等细胞外基质沉积，参与纤维化。",
        "nodes": [
            _n("tg_01", "TGF-β1", "LIGAND", "TGFB1"),
            _n("tg_02", "TGFBR1/2", "RECEPTOR", "TGFBR1"),
            _n("tg_03", "Smad2/3", "TF", "SMAD3"),
            _n("tg_04", "Smad4", "TF", "SMAD4"),
            _n("tg_05", "胶原 / COL1A1", "EFFECT", "COL1A1"),
            _n("tg_06", "MMP2 / MMP9", "EFFECT", "MMP9"),
            _n("tg_07", "组织纤维化", "EFFECT"),
        ],
        "edges": [("tg_01", "tg_02", "结合"), ("tg_02", "tg_03", "磷酸化"),
                  ("tg_03", "tg_04", "复合"), ("tg_04", "tg_05", "转录"),
                  ("tg_04", "tg_06", "转录"), ("tg_05", "tg_07", "沉积")],
    },
    "Wnt/β-catenin": {
        "desc": "Wnt 结合 Frizzled 抑制降解复合体，β-catenin 累积入核，与 TCF/LEF 启动 MYC、CCND1 转录。",
        "nodes": [
            _n("wn_01", "Wnt 配体", "LIGAND", "WNT3A"),
            _n("wn_02", "Frizzled / LRP", "RECEPTOR", "FZD1"),
            _n("wn_03", "β-catenin", "TF", "CTNNB1"),
            _n("wn_04", "GSK3β 降解复合体", "KINASE", "GSK3B"),
            _n("wn_05", "TCF / LEF", "TF", "TCF7L2"),
            _n("wn_06", "MYC / CCND1", "EFFECT", "MYC"),
        ],
        "edges": [("wn_01", "wn_02", "结合"), ("wn_02", "wn_04", "抑制"),
                  ("wn_04", "wn_03", "促进降解"), ("wn_03", "wn_05", "结合"),
                  ("wn_05", "wn_06", "转录")],
    },
    "HIF-1/VEGF 血管新生": {
        "desc": "缺氧稳定 HIF-1α，诱导 VEGF 与糖酵解基因，促进血管新生与代谢重编程。",
        "nodes": [
            _n("hi_01", "缺氧", "LIGAND"),
            _n("hi_02", "HIF-1α", "TF", "HIF1A"),
            _n("hi_03", "VEGF", "LIGAND", "VEGFA"),
            _n("hi_04", "VEGFR2", "RECEPTOR", "KDR"),
            _n("hi_05", "血管新生", "EFFECT"),
            _n("hi_06", "糖酵解重编程", "EFFECT", "SLC2A1"),
        ],
        "edges": [("hi_01", "hi_02", "稳定"), ("hi_02", "hi_03", "转录"),
                  ("hi_03", "hi_04", "结合"), ("hi_04", "hi_05", "促进"),
                  ("hi_02", "hi_06", "转录")],
    },
    "线粒体凋亡": {
        "desc": "应激使 BAX 上调 / BCL2 下调，线粒体释放细胞色素 c，启动 CASP9 → CASP3 级联，执行凋亡。",
        "nodes": [
            _n("ap_01", "细胞应激 / 损伤", "LIGAND"),
            _n("ap_02", "TP53", "TF", "TP53"),
            _n("ap_03", "BAX ↑", "EFFECT", "BAX"),
            _n("ap_04", "BCL2 ↓", "EFFECT", "BCL2"),
            _n("ap_05", "线粒体通透化", "EFFECT"),
            _n("ap_06", "CASP9", "KINASE", "CASP9"),
            _n("ap_07", "CASP3", "KINASE", "CASP3"),
            _n("ap_08", "细胞凋亡", "EFFECT"),
        ],
        "edges": [("ap_01", "ap_02", "激活"), ("ap_02", "ap_03", "转录"),
                  ("ap_02", "ap_04", "抑制"), ("ap_03", "ap_05", "促发"),
                  ("ap_04", "ap_05", "抑制"), ("ap_05", "ap_06", "激活"),
                  ("ap_06", "ap_07", "激活"), ("ap_07", "ap_08", "执行")],
    },
    "铁死亡": {
        "desc": "GPX4 失活或 System Xc⁻ 受抑 → 脂质过氧化物累积 → 铁依赖的细胞死亡；Nrf2 通路可拮抗。",
        "nodes": [
            _n("fe_01", "System Xc⁻ (SLC7A11)", "RECEPTOR", "SLC7A11"),
            _n("fe_02", "GPX4", "EFFECT", "GPX4"),
            _n("fe_03", "脂质过氧化物 ↑", "EFFECT", "ACSL4"),
            _n("fe_04", "铁离子蓄积", "EFFECT"),
            _n("fe_05", "铁死亡", "EFFECT"),
            _n("fe_06", "Nrf2 (拮抗)", "TF", "NFE2L2"),
        ],
        "edges": [("fe_01", "fe_02", "支持"), ("fe_02", "fe_03", "清除"),
                  ("fe_03", "fe_05", "驱动"), ("fe_04", "fe_05", "驱动"),
                  ("fe_06", "fe_02", "上调"), ("fe_06", "fe_03", "抑制")],
    },
}

KIND_STYLE = {
    "LIGAND": ("配体", "#FF7F3E"),
    "RECEPTOR": ("受体", "#E4572E"),
    "KINASE": ("激酶", "#2A629A"),
    "ADAPTOR": ("接头/调控", "#7D5BA6"),
    "TF": ("转录因子", "#17A398"),
    "EFFECT": ("生物学效应", "#3CB371"),
}


def pathways_of_target(gene: str) -> List[str]:
    """该基因符号参与的通路名列表。"""
    gene = gene.upper()
    return [name for name, pw in PATHWAYS.items()
            if any(n.get("gene", "").upper() == gene for n in pw["nodes"])]


def targets_in_pathway(name: str) -> List[str]:
    """通路里的节点对应的本方靶点基因符号。"""
    return [n["gene"] for n in PATHWAYS[name]["nodes"] if n.get("gene")]


def pathway_elements(name: str) -> dict:
    """把通路转成 st-link-analysis（Cytoscape.js）可用的 elements。"""
    pw = PATHWAYS[name]
    nodes = [{"data": {"id": n["id"], "label": "PATH_NODE", "name": n["label"],
                       "kind": n["kind"]}} for n in pw["nodes"]]
    edges = [{"data": {"id": f"{s}-{t}", "label": "PATH_STEP", "source": s,
                       "target": t, "name": label}}
             for s, t, label in pw["edges"]]
    return {"nodes": nodes, "edges": edges}


# --------------------------------------------------------------- 靶点功能归类
RECEPTOR_HINTS = [
    ("膜受体", ["receptor", "R ", "IL6R", "TLR", "TNFR", "EGFR", "VEGFR"]),
    ("核受体/转录因子", ["transcription factor", "nuclear receptor", "activator"]),
    ("蛋白激酶", ["kinase"]),
    ("分泌型配体/细胞因子", ["interleukin", "factor", "cytokine", "hormone"]),
    ("蛋白酶", ["protease", "caspase", "matrix metalloproteinase", "peptidase"]),
    ("转运蛋白/通道", ["channel", "transporter", "carrier"]),
    ("抗氧化酶/代谢酶", ["reductase", "oxidase", "synthase", "oxygenase",
                        "dehydrogenase", "peroxidase", "catalase"]),
]


def target_class(gene: str, info: Optional[dict]) -> str:
    """按 UniProt 蛋白名与亚细胞定位推断靶点类别（用于受体点位展示）。"""
    protein = (info or {}).get("protein", "") or ""
    locations = " ".join((info or {}).get("locations", []) or [])
    low = protein.lower()

    for label, hints in RECEPTOR_HINTS:
        if any(h.lower() in low for h in hints):
            return label
    if "Nucleus" in locations:
        return "核内调控因子"
    if "Cell membrane" in locations or "Membrane" in locations:
        return "膜定位蛋白"
    if "Secreted" in locations or "Extracellular" in locations:
        return "分泌型蛋白"
    if "Mitochondrion" in locations:
        return "线粒体蛋白"
    return "其他 / 未分类"
