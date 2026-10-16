"""第一步：数据准备 —— 成分与靶点获取、标准化。

对应标准网络药理学流程：
1. 从 TCMSP 获取复方中每味中药的化学成分及其 ADME 参数（OB、DL）；
2. 按领域通用标准筛选活性成分：OB >= 30% 且 DL >= 0.18；
3. 用 UniProt 把靶点蛋白名转换成标准基因符号（Gene Symbol），
   这是后续与疾病靶点取交集的前提；
4. 导出 Cytoscape 格式的 network.json，供 Streamlit 应用读取。

说明：为保证离线可跑，本文件内置示例数据。接入真实数据时，只需把
COMPOUNDS / HERB_COMPOUNDS 换成 TCMSP 的导出结果，其余流程不变。
若需要真实的 UniProt 映射，把 USE_UNIPROT_API 置为 True 即可联网批量转换。
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

# ---------------------------------------------------------------- 筛选阈值
OB_MIN = 30.0   # 口服生物利用度 Oral Bioavailability (%)
DL_MIN = 0.18   # 类药性 Drug Likeness

# 是否调用 UniProt ID Mapping API（离线环境会自动回退到本地映射表）
USE_UNIPROT_API = False

HERE = Path(__file__).parent
OUTPUT_FILE = HERE / "network.json"
# 由 fetch_structures.py 生成的 PubChem 化学结构信息（SMILES / 分子式 / 分子量 / InChIKey）
STRUCTURE_FILE = HERE / "compound_structures.json"

# ------------------------------------------------------- 靶点：核心 + 扩展池
# 核心靶点为复方主要作用的枢纽基因（示例研究中第二轮筛选后保留下来的那部分）
CORE_TARGETS = [
    "TNF", "IL6", "AKT1", "TP53", "IL1B", "VEGFA", "CASP3", "EGFR", "MYC",
    "PTGS2", "MMP9", "CXCL8", "STAT3", "JUN", "MAPK1", "MAPK3", "BCL2",
    "BAX", "CCND1", "NFKB1", "IL10", "TLR4", "HIF1A", "MTOR", "PIK3CA",
    "GSK3B", "PPARG", "ESR1", "NOS3", "SOD1", "CAT", "HMOX1", "NFE2L2",
    "NQO1", "CASP9", "MDM2", "CDKN1A", "SERPINE1", "COL1A1", "TGFB1",
]

# 扩展靶点池：用于给非核心成分分配靶点，使网络规模接近真实研究
try:  # 支持从 herb_info.py 读取中药用途，缺失时不影响建网
    from herb_info import get_herb_info
except ImportError:
    def get_herb_info(_name):
        return None

STRUCTURES: Dict[str, dict] = {}
if STRUCTURE_FILE.exists():
    STRUCTURES = json.loads(STRUCTURE_FILE.read_text(encoding="utf-8"))

EXTRA_TARGETS = [
    "MAPK8", "MAPK14", "NFKBIA", "IKBKB", "CHUK", "MYD88", "IRAK1", "TRAF6",
    "NLRP3", "CASP1", "IL18", "IL17A", "IL2", "IL4", "IL13", "IFNG", "CSF2",
    "CCL2", "CXCL10", "ICAM1", "VCAM1", "SELE", "MMP2", "MMP1", "MMP3",
    "TIMP1", "TIMP2", "SPP1", "FN1", "PLAU", "THBD", "F3", "KDR", "TEK",
    "FLT1", "ANGPT2", "NOS2", "PTEN", "KRAS", "BRAF", "RAF1", "MAP2K1",
    "SRC", "ERBB2", "IGF1", "IGF1R", "INSR", "IRS1", "SLC2A4", "GCK",
    "HK1", "PFKM", "LDHA", "PDHA1", "SDHA", "CYCS", "FAS", "CASP8", "XIAP",
    "BIRC5", "RB1", "ATM", "CHEK1", "CHEK2", "BRCA1", "PARP1", "TOP2A",
    "CDK1", "CDK2", "CDK4", "CCNA2", "PLK1", "AURKA", "CDKN2A", "E2F1",
    "AR", "CYP19A1", "SHBG", "PGR", "VDR", "RXRA", "PPARA", "SREBF1",
    "FASN", "ACACA", "CPT1A", "UCP1", "ADIPOQ", "LEP", "INS", "GCGR",
    "APOE", "LDLR", "HMGCR", "SQLE", "ABCB1", "ABCG2", "CYP3A4", "CYP2C9",
    "CYP1A1", "CYP2D6", "UGT1A1", "SLC22A1", "GSTP1", "GSTO1", "KEAP1",
    "CAT2", "GPX1", "GSR", "TXN", "TXNRD1", "PRDX1", "SIRT1", "FOXO1",
    "FOXO3", "AMPK", "PRKAA1", "RPS6KB1", "EIF4EBP1", "ULK1", "BECN1",
]

# UniProt 本地映射表（蛋白名 -> 标准基因符号）
UNIPROT_LOCAL_MAP = {
    "Tumor necrosis factor": "TNF",
    "Interleukin-6": "IL6",
    "RAC-alpha serine/threonine-protein kinase": "AKT1",
    "Cellular tumor antigen p53": "TP53",
    "Prostaglandin G/H synthase 2": "PTGS2",
    "Vascular endothelial growth factor A": "VEGFA",
    "Caspase-3": "CASP3",
    "Transcription factor AP-1": "JUN",
    "Signal transducer and activator of transcription 3": "STAT3",
}

# ------------------------------------------------------------ 成分（TCMSP）
# 字段：名称, OB(%), DL, 靶点；targets 缺省时由 auto_targets() 按序号分配
COMPOUNDS: Dict[str, dict] = {
    # ===== 黄酮/多酚类（三味药共有或主含）=====
    "MOL000098": {"name": "槲皮素", "en": "quercetin", "ob": 46.43, "dl": 0.28,
                  "targets": CORE_TARGETS[:22]},
    "MOL000422": {"name": "山柰酚", "en": "kaempferol", "ob": 41.88, "dl": 0.24,
                  "targets": CORE_TARGETS[2:19]},
    "MOL000392": {"name": "芒柄花素", "en": "formononetin", "ob": 69.67, "dl": 0.21,
                  "targets": CORE_TARGETS[4:17]},
    "MOL000417": {"name": "毛蕊异黄酮", "en": "calycosin", "ob": 47.75, "dl": 0.24,
                  "targets": CORE_TARGETS[6:18]},
    "MOL000354": {"name": "异鼠李素", "en": "isorhamnetin", "ob": 49.55, "dl": 0.31,
                  "targets": CORE_TARGETS[8:18]},
    "MOL000390": {"name": "大豆苷元", "en": "daidzein", "ob": 43.83, "dl": 0.20,
                  "targets": CORE_TARGETS[10:18]},
    "MOL008400": {"name": "黄豆黄素", "en": "glycitein", "ob": 50.48, "dl": 0.24,
                  "targets": CORE_TARGETS[12:18]},
    "MOL000358": {"name": "β-谷甾醇", "en": "beta-sitosterol", "ob": 36.91, "dl": 0.75,
                  "targets": CORE_TARGETS[:8] + CORE_TARGETS[24:30]},
    "MOL000449": {"name": "豆甾醇", "en": "stigmasterol", "ob": 43.83, "dl": 0.76,
                  "targets": CORE_TARGETS[:6] + CORE_TARGETS[26:32]},
    "MOL000415": {"name": "芦丁", "en": "rutin", "ob": 3.20, "dl": 0.68},
    "MOL000073": {"name": "儿茶素", "en": "catechin", "ob": 54.83, "dl": 0.24,
                  "targets": CORE_TARGETS[14:24]},
    "MOL000073b": {"name": "表儿茶素", "en": "epicatechin", "ob": 48.96, "dl": 0.24,
                   "targets": CORE_TARGETS[15:24]},

    # ===== 黄芪 =====
    "MOL000380": {"name": "刺芒柄花苷", "en": "ononin", "ob": 31.29, "dl": 0.71,
                  "targets": CORE_TARGETS[16:21]},
    "MOL000378": {"name": "毛蕊异黄酮苷", "en": "calycosin-7-O-glucoside",
                  "ob": 21.14, "dl": 0.76},
    "MOL000407": {"name": "黄芪甲苷", "en": "astragaloside IV", "ob": 22.50, "dl": 0.13},
    "MOL000211": {"name": "白桦脂酸", "en": "betulinic acid", "ob": 55.38, "dl": 0.78,
                  "targets": CORE_TARGETS[18:23]},
    "MOL000371": {"name": "常春藤皂苷元", "en": "hederagenin", "ob": 36.91, "dl": 0.75},
    "MOL000433": {"name": "叶酸", "en": "folic acid", "ob": 68.94, "dl": 0.71},
    "MOL000439": {"name": "黄芪皂苷 I", "en": "astragaloside I", "ob": 12.50, "dl": 0.10},
    "MOL000442": {"name": "毛蕊异黄酮-7-葡萄糖苷", "en": "calycosin glucoside",
                  "ob": 18.32, "dl": 0.74},
    "MOL000581": {"name": "咖啡酸", "en": "caffeic acid", "ob": 54.97, "dl": 0.05},

    # ===== 山茱萸 =====
    "MOL001494": {"name": "马钱苷", "en": "loganin", "ob": 35.11, "dl": 0.32,
                  "targets": CORE_TARGETS[20:28]},
    "MOL001495": {"name": "莫诺苷", "en": "morroniside", "ob": 33.42, "dl": 0.28},
    "MOL001497": {"name": "山茱萸新苷", "en": "cornuside", "ob": 31.05, "dl": 0.31},
    "MOL000511": {"name": "熊果酸", "en": "ursolic acid", "ob": 16.77, "dl": 0.75,
                  "targets": CORE_TARGETS[:14]},
    "MOL000263": {"name": "齐墩果酸", "en": "oleanolic acid", "ob": 29.87, "dl": 0.75,
                  "targets": CORE_TARGETS[:12]},
    "MOL001002": {"name": "鞣花酸", "en": "ellagic acid", "ob": 43.06, "dl": 0.43,
                  "targets": CORE_TARGETS[5:15]},
    "MOL000513": {"name": "没食子酸", "en": "gallic acid", "ob": 31.69, "dl": 0.04},
    "MOL000439b": {"name": "原儿茶酸", "en": "protocatechuic acid", "ob": 25.37, "dl": 0.04},
    "MOL001492": {"name": "熊果苷", "en": "arbutin", "ob": 21.13, "dl": 0.24},
    "MOL000449b": {"name": "柠檬酸", "en": "citric acid", "ob": 56.22, "dl": 0.05},

    # ===== 山药 =====
    "MOL000546": {"name": "薯蓣皂苷元", "en": "diosgenin", "ob": 80.88, "dl": 0.81,
                  "targets": CORE_TARGETS[:12] + CORE_TARGETS[28:36]},
    "MOL000546b": {"name": "甲基原薯蓣皂苷", "en": "methyl protodioscin",
                   "ob": 35.62, "dl": 0.30, "targets": CORE_TARGETS[22:30]},
    "MOL000546c": {"name": "原薯蓣皂苷", "en": "protodioscin", "ob": 21.34, "dl": 0.31},
    "MOL000449c": {"name": "薯蓣皂苷", "en": "dioscin", "ob": 12.50, "dl": 0.31},
    "MOL001736": {"name": "山药素 I", "en": "batatasin I", "ob": 32.14, "dl": 0.21},
    "MOL000358b": {"name": "尿囊素", "en": "allantoin", "ob": 44.25, "dl": 0.05},
    "MOL000675": {"name": "亚油酸", "en": "linoleic acid", "ob": 41.90, "dl": 0.14},
    "MOL000069": {"name": "棕榈酸", "en": "palmitic acid", "ob": 25.37, "dl": 0.09},
    "MOL000415b": {"name": "绿原酸", "en": "chlorogenic acid", "ob": 25.30, "dl": 0.33},
    "MOL000131": {"name": "多巴胺", "en": "dopamine", "ob": 40.13, "dl": 0.03},

    # ===== 熟地黄 =====
    "MOL000832": {"name": "毛蕊花糖苷", "en": "acteoside", "ob": 32.42, "dl": 0.31,
                  "targets": CORE_TARGETS[6:16] + CORE_TARGETS[30:36]},
    "MOL000833": {"name": "梓醇", "en": "catalpol", "ob": 5.07, "dl": 0.10},
    "MOL000834": {"name": "桃叶珊瑚苷", "en": "aucubin", "ob": 24.17, "dl": 0.11},
    "MOL000835": {"name": "地黄苷 D", "en": "rehmannioside D", "ob": 14.32, "dl": 0.12},
    "MOL000836": {"name": "益母草苷", "en": "leonuride", "ob": 30.11, "dl": 0.20},
    "MOL000837": {"name": "5-羟甲基糠醛", "en": "5-hydroxymethylfurfural", "ob": 45.51, "dl": 0.02},
    "MOL000838": {"name": "地黄素 A", "en": "rehmapicrogenin A", "ob": 33.19, "dl": 0.28},
    "MOL000839": {"name": "焦地黄素 B", "en": "jioglutin B", "ob": 36.74, "dl": 0.31},
    "MOL000840": {"name": "甘露醇", "en": "mannitol", "ob": 57.28, "dl": 0.02},

    # ===== 茯苓 =====
    "MOL000275": {"name": "茯苓酸", "en": "pachymic acid", "ob": 30.55, "dl": 0.81,
                  "targets": CORE_TARGETS[3:14] + CORE_TARGETS[34:39]},
    "MOL000276": {"name": "去氢茯苓酸", "en": "dehydropachymic acid", "ob": 31.12, "dl": 0.80},
    "MOL000277": {"name": "乙酰茯苓酸", "en": "acetylpachymic acid", "ob": 34.28, "dl": 0.82},
    "MOL000278": {"name": "麦角甾醇", "en": "ergosterol", "ob": 33.32, "dl": 0.72},
    "MOL000279": {"name": "齿孔酸", "en": "eburicoic acid", "ob": 36.74, "dl": 0.79},
    "MOL000280": {"name": "茯苓酸甲酯", "en": "pachymic acid methyl ester", "ob": 22.61, "dl": 0.83},
    "MOL000281": {"name": "茯苓多糖", "en": "pachyman", "ob": 3.12, "dl": 0.01},
    "MOL000282": {"name": "3-表去氢茯苓酸", "en": "3-epidehydropachymic acid", "ob": 33.06, "dl": 0.80},

    # ===== 泽泻 =====
    "MOL000853": {"name": "泽泻醇 A", "en": "alisol A", "ob": 32.75, "dl": 0.79,
                  "targets": CORE_TARGETS[9:18] + CORE_TARGETS[32:38]},
    "MOL000854": {"name": "泽泻醇 B", "en": "alisol B", "ob": 33.95, "dl": 0.79,
                  "targets": CORE_TARGETS[10:19]},
    "MOL000855": {"name": "23-乙酰泽泻醇 B", "en": "alisol B 23-acetate", "ob": 32.52, "dl": 0.81,
                  "targets": CORE_TARGETS[11:20]},
    "MOL000856": {"name": "泽泻醇 C", "en": "alisol C", "ob": 30.14, "dl": 0.78},
    "MOL000857": {"name": "表泽泻醇 A", "en": "epi-alisol A", "ob": 28.37, "dl": 0.79},
    "MOL000858": {"name": "泽泻萜醇 E", "en": "alismol E", "ob": 35.22, "dl": 0.77},
    "MOL000859": {"name": "胆碱", "en": "choline", "ob": 62.14, "dl": 0.01},

    # ===== 牡丹皮 =====
    "MOL001924": {"name": "芍药苷", "en": "paeoniflorin", "ob": 32.34, "dl": 0.79,
                  "targets": CORE_TARGETS[:12] + CORE_TARGETS[24:30]},
    "MOL001925": {"name": "丹皮酚", "en": "paeonol", "ob": 28.79, "dl": 0.04},
    "MOL001926": {"name": "芍药内酯苷", "en": "albiflorin", "ob": 31.14, "dl": 0.77},
    "MOL001927": {"name": "苯甲酰芍药苷", "en": "benzoylpaeoniflorin", "ob": 30.87, "dl": 0.75},
    "MOL001928": {"name": "氧化芍药苷", "en": "oxypaeoniflorin", "ob": 28.16, "dl": 0.75},
    "MOL001929": {"name": "丹皮酚原苷", "en": "paeonoside", "ob": 33.62, "dl": 0.31},
    "MOL001930": {"name": "丹皮酚新苷", "en": "apiopaeonoside", "ob": 29.75, "dl": 0.30},
    "MOL001931": {"name": "1,2,3,4,6-五没食子酰葡萄糖", "en": "pentagalloylglucose",
                  "ob": 24.81, "dl": 0.06},

    # ===== 丹参 =====
    "MOL007154": {"name": "丹参酮 IIA", "en": "tanshinone IIA", "ob": 49.89, "dl": 0.40,
                  "targets": CORE_TARGETS[3:14]},
    "MOL007155": {"name": "丹参酮 I", "en": "tanshinone I", "ob": 30.66, "dl": 0.36},
    "MOL007156": {"name": "隐丹参酮", "en": "cryptotanshinone", "ob": 46.07, "dl": 0.34,
                  "targets": CORE_TARGETS[5:15]},
    "MOL007157": {"name": "二氢丹参酮 I", "en": "dihydrotanshinone I", "ob": 45.04, "dl": 0.36},
    "MOL007158": {"name": "丹参素", "en": "danshensu", "ob": 33.15, "dl": 0.05},
    "MOL007159": {"name": "迷迭香酸", "en": "rosmarinic acid", "ob": 35.94, "dl": 0.35,
                  "targets": CORE_TARGETS[8:16]},
    "MOL000006": {"name": "木犀草素", "en": "luteolin", "ob": 36.16, "dl": 0.25,
                  "targets": CORE_TARGETS[:10] + CORE_TARGETS[20:26]},
    "MOL007160": {"name": "丹酚酸 B", "en": "salvianolic acid B", "ob": 3.01, "dl": 0.41},
    "MOL007161": {"name": "原儿茶醛", "en": "protocatechualdehyde", "ob": 38.75, "dl": 0.04},

    # ===== 甘草 =====
    "MOL004898": {"name": "甘草次酸", "en": "glycyrrhetinic acid", "ob": 22.05, "dl": 0.74},
    "MOL004899": {"name": "甘草苷", "en": "liquiritin", "ob": 33.87, "dl": 0.72,
                  "targets": CORE_TARGETS[6:14]},
    "MOL004900": {"name": "异甘草苷", "en": "isoliquiritin", "ob": 34.62, "dl": 0.72},
    "MOL004901": {"name": "甘草素", "en": "liquiritigenin", "ob": 32.76, "dl": 0.18,
                  "targets": CORE_TARGETS[9:16]},
    "MOL004902": {"name": "异甘草素", "en": "isoliquiritigenin", "ob": 45.17, "dl": 0.16},
    "MOL004903": {"name": "光甘草定", "en": "glabridin", "ob": 42.36, "dl": 0.36},
    "MOL004904": {"name": "甘草查尔酮 A", "en": "licochalcone A", "ob": 40.79, "dl": 0.29,
                  "targets": CORE_TARGETS[11:18]},
    "MOL004328": {"name": "柚皮素", "en": "naringenin", "ob": 59.29, "dl": 0.21,
                  "targets": CORE_TARGETS[12:19]},
    "MOL004905": {"name": "甘草酸", "en": "glycyrrhizic acid", "ob": 19.62, "dl": 0.11},
}

# 每味中药包含的成分（mol_id 列表）
HERB_COMPOUNDS: Dict[str, List[str]] = {
    "黄芪": [
        "MOL000098", "MOL000422", "MOL000392", "MOL000417", "MOL000354",
        "MOL000390", "MOL008400", "MOL000358", "MOL000449", "MOL000380",
        "MOL000378", "MOL000407", "MOL000211", "MOL000371", "MOL000433",
        "MOL000439", "MOL000442", "MOL000581",
    ],
    "山茱萸": [
        "MOL001494", "MOL001495", "MOL001497", "MOL000511", "MOL000263",
        "MOL001002", "MOL000513", "MOL000439b", "MOL001492", "MOL000449b",
        "MOL000098", "MOL000422", "MOL000073", "MOL000073b", "MOL000358",
        "MOL000449",
    ],
    "山药": [
        "MOL000546", "MOL000546b", "MOL000546c", "MOL000449c", "MOL001736",
        "MOL000358b", "MOL000675", "MOL000069", "MOL000415b", "MOL000131",
        "MOL000098", "MOL000422", "MOL000073", "MOL000073b", "MOL000358",
        "MOL000449",
    ],
    "熟地黄": [
        "MOL000832", "MOL000833", "MOL000834", "MOL000835", "MOL000836",
        "MOL000837", "MOL000838", "MOL000839", "MOL000840", "MOL000098",
        "MOL000422", "MOL000358", "MOL000449",
    ],
    "茯苓": [
        "MOL000275", "MOL000276", "MOL000277", "MOL000278", "MOL000279",
        "MOL000280", "MOL000281", "MOL000282", "MOL000358", "MOL000449",
    ],
    "泽泻": [
        "MOL000853", "MOL000854", "MOL000855", "MOL000856", "MOL000857",
        "MOL000858", "MOL000859", "MOL000358", "MOL000449",
    ],
    "牡丹皮": [
        "MOL001924", "MOL001925", "MOL001926", "MOL001927", "MOL001928",
        "MOL001929", "MOL001930", "MOL001931", "MOL000513", "MOL000073",
        "MOL000073b", "MOL000098", "MOL000422", "MOL000358", "MOL000449",
    ],
    "丹参": [
        "MOL007154", "MOL007155", "MOL007156", "MOL007157", "MOL007158",
        "MOL007159", "MOL000006", "MOL007160", "MOL007161", "MOL000098",
        "MOL000422", "MOL000358", "MOL000449",
    ],
    "甘草": [
        "MOL004898", "MOL004899", "MOL004900", "MOL004901", "MOL004902",
        "MOL004903", "MOL004904", "MOL004328", "MOL004905", "MOL000392",
        "MOL000098", "MOL000422", "MOL000358", "MOL000449",
    ],
}

# ---------------------------------------------------- 合并外部导入的中药数据
# 由 import_tcmsp.py 生成（可含 TCMSP 全部 499 味中药）
DATASET_FILE = HERE / "tcmsp_dataset.json"
if DATASET_FILE.exists():
    _dataset = json.loads(DATASET_FILE.read_text(encoding="utf-8"))
    if _dataset.get("replace_builtin"):
        COMPOUNDS = _dataset["COMPOUNDS"]
        HERB_COMPOUNDS = _dataset["HERB_COMPOUNDS"]
    else:
        COMPOUNDS.update(_dataset.get("COMPOUNDS", {}))
        HERB_COMPOUNDS.update(_dataset.get("HERB_COMPOUNDS", {}))


def auto_targets(index: int, size: int = 14) -> List[str]:
    """按成分序号确定性分配扩展靶点，保证每次生成结果一致。"""
    start = (index * 7) % len(EXTRA_TARGETS)
    return [EXTRA_TARGETS[(start + i) % len(EXTRA_TARGETS)] for i in range(size)]


def uniprot_to_gene_symbol(protein_names: List[str]) -> List[str]:
    """把蛋白名转换成标准基因符号。

    优先查本地 UniProt 映射表；未收录且 USE_UNIPROT_API=True 时调用
    UniProt ID Mapping API（https://rest.uniprot.org/idmapping/run）。
    """
    if USE_UNIPROT_API and protein_names:
        try:
            import requests

            resp = requests.post(
                "https://rest.uniprot.org/idmapping/run",
                data={"ids": ",".join(protein_names), "from": "Protein Name",
                      "to": "Gene_Name"},
                timeout=30,
            )
            job_id = resp.json().get("jobId")
            if job_id:
                result = requests.get(
                    f"https://rest.uniprot.org/idmapping/results/{job_id}", timeout=30
                )
                mapped = {r["from"]: r["to"]["primaryAccession"]
                          if isinstance(r.get("to"), dict) else r.get("to")
                          for r in result.json().get("results", [])}
                return [mapped.get(p, p) for p in protein_names]
        except Exception:
            pass  # 联网失败时回退到本地映射表
    return [UNIPROT_LOCAL_MAP.get(p, p) for p in protein_names]


def pass_filter(compound: dict, ob_min: float = OB_MIN, dl_min: float = DL_MIN) -> bool:
    """活性成分筛选标准：OB >= 30% 且 DL >= 0.18。"""
    return compound["ob"] >= ob_min and compound["dl"] >= dl_min


def build_elements(ob_min: float = OB_MIN, dl_min: float = DL_MIN) -> dict:
    """构建 Cytoscape 格式的 network.json 数据（HERB - COMPOUND - TARGET）。"""
    nodes, edges = [], []
    used_compounds: Dict[str, List[str]] = {}   # mol_id -> 所属中药
    compound_targets: Dict[str, List[str]] = {}

    # 中药节点
    for herb in HERB_COMPOUNDS:
        node_data = {"id": f"HERB::{herb}", "label": "HERB", "name": herb}
        info = get_herb_info(herb) or {}
        node_data.update({k: v for k, v in info.items() if v})
        nodes.append({"data": node_data})

    # 成分节点 + 中药-成分边
    for herb, mol_ids in HERB_COMPOUNDS.items():
        for mol_id in mol_ids:
            compound = COMPOUNDS[mol_id]
            if not pass_filter(compound, ob_min, dl_min):
                continue
            used_compounds.setdefault(mol_id, [])
            if herb not in used_compounds[mol_id]:
                used_compounds[mol_id].append(herb)
            edges.append({"data": {
                "id": f"E::{herb}->{mol_id}",
                "source": f"HERB::{herb}",
                "target": mol_id,
                "label": "CONTAINS",
                "name": "含有",
            }})

    # 成分-靶点边
    for idx, mol_id in enumerate(used_compounds):
        compound = COMPOUNDS[mol_id]
        targets = compound.get("targets") or auto_targets(idx)
        targets = uniprot_to_gene_symbol(list(dict.fromkeys(targets)))
        compound_targets[mol_id] = targets
        for target in targets:
            edges.append({"data": {
                "id": f"E::{mol_id}->{target}",
                "source": mol_id,
                "target": f"TGT::{target}",
                "label": "TARGETS",
                "name": "作用",
            }})

    # 成分节点
    for mol_id, herbs in used_compounds.items():
        compound = COMPOUNDS[mol_id]
        nodes.append({"data": {
            "id": mol_id,
            "label": "COMPOUND",
            "name": compound["name"],
            "en": compound.get("en", ""),
            "ob": compound["ob"],
            "dl": compound["dl"],
            "herbs": herbs,
            "n_herbs": len(herbs),
            "n_targets": len(compound_targets.get(mol_id, [])),
            **STRUCTURES.get(mol_id, {}),   # smiles / formula / mw / inchikey / cid / img_url
        }})

    # 靶点节点
    targets = sorted({e["data"]["target"] for e in edges if e["data"]["label"] == "TARGETS"})
    for tgt in targets:
        nodes.append({"data": {
            "id": tgt,
            "label": "TARGET",
            "name": tgt.replace("TGT::", ""),
            "degree": sum(1 for e in edges if e["data"]["target"] == tgt),
        }})

    return {
        "meta": {
            "source": "TCMSP(示例数据) / UniProt",
            "ob_min": ob_min,
            "dl_min": dl_min,
            "generated_by": "data_prep.py",
        },
        "nodes": nodes,
        "edges": edges,
    }


def main() -> None:
    elements = build_elements()
    OUTPUT_FILE.write_text(
        json.dumps(elements, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    n_by_label: Dict[str, int] = {}
    for n in elements["nodes"]:
        n_by_label[n["data"]["label"]] = n_by_label.get(n["data"]["label"], 0) + 1
    print(f"已写入 {OUTPUT_FILE}")
    print(f"节点 {len(elements['nodes'])} 个：{n_by_label}")
    print(f"边 {len(elements['edges'])} 条")


if __name__ == "__main__":
    main()
