"""为《中国药典》/《本草纲目》目录药材生成演示用的成分-靶点关联数据。

重要说明（请务必保留此声明）：
- 本模块生成的"药材-成分"关联属于**演示数据**，不是 TCMSP 实测结果；
- 所用成分均为真实存在的天然产物，OB / DL 取文献报道量级；
- 某味药材"含有"哪些成分由确定性算法分配（同一药材每次结果一致），
  目的是让 302 味目录药材能完整走通活性成分 / 靶点 / 配伍协同 /
  核心靶点 / 靶点详情 / 药动学 / 动物模拟等分析流程；
- 正式科研分析请通过 ``import_tcmsp.py`` 导入 TCMSP 导出表（--replace），
  再运行 ``data_prep.py``，届时以真实数据覆盖本演示数据。

被 data_prep.py 调用：expand_catalog_herbs(COMPOUNDS, HERB_COMPOUNDS)
"""

from __future__ import annotations

import hashlib
from typing import Dict, List

from pharmacopoeia import HERB_CATALOG

# 补充的常见天然产物（真实化合物，OB ≥ 30 且 DL ≥ 0.18，取文献报道量级）
DEMO_COMPOUNDS: Dict[str, dict] = {
    "SIM001": {"name": "小檗碱", "en": "berberine", "ob": 36.86, "dl": 0.78},
    "SIM002": {"name": "黄芩苷", "en": "baicalin", "ob": 33.67, "dl": 0.73},
    "SIM003": {"name": "黄芩素", "en": "baicalein", "ob": 33.52, "dl": 0.21},
    "SIM004": {"name": "汉黄芩素", "en": "wogonin", "ob": 30.68, "dl": 0.37},
    "SIM005": {"name": "姜黄素", "en": "curcumin", "ob": 46.36, "dl": 0.26},
    "SIM006": {"name": "胡椒碱", "en": "piperine", "ob": 42.52, "dl": 0.35},
    "SIM007": {"name": "欧前胡素", "en": "imperatorin", "ob": 40.10, "dl": 0.30},
    "SIM008": {"name": "异欧前胡素", "en": "isoimperatorin", "ob": 33.55, "dl": 0.30},
    "SIM009": {"name": "花椒毒素", "en": "xanthotoxin", "ob": 43.79, "dl": 0.42},
    "SIM010": {"name": "五味子醇甲", "en": "schisandrin", "ob": 55.46, "dl": 0.66},
    "SIM011": {"name": "五味子乙素", "en": "schisandrin B", "ob": 30.62, "dl": 0.66},
    "SIM012": {"name": "苦参碱", "en": "matrine", "ob": 63.77, "dl": 0.25},
    "SIM013": {"name": "氧化苦参碱", "en": "oxymatrine", "ob": 30.12, "dl": 0.18},
    "SIM014": {"name": "厚朴酚", "en": "magnolol", "ob": 30.05, "dl": 0.29},
    "SIM015": {"name": "吴茱萸碱", "en": "evodiamine", "ob": 45.10, "dl": 0.40},
    "SIM016": {"name": "吴茱萸次碱", "en": "rutaecarpine", "ob": 41.66, "dl": 0.32},
    "SIM017": {"name": "藁本内酯", "en": "ligustilide", "ob": 45.70, "dl": 0.30},
    "SIM018": {"name": "淫羊藿苷", "en": "icariin", "ob": 41.58, "dl": 0.61},
    "SIM019": {"name": "香叶木素", "en": "diosmetin", "ob": 31.14, "dl": 0.27},
}

MIN_COMPOUNDS = 6   # 每味药材最少分配成分数
MAX_COMPOUNDS = 10  # 每味药材最多分配成分数


def _seed(name: str) -> int:
    """由药材名生成稳定的整数种子（同一药材结果恒定）。"""
    return int(hashlib.md5(name.encode("utf-8")).hexdigest()[:8], 16)


def expand_catalog_herbs(compounds: Dict[str, dict],
                         herb_compounds: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """给尚无成分数据的目录药材分配演示成分，原地更新两个字典。

    可分配成分池 = 所有满足 OB ≥ 30%、DL ≥ 0.18 的成分（含 TCMSP 内置成分
    与新补充的天然产物）。每味药材取 6–10 个，按药材名哈希决定起点，
    步长 13 与池长度互质，保证同一味药内不重复且不同药材间有差异。

    返回本次新增的 {药材: [成分ID, ...]}。
    """
    for mol_id, comp in DEMO_COMPOUNDS.items():
        compounds.setdefault(mol_id, dict(comp))

    pool = sorted(mid for mid, c in compounds.items()
                  if c["ob"] >= 30.0 and c["dl"] >= 0.18)
    if not pool:
        return {}

    added: Dict[str, List[str]] = {}
    for herb in HERB_CATALOG:
        if herb in herb_compounds:
            continue  # 已有真实 TCMSP 数据，不覆盖
        seed = _seed(herb)
        n = MIN_COMPOUNDS + seed % (MAX_COMPOUNDS - MIN_COMPOUNDS + 1)
        start = seed % len(pool)
        picks: List[str] = []
        for j in range(n):
            mol_id = pool[(start + j * 13) % len(pool)]
            if mol_id not in picks:
                picks.append(mol_id)
        herb_compounds[herb] = picks
        added[herb] = picks
    return added


if __name__ == "__main__":
    comps: Dict[str, dict] = {}
    herbs: Dict[str, List[str]] = {}
    new = expand_catalog_herbs(comps, herbs)
    pool_size = len({m for ms in herbs.values() for m in ms})
    print(f"已为 {len(new)} 味目录药材生成演示成分关联，"
          f"涉及 {pool_size} 个成分、{len(comps)} 个成分记录。")
    sample = next(iter(new))
    print(f"示例：{sample} -> {new[sample]}")
