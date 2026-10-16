"""动物模型与细胞示意图的 SVG 可视化。

提供三张图：
1. 动物（小鼠 / 大鼠 / 兔）轮廓 + 器官点位：按靶点分布数量着色，可高亮所选器官；
2. 细胞分区图：把靶点按 UniProt 亚细胞定位放到 细胞外 / 细胞膜 / 细胞质 / 细胞核 /
   线粒体 / 内质网 等区室，直观点出"受体点位"；
3. 组织来源判定：用关键词把 UniProt 组织特异性文本映射到器官。

注意：器官位置是示意性布局（不是解剖学精确作图），用于表达"靶点-器官"关系。
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional, Tuple

# ------------------------------------------------------- 器官关键词（组织来源）
ORGAN_PATTERNS: Dict[str, str] = {
    "脑": r"brain|cerebr|neuron|hippocamp|astrocyte|microglia|cortex|中枢神经",
    "心": r"heart|cardiac|myocard|cardiomyocyte|心肌",
    "肺": r"lung|pulmonary|alveolar|bronch|肺",
    "肝": r"liver|hepatic|hepatocyte|肝",
    "脾": r"spleen|splenic|脾",
    "肾": r"kidney|renal|nephron|glomerul|肾",
    "胃": r"stomach|gastric|胃",
    "肠": r"intestin|colon|gut|ileum|jejunum|cecum|rectum|肠道",
    "胰腺": r"pancrea|islet|胰岛",
    "骨骼肌": r"skeletal muscle|myocyte|myoblast|muscle|骨骼肌",
    "脂肪": r"adipose|adipocyte|fat tissue|脂肪",
    "血管": r"endotheli|vascular|vessel|aorta|arter|血管",
    "生殖系统": r"testis|testicular|ovary|ovarian|uterus|endometrium|prostate|placenta",
    "骨与关节": r"bone|osteoblast|osteoclast|cartilage|chondrocyte|joint|synovi",
    "皮肤": r"skin|keratinocyte|epiderm|dermal|fibroblast|皮肤",
    "免疫细胞": r"macrophage|monocyte|lymphocyte|dendritic|neutrophil|t-cell|b-cell|"
                r"thymus|bone marrow|巨噬|免疫细胞",
    "血液": r"erythrocyte|platelet|blood cell|serum|plasma|leukocyte",
}

# -------------------------------------------------------------- 动物轮廓定义
ANIMALS: Dict[str, dict] = {
    "小鼠（C57BL/6）": {
        "body": """
<ellipse cx="150" cy="115" rx="72" ry="40" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<circle cx="245" cy="92" r="25" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<circle cx="252" cy="66" r="12" fill="#efe6dd" stroke="#8a8a8a" stroke-width="2"/>
<circle cx="258" cy="88" r="3" fill="#333"/>
<circle cx="266" cy="100" r="4" fill="#e0a3a3" stroke="#8a8a8a" stroke-width="1"/>
<path d="M80 112 q -34 -10 -56 20" fill="none" stroke="#8a8a8a" stroke-width="3"/>
<ellipse cx="112" cy="152" rx="12" ry="8" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="140" cy="156" rx="12" ry="8" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="192" cy="152" rx="12" ry="8" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="216" cy="154" rx="12" ry="8" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
""",
        "organs": {
            "脑": (245, 80), "心": (150, 100), "肺": (174, 98), "肝": (130, 122),
            "脾": (182, 124), "胃": (150, 134), "肾": (106, 130), "肠": (196, 140),
            "胰腺": (162, 144), "血管": (150, 84), "骨骼肌": (214, 152),
            "脂肪": (92, 104), "生殖系统": (214, 128), "骨与关节": (96, 150),
            "皮肤": (120, 88), "免疫细胞": (168, 110),
        },
    },
    "大鼠（SD）": {
        "body": """
<ellipse cx="145" cy="112" rx="80" ry="42" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<circle cx="250" cy="95" r="24" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<circle cx="256" cy="73" r="9" fill="#efe6dd" stroke="#8a8a8a" stroke-width="2"/>
<circle cx="260" cy="92" r="3" fill="#333"/>
<circle cx="268" cy="104" r="4" fill="#e0a3a3" stroke="#8a8a8a" stroke-width="1"/>
<path d="M66 110 q -46 -6 -62 26" fill="none" stroke="#8a8a8a" stroke-width="5"/>
<ellipse cx="105" cy="150" rx="13" ry="9" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="135" cy="154" rx="13" ry="9" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="188" cy="150" rx="13" ry="9" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="214" cy="152" rx="13" ry="9" fill="#f4f1ea" stroke="#8a8a8a" stroke-width="2"/>
""",
        "organs": {
            "脑": (250, 74), "心": (145, 96), "肺": (170, 94), "肝": (124, 118),
            "脾": (178, 120), "胃": (145, 130), "肾": (100, 126), "肠": (192, 136),
            "胰腺": (158, 140), "血管": (145, 80), "骨骼肌": (212, 150),
            "脂肪": (86, 100), "生殖系统": (210, 124), "骨与关节": (90, 148),
            "皮肤": (115, 84), "免疫细胞": (164, 106),
        },
    },
    "兔（新西兰白兔）": {
        "body": """
<ellipse cx="158" cy="122" rx="70" ry="46" fill="#f7f5f1" stroke="#8a8a8a" stroke-width="2"/>
<circle cx="248" cy="100" r="23" fill="#f7f5f1" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="232" cy="56" rx="9" ry="28" fill="#f0e8e0" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="258" cy="58" rx="9" ry="26" fill="#f0e8e0" stroke="#8a8a8a" stroke-width="2"/>
<circle cx="262" cy="98" r="3" fill="#333"/>
<circle cx="270" cy="108" r="4" fill="#e0a3a3" stroke="#8a8a8a" stroke-width="1"/>
<circle cx="88" cy="118" r="10" fill="#f7f5f1" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="118" cy="164" rx="14" ry="9" fill="#f7f5f1" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="150" cy="168" rx="14" ry="9" fill="#f7f5f1" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="196" cy="164" rx="14" ry="9" fill="#f7f5f1" stroke="#8a8a8a" stroke-width="2"/>
<ellipse cx="220" cy="166" rx="14" ry="9" fill="#f7f5f1" stroke="#8a8a8a" stroke-width="2"/>
""",
        "organs": {
            "脑": (248, 92), "心": (158, 106), "肺": (182, 104), "肝": (138, 128),
            "脾": (190, 130), "胃": (158, 140), "肾": (114, 136), "肠": (204, 146),
            "胰腺": (170, 150), "血管": (158, 90), "骨骼肌": (218, 164),
            "脂肪": (100, 108), "生殖系统": (216, 134), "骨与关节": (104, 158),
            "皮肤": (128, 94), "免疫细胞": (176, 116),
        },
    },
}


def organs_from_text(text: str) -> List[str]:
    """从 UniProt 组织特异性文本中识别器官。"""
    if not text:
        return []
    return [organ for organ, pattern in ORGAN_PATTERNS.items()
            if re.search(pattern, text, re.I)]


def organ_target_map(target_info: Dict[str, dict]) -> Dict[str, List[str]]:
    """器官 -> 在该器官表达的靶点列表（依据 UniProt 组织特异性注释）。"""
    mapping: Dict[str, List[str]] = {}
    for gene, item in target_info.items():
        for organ in organs_from_text(item.get("tissue", "")):
            mapping.setdefault(organ, []).append(gene)
    return {k: sorted(set(v)) for k, v in mapping.items()}


def render_animal(animal: str, counts: Dict[str, int],
                  highlight: Optional[str] = None) -> str:
    """生成动物 + 器官点位的 SVG 网页（器官越大越红表示该器官相关靶点越多）。"""
    spec = ANIMALS.get(animal, ANIMALS["小鼠（C57BL/6）"])
    max_count = max(counts.values()) if counts else 1
    markers = []

    for organ, (x, y) in spec["organs"].items():
        count = counts.get(organ, 0)
        ratio = count / max_count if max_count else 0
        radius = 7 + 7 * ratio
        red = int(90 + 150 * ratio)
        green = int(150 - 90 * ratio)
        fill = f"rgb({red},{green},70)"
        stroke = "#E4572E" if organ == highlight else "#5a5a5a"
        width = 3 if organ == highlight else 1
        opacity = 0.35 + 0.6 * ratio if count else 0.18
        markers.append(
            f'<circle cx="{x}" cy="{y}" r="{radius:.1f}" fill="{fill}" '
            f'fill-opacity="{opacity:.2f}" stroke="{stroke}" stroke-width="{width}"/>'
        )
        if count:
            markers.append(
                f'<text x="{x}" y="{y + 4}" text-anchor="middle" '
                f'font-size="11" font-weight="bold" fill="#1a1a1a">{count}</text>')
        label_y = y + radius + 13
        markers.append(
            f'<text x="{x}" y="{label_y:.0f}" text-anchor="middle" font-size="10" '
            f'fill="{"#E4572E" if organ == highlight else "#555"}">{organ}</text>')

    svg = (f'<svg viewBox="0 0 320 210" width="100%" height="100%" '
           f'xmlns="http://www.w3.org/2000/svg">'
           f'<rect width="320" height="210" fill="#ffffff"/>'
           f'{spec["body"]}' + "".join(markers) + "</svg>")
    return (f'<html><body style="margin:0;background:#ffffff;font-family:'
            f'Microsoft YaHei,sans-serif">{svg}</body></html>')


# ---------------------------------------------------------------- 细胞示意图
COMPARTMENTS = {
    "细胞外（分泌/配体）": (60, 34),
    "细胞膜（受体）": (60, 78),
    "细胞质（激酶/信号）": (60, 148),
    "细胞核（转录因子）": (250, 150),
    "线粒体": (250, 78),
    "内质网/高尔基体": (60, 232),
}

COMPARTMENT_RULES = [
    ("细胞外（分泌/配体）", r"secreted|extracellular|extracellular space"),
    ("细胞膜（受体）", r"cell membrane|membrane|cell surface|cell junction|apical"),
    ("细胞核（转录因子）", r"nucleus|nucleolus|nucleoplasm|chromosome|pml body"),
    ("线粒体", r"mitochondri|mitochondrial"),
    ("内质网/高尔基体", r"endoplasmic reticulum|golgi|er\b|lysosome|peroxisome|"
                        r"cytoplasmic vesicle"),
    ("细胞质（激酶/信号）", r"cytoplasm|cytosol|cytoskeleton|centrosome"),
]


def compartment_of(location: str) -> str:
    """把 UniProt 亚细胞定位文本映射到细胞区室。"""
    for compartment, pattern in COMPARTMENT_RULES:
        if re.search(pattern, location, re.I):
            return compartment
    return "细胞质（激酶/信号）"


def render_cell(targets_with_location: Iterable[Tuple[str, str]]) -> str:
    """生成细胞分区示意图，把靶点放到对应区室（受体点位可视化）。"""
    buckets: Dict[str, List[str]] = {}
    for gene, location in targets_with_location:
        buckets.setdefault(compartment_of(location), []).append(gene)

    parts = [
        '<svg viewBox="0 0 420 300" width="100%" height="100%" '
        'xmlns="http://www.w3.org/2000/svg">',
        '<rect width="420" height="300" fill="#ffffff"/>',
        # 细胞外形
        '<ellipse cx="210" cy="160" rx="195" ry="120" fill="#fdfaf3" '
        'stroke="#8a8a8a" stroke-width="2"/>',
        # 细胞核
        '<circle cx="250" cy="160" r="58" fill="#e8f0fb" stroke="#5677a8" '
        'stroke-width="2"/>',
        '<text x="250" y="163" text-anchor="middle" font-size="12" '
        'fill="#5677a8">细胞核</text>',
        # 细胞膜（下方横线 + 标注）
        '<path d="M18 100 A195 120 0 0 0 402 100" fill="none" '
        'stroke="#E4572E" stroke-width="4"/>',
        # 线粒体
        '<ellipse cx="330" cy="95" rx="34" ry="18" fill="#ffe9d6" '
        'stroke="#d98324" stroke-width="2"/>',
        '<text x="330" y="99" text-anchor="middle" font-size="10" '
        'fill="#d98324">线粒体</text>',
        # 内质网
        '<path d="M60 236 q 30 -18 60 0 q 30 18 60 0" fill="none" '
        'stroke="#7D5BA6" stroke-width="2"/>',
        '<text x="150" y="258" text-anchor="middle" font-size="10" '
        'fill="#7D5BA6">内质网 / 高尔基体</text>',
        '<text x="30" y="26" font-size="12" fill="#888">细胞外（配体 / 分泌蛋白）</text>',
        '<text x="30" y="118" font-size="12" fill="#E4572E">细胞膜（受体）</text>',
        '<text x="30" y="200" font-size="12" fill="#555">细胞质（激酶 / 信号蛋白）</text>',
    ]

    for compartment, genes in buckets.items():
        x, y = COMPARTMENTS.get(compartment, (60, 148))
        for i, gene in enumerate(sorted(genes)[:8]):
            rows, cols = divmod(i, 3)
            cx, cy = x + cols * 62, y + rows * 20
            parts.append(
                f'<rect x="{cx}" y="{cy - 12}" rx="8" ry="8" '
                f'width="56" height="17" fill="#2A629A" fill-opacity="0.85"/>')
            parts.append(
                f'<text x="{cx + 28}" y="{cy + 1}" text-anchor="middle" '
                f'font-size="10" fill="#ffffff">{gene}</text>')
        if len(genes) > 8:
            parts.append(
                f'<text x="{x}" y="{y + 68}" font-size="10" fill="#888">'
                f'+{len(genes) - 8} 个</text>')

    parts.append("</svg>")
    return ('<html><body style="margin:0;background:#ffffff;font-family:'
            'Microsoft YaHei,sans-serif">' + "".join(parts) + "</body></html>")
