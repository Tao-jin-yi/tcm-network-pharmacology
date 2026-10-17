"""解剖层级动物模型 SVG（小鼠 / 大鼠 / 兔，侧位观）。

与 animal_view.py 的示意圆点图不同，本模块按解剖比例绘制：
- 体壁轮廓、头部、耳、四肢、尾按物种形态分别绘制；
- 内脏按器官形态绘制（心脏四腔、肺叶与支气管、肝分叶、胃、脾、胰腺、
  小肠盘曲、盲肠、结肠、肾（豆形+肾门+输尿管）、膀胱、生殖腺、胸腺、淋巴结）；
- 叠加骨骼层（脊柱、肋骨、颅骨、骨盆、股骨/胫骨与关节）与大血管（主动脉/腔静脉）；
- 器官按"该器官相关靶点数"着色（热图），选中器官高亮描边并引出标注线。

说明：这是依据解剖学比例绘制的矢量示意图，用于表达器官层级关系，
不等同于影像学或组织切片图像。
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

# ------------------------------------------------------------ 器官局部形状库
# 每个形状以自身中心为 (0,0)，尺寸约 ±40，使用 {fill}/{op}/{st} 占位
ORGAN_SHAPES: Dict[str, str] = {
    "脑": """
<path d="M0,-16 C18,-26 36,-16 34,0 C32,12 16,18 0,18 C-16,18 -34,12 -34,0
         C-36,-16 -18,-26 0,-16 Z" fill="{fill}" fill-opacity="{op}"
      stroke="{st}" stroke-width="1.6"/>
<path d="M-4,-18 C10,-10 14,0 8,8" fill="none" stroke="{st}" stroke-width="1"/>
<circle cx="26" cy="14" r="9" fill="{fill}" fill-opacity="{op}"
        stroke="{st}" stroke-width="1.4"/>
<path d="M2,18 L6,30" stroke="{st}" stroke-width="4" stroke-linecap="round"/>
""",
    "心": """
<path d="M0,-20 C16,-30 32,-22 30,-6 C28,8 14,22 0,32 C-14,22 -28,8 -30,-6
         C-32,-22 -16,-30 0,-20 Z" fill="{fill}" fill-opacity="{op}"
      stroke="{st}" stroke-width="1.8"/>
<path d="M0,-20 L0,32" stroke="{st}" stroke-width="1.4"/>
<path d="M-14,-14 C-8,-24 8,-24 14,-14" fill="none" stroke="{st}" stroke-width="1.4"/>
<path d="M-6,-22 C-6,-34 10,-38 16,-28 C20,-20 14,-12 6,-12"
      fill="none" stroke="#c0392b" stroke-width="4"/>
<path d="M-2,-24 C-4,-36 -18,-34 -22,-24" fill="none" stroke="#2d6fb3"
      stroke-width="3.5"/>
""",
    "肺": """
<path d="M-4,-32 L4,-32 L8,-8 L-6,-6 Z" fill="#cfd8e3" stroke="{st}" stroke-width="1"/>
<path d="M-6,-6 C-26,-10 -36,4 -32,16 C-28,28 -12,30 -6,20 Z" fill="{fill}"
      fill-opacity="{op}" stroke="{st}" stroke-width="1.6"/>
<path d="M8,-6 C28,-10 38,4 34,16 C30,28 14,30 8,20 Z" fill="{fill}"
      fill-opacity="{op}" stroke="{st}" stroke-width="1.6"/>
<path d="M2,-8 L-16,8 M2,-8 L20,8" stroke="{st}" stroke-width="1.2"/>
""",
    "肝": """
<path d="M-44,-12 C-24,-26 20,-28 42,-12 C48,2 36,18 18,20 C-2,22 -22,18 -34,10
         C-44,2 -50,-4 -44,-12 Z" fill="{fill}" fill-opacity="{op}"
      stroke="{st}" stroke-width="1.8"/>
<path d="M6,-24 L2,20" stroke="{st}" stroke-width="1.2"/>
<path d="M-24,-20 C-18,-8 -20,4 -28,12" fill="none" stroke="{st}" stroke-width="1"/>
""",
    "胃": """
<path d="M-20,-18 C-4,-26 16,-20 20,-6 C24,8 12,22 -4,24 C-16,26 -24,18 -22,8
         C-20,-2 -24,-12 -20,-18 Z" fill="{fill}" fill-opacity="{op}"
      stroke="{st}" stroke-width="1.6"/>
<path d="M20,-6 C30,-4 34,6 30,14" fill="none" stroke="{st}" stroke-width="2.4"/>
<path d="M-20,-18 C-28,-24 -34,-20 -32,-12" fill="none" stroke="{st}" stroke-width="2.4"/>
""",
    "脾": """
<path d="M0,-22 C16,-14 18,8 0,24 C-16,8 -14,-14 0,-22 Z" fill="{fill}"
      fill-opacity="{op}" stroke="{st}" stroke-width="1.6"/>
<path d="M0,-16 L0,18" stroke="{st}" stroke-width="1"/>
""",
    "胰腺": """
<path d="M-28,-2 C-12,-14 14,-16 30,-6 C16,8 -10,10 -28,-2 Z" fill="{fill}"
      fill-opacity="{op}" stroke="{st}" stroke-width="1.6"/>
<circle cx="-14" cy="-4" r="4" fill="{fill}" stroke="{st}" stroke-width="1"/>
<circle cx="4" cy="-6" r="4" fill="{fill}" stroke="{st}" stroke-width="1"/>
<circle cx="20" cy="-2" r="4" fill="{fill}" stroke="{st}" stroke-width="1"/>
""",
    "肠": """
<path d="M-30,-10 q14,-16 28,-2 q14,14 28,-2 q14,-16 28,-2" fill="none"
      stroke="{fill}" stroke-opacity="{op}" stroke-width="9" stroke-linecap="round"/>
<path d="M-30,-10 q14,-16 28,-2 q14,14 28,-2 q14,-16 28,-2" fill="none"
      stroke="{st}" stroke-width="1.2"/>
<path d="M54,-14 q16,-14 26,4" fill="none" stroke="{fill}" stroke-opacity="{op}"
      stroke-width="8" stroke-linecap="round"/>
""",
    "盲肠": """
<path d="M0,-16 C18,-16 26,-2 20,10 C14,22 -8,26 -16,14 C-24,2 -18,-16 0,-16 Z"
      fill="{fill}" fill-opacity="{op}" stroke="{st}" stroke-width="1.6"/>
<path d="M-10,-14 L-16,14" stroke="{st}" stroke-width="1"/>
""",
    "结肠": """
<path d="M-34,-14 L-34,18 L34,18 L34,-6" fill="none" stroke="{fill}"
      stroke-opacity="{op}" stroke-width="8" stroke-linecap="round"/>
<path d="M-34,-14 L-34,18 L34,18 L34,-6" fill="none" stroke="{st}" stroke-width="1.2"/>
""",
    "肾": """
<path d="M0,-20 C18,-20 26,-6 22,8 C18,20 6,26 -6,20 C-14,14 -10,8 -8,2
         C-6,-4 -10,-8 -16,-10 C-22,-12 -18,-20 0,-20 Z" fill="{fill}"
      fill-opacity="{op}" stroke="{st}" stroke-width="1.7"/>
<circle cx="-4" cy="2" r="4" fill="#ffffff" stroke="{st}" stroke-width="1.2"/>
<path d="M-4,6 C-10,16 -14,26 -12,36" fill="none" stroke="{st}" stroke-width="2"/>
""",
    "膀胱": """
<path d="M0,-12 C14,-12 20,2 12,12 C4,20 -10,18 -14,8 C-18,-2 -12,-12 0,-12 Z"
      fill="{fill}" fill-opacity="{op}" stroke="{st}" stroke-width="1.6"/>
<path d="M-2,-12 C-2,-22 6,-26 10,-22" fill="none" stroke="{st}" stroke-width="2"/>
""",
    "生殖系统": """
<ellipse cx="0" cy="0" rx="16" ry="11" fill="{fill}" fill-opacity="{op}"
         stroke="{st}" stroke-width="1.6"/>
<path d="M16,-4 C26,-8 30,2 24,8" fill="none" stroke="{st}" stroke-width="2"/>
<path d="M-2,11 C2,20 10,22 14,18" fill="none" stroke="{st}" stroke-width="1.6"/>
""",
    "胸腺": """
<ellipse cx="-10" cy="0" rx="8" ry="12" fill="{fill}" fill-opacity="{op}"
         stroke="{st}" stroke-width="1.4"/>
<ellipse cx="10" cy="4" rx="7" ry="10" fill="{fill}" fill-opacity="{op}"
         stroke="{st}" stroke-width="1.4"/>
""",
    "免疫细胞": """
<circle cx="-14" cy="-6" r="7" fill="{fill}" fill-opacity="{op}" stroke="{st}"
        stroke-width="1.4"/>
<circle cx="6" cy="6" r="7" fill="{fill}" fill-opacity="{op}" stroke="{st}"
        stroke-width="1.4"/>
<circle cx="20" cy="-10" r="6" fill="{fill}" fill-opacity="{op}" stroke="{st}"
        stroke-width="1.4"/>
<path d="M-14,-6 L6,6 L20,-10" fill="none" stroke="{st}" stroke-width="1"/>
""",
    "血管": """
<path d="M-30,-20 C-6,-24 18,-16 30,0" fill="none" stroke="#c0392b" stroke-width="5"
      stroke-linecap="round"/>
<path d="M-30,-8 C-4,-10 16,-4 28,10" fill="none" stroke="#2d6fb3" stroke-width="4"
      stroke-linecap="round"/>
<path d="M-30,-14 C-8,-16 14,-10 26,4" fill="none" stroke="{fill}"
      stroke-opacity="{op}" stroke-width="2"/>
""",
    "骨骼肌": """
<path d="M-26,-14 C-6,-24 20,-18 26,-2 C30,12 16,22 -2,22 C-18,22 -30,10 -26,-14 Z"
      fill="{fill}" fill-opacity="{op}" stroke="{st}" stroke-width="1.6"/>
<path d="M-18,-8 L14,-2 M-16,4 L16,8 M-12,14 L10,16" stroke="{st}" stroke-width="1"/>
""",
    "骨与关节": """
<path d="M-18,-26 L-6,-4 L-2,18" stroke="#b9b9b9" stroke-width="9"
      stroke-linecap="round" fill="none"/>
<path d="M-2,18 L14,34" stroke="#b9b9b9" stroke-width="8" stroke-linecap="round"
      fill="none"/>
<circle cx="-2" cy="18" r="7" fill="{fill}" fill-opacity="{op}" stroke="{st}"
        stroke-width="1.6"/>
<circle cx="-18" cy="-26" r="6" fill="#dcdcdc" stroke="{st}" stroke-width="1.2"/>
""",
    "脂肪": """
<path d="M-30,-8 C-10,-16 16,-14 30,-6 C34,4 22,14 0,14 C-20,14 -34,4 -30,-8 Z"
      fill="{fill}" fill-opacity="{op}" stroke="{st}" stroke-width="1.4"/>
<circle cx="-12" cy="0" r="5" fill="#ffffff" stroke="{st}" stroke-width="0.8"/>
<circle cx="8" cy="4" r="5" fill="#ffffff" stroke="{st}" stroke-width="0.8"/>
""",
    "皮肤": """
<path d="M-30,-10 C-8,-18 18,-16 30,-8" fill="none" stroke="{fill}"
      stroke-opacity="{op}" stroke-width="6" stroke-linecap="round"/>
<path d="M-30,-2 C-8,-10 18,-8 30,0" fill="none" stroke="{st}" stroke-width="1.2"/>
""",
}

# 器官在躯干中的相对位置（fx = 躯干宽度比例, fy = 躯干高度比例），s 为缩放系数
ORGAN_LAYOUT: Dict[str, Tuple[float, float, float]] = {
    "胸腺": (0.22, 0.30, 0.55),
    "心": (0.38, 0.46, 0.62),
    "肺": (0.44, 0.26, 0.66),
    "血管": (0.30, 0.62, 0.55),
    "肝": (0.60, 0.44, 0.78),
    "胃": (0.70, 0.34, 0.58),
    "胰腺": (0.64, 0.62, 0.58),
    "脾": (0.80, 0.34, 0.50),
    "肠": (0.78, 0.62, 0.72),
    "结肠": (0.86, 0.40, 0.60),
    "盲肠": (0.90, 0.72, 0.62),
    "肾": (0.62, 0.80, 0.62),
    "膀胱": (0.94, 0.86, 0.48),
    "生殖系统": (0.96, 0.66, 0.52),
    "免疫细胞": (0.18, 0.76, 0.60),
    "骨骼肌": (0.88, 1.06, 0.72),
    "骨与关节": (0.80, 1.16, 0.70),
    "脂肪": (0.60, 0.96, 0.50),
    "皮肤": (0.30, 1.02, 0.50),
}

# 头部器官（相对头心半径）
HEAD_ORGANS = {"脑": (-0.05, -0.12, 0.95)}

SPECIES: Dict[str, dict] = {
    "小鼠（C57BL/6）": {"x0": 150, "y0": 148, "w": 250, "h": 122, "head": 34,
                        "torso_top": 22, "torso_bottom": 16, "tail": "thin-long",
                        "ear": "round", "ear_r": 17, "cecum": 0.55},
    "大鼠（SD）": {"x0": 140, "y0": 156, "w": 292, "h": 146, "head": 38,
                   "torso_top": 26, "torso_bottom": 20, "tail": "thick-long",
                   "ear": "round", "ear_r": 14, "cecum": 0.7},
    "兔（新西兰白兔）": {"x0": 150, "y0": 164, "w": 278, "h": 158, "head": 34,
                         "torso_top": 30, "torso_bottom": 22, "tail": "short",
                         "ear": "long", "ear_r": 0, "cecum": 1.15},
}


def _heat(value: int, max_value: int) -> Tuple[str, float]:
    """靶点数 → (填充色, 不透明度)。"""
    ratio = value / max_value if max_value else 0.0
    ratio = max(0.0, min(1.0, ratio))
    red = int(250 - 40 * ratio)
    green = int(238 - 190 * ratio)
    blue = int(196 - 160 * ratio)
    return f"rgb({red},{green},{blue})", round(0.35 + 0.6 * ratio, 2)


def _torso_outline(s: dict) -> str:
    x0, y0, w, h = s["x0"], s["y0"], s["w"], s["h"]
    top, bottom = s["torso_top"], s["torso_bottom"]
    return (f'<path d="M{x0},{y0 + h / 2:.0f} '
            f'C{x0 - 14},{y0 - top} {x0 + w * 0.28:.0f},{y0 - top - 6} '
            f'{x0 + w * 0.55:.0f},{y0 - top + 2} '
            f'C{x0 + w * 0.82:.0f},{y0 - top + 4} {x0 + w + 16},{y0 + 10} '
            f'{x0 + w},{y0 + h * 0.62:.0f} '
            f'C{x0 + w - 8},{y0 + h + bottom} {x0 + w * 0.42:.0f},{y0 + h + bottom + 6} '
            f'{x0 + w * 0.24:.0f},{y0 + h + bottom} '
            f'C{x0 + 4},{y0 + h + bottom - 4} {x0 - 10},{y0 + h * 0.72:.0f} '
            f'{x0},{y0 + h / 2:.0f} Z" fill="#f6f2ea" stroke="#8a8a8a" stroke-width="2"/>')


def _head(s: dict) -> Tuple[str, Tuple[float, float, float]]:
    x0, y0, w, h, r = s["x0"], s["y0"], s["w"], s["h"], s["head"]
    cx, cy = x0 + w + 18, y0 + h * 0.30
    parts = [f'<circle cx="{cx:.0f}" cy="{cy:.0f}" r="{r}" fill="#f6f2ea" '
             f'stroke="#8a8a8a" stroke-width="2"/>',
             f'<path d="M{cx + r * 0.6:.0f},{cy + r * 0.1:.0f} '
             f'L{cx + r * 1.45:.0f},{cy + r * 0.42:.0f} '
             f'L{cx + r * 0.55:.0f},{cy + r * 0.72:.0f} Z" fill="#f6f2ea" '
             f'stroke="#8a8a8a" stroke-width="2"/>',
             f'<circle cx="{cx + r * 1.28:.0f}" cy="{cy + r * 0.42:.0f}" r="2.6" '
             f'fill="#333"/>',
             f'<circle cx="{cx + r * 0.35:.0f}" cy="{cy - r * 0.15:.0f}" r="3" '
             f'fill="#222"/>']
    if s["ear"] == "round":
        parts.append(f'<circle cx="{cx + r * 0.42:.0f}" cy="{cy - r * 0.92:.0f}" '
                     f'r="{s["ear_r"]}" fill="#efe4d8" stroke="#8a8a8a" stroke-width="2"/>')
    else:
        for dx in (-0.35, 0.35):
            parts.append(
                f'<ellipse cx="{cx + r * dx:.0f}" cy="{cy - r * 1.9:.0f}" rx="9" ry="34" '
                f'fill="#efe4d8" stroke="#8a8a8a" stroke-width="2"/>')
    return "".join(parts), (cx, cy, r)


def _tail_legs(s: dict) -> str:
    x0, y0, w, h = s["x0"], s["y0"], s["w"], s["h"]
    parts = []
    if s["tail"] == "thin-long":
        parts.append(f'<path d="M{x0},{y0 + h * 0.36:.0f} q -70 -6 -124 34" fill="none" '
                     f'stroke="#8a8a8a" stroke-width="3.5"/>')
    elif s["tail"] == "thick-long":
        parts.append(f'<path d="M{x0},{y0 + h * 0.34:.0f} q -84 -4 -138 44" fill="none" '
                     f'stroke="#8a8a8a" stroke-width="6"/>')
    else:
        parts.append(f'<circle cx="{x0 - 6}" cy="{y0 + h * 0.42:.0f}" r="12" '
                     f'fill="#f6f2ea" stroke="#8a8a8a" stroke-width="2"/>')
    for fx in (0.16, 0.86):
        lx = x0 + w * fx
        parts.append(f'<path d="M{lx:.0f},{y0 + h * 0.92:.0f} l -6,44 l 14,10" '
                     f'fill="none" stroke="#8a8a8a" stroke-width="7" '
                     f'stroke-linecap="round"/>')
    return "".join(parts)


def _skeleton(s: dict) -> str:
    """脊柱 / 肋骨 / 骨盆 / 颅骨的浅灰底衬。"""
    x0, y0, w, h = s["x0"], s["y0"], s["w"], s["h"]
    parts = [f'<path d="M{x0 + 6},{y0 + 6} C{x0 + w * 0.35:.0f},{y0 - 8} '
             f'{x0 + w * 0.75:.0f},{y0 - 2} {x0 + w - 10},{y0 + h * 0.22:.0f}" '
             f'fill="none" stroke="#cfcfcf" stroke-width="7" stroke-linecap="round"/>']
    for i in range(9):  # 椎骨
        fx = 0.10 + i * 0.09
        px = x0 + w * fx
        py = y0 + 4 + (i * 1.5)
        parts.append(f'<rect x="{px:.0f}" y="{py:.0f}" width="9" height="6" rx="2" '
                     f'fill="#e2e2e2" stroke="#bdbdbd" stroke-width="0.8"/>')
    for i in range(7):  # 肋骨
        fx = 0.30 + i * 0.075
        px = x0 + w * fx
        parts.append(f'<path d="M{px:.0f},{y0 + 8} q 6,{h * 0.42:.0f} '
                     f'-14,{h * 0.62:.0f}" fill="none" stroke="#dcdcdc" '
                     f'stroke-width="3"/>')
    parts.append(f'<ellipse cx="{x0 + w * 0.86:.0f}" cy="{y0 + h * 0.92:.0f}" rx="30" '
                 f'ry="16" fill="none" stroke="#d8d8d8" stroke-width="4"/>')
    return "".join(parts)


def render_anatomy(animal: str, counts: Dict[str, int],
                   highlight: Optional[str] = None,
                   width: int = 640, height: int = 430) -> str:
    """生成解剖层级动物 SVG（含骨骼层、器官层、标注层）。"""
    s = SPECIES.get(animal, SPECIES["小鼠（C57BL/6）"])
    x0, y0, w, h = s["x0"], s["y0"], s["w"], s["h"]
    max_count = max(counts.values()) if counts else 1
    base_scale = h / 130.0

    layers = [
        f'<svg viewBox="0 0 {width} {height}" width="100%" height="100%" '
        f'xmlns="http://www.w3.org/2000/svg">',
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
        f'<text x="14" y="24" font-size="15" font-weight="bold" fill="#333">{animal}'
        f'（侧位观 · 解剖示意）</text>',
        _tail_legs(s), _torso_outline(s), _skeleton(s),
    ]

    head_parts, (hx, hy, hr) = _head(s)
    layers.append(head_parts)

    labels: List[str] = []
    drawn: Dict[str, Tuple[float, float]] = {}

    def draw_organ(name: str, cx: float, cy: float, scale: float) -> None:
        shape = ORGAN_SHAPES.get(name)
        if not shape:
            return
        value = counts.get(name, 0)
        fill, opacity = _heat(value, max_count)
        stroke = "#E4572E" if name == highlight else "#5a5a5a"
        extra = ' stroke-dasharray=""'
        shape = shape.format(fill=fill, op=opacity, st=stroke)
        thick = f'stroke-width="{3 if name == highlight else 1}"'
        layers.append(f'<g transform="translate({cx:.0f},{cy:.0f}) '
                      f'scale({scale:.2f})">{shape}</g>')
        if name == highlight:
            layers.append(f'<circle cx="{cx:.0f}" cy="{cy:.0f}" r="{34 * scale:.0f}" '
                          f'fill="none" stroke="#E4572E" stroke-width="2.5" '
                          f'stroke-dasharray="5 4"/>')
        drawn[name] = (cx, cy)
        if value:
            layers.append(f'<text x="{cx:.0f}" y="{cy + 4:.0f}" text-anchor="middle" '
                          f'font-size="11" font-weight="bold" fill="#1a1a1a">{value}</text>')

    # 头部器官
    for name, (fx, fy, fs) in HEAD_ORGANS.items():
        draw_organ(name, hx + hr * fx, hy + hr * fy, (hr / 34) * fs)

    # 躯干器官
    for name, (fx, fy, fs) in ORGAN_LAYOUT.items():
        cx, cy = x0 + w * fx, y0 + h * fy
        scale = base_scale * fs
        if name == "盲肠":
            scale *= s["cecum"]
        draw_organ(name, cx, cy, scale)

    # 标注：上下两排引出线
    top_labels = [n for n in ORGAN_LAYOUT if ORGAN_LAYOUT[n][1] <= 0.62
                  and n in drawn] + ["脑"]
    bottom_labels = [n for n in ORGAN_LAYOUT if ORGAN_LAYOUT[n][1] > 0.62 and n in drawn]
    for group, base_y, direction in ((top_labels, 44, -1), (bottom_labels, height - 34, 1)):
        for i, name in enumerate(group):
            cx, cy = drawn[name]
            ly = base_y + (i % 2) * 18
            lx = 24 + (i // 2) * 78
            color = "#E4572E" if name == highlight else "#666"
            layers.append(f'<path d="M{cx:.0f},{cy:.0f} L{lx},{ly}" stroke="{color}" '
                          f'stroke-width="0.9" fill="none" stroke-dasharray="3 3"/>')
            layers.append(f'<text x="{lx}" y="{ly - 3}" font-size="11" fill="{color}">'
                          f'{name}（{counts.get(name, 0)}）</text>')

    layers.append("</svg>")
    svg = "".join(layers)
    return ('<html><body style="margin:0;background:#ffffff;font-family:'
            'Microsoft YaHei,sans-serif">' + svg + "</body></html>")
