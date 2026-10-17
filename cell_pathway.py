"""细胞水平的通路可视化：把通路节点放进细胞结构里。

统一坐标系：x ∈ [0,10]，y ∈ [0,7]（SVG 与视频共用同一套坐标，保证两处图形一致）

    0 ── 2.5          2.5 ── 5.6            5.6 ── 10
    细胞外（配体） │ 细胞质（接头 / 激酶级联） │ 细胞核（转录因子 / 效应）
    膜受体画在 x≈2.55 的细胞膜上

节点按类型自动分配到对应区室：
    LIGAND → 细胞外，RECEPTOR → 细胞膜，ADAPTOR/KINASE → 细胞质，
    TF → 核膜内侧，EFFECT → 核内（基因转录 → 生物学效应）
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from pathways import KIND_STYLE, PATHWAYS

# 区室中心 x（数据坐标）
X_LIGAND, X_RECEPTOR, X_ADAPTOR, X_KINASE, X_TF, X_EFFECT = 0.95, 2.55, 3.65, 4.85, 6.55, 8.35
MEMBRANE_X = 2.5
NUCLEUS = (7.6, 3.3, 2.05, 2.35)  # cx, cy, rx, ry（数据坐标）
Y_RANGE = (0.95, 6.05)            # 细胞质内节点可放置的纵向范围
NUC_Y_RANGE = (1.35, 5.25)        # 核内节点纵向范围


def _spread(count: int, y_range: Tuple[float, float], max_gap: float = 1.15) -> List[float]:
    """把 count 个节点在纵向范围内均匀铺开（居中）。"""
    if count <= 0:
        return []
    lo, hi = y_range
    gap = min(max_gap, (hi - lo) / max(count, 1))
    center = (lo + hi) / 2
    return [center + (i - (count - 1) / 2) * gap for i in range(count)]


def cell_layout(name: str) -> Tuple[Dict[str, Tuple[float, float]], List[str]]:
    """返回 (节点坐标, 级联顺序)。级联顺序按拓扑深度排序，保证动画按信号流向推进。"""
    pw = PATHWAYS[name]
    nodes = {n["id"]: n for n in pw["nodes"]}
    edges = pw["edges"]

    depth: Dict[str, int] = {nid: 0 for nid in nodes}
    for _ in range(len(nodes) + 1):
        changed = False
        for s, t, _ in edges:
            if depth[t] < depth[s] + 1:
                depth[t] = depth[s] + 1
                changed = True
        if not changed:
            break

    by_kind: Dict[str, List[str]] = {}
    for nid, node in nodes.items():
        by_kind.setdefault(node["kind"], []).append(nid)

    pos: Dict[str, Tuple[float, float]] = {}
    for kind, x in (("LIGAND", X_LIGAND), ("RECEPTOR", X_RECEPTOR),
                    ("ADAPTOR", X_ADAPTOR), ("KINASE", X_KINASE)):
        members = sorted(by_kind.get(kind, []), key=lambda k: (depth[k], k))
        for nid, y in zip(members, _spread(len(members), Y_RANGE)):
            pos[nid] = (x, y)
    for kind, x, y_range in (("TF", X_TF, NUC_Y_RANGE), ("EFFECT", X_EFFECT, NUC_Y_RANGE)):
        members = sorted(by_kind.get(kind, []), key=lambda k: (depth[k], k))
        for nid, y in zip(members, _spread(len(members), y_range)):
            pos[nid] = (x, y)

    order = [nid for nid in sorted(nodes, key=lambda k: (depth[k], k))]
    return pos, order


# ------------------------------------------------------------------ SVG 渲染
SCALE = 100  # 数据坐标 → SVG 像素


def _wrap(text: str, width: int = 6) -> List[str]:
    if " / " in text:
        parts = text.split(" / ")
        return ["/".join(parts[:2]), "/".join(parts[2:])] if len(parts) > 2 else parts
    return [text[i:i + width] for i in range(0, len(text), width)] or [text]


def _svg_pos(x: float, y: float) -> Tuple[float, float]:
    return x * SCALE, y * SCALE


def render_cell_pathway(name: str, highlight: str = "", targets: set = frozenset(),
                        show_legend: bool = True) -> str:
    """生成细胞水平通路 SVG（HTML 片段，供 components.html 使用）。"""
    pw = PATHWAYS[name]
    nodes = {n["id"]: n for n in pw["nodes"]}
    pos, _ = cell_layout(name)

    w, h = 10 * SCALE, 7 * SCALE
    parts: List[str] = [
        f'<svg viewBox="0 0 {w} {h}" width="100%" height="100%" '
        f'xmlns="http://www.w3.org/2000/svg">',
        '<defs>'
        '<marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
        'markerHeight="7" orient="auto-start-reverse">'
        '<path d="M0,0 L10,5 L0,10 z" fill="#7a7a7a"/></marker>'
        '<marker id="arrowHot" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
        'markerHeight="7" orient="auto-start-reverse">'
        '<path d="M0,0 L10,5 L0,10 z" fill="#E4572E"/></marker>'
        '</defs>',
        # 细胞外
        f'<rect x="0" y="0" width="{MEMBRANE_X * SCALE}" height="{h}" fill="#eaf2fb"/>',
        '<text x="14" y="30" font-size="20" fill="#4a76a8" font-weight="bold">细胞外</text>',
        '<text x="14" y="54" font-size="15" fill="#6b8fb5">配体 / 细胞因子</text>',
        # 细胞体
        f'<rect x="{MEMBRANE_X * SCALE}" y="20" width="{(10 - MEMBRANE_X) * SCALE - 20}" '
        f'height="{h - 40}" rx="60" fill="#fffdf5" stroke="#8a8a8a" stroke-width="2"/>',
        # 细胞膜（磷脂双层示意）
        f'<rect x="{MEMBRANE_X * SCALE - 8}" y="20" width="16" height="{h - 40}" rx="8" '
        f'fill="#f3d7a4" stroke="#c9a45c" stroke-width="1.5"/>',
        f'<text x="{MEMBRANE_X * SCALE + 16}" y="{h - 26}" font-size="15" '
        f'fill="#b08a3c">细胞膜（受体）</text>',
        # 细胞核
        f'<ellipse cx="{NUCLEUS[0] * SCALE}" cy="{NUCLEUS[1] * SCALE}" '
        f'rx="{NUCLEUS[2] * SCALE}" ry="{NUCLEUS[3] * SCALE}" fill="#e8f0fb" '
        f'stroke="#5677a8" stroke-width="2.5"/>',
        f'<ellipse cx="{NUCLEUS[0] * SCALE}" cy="{NUCLEUS[1] * SCALE}" '
        f'rx="{NUCLEUS[2] * SCALE - 9}" ry="{NUCLEUS[3] * SCALE - 9}" fill="none" '
        f'stroke="#9db8d8" stroke-width="1.5" stroke-dasharray="6 5"/>',
        f'<circle cx="{NUCLEUS[0] * SCALE + 55}" cy="{NUCLEUS[1] * SCALE + 70}" r="34" '
        f'fill="#d5e3f5" stroke="#9db8d8" stroke-width="1.2"/>',
        f'<text x="{NUCLEUS[0] * SCALE - 60}" y="{NUCLEUS[1] * SCALE - NUCLEUS[3] * SCALE + 34}" '
        f'font-size="20" fill="#5677a8" font-weight="bold">细胞核</text>',
        '<text x="290" y="655" font-size="16" fill="#8a8a8a">细胞质（激酶级联）</text>',
        # 线粒体（示意）
        '<ellipse cx="330" cy="560" rx="42" ry="22" fill="#ffe9d6" stroke="#d98324" '
        'stroke-width="1.5"/>'
        '<path d="M306,552 q12,10 24,0 q12,-10 24,0" fill="none" stroke="#d98324" '
        'stroke-width="1.2"/>'
        '<text x="330" y="600" font-size="13" fill="#d98324" text-anchor="middle">线粒体</text>',
        '<ellipse cx="380" cy="120" rx="36" ry="19" fill="#ffe9d6" stroke="#d98324" '
        'stroke-width="1.5"/>',
    ]

    # 边
    for s, t, label in pw["edges"]:
        x1, y1 = _svg_pos(*pos[s])
        x2, y2 = _svg_pos(*pos[t])
        hot = (highlight and (highlight in nodes[s].get("gene", "") or
                              highlight in nodes[t].get("gene", "")))
        marker = "arrowHot" if hot else "arrow"
        color = "#E4572E" if hot else "#7a7a7a"
        # 从源节点框边缘出发，避免箭头被框遮住
        import math
        dx, dy = x2 - x1, y2 - y1
        dist = math.hypot(dx, dy) or 1
        sx = x1 + dx / dist * 62
        sy = y1 + dy / dist * 34
        tx = x2 - dx / dist * 62
        ty = y2 - dy / dist * 34
        parts.append(f'<line x1="{sx:.0f}" y1="{sy:.0f}" x2="{tx:.0f}" y2="{ty:.0f}" '
                     f'stroke="{color}" stroke-width="{3 if hot else 2}" '
                     f'marker-end="url(#{marker})"/>')
        if label:
            parts.append(f'<text x="{(x1 + x2) / 2:.0f}" y="{(y1 + y2) / 2 - 12:.0f}" '
                         f'font-size="13" fill="{color}" text-anchor="middle">{label}</text>')

    # 节点
    for nid, node in nodes.items():
        x, y = _svg_pos(*pos[nid])
        label_text, color = KIND_STYLE[node["kind"]]
        lines = _wrap(node["label"])
        box_w = 148 if node["kind"] != "RECEPTOR" else 132
        box_h = 44 + (len(lines) - 1) * 20
        is_hot = bool(highlight) and highlight.upper() == node.get("gene", "").upper()
        is_target = node.get("gene", "").upper() in targets
        parts.append(
            f'<rect x="{x - box_w / 2:.0f}" y="{y - box_h / 2:.0f}" width="{box_w}" '
            f'height="{box_h}" rx="12" fill="{color}" fill-opacity="0.92" '
            f'stroke="{"#E4572E" if is_hot else "#ffffff"}" stroke-width="'
            f'{3.5 if is_hot else 2}"/>')
        for i, line in enumerate(lines):
            parts.append(
                f'<text x="{x:.0f}" y="{y - (len(lines) - 1) * 10 + i * 20 + 4:.0f}" '
                f'text-anchor="middle" font-size="15" fill="#ffffff" '
                f'font-weight="bold">{"★ " if is_target else ""}{line}</text>')
        if is_target or is_hot:
            parts.append(f'<text x="{x:.0f}" y="{y + box_h / 2 + 17:.0f}" '
                         f'text-anchor="middle" font-size="12" fill="#E4572E">'
                         f'{label_text}'
                         f'{"｜本方靶点" if is_target else ""}</text>')

    parts.append(f'<text x="14" y="{h - 16}" font-size="15" fill="#666">'
                 f'{name}：{pw["desc"][:40]}</text>')
    parts.append("</svg>")

    legend = ""
    if show_legend:
        legend = ("<div style='padding:6px 4px;font-size:13px;color:#555'>"
                  + "　".join(
                      f"<span style='color:{c}'>■</span>{l}"
                      for l, c in KIND_STYLE.values())
                  + "　<span style='color:#E4572E'>★ 本方成分作用靶点</span></div>")

    return ('<html><body style="margin:0;background:#ffffff;font-family:'
            'Microsoft YaHei,sans-serif">' + legend + "".join(parts)
            + "</body></html>")
