"""把作用机制通路生成动画视频（MP4，失败时自动回退 GIF），便于学习观看。

动画逻辑：按信号级联顺序（配体 → 受体 → 接头 → 激酶 → 转录因子 → 效应）
逐节点揭示，新节点带缩放淡入，末帧停留并标注"★ 本方成分作用靶点"。

用法：
    python make_video.py --pathway "NF-κB 通路"     # 生成单条通路
    python make_video.py --all                      # 生成全部通路
    python make_video.py --overview                 # 生成总览视频（列出全部通路）
    python make_video.py --all --fps 24             # 指定帧率

依赖：matplotlib、imageio-ffmpeg（提供 ffmpeg 可执行文件）。
输出目录：videos/
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from matplotlib import font_manager
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter

from pathways import KIND_STYLE, PATHWAYS

HERE = Path(__file__).parent
NETWORK_FILE = HERE / "network.json"
VIDEO_DIR = HERE / "videos"

FONT_CANDIDATES = ["C:/Windows/Fonts/msyh.ttc", "C:/Windows/Fonts/simhei.ttf"]


def setup_font() -> None:
    for path in FONT_CANDIDATES:
        if Path(path).exists():
            font_manager.fontManager.addfont(path)
            name = font_manager.FontProperties(fname=path).get_name()
            plt.rcParams["font.family"] = name
            break
    plt.rcParams["axes.unicode_minus"] = False


def configure_ffmpeg() -> None:
    """用 imageio-ffmpeg 自带的 ffmpeg（matplotlib 3.11 的键名是 animation.ffmpeg_path）。"""
    try:
        import imageio_ffmpeg
        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        pass


def pathway_targets() -> set:
    """本方靶点基因符号集合（用于标注 ★）。"""
    if not NETWORK_FILE.exists():
        return set()
    elements = NETWORK_FILE.read_text(encoding="utf-8")
    import json
    nodes = json.loads(elements)["nodes"]
    return {n["data"]["name"].upper() for n in nodes if n["data"]["label"] == "TARGET"}


# ------------------------------------------------------------------ 布局计算
def layout(name: str) -> Tuple[Dict[str, Tuple[float, float]], List[str]]:
    """按级联深度分层布局，返回 (节点坐标, 揭示顺序)。"""
    pw = PATHWAYS[name]
    nodes = {n["id"]: n for n in pw["nodes"]}
    edges = pw["edges"]

    indeg = {nid: 0 for nid in nodes}
    for s, t, _ in edges:
        indeg[t] += 1

    depth: Dict[str, int] = {nid: 0 for nid in nodes}
    for _ in range(len(nodes) + 1):  # 简单松弛，取最长路径深度
        changed = False
        for s, t, _ in edges:
            if depth[t] < depth[s] + 1:
                depth[t] = depth[s] + 1
                changed = True
        if not changed:
            break

    layers: Dict[int, List[str]] = {}
    for nid, d in depth.items():
        layers.setdefault(d, []).append(nid)
    max_depth = max(layers) if layers else 0

    pos: Dict[str, Tuple[float, float]] = {}
    for d, members in layers.items():
        x = 1.1 + d * (8.4 / max(max_depth, 1))
        count = len(members)
        for i, nid in enumerate(sorted(members)):
            y = 5.6 - (i + 1) * (5.0 / (count + 1))
            pos[nid] = (x, y)

    order = [nid for nid in sorted(nodes, key=lambda k: (depth[k], k))]
    return pos, order


# -------------------------------------------------------------------- 绘图
def _wrap(text: str, width: int = 8) -> str:
    """节点名换行：优先按 ' / ' 分行，否则按宽度截断。"""
    if " / " in text:
        return "\n".join(text.split(" / "))
    if len(text) <= width:
        return text
    return "\n".join(text[i:i + width] for i in range(0, len(text), width))


def draw_frame(fig, ax, name: str, pos: Dict[str, Tuple[float, float]],
               order: List[str], step: int, progress: float,
               targets: set, subtitle: str = "") -> None:
    pw = PATHWAYS[name]
    nodes = {n["id"]: n for n in pw["nodes"]}
    visible = set(order[:step + 1])

    ax.clear()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7.2)
    ax.axis("off")

    stage = nodes[order[step]]
    ax.text(0.12, 6.85, f"{name}", fontsize=17, fontweight="bold", color="#222")
    ax.text(0.12, 6.45, subtitle or pw["desc"][:52], fontsize=10.5, color="#666")

    # 边
    for s, t, label in pw["edges"]:
        if t not in visible:
            continue
        x1, y1 = pos[s]
        x2, y2 = pos[t]
        is_new = (t == order[step])
        alpha = 0.25 + 0.7 * progress if is_new else 0.85
        ax.add_patch(mp.FancyArrowPatch(
            (x1 + 0.78, y1), (x2 - 0.82, y2),
            arrowstyle="-|>", mutation_scale=13, linewidth=1.8,
            color="#E4572E" if is_new else "#7a7a7a", alpha=alpha,
            connectionstyle="arc3,rad=0.08"))
        if label:
            ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.18, label, fontsize=8,
                    color="#E4572E" if is_new else "#888", ha="center",
                    alpha=alpha)

    # 节点
    for nid in visible:
        node = nodes[nid]
        x, y = pos[nid]
        label_text, color = KIND_STYLE[node["kind"]]
        is_new = (nid == order[step])
        scale = 0.72 + 0.28 * progress if is_new else 1.0
        alpha = 0.25 + 0.75 * progress if is_new else 1.0
        star = "★ " if node.get("gene", "").upper() in targets else ""
        wrapped = _wrap(node["label"])

        width, height = 1.7 * scale, (0.9 if "\n" in wrapped else 0.74) * scale
        ax.add_patch(mp.FancyBboxPatch(
            (x - width / 2, y - height / 2), width, height,
            boxstyle="round,pad=0.06,rounding_size=0.16",
            facecolor=color, edgecolor="#ffffff", linewidth=2.2, alpha=alpha,
            zorder=3))
        if is_new:
            ax.add_patch(mp.FancyBboxPatch(
                (x - width / 2 - 0.09, y - height / 2 - 0.09), width + 0.18,
                height + 0.18, boxstyle="round,pad=0.06,rounding_size=0.2",
                facecolor="none", edgecolor="#E4572E", linewidth=2.4,
                alpha=0.35 + 0.6 * progress, zorder=4))
        ax.text(x, y + 0.02, f"{star}{wrapped}", ha="center", va="center",
                fontsize=9.5, color="#ffffff", fontweight="bold", alpha=alpha,
                zorder=5)
        ax.text(x, y - height / 2 - 0.14, label_text, ha="center", fontsize=7.5,
                color="#777", alpha=alpha)

    # 图例
    handles = [mp.Patch(facecolor=color, edgecolor="none", label=label)
               for label, color in KIND_STYLE.values()]
    ax.legend(handles=handles, loc="lower left", fontsize=8.5, ncol=6,
              frameon=False, bbox_to_anchor=(0.0, -0.02))
    ax.text(9.9, 0.12, "★ = 本方成分作用靶点", fontsize=8.5, color="#E4572E",
            ha="right")
    ax.text(9.9, 6.85, f"{step + 1}/{len(order)}", fontsize=10, color="#999",
            ha="right")


def safe_name(name: str) -> str:
    return re.sub(r"[^\w\u4e00-\u9fa5-]", "_", name)


def make_video(name: str, fps: int = 20, out_dir: Path = VIDEO_DIR,
               targets: Optional[set] = None) -> Optional[Path]:
    """生成单条通路的动画视频。"""
    setup_font()
    out_dir.mkdir(exist_ok=True)
    targets = targets if targets is not None else pathway_targets()
    pos, order = layout(name)

    fig, ax = plt.subplots(figsize=(12, 7.2), dpi=110)
    fig.patch.set_facecolor("#ffffff")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.06)

    grow, hold = 10, 12  # 每个节点的淡入帧数与停留帧数
    frames_per_step = grow + hold
    total = len(order) * frames_per_step + 40  # 末尾额外停留（全图）

    def update(i: int) -> None:
        if i < len(order) * frames_per_step:
            step = i // frames_per_step
            progress = min(1.0, (i % frames_per_step) / grow)
        else:
            step = len(order) - 1
            progress = 1.0
        draw_frame(fig, ax, name, pos, order, step, progress, targets)

    draw_frame(fig, ax, name, pos, order, 0, 0.0, targets)
    anim = FuncAnimation(fig, update, frames=total, interval=1000 / fps)

    mp4 = out_dir / f"{safe_name(name)}.mp4"
    try:
        configure_ffmpeg()
        anim.save(mp4, writer=FFMpegWriter(fps=fps, codec="libx264",
                                           extra_args=["-pix_fmt", "yuv420p"]))
    except Exception:
        try:
            gif = out_dir / f"{safe_name(name)}.gif"
            anim.save(gif, writer=PillowWriter(fps=max(6, fps // 3)))
            plt.close(fig)
            return gif
        except Exception as exc:
            print(f"  {name} 视频生成失败：{exc}")
            plt.close(fig)
            return None
    plt.close(fig)
    return mp4


def make_overview(fps: int = 20, out_dir: Path = VIDEO_DIR) -> Optional[Path]:
    """总览视频：依次列出全部通路名与其关键节点。"""
    setup_font()
    configure_ffmpeg()
    out_dir.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 7.2), dpi=110)
    fig.patch.set_facecolor("#ffffff")
    names = list(PATHWAYS)
    per = 30

    def update(i: int) -> None:
        idx = min(i // per, len(names) - 1)
        ax.clear()
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 7.2)
        ax.axis("off")
        ax.text(0.5, 6.5, "中药作用机制通路总览", fontsize=22, fontweight="bold",
                color="#222")
        ax.text(0.5, 6.05, f"共 {len(names)} 条通路 · 配体→受体→激酶→转录因子→效应",
                fontsize=12, color="#666")
        for j, nm in enumerate(names):
            y = 5.4 - j * 0.42
            active = j == idx
            pw = PATHWAYS[nm]
            genes = [n["gene"] for n in pw["nodes"] if n.get("gene")]
            ax.add_patch(mp.FancyBboxPatch(
                (0.4, y - 0.17), 9.2, 0.34,
                boxstyle="round,pad=0.03,rounding_size=0.08",
                facecolor="#2A629A" if active else "#dfe6ef",
                alpha=0.95 if active else 0.6, edgecolor="none"))
            ax.text(0.6, y, f"{j + 1}. {nm}", fontsize=12, color="#ffffff" if active
                    else "#333", va="center", fontweight="bold" if active else "normal")
            ax.text(6.0, y, "关键节点：" + "、".join(genes[:6]), fontsize=9.5,
                    color="#ffffff" if active else "#555", va="center")

    update(0)
    anim = FuncAnimation(fig, update, frames=len(names) * per + 40,
                         interval=1000 / fps)
    out = out_dir / "00_通路总览.mp4"
    try:
        anim.save(out, writer=FFMpegWriter(fps=fps, codec="libx264",
                                           extra_args=["-pix_fmt", "yuv420p"]))
    except Exception:
        out = out_dir / "00_通路总览.gif"
        anim.save(out, writer=PillowWriter(fps=max(6, fps // 3)))
    plt.close(fig)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pathway", help="通路名（与 pathways.py 中的键一致）")
    parser.add_argument("--all", action="store_true", help="生成全部通路")
    parser.add_argument("--overview", action="store_true", help="生成总览视频")
    parser.add_argument("--fps", type=int, default=20)
    args = parser.parse_args()

    if not (args.pathway or args.all or args.overview):
        print("可用通路：")
        for name in PATHWAYS:
            print("  -", name)
        print("示例：python make_video.py --pathway \"NF-κB 通路\"")
        return

    made: List[Path] = []
    if args.overview:
        out = make_overview(args.fps)
        if out:
            made.append(out)
            print("已生成", out)

    names = list(PATHWAYS) if args.all else ([args.pathway] if args.pathway else [])
    for name in names:
        if name not in PATHWAYS:
            print("未找到通路：", name)
            continue
        out = make_video(name, args.fps)
        if out:
            made.append(out)
            print("已生成", out)

    print(f"完成：{len(made)} 个视频，输出目录 {VIDEO_DIR}")


if __name__ == "__main__":
    main()
