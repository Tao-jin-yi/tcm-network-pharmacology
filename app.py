"""第三步：Streamlit 应用 —— 复方网络药理学分析平台。

运行方式（不要用 PyCharm 的 Run 按钮）：
    streamlit run app.py

界面由四部分组成：
1. 网络图：st-link-analysis（Cytoscape.js）渲染“中药-成分-靶点”网络；
2. 配伍协同分析：集合运算比较各味药的靶点覆盖、共享与独有靶点；
3. 核心靶点：CytoNCA 六指标迭代拓扑筛选（见 topology.py）；
4. 成分明细：活性成分的 OB / DL 与所属中药。
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd
import streamlit as st
from st_link_analysis import EdgeStyle, NodeStyle, st_link_analysis

from topology import cyto_nca_filter, elements_to_graph

HERE = Path(__file__).parent
DATA_FILE = HERE / "network.json"

LAYOUTS = ["cose", "fcose", "breadthfirst", "circle", "concentric", "grid", "cola", "dagre"]

st.set_page_config(page_title="复方网络药理学分析平台", page_icon=":pill:", layout="wide")
st.title("复方网络药理学分析平台")
st.caption("TCMSP 活性成分筛选（OB ≥ 30%，DL ≥ 0.18） → UniProt 基因符号标准化 → "
           "Cytoscape.js 网络可视化 → CytoNCA 拓扑筛选 → 配伍协同分析")


# ------------------------------------------------------------------ 数据读取
@st.cache_data(show_spinner=False)
def load_network(path_str: str) -> dict:
    """读取 Cytoscape 格式的网络数据（network.json）。"""
    path = Path(path_str)
    if not path.exists():
        st.error(f"找不到数据文件 {path}，请先运行 `python data_prep.py` 生成。")
        st.stop()
    return json.loads(path.read_text(encoding="utf-8"))


def filter_elements(elements: dict, ob_min: float, dl_min: float) -> dict:
    """按 OB / DL 阈值过滤成分，并删除因此孤立的中药与靶点。"""
    compounds = {n["data"]["id"]: n["data"] for n in elements["nodes"]
                 if n["data"]["label"] == "COMPOUND"}
    keep_compounds = {cid for cid, d in compounds.items()
                      if d.get("ob", 0) >= ob_min and d.get("dl", 0) >= dl_min}

    edges = [e for e in elements["edges"]
             if e["data"]["label"] == "CONTAINS" and e["data"]["target"] in keep_compounds]
    used_herbs = {e["data"]["source"] for e in edges}
    used_targets = {e["data"]["target"] for e in elements["edges"]
                    if e["data"]["label"] == "TARGETS" and e["data"]["source"] in keep_compounds}

    nodes = [n for n in elements["nodes"]
             if n["data"]["id"] in used_herbs
             or n["data"]["id"] in keep_compounds
             or n["data"]["id"] in used_targets]
    edges += [e for e in elements["edges"]
              if e["data"]["label"] == "TARGETS" and e["data"]["source"] in keep_compounds]
    return {"meta": elements.get("meta", {}), "nodes": nodes, "edges": edges}


def herb_target_map(elements: dict) -> Dict[str, Set[str]]:
    """每味中药 -> 其成分作用到的靶点集合。"""
    compound_herbs: Dict[str, Set[str]] = {}
    for edge in elements["edges"]:
        if edge["data"]["label"] == "CONTAINS":
            compound_herbs.setdefault(edge["data"]["target"], set()).add(
                edge["data"]["source"].replace("HERB::", ""))

    result: Dict[str, Set[str]] = {}
    for edge in elements["edges"]:
        if edge["data"]["label"] == "TARGETS":
            target = edge["data"]["target"].replace("TGT::", "")
            for herb in compound_herbs.get(edge["data"]["source"], set()):
                result.setdefault(herb, set()).add(target)
    return result


# ------------------------------------------------------------------ 侧边栏
raw_elements = load_network(str(DATA_FILE))
meta = raw_elements.get("meta", {})

with st.sidebar:
    st.header("参数")
    st.subheader("活性成分筛选")
    ob_min = st.slider("OB 口服生物利用度 ≥ (%)", 0.0, 100.0,
                       float(meta.get("ob_min", 30.0)), 1.0)
    dl_min = st.slider("DL 类药性 ≥", 0.0, 1.0,
                       float(meta.get("dl_min", 0.18)), 0.01)

    st.subheader("网络显示")
    layout = st.selectbox("布局算法", LAYOUTS, index=0)
    height = st.slider("画布高度 (px)", 400, 1000, 620, 20)

    st.subheader("CytoNCA 拓扑筛选")
    n_rounds = st.slider("迭代轮次", 1, 4, 1,
                         help="每轮保留 DC/BC/CC/EC/LAC/NC 六个指标均高于中位数的节点")
    st.markdown("数据来源：`network.json`（由 `data_prep.py` 生成）")

elements = filter_elements(raw_elements, ob_min, dl_min)

herbs = sorted({n["data"]["name"] for n in elements["nodes"]
                if n["data"]["label"] == "HERB"})
compounds_df = pd.DataFrame([n["data"] for n in elements["nodes"]
                             if n["data"]["label"] == "COMPOUND"])
n_targets = sum(1 for n in elements["nodes"] if n["data"]["label"] == "TARGET")

col1, col2, col3, col4 = st.columns(4)
col1.metric("中药", len(herbs))
col2.metric("活性成分", len(compounds_df))
col3.metric("靶点", n_targets)
col4.metric("关系（边）", len(elements["edges"]))

# ------------------------------------------------------------------ 主体
tab_graph, tab_synergy, tab_core, tab_compound = st.tabs(
    ["网络图", "配伍协同分析", "核心靶点（CytoNCA）", "成分明细"]
)

# ---------------------------------------------------------------- 1. 网络图
with tab_graph:
    if not elements["nodes"]:
        st.warning("当前阈值下没有符合条件的成分，请调低 OB / DL。")
    else:
        node_styles = [
            NodeStyle("HERB", "#FF7F3E", "name", "store"),
            NodeStyle("COMPOUND", "#2A629A", "name", "science"),
            NodeStyle("TARGET", "#3CB371", "name", "medical_services"),
        ]
        edge_styles = [
            EdgeStyle("CONTAINS", caption="name", directed=True),
            EdgeStyle("TARGETS", caption="name", directed=True),
        ]
        st_link_analysis(elements, layout, node_styles, edge_styles,
                         height=height, key="tcm_network")
        st.markdown(
            "**图例**："
            "<span style='color:#FF7F3E'>■ 中药</span>　"
            "<span style='color:#2A629A'>■ 成分</span>　"
            "<span style='color:#3CB371'>■ 靶点</span>　　"
            "点击节点高亮邻居，可拖拽调整布局。",
            unsafe_allow_html=True,
        )

# -------------------------------------------------------- 2. 配伍协同分析
with tab_synergy:
    herb_targets = herb_target_map(elements)
    if not herb_targets:
        st.warning("当前阈值下没有可分析的数据。")
    else:
        all_targets: Set[str] = set().union(*herb_targets.values())
        selected = st.selectbox("选择一味药", sorted(herb_targets))
        herb_set = herb_targets[selected]

        c1, c2, c3 = st.columns(3)
        c1.metric("覆盖靶点数", len(herb_set))
        c2.metric("靶点覆盖率", f"{len(herb_set & all_targets) / len(all_targets):.1%}")
        c3.metric("独有靶点数",
                  len(herb_set - set().union(*[v for k, v in herb_targets.items()
                                               if k != selected])))

        others = [h for h in sorted(herb_targets) if h != selected]
        compare = st.multiselect("与其他药对比", others, default=others[:1])
        if compare:
            rows = []
            for other in compare:
                shared = herb_set & herb_targets[other]
                union = herb_set | herb_targets[other]
                rows.append({
                    "对比药": other,
                    "共享靶点数": len(shared),
                    "Jaccard 相似度": round(len(shared) / len(union), 3) if union else 0.0,
                    "共享靶点": ", ".join(sorted(shared)[:30]),
                })
            st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)

            st.markdown("**共享靶点（交集）**")
            shared_all: Set[str] = set(herb_set)
            for other in compare:
                shared_all &= herb_targets[other]
            st.write(", ".join(sorted(shared_all)) if shared_all else "无")

        st.markdown("**各味药靶点数**")
        chart_df = pd.DataFrame(
            {"靶点数": [len(herb_targets[h]) for h in sorted(herb_targets)]},
            index=sorted(herb_targets),
        )
        st.bar_chart(chart_df)

        st.markdown("**该味药的靶点列表**")
        st.dataframe(pd.DataFrame({"靶点": sorted(herb_set)}),
                     hide_index=True, use_container_width=True)

# ------------------------------------------------------------- 3. 核心靶点
with tab_core:
    if not elements["nodes"]:
        st.warning("当前阈值下没有可分析的网络。")
    else:
        graph = elements_to_graph(elements)
        with st.spinner("正在计算 DC / BC / CC / EC / LAC / NC ..."):
            core_graph, rounds_df, metrics_df = cyto_nca_filter(graph, max_rounds=n_rounds)

        st.markdown("**迭代筛选过程**")
        st.dataframe(rounds_df, hide_index=True, use_container_width=True)

        core_targets = sorted(n.replace("TGT::", "") for n in core_graph.nodes
                              if str(n).startswith("TGT::"))
        st.metric("核心靶点数", len(core_targets))
        if core_targets:
            st.write(", ".join(core_targets))
            table = metrics_df.loc[[f"TGT::{t}" for t in core_targets]].copy()
            table.index = core_targets
            table = table.sort_values("DC", ascending=False).round(4)
            st.markdown("**核心靶点的拓扑指标**")
            st.dataframe(table, use_container_width=True)
            st.download_button(
                "下载核心靶点 CSV",
                data=table.reset_index().rename(columns={"index": "靶点"}).to_csv(index=False)
                .encode("utf-8-sig"),
                file_name="core_targets.csv",
                mime="text/csv",
            )

# ------------------------------------------------------------- 4. 成分明细
with tab_compound:
    if compounds_df.empty:
        st.warning("当前阈值下没有成分。")
    else:
        show = compounds_df[["name", "en", "ob", "dl", "herbs", "n_targets"]].copy()
        show.columns = ["成分", "英文名", "OB(%)", "DL", "所属中药", "靶点数"]
        show["所属中药"] = show["所属中药"].apply(lambda x: "、".join(x))
        show = show.sort_values(["靶点数", "OB(%)"], ascending=[False, False])
        st.dataframe(show, hide_index=True, use_container_width=True)
        st.download_button(
            "下载成分表 CSV",
            data=show.to_csv(index=False).encode("utf-8-sig"),
            file_name="compounds.csv",
            mime="text/csv",
        )
