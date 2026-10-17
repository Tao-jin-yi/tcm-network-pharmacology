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

import altair as alt
import numpy as np
import pandas as pd
import subprocess
import sys
from urllib.parse import quote
import requests
import streamlit as st
import streamlit.components.v1 as components
from st_link_analysis import EdgeStyle, NodeStyle, st_link_analysis

from anatomy import SPECIES, render_anatomy
from animal_view import (compartment_of, organ_target_map, render_cell)
from cell_pathway import render_cell_pathway
from formulas import (DEFAULT_FORMULA, NOTES, REFERENCES, formula_coverage,
                      formula_herbs, formula_names, formula_rows, get_formula)
from herb_info import HERB_INFO, get_herb_info
from pathways import (KIND_STYLE, PATHWAYS, pathway_elements, pathways_of_target,
                      target_class, targets_in_pathway)
from pkpd import (fit_dose_response, fit_pk, hill_emax, pk_metrics,
                  simulate_dose_response, simulate_effect_time, simulate_pk)
from topology import cyto_nca_filter, elements_to_graph

HERE = Path(__file__).parent
DATA_FILE = HERE / "network.json"
TARGET_INFO_FILE = HERE / "target_info.json"  # 由 fetch_target_info.py 生成（UniProt）
LITERATURE_FILE = HERE / "literature.json"    # 由 fetch_literature.py 生成（PubMed）

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


@st.cache_data(show_spinner=False)
def load_target_info() -> dict:
    """读取 UniProt 蛋白信息（蛋白名、登录号、亚细胞定位、功能）。"""
    if not TARGET_INFO_FILE.exists():
        return {}
    return json.loads(TARGET_INFO_FILE.read_text(encoding="utf-8"))


@st.cache_data(show_spinner="正在从 AlphaFold 下载蛋白结构 ...")
def load_alphafold_pdb(accession: str) -> str:
    """按 UniProt 登录号取 AlphaFold 预测的 PDB 文本，失败返回空串。"""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb"
    try:
        resp = requests.get(url, timeout=45)
    except Exception:
        return ""
    return resp.text if resp.status_code == 200 and "ATOM" in resp.text else ""


def render_protein_3d(pdb_text: str, height: int = 430) -> None:
    """用 3Dmol.js 在页面内渲染蛋白三维结构（PDB 内容直接内联，无需浏览器联网）。"""
    html = f"""
<div id="protein-viewer" style="height:{height}px; width:100%; background:#ffffff;"></div>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<script>
  const viewer = $3Dmol.createViewer(document.getElementById("protein-viewer"),
                                     {{backgroundColor: "white"}});
  viewer.addModel({json.dumps(pdb_text)}, "pdb");
  viewer.setStyle({{}}, {{cartoon: {{color: "spectrum"}}}});
  viewer.zoomTo();
  viewer.render();
</script>
"""
    components.html(html, height=height + 25)


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

    st.subheader("分析方剂")
    formula_choice = st.selectbox(
        "选择复方", ["全方（全部药材）"] + formula_names(),
        help="切换到某首方剂时，网络与后续分析会只保留该方剂所含药材相关的"
             "成分与靶点；暂无成分数据的方剂只展示方解。")
    st.caption(f"复方库共 {len(formula_names())} 首方剂（formulas.py）")

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

def filter_by_formula(elements: dict, herb_names: List[str]) -> dict:
    """只保留指定方剂药材相关的成分与靶点（CONTAINS / TARGETS 边）。"""
    if not herb_names:
        return elements
    herb_ids = {f"HERB::{h}" for h in herb_names}
    keep = set(herb_ids)
    for edge in elements["edges"]:
        data = edge["data"]
        src, tgt = data.get("source", ""), data.get("target", "")
        if data.get("label") == "CONTAINS":
            if src in herb_ids:
                keep.add(tgt)
            elif tgt in herb_ids:
                keep.add(src)
    compounds = {n for n in keep if n.startswith("COMP::")}
    for edge in elements["edges"]:
        data = edge["data"]
        if data.get("label") == "TARGETS" and data.get("source", "") in compounds:
            keep.add(data["target"])
    return {
        "meta": elements.get("meta", {}),
        "nodes": [n for n in elements["nodes"] if n["data"]["id"] in keep],
        "edges": [e for e in elements["edges"]
                  if e["data"]["source"] in keep and e["data"]["target"] in keep],
    }


elements = filter_elements(raw_elements, ob_min, dl_min)
selected_formula = None if formula_choice.startswith("全方") else formula_choice
if selected_formula:
    elements = filter_by_formula(elements, formula_herbs(selected_formula))

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
tab_formula, tab_graph, tab_synergy, tab_core, tab_target, tab_lit, tab_pk, tab_animal, tab_compound = st.tabs(
    ["复方说明", "网络图", "配伍协同分析", "核心靶点（CytoNCA）",
     "靶点详情（蛋白 / 定位 / 结构）",
     "文献机制（PubMed）", "药动学 / 量-效-时", "动物模型模拟实验", "成分明细"]
)
core_targets: List[str] = []

# -------------------------------------------------------------- 0. 复方说明
with tab_formula:
    shown = selected_formula or DEFAULT_FORMULA
    formula = get_formula(shown)
    st.subheader(shown)
    st.caption(f"出处：{formula['source']}　|　"
               f"复方库共 {len(formula_names())} 首，可在左侧「分析方剂」切换")

    hit, total, hit_herbs = formula_coverage(shown, herbs)
    f1, f2, f3, f4 = st.columns(4)
    f1.metric("方中药材", total)
    f2.metric("有成分数据的药材", f"{hit}/{total}")
    f3.metric("预测靶点（当前范围）", n_targets)
    f4.metric("关系（边）", len(elements["edges"]))
    if hit == 0:
        st.warning("该方剂的药材在数据库中暂无成分-靶点数据，本页仅展示方剂说明；"
                   "如需纳入分析，请用 `import_tcmsp.py` 导入成分与靶点后重建网络。")
    elif hit < total:
        st.info(f"该方剂中 {hit}/{total} 味药材有成分数据：{'、'.join(hit_herbs)}；"
                f"缺：{'、'.join([h for h in formula_herbs(shown) if h not in hit_herbs])}。"
                f"网络分析仅覆盖有数据的药材。")

    st.markdown(f"**组成**：" + "、".join(
        f"{h} {d} g" for h, d in formula["dosage"].items()))
    st.markdown(f"**功用**：{formula['action']}")
    st.markdown(f"**主治**：{formula['indication']}")
    st.markdown(f"**辨证要点**：{formula['syndrome']}")
    st.markdown(f"**用法**：{formula['usage']}")

    st.markdown("### 组成与方解（君臣佐使）")
    st.dataframe(pd.DataFrame(formula_rows(shown)), hide_index=True,
                 use_container_width=True)

    for role, members, explain in formula["roles"]:
        with st.expander(f"{role}：{'、'.join(members)}", expanded=(role == "君药")):
            st.write(explain)
            st.caption("剂量：" + "；".join(
                f"{m} {formula['dosage'].get(m, '—')} g" for m in members))

    st.markdown("### 配伍特点")
    for title, text in formula["compatibility"]:
        st.markdown(f"- **{title}**：{text}")

    st.markdown("### 现代药理研究方向（文献报道）")
    for title, text in formula["pharmacology"]:
        st.markdown(f"- **{title}**：{text}")
    st.caption("以上为单味药 / 同类方已有研究报道的方向，"
               "本平台的网络分析结果为预测，不等同于实验或临床结论。")

    st.markdown("### 使用注意")
    for item in formula["precautions"]:
        st.markdown(f"- {item}")

    st.markdown("### 说明与数据来源")
    for item in NOTES:
        st.caption(f"- {item}")
    with st.expander("参考来源"):
        for ref in REFERENCES:
            st.markdown(f"- {ref}")

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

        info = get_herb_info(selected)
    if info:
        st.markdown(f"**{selected}**　{info['nature']}　{info['meridian']}　"
                    f"方中角色：{info['role']}　分类：{info['category']}")
        st.caption(f"功效主治：{info['action']}")
    else:
        st.caption("该药尚未录入功效主治，可在 herb_info.py 中补充。")


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

        with st.expander("全方药材用途一览（性味归经 / 功效主治 / 君臣佐使）"):
            rows = []
            for herb_name in sorted(herb_targets):
                item = get_herb_info(herb_name)
                if item:
                    rows.append({"中药": herb_name, "性味": item["nature"],
                                 "归经": item["meridian"], "分类": item["category"],
                                 "角色": item["role"], "功效主治": item["action"]})
            if rows:
                st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
            else:
                st.caption("尚未录入药材用途，可在 herb_info.py 中补充。")

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
            info = load_target_info()
            if info:
                table.insert(0, "蛋白名",
                             [info.get(t, {}).get("protein", "") for t in table.index])
                table.insert(1, "UniProt",
                             [info.get(t, {}).get("accession", "") for t in table.index])
                table.insert(2, "亚细胞定位",
                             ["、".join(info.get(t, {}).get("locations", [])[:2])
                              for t in table.index])
            st.markdown("**核心靶点的拓扑指标**")
            st.dataframe(table, use_container_width=True)
            st.download_button(
                "下载核心靶点 CSV",
                data=table.reset_index().rename(columns={"index": "靶点"}).to_csv(index=False)
                .encode("utf-8-sig"),
                file_name="core_targets.csv",
                mime="text/csv",
            )

# --------------------------------- 4. 靶点详情：蛋白信息 / 亚细胞定位 / 三维结构
with tab_target:
    info = load_target_info()
    if not info:
        st.info("未找到 target_info.json，请先运行 `python fetch_target_info.py` "
                "从 UniProt 抓取蛋白信息。")
    else:
        all_t = sorted({n["data"]["name"] for n in elements["nodes"]
                        if n["data"]["label"] == "TARGET"})
        options = core_targets + [t for t in all_t if t not in core_targets]
        picked = st.selectbox("选择靶点（核心靶点排在前面）", options)

        item = info.get(picked.upper())
        if not item:
            st.warning(f"UniProt 中暂无 {picked} 的人源已审阅记录。")
        else:
            st.markdown(f"### {picked}　{item['protein']}")
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("UniProt 登录号", item["accession"] or "-")
            c2.metric("物种", item["organism"] or "Homo sapiens")
            c3.metric("氨基酸长度", item["length"] or "-")
            c4.metric("是否核心靶点", "是" if picked in core_targets else "否")
            st.markdown(f"[在 UniProt 打开 {item['accession']}]({item['uniprot_url']})")

            if item["locations"]:
                st.markdown("**亚细胞定位（UniProt 注释）**")
                st.markdown("　".join(f"`{loc}`" for loc in item["locations"]))
            if item["function"]:
                st.markdown("**功能描述**")
                st.caption(item["function"])

            st.markdown("**蛋白三维结构**（AlphaFold DB 预测模型，可旋转缩放）")
            pdb_text = load_alphafold_pdb(item["accession"])
            if pdb_text:
                render_protein_3d(pdb_text)
                st.download_button(
                    "下载 PDB 结构文件",
                    data=pdb_text,
                    file_name=f"AF-{item['accession']}-F1-model_v6.pdb",
                    mime="chemical/x-pdb",
                )
            else:
                st.info("AlphaFold 暂无该蛋白的预测结构。")

        if core_targets:
            st.markdown("**核心靶点的亚细胞定位分布**")
            rows = []
            for t in core_targets:
                locs = info.get(t.upper(), {}).get("locations") or ["未标注"]
                for loc in locs:
                    rows.append({"靶点": t, "定位": loc.split(",")[0]})
            if rows:
                dist = (pd.DataFrame(rows)
                        .groupby("定位")["靶点"].nunique()
                        .sort_values(ascending=False))
                st.bar_chart(dist)
                st.caption("同一蛋白可定位于多个区室，故计数之和大于核心靶点数。")

# ---------------------------------------------- 5. 文献机制（PubMed 抓取与抽取）
with tab_lit:
    if not LITERATURE_FILE.exists():
        st.info("未找到 literature.json，请先运行 `python fetch_literature.py` "
                "从 PubMed 抓取文献。")
    else:
        lit = json.loads(LITERATURE_FILE.read_text(encoding="utf-8"))
        herb_options = ["全部中药"] + sorted(lit)
        picked_herb = st.selectbox("选择中药", herb_options)
        scope = (lit if picked_herb == "全部中药"
                 else {picked_herb: lit.get(picked_herb, [])})
        articles = [dict(a, herb=h) for h, arts in scope.items() for a in arts]

        if not articles:
            st.warning("该药暂无抓取到的文献。")
        else:
            c1, c2, c3 = st.columns(3)
            c1.metric("文献数", len(articles))
            c2.metric("含机制标注", sum(1 for a in articles if a.get("mechanisms")))
            c3.metric("涉及靶点数",
                      len({t for a in articles for t in a.get("targets", [])}))

            mech_rows = [{"机制": m} for a in articles for m in a.get("mechanisms", [])]
            if mech_rows:
                st.markdown("**机制 / 通路分布**")
                mech_count = (pd.DataFrame(mech_rows).value_counts("机制")
                              .rename("文献数").reset_index()
                              .sort_values("文献数", ascending=False))
                st.altair_chart(
                    alt.Chart(mech_count).mark_bar(cornerRadiusEnd=3).encode(
                        x=alt.X("文献数:Q"), y=alt.Y("机制:N", sort="-x"),
                        color=alt.Color("文献数:Q", scale=alt.Scale(scheme="oranges")),
                        tooltip=["机制", "文献数"],
                    ).properties(height=max(160, 26 * len(mech_count))),
                    use_container_width=True)

            target_rows = [{"靶点": t} for a in articles for t in a.get("targets", [])]
            if target_rows:
                st.markdown("**文献中出现的靶点（与本方靶点交集）**")
                t_count = (pd.DataFrame(target_rows).value_counts("靶点")
                           .rename("文献数").reset_index()
                           .sort_values("文献数", ascending=False).head(20))
                st.altair_chart(
                    alt.Chart(t_count).mark_bar(cornerRadiusEnd=3).encode(
                        x=alt.X("文献数:Q"), y=alt.Y("靶点:N", sort="-x"),
                        color=alt.Color("文献数:Q", scale=alt.Scale(scheme="teals")),
                        tooltip=["靶点", "文献数"],
                    ).properties(height=max(160, 24 * len(t_count))),
                    use_container_width=True)

            st.markdown("**文献列表**")
            table = pd.DataFrame([{
                "中药": a["herb"], "年份": a.get("year"), "期刊": a.get("journal"),
                "题目": a.get("title"), "机制": "、".join(a.get("mechanisms", [])),
                "靶点": "、".join(a.get("targets", [])),
                "PMID": a.get("pmid"),
            } for a in articles]).sort_values("年份", ascending=False)
            st.dataframe(table, hide_index=True, use_container_width=True)

            # -------------------------------------------- 一键直达原文链接
            st.markdown("### 打开原文（一键直达）")
            st.caption("链接由 Europe PMC 解析（DOI / PMCID / 开放获取全文）；"
                       "若期刊需订阅，请用学术检索找免费版本。")
            titles = [(f"{a.get('year')}　{a.get('title')[:60]}", a) for a in articles]
            picked_title = st.selectbox("选择文献", [t for t, _ in titles],
                                        key="lit_pick")
            picked_article = next(a for t, a in titles if t == picked_title)
            pmid = str(picked_article.get("pmid") or "")
            links = [
                ("PubMed 摘要",
                 picked_article.get("url") or
                 (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None)),
                ("Europe PMC",
                 picked_article.get("epmc_url") or
                 (f"https://europepmc.org/article/MED/{pmid}" if pmid else None)),
                ("DOI 原文页", picked_article.get("doi_url")),
                ("PMC 全文页", picked_article.get("pmc_url")),
                ("开放获取 PDF", picked_article.get("fulltext_url")),
                ("学术检索（找免费版）",
                 "https://scholar.google.com/scholar?q=" +
                 quote(str(picked_article.get("title") or ""))),
            ]
            available = [(label, url) for label, url in links if url]
            if picked_article.get("is_oa"):
                st.success("该文献为开放获取（Open Access），可直接打开全文。")
            cols = st.columns(min(3, max(1, len(available))))
            for i, (label, url) in enumerate(available):
                cols[i % len(cols)].link_button(label, url)
            st.caption(f"PMID {pmid}"
                       + (f"　DOI {picked_article['doi']}"
                          if picked_article.get("doi") else "")
                       + (f"　PMCID {picked_article['pmcid']}"
                          if picked_article.get("pmcid") else ""))

            with st.expander("查看摘要（前 5 篇）"):
                for a in articles[:5]:
                    st.markdown(f"**{a.get('title')}**　{a.get('year')}　{a.get('journal')}")
                    st.caption(f"机制：{'、'.join(a.get('mechanisms', [])) or '未标注'}　"
                               f"靶点：{'、'.join(a.get('targets', [])) or '未命中'}")
                    st.write(a.get("abstract", "")[:600] + " ...")
                    st.markdown(f"[PubMed {a.get('pmid')}]({a.get('url')})")

# --------------------------------------- 6. 药动学（PK）与量-效-时（PD）曲线
with tab_pk:
    st.caption("一室一级吸收 PK 模型 + Hill 量-效方程；可用参数正向模拟，"
               "也可上传实验测得数据反推参数。")
    mode = st.radio("模式", ["参数模拟", "上传实验数据拟合"], horizontal=True)

    if mode == "参数模拟":
        st.markdown("**药动学参数**")
        p1, p2, p3, p4, p5 = st.columns(5)
        dose = p1.number_input("剂量 D (mg/kg)", 1.0, 1000.0, 100.0, 10.0)
        f_bio = p2.number_input("生物利用度 F", 0.01, 1.0, 0.40, 0.05)
        vd = p3.number_input("表观分布容积 V (L/kg)", 0.1, 100.0, 5.0, 0.5)
        ka = p4.number_input("吸收速率 ka (1/h)", 0.05, 5.0, 1.20, 0.05)
        ke = p5.number_input("消除速率 ke (1/h)", 0.01, 2.0, 0.25, 0.01)

        st.markdown("**量-效参数**")
        q1, q2, q3, q4 = st.columns(4)
        e0 = q1.number_input("基线效应 E0 (%)", -100.0, 100.0, 5.0, 1.0)
        emax = q2.number_input("最大效应 Emax (%)", 0.0, 200.0, 95.0, 5.0)
        ec50 = q3.number_input("半数有效量 EC50 (mg/L)", 0.1, 500.0, 15.0, 1.0)
        hill = q4.number_input("Hill 系数 h", 0.2, 8.0, 1.2, 0.1)

        pk_df, cmax, tmax = simulate_pk(dose, f_bio, vd, ka, ke)
        metrics = pk_metrics(ka, ke, dose, f_bio, vd)
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Cmax (mg/L)", f"{metrics['Cmax']:.2f}")
        m2.metric("Tmax (h)", f"{metrics['Tmax']:.2f}")
        m3.metric("AUC₀₋t (mg·h/L)", f"{metrics['AUC0-t']:.1f}")
        m4.metric("消除半衰期 (h)", f"{metrics['t1/2(ke)']:.2f}")

        col_a, col_b = st.columns(2)
        with col_a:
            st.markdown("**血药浓度 - 时间曲线（PK）**")
            st.altair_chart(
                alt.Chart(pk_df).mark_line(color="#2A629A").encode(
                    x=alt.X("时间 (h):Q"), y=alt.Y("血药浓度 (mg/L):Q")),
                use_container_width=True)
            st.markdown("**剂量 - 效应曲线（量-效）**")
            dr_df = simulate_dose_response(e0, emax, ec50, hill,
                                           d_max=max(200.0, dose * 3))
            st.altair_chart(
                alt.Chart(dr_df).mark_line(color="#FF7F3E").encode(
                    x=alt.X("剂量 (mg/kg):Q"), y=alt.Y("效应 (%):Q")),
                use_container_width=True)
        with col_b:
            st.markdown("**效应 - 时间曲线（时-效，PK 驱动）**")
            et_df = simulate_effect_time(dose, f_bio, vd, ka, ke, e0, emax, ec50, hill)
            st.altair_chart(
                alt.Chart(et_df).mark_line(color="#3CB371").encode(
                    x=alt.X("时间 (h):Q"), y=alt.Y("效应 (%):Q")),
                use_container_width=True)
            st.markdown("**多剂量 PK 比较**")
            multi = []
            for d in [dose * 0.5, dose, dose * 2]:
                frame, _, _ = simulate_pk(d, f_bio, vd, ka, ke)
                frame["剂量"] = f"{d:.0f} mg/kg"
                multi.append(frame)
            st.altair_chart(
                alt.Chart(pd.concat(multi)).mark_line().encode(
                    x=alt.X("时间 (h):Q"), y=alt.Y("血药浓度 (mg/L):Q"),
                    color=alt.Color("剂量:N")),
                use_container_width=True)
    else:
        st.markdown("上传 CSV：一列自变量（剂量或时间），一列因变量（效应或浓度）。")
        up = st.file_uploader("实验数据 CSV", type=["csv"])
        if up is not None:
            data = pd.read_csv(up)
            st.dataframe(data, hide_index=True, use_container_width=True)
            kind = st.radio("数据类型", ["量-效（剂量 vs 效应）", "药-时（时间 vs 浓度）"],
                            horizontal=True)
            xcol = st.selectbox("自变量列", data.columns)
            ycol = st.selectbox("因变量列", data.columns)
            if st.button("拟合曲线"):
                if kind.startswith("量"):
                    params, r2 = fit_dose_response(data[xcol], data[ycol])
                else:
                    params, r2 = fit_pk(data[xcol], data[ycol])
                if not params:
                    st.error("拟合失败：数据点太少或噪声过大（至少需要 4 个点）。")
                else:
                    st.success(f"拟合完成，R² = {r2:.4f}")
                    st.write({k: round(float(v), 4) for k, v in params.items()})
                    xs = np.linspace(float(data[xcol].min()), float(data[xcol].max()), 200)
                    if kind.startswith("量"):
                        ys = hill_emax(xs, params["E0"], params["Emax"],
                                       params["EC50"], params["h"])
                    else:
                        ys = params["A"] * (np.exp(-params["ke"] * xs)
                                            - np.exp(-params["ka"] * xs))
                    fit_df = pd.DataFrame({xcol: xs, "拟合值": ys})
                    base = alt.Chart(data).mark_circle(size=70, color="#2A629A").encode(
                        x=f"{xcol}:Q", y=f"{ycol}:Q")
                    line = alt.Chart(fit_df).mark_line(color="#FF7F3E").encode(
                        x=f"{xcol}:Q", y="拟合值:Q")
                    st.altair_chart(base + line, use_container_width=True)

    st.markdown("---")
    st.markdown("### 实验动物模型")
    st.caption("默认示例为可编辑模板，请按真实实验方案修改；时间轴与分组可直接导出。")
    a1, a2, a3, a4 = st.columns(4)
    strain = a1.text_input("动物品系", "C57BL/6J 小鼠（雄性，8 周龄）")
    model = a2.text_input("造模方法", "高脂饮食（HFD）诱导 8 周")
    route = a3.text_input("给药途径", "灌胃（ig），每日 1 次")
    endpoint = a4.text_input("主要检测指标", "体重、血脂四项、肝脏 TNF/IL6、HE 染色")

    st.markdown(
        f"<div style='display:flex;gap:10px;flex-wrap:wrap'>"
        f"<div style='padding:8px 12px;border:1px solid #ddd;border-radius:6px'>"
        f"品系<br><b>{strain}</b></div>"
        f"<div style='padding:8px 12px;border:1px solid #ddd;border-radius:6px'>"
        f"造模<br><b>{model}</b></div>"
        f"<div style='padding:8px 12px;border:1px solid #ddd;border-radius:6px'>"
        f"给药<br><b>{route}</b></div>"
        f"<div style='padding:8px 12px;border:1px solid #ddd;border-radius:6px'>"
        f"终点指标<br><b>{endpoint}</b></div></div>",
        unsafe_allow_html=True,
    )

    stages = pd.DataFrame([
        {"阶段": "适应期", "开始": -7, "结束": 0},
        {"阶段": "造模期", "开始": 0, "结束": 56},
        {"阶段": "给药干预", "开始": 28, "结束": 56},
        {"阶段": "取材与检测", "开始": 56, "结束": 60},
    ])
    st.altair_chart(
        alt.Chart(stages).mark_bar(size=26, cornerRadius=4).encode(
            x=alt.X("开始:Q", title="实验天数（天）"), x2="结束:Q",
            y=alt.Y("阶段:N", sort=None, title=None),
            color=alt.Color("阶段:N", legend=None),
            tooltip=["阶段", "开始", "结束"],
        ).properties(height=180),
        use_container_width=True)

    default_groups = pd.DataFrame([
        {"组别": "空白对照组", "动物数": 12, "造模": "普通饲料", "剂量 (mg/kg)": 0,
         "给药": "溶剂灌胃"},
        {"组别": "模型组", "动物数": 12, "造模": "高脂饮食", "剂量 (mg/kg)": 0,
         "给药": "溶剂灌胃"},
        {"组别": "受试药低剂量", "动物数": 12, "造模": "高脂饮食", "剂量 (mg/kg)": 50,
         "给药": "灌胃"},
        {"组别": "受试药中剂量", "动物数": 12, "造模": "高脂饮食", "剂量 (mg/kg)": 100,
         "给药": "灌胃"},
        {"组别": "受试药高剂量", "动物数": 12, "造模": "高脂饮食", "剂量 (mg/kg)": 200,
         "给药": "灌胃"},
        {"组别": "阳性药组", "动物数": 12, "造模": "高脂饮食", "剂量 (mg/kg)": 10,
         "给药": "灌胃"},
    ])
    groups = st.data_editor(default_groups, num_rows="dynamic",
                            use_container_width=True)
    st.altair_chart(
        alt.Chart(groups).mark_bar(cornerRadiusEnd=3).encode(
            x=alt.X("剂量 (mg/kg):Q"), y=alt.Y("组别:N", sort=None),
            color=alt.Color("动物数:Q", scale=alt.Scale(scheme="blues")),
            tooltip=list(groups.columns),
        ).properties(height=220),
        use_container_width=True)
    st.download_button("下载分组方案 CSV",
                       data=groups.to_csv(index=False).encode("utf-8-sig"),
                       file_name="animal_groups.csv", mime="text/csv")

# --------------------------------------- 7. 动物模型模拟实验（器官 / 受体 / 通路）
with tab_animal:
    info = load_target_info()
    if not info:
        st.info("未找到 target_info.json，请先运行 `python fetch_target_info.py`。")
    else:
        st.caption("器官点位依据 UniProt 组织特异性注释，受体点位依据亚细胞定位注释，"
                   "通路图依据 pathways.py 中整理的信号通路拓扑。")
        all_genes = sorted({n["data"]["name"] for n in elements["nodes"]
                            if n["data"]["label"] == "TARGET"})
        scope = st.radio("靶点范围", ["核心靶点", "全方靶点"], horizontal=True)
        genes = core_targets if (scope == "核心靶点" and core_targets) else all_genes

        subset = {g: info[g.upper()] for g in genes if g.upper() in info}
        o_map = organ_target_map(subset)
        counts = {organ: len(v) for organ, v in o_map.items()}

        left, right = st.columns([1.2, 1])
        animal = left.selectbox("实验动物", list(SPECIES))
        organs_sorted = sorted(counts, key=lambda o: -counts[o])
        picked_organ = right.selectbox("选择器官 / 组织", organs_sorted or ["（暂无）"])
        with left:
            components.html(render_anatomy(animal, counts, picked_organ), height=430)
            st.caption("解剖层级示意：心脏四腔、肺叶、肝分叶、豆形肾、小肠盘曲、盲肠、"
                       "骨骼（脊柱/肋骨/骨盆）与大血管均按解剖比例绘制；"
                       "器官颜色深浅 = 该器官相关靶点数，红框为当前选中器官。")

        with right:
            organ_genes = o_map.get(picked_organ, [])
            st.markdown(f"**{picked_organ}** 相关靶点：{len(organ_genes)} 个")
            if organ_genes:
                st.dataframe(pd.DataFrame([{
                    "靶点": g,
                    "蛋白": subset[g].get("protein", ""),
                    "类别": target_class(g, subset[g]),
                    "亚细胞定位": "、".join(subset[g].get("locations", [])[:2]),
                    "组织特异性": subset[g].get("tissue", "")[:60],
                } for g in organ_genes]), hide_index=True, use_container_width=True)
            else:
                st.caption("该器官在当前靶点范围内没有 UniProt 组织特异性注释支持。")

        st.markdown("---")
        st.markdown("### 器官 → 细胞 → 通路（两级联动）")
        st.caption("选中的器官 → 该器官表达的本方靶点 → 这些靶点参与的信号通路；"
                   "★ 为该器官表达的靶点。")
        organ_upper = {g.upper() for g in organ_genes}
        if not organ_upper:
            st.caption("该器官暂无可映射靶点（UniProt 组织特异性未标注），无法展示通路。")
        else:
            pw_hits = []
            for pw_name in PATHWAYS:
                hits = [g for g in targets_in_pathway(pw_name) if g.upper() in organ_upper]
                if hits:
                    pw_hits.append((pw_name, hits))
            pw_hits.sort(key=lambda x: -len(x[1]))
            if not pw_hits:
                st.caption("本方 12 条通路中没有靶点在该器官表达。")
            else:
                organ_options = [f"{n}（{len(h)} 个靶点：{'、'.join(h[:5])}"
                                 f"{'…' if len(h) > 5 else ''}）" for n, h in pw_hits]
                organ_pick = st.selectbox("该器官涉及的通路", organ_options)
                pw_name, hits = pw_hits[organ_options.index(organ_pick)]
                st.caption(f"{pw_name} 在 {picked_organ} 表达的靶点：{'、'.join(hits)}")
                components.html(
                    render_cell_pathway(pw_name, highlight=hits[0].upper(),
                                        targets=organ_upper),
                    height=600)

        st.markdown("---")
        st.markdown("### 受体点位（细胞分区）")
        pool = organ_genes or genes
        picked_gene = st.selectbox("选择靶点", pool)
        gene_item = info.get(picked_gene.upper(), {})
        show_all = st.checkbox("显示当前靶点范围全部靶点的细胞分布", value=False)

        if show_all:
            pairs = [(g, loc) for g in genes
                     for loc in (info.get(g.upper(), {}).get("locations") or [])]
        else:
            pairs = [(picked_gene, loc)
                     for loc in (gene_item.get("locations") or ["未标注"])]

        cell_left, cell_right = st.columns([1.3, 1])
        with cell_left:
            components.html(render_cell(pairs), height=340)
        with cell_right:
            st.markdown(f"**{picked_gene}**　{gene_item.get('protein', '-')}")
            st.metric("功能类别", target_class(picked_gene, gene_item))
            st.metric("UniProt", gene_item.get("accession", "-"))
            locs = gene_item.get("locations") or []
            if locs:
                st.markdown("**亚细胞定位**")
                for loc in locs:
                    st.markdown(f"- `{loc}` → **{compartment_of(loc)}**")
            else:
                st.caption("UniProt 未提供亚细胞定位注释。")

        st.markdown("---")
        st.markdown("### 作用机制通路图")
        pw_names = pathways_of_target(picked_gene)
        options = list(PATHWAYS)
        default_index = options.index(pw_names[0]) if pw_names else 0
        picked_pw = st.selectbox("选择通路", options, index=default_index)
        if pw_names:
            st.caption(f"{picked_gene} 直接参与的通路：{'、'.join(pw_names)}；"
                       f"下方展示的是所选通路拓扑。")
        st.write(PATHWAYS[picked_pw]["desc"])

        view_pw = st.radio("通路视图", ["细胞水平", "网络拓扑"], horizontal=True)
        gene_set = {g.upper() for g in genes}
        if view_pw == "细胞水平":
            st.caption("信号流向：细胞外配体 → 膜受体 → 胞质激酶级联 → 入核转录 → 效应；"
                       "红框为当前选中靶点，★ 为本方成分作用靶点。")
            components.html(
                render_cell_pathway(picked_pw, highlight=picked_gene.upper(),
                                    targets=gene_set),
                height=600)
        else:
            layout_pw = st.selectbox("通路布局",
                                     ["breadthfirst", "dagre", "cose", "grid"], index=0)
            node_styles_pw = [NodeStyle(kind, color, "name") for kind, (_, color)
                              in KIND_STYLE.items()]
            st_link_analysis(pathway_elements(picked_pw), layout_pw, node_styles_pw,
                             [EdgeStyle("PATH_STEP", "#8a8a8a", "name", directed=True)],
                             height=430, key="pathway_network")
            st.markdown("　".join(f"<span style='color:{color}'>■ {label}</span>"
                                  for label, color in KIND_STYLE.values()),
                        unsafe_allow_html=True)

        # ------------------------------------------------------ 教学视频
        st.markdown("---")
        st.markdown("### 教学视频（成分 → 靶点 → 通路 三层联动动画）")
        st.caption("默认形态：本方成分（顶部紫色条带）→ 细胞内靶点（膜受体 / 胞质激酶 / "
                   "核内因子）→ 信号级联推进至核内效应，★ 标注本方成分作用靶点；"
                   "由 `make_video.py` 渲染。")
        video_dir = HERE / "videos"
        videos = ({v.stem: v for v in video_dir.glob("*.mp4")}
                  if video_dir.exists() else {})
        if videos:
            pick_video = st.selectbox("选择视频", sorted(videos))
            st.video(str(videos[pick_video]))
            st.download_button("下载该视频", data=videos[pick_video].read_bytes(),
                               file_name=videos[pick_video].name, mime="video/mp4")
        else:
            st.info("还没有生成视频，点击下方按钮生成后即可在此播放。")

        b1, b2, b3 = st.columns(3)
        if b1.button(f"生成「{picked_pw}」三层联动动画"):
            with st.spinner("正在渲染三层联动动画（约 3 分钟）…"):
                subprocess.run([sys.executable, "make_video.py", "--pathway", picked_pw],
                               cwd=str(HERE), capture_output=True)
            st.rerun()
        if b2.button(f"生成「{picked_pw}」细胞水平动画"):
            with st.spinner("正在渲染细胞水平动画（约 1–2 分钟）…"):
                subprocess.run([sys.executable, "make_video.py", "--pathway", picked_pw,
                                "--cell"], cwd=str(HERE), capture_output=True)
            st.rerun()
        if b3.button("重新生成全部 12 条三层动画 + 总览"):
            with st.spinner("正在渲染全部视频（约 35 分钟）…"):
                subprocess.run([sys.executable, "make_video.py", "--all", "--overview"],
                               cwd=str(HERE), capture_output=True)
            st.rerun()

# ------------------------------------------------------------- 8. 成分明细
with tab_compound:
    if compounds_df.empty:
        st.warning("当前阈值下没有成分。")
    else:
        for extra in ["formula", "mw", "smiles", "cid", "img_url"]:
            if extra not in compounds_df.columns:
                compounds_df[extra] = None
        show = compounds_df[["name", "en", "ob", "dl", "herbs", "n_targets",
                             "formula", "mw", "smiles"]].copy()
        show.columns = ["成分", "英文名", "OB(%)", "DL", "所属中药", "靶点数",
                        "分子式", "分子量", "SMILES"]
        show["所属中药"] = show["所属中药"].apply(lambda x: "、".join(x))
        show = show.sort_values(["靶点数", "OB(%)"], ascending=[False, False])
        st.dataframe(show, hide_index=True, use_container_width=True)

        st.markdown("**化学结构式**")
        pick = st.selectbox("选择成分查看结构", show["成分"].tolist())
        row = compounds_df[compounds_df["name"] == pick].iloc[0]
        left, right = st.columns([1, 2])
        if row.get("img_url"):
            left.image(row["img_url"], caption=f"{pick}（PubChem CID {row.get('cid')}）",
                       width=320)
        else:
            left.info("该成分暂无结构图，可运行 `python fetch_structures.py` 从 PubChem 补全。")
        if row.get("smiles"):
            right.markdown("**SMILES（可复制）**")
            right.code(row["smiles"])
            right.caption(f"分子式 {row.get('formula')}　分子量 {row.get('mw')}")
        st.download_button(
            "下载成分表 CSV",
            data=show.to_csv(index=False).encode("utf-8-sig"),
            file_name="compounds.csv",
            mime="text/csv",
        )
