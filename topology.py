"""第二步：网络拓扑筛选 —— 复刻 CytoNCA 插件的六指标迭代筛选。

CytoNCA 计算 6 个中心性指标：
    DC  Degree Centrality        度中心性
    BC  Betweenness Centrality   介数中心性
    CC  Closeness Centrality     紧密中心性
    EC  Eigenvector Centrality   特征向量中心性
    LAC Local Average Connectivity 局部平均连通性
    NC  Network Centrality       网络中心性（邻居度之和）

筛选规则：保留 6 个指标均高于各自中位数的节点，重复迭代直到网络收敛
（节点数不再变化或达到最大轮次）。最终留下的 TARGET 节点即核心靶点。
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import networkx as nx
import pandas as pd

METRICS = ["DC", "BC", "CC", "EC", "LAC", "NC"]


def local_average_connectivity(graph: nx.Graph, node) -> float:
    """LAC(v)：v 的邻居与 v 的共同邻居数的平均值。"""
    neighbors = set(graph[node])
    if not neighbors:
        return 0.0
    common = sum(len(set(graph[n]) & neighbors) for n in neighbors)
    return common / len(neighbors)


def network_centrality(graph: nx.Graph, node) -> float:
    """NC(v)：v 的所有邻居的度之和。"""
    return float(sum(graph.degree(n) for n in graph[node]))


def metric_table(graph: nx.Graph) -> pd.DataFrame:
    """计算当前子图的 6 个中心性指标。"""
    n_nodes = graph.number_of_nodes()
    if n_nodes == 0:
        return pd.DataFrame(columns=METRICS)

    dc = nx.degree_centrality(graph)
    bc = nx.betweenness_centrality(graph)
    cc = nx.closeness_centrality(graph)
    try:
        ec = nx.eigenvector_centrality_numpy(graph)
    except Exception:  # 不收敛时退化为 PageRank
        ec = nx.pagerank(graph)

    rows: Dict[str, List[float]] = {m: [] for m in METRICS}
    index: List[str] = []
    for node in graph.nodes():
        index.append(node)
        rows["DC"].append(dc[node])
        rows["BC"].append(bc[node])
        rows["CC"].append(cc[node])
        rows["EC"].append(ec[node])
        rows["LAC"].append(local_average_connectivity(graph, node))
        rows["NC"].append(network_centrality(graph, node))

    df = pd.DataFrame(rows, index=index)
    # NC 保持 CytoNCA 原始定义（邻居度之和，未归一化），跨轮比较时只看相对大小
    return df


def _filter_round(graph: nx.Graph, table: pd.DataFrame, strict: bool = True) -> List:
    """保留 6 个指标均高于（strict）或不低于各自中位数的节点。"""
    keep = pd.Series(True, index=table.index)
    for metric in METRICS:
        median = table[metric].median()
        keep &= table[metric] > median if strict else table[metric] >= median
    kept = [n for n in table.index if keep[n]]
    return kept


def cyto_nca_filter(
    graph: nx.Graph, max_rounds: int = 10, min_nodes: int = 1
) -> Tuple[nx.Graph, pd.DataFrame, pd.DataFrame]:
    """迭代筛选。

    返回
    ----
    core_graph : 收敛后的子图
    rounds     : 每轮的节点数 / 靶点数统计
    metrics    : 最终子图的 6 指标表
    """
    current = graph.copy()
    rounds = []

    for round_id in range(1, max_rounds + 1):
        if current.number_of_nodes() <= min_nodes:
            break
        table = metric_table(current)
        kept = _filter_round(current, table, strict=True)
        if not kept:  # 严格“高于中位数”无解时退化为“不低于中位数”
            kept = _filter_round(current, table, strict=False)
        if not kept or len(kept) == current.number_of_nodes():
            break

        n_targets = sum(1 for n in kept if str(n).startswith("TGT::"))
        rounds.append({"轮次": f"第 {round_id} 轮", "节点数": len(kept), "靶点数": n_targets})
        current = current.subgraph(kept).copy()

    rounds_df = pd.DataFrame(rounds) if rounds else pd.DataFrame(
        columns=["轮次", "节点数", "靶点数"]
    )
    return current, rounds_df, metric_table(current)


def elements_to_graph(elements: dict) -> nx.Graph:
    """把 Cytoscape 格式的 network.json 转换成 networkx 图。"""
    graph = nx.Graph()
    for node in elements.get("nodes", []):
        data = node["data"]
        graph.add_node(data["id"], **{k: v for k, v in data.items() if k != "id"})
    for edge in elements.get("edges", []):
        data = edge["data"]
        graph.add_edge(data["source"], data["target"], label=data.get("label", ""))
    return graph
