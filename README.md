# 复方网络药理学分析平台

用 Streamlit + Cytoscape.js（`st-link-analysis`）把网络药理学的标准研究流程串成一个可交互网页应用。

## 项目结构

| 文件 | 作用 | 对应流程 |
| --- | --- | --- |
| `data_prep.py` | 成分/靶点获取、OB/DL 筛选、UniProt 基因符号标准化，输出 `network.json` | 第一步：数据准备 |
| `network.json` | Cytoscape 格式的网络数据（`HERB` / `COMPOUND` / `TARGET`） | 数据交换文件 |
| `topology.py` | CytoNCA 六指标（DC、BC、CC、EC、LAC、NC）迭代筛选 | 第二步：拓扑筛选 |
| `app.py` | Streamlit 界面：网络图、配伍协同分析、核心靶点、成分明细 | 第三步：应用开发 |
| `requirements.txt` | 依赖清单 | 部署 |

## 安装

```bash
python -m venv .venv
.venv\Scripts\activate          # macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
```

国内网络可加镜像：`pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple`

## 运行

**不要用 PyCharm 的 Run 按钮**，必须在终端执行：

```bash
streamlit run app.py
```

浏览器自动打开 `http://localhost:8501`，改完代码保存即自动刷新。

## 界面功能

- **网络图**：`st-link-analysis` 渲染“中药-成分-靶点”网络，点击节点高亮邻居、可拖拽布局，
  支持切换 `cose / fcose / breadthfirst / circle / concentric / grid / cola / dagre` 布局。
- **配伍协同分析**：选一味药看靶点覆盖率、独有靶点，与其他药做交集与 Jaccard 相似度比较。
- **核心靶点**：CytoNCA 六指标迭代筛选，输出每轮节点数/靶点数与核心靶点的拓扑指标表（可下载 CSV）。
- **成分明细**：OB / DL、所属中药、靶点数（侧边栏可实时调整 OB / DL 阈值重新构图）。

## 换成真实数据

1. 从 TCMSP 导出每味中药的成分表（含 OB、DL、靶点蛋白名）；
2. 替换 `data_prep.py` 中的 `COMPOUNDS` 与 `HERB_COMPOUNDS`
   （或改写成读取自己的 CSV，保持 `build_elements()` 的输出格式不变）；
3. 靶点标准化用 `uniprot_to_gene_symbol()`，把 `USE_UNIPROT_API` 置为 `True` 即调用
   UniProt ID Mapping API，离线时自动回退到本地映射表；
4. 重新生成数据：`python data_prep.py`。

当前内置示例数据规模：6 味中药、37 个活性成分、153 个靶点、519 条关系；
第一轮 CytoNCA 筛选得到 38 个核心靶点（含 TNF、IL6、AKT1、TP53），与示例研究量级一致。

## 部署到 Streamlit Community Cloud

1. 推送到 GitHub（需包含 `app.py`、`data_prep.py`、`topology.py`、`network.json`、`requirements.txt`）；
2. 打开 <https://share.streamlit.io> 连接仓库，入口文件选 `app.py`；
3. 部署完成后获得公开链接。
