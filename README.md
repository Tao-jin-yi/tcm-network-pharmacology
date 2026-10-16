# 复方网络药理学分析平台

用 Streamlit + Cytoscape.js（`st-link-analysis`）把网络药理学的标准研究流程串成一个可交互网页应用。

## 项目结构

| 文件 | 作用 | 对应流程 |
| --- | --- | --- |
| `data_prep.py` | 成分/靶点获取、OB/DL 筛选、UniProt 基因符号标准化，输出 `network.json` | 第一步：数据准备 |
| `network.json` | Cytoscape 格式的网络数据（`HERB` / `COMPOUND` / `TARGET`） | 数据交换文件 |
| `topology.py` | CytoNCA 六指标（DC、BC、CC、EC、LAC、NC）迭代筛选 | 第二步：拓扑筛选 |
| `app.py` | Streamlit 界面：网络图、配伍协同分析、核心靶点、成分明细 | 第三步：应用开发 |
| `fetch_structures.py` | 用 PubChem 补全成分化学结构（SMILES、分子式、分子量、InChIKey、2D 图） | 结构式数据 |
| `compound_structures.json` | 上一步的产物，被 `data_prep.py` 自动合并进网络 | 结构式缓存 |
| `herb_info.py` | 中药用途：性味归经、功效主治、君臣佐使、分类 | 用途数据 |
| `import_tcmsp.py` | 导入 TCMSP 导出的 CSV/Excel（支持录入全部 499 味中药） | 批量入库 |
| `fetch_target_info.py` | 用 UniProt 补全靶点蛋白信息（蛋白名、登录号、亚细胞定位、功能） | 蛋白数据 |
| `target_info.json` | 上一步的产物，被 `app.py` 直接读取 | 蛋白数据缓存 |
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

- **网络图**：渲染“中药-成分-靶点”网络，点击节点高亮邻居、可拖拽布局，
  支持切换 `cose / fcose / breadthfirst / circle / concentric / grid / cola / dagre`。
- **配伍协同分析**：选一味药查看靶点覆盖率、独有靶点，与其他药做交集与 Jaccard 相似度比较；
  同时显示该药性味归经、功效主治与方中角色（君/臣/佐/使）。
- **核心靶点**：CytoNCA 六指标迭代筛选，输出每轮节点数/靶点数与核心靶点的拓扑指标表（可下载 CSV）。
- **成分明细**：OB / DL、所属中药、靶点数、**分子式、分子量、SMILES**，并展示 PubChem 的 2D 结构式图片。

## 化学结构式

```bash
python fetch_structures.py          # 只补缺失的成分
python fetch_structures.py --all    # 重新抓取全部
python data_prep.py                 # 把结构信息写进 network.json
```

数据来自 PubChem PUG REST（免费、无需 key），字段为 SMILES、分子式、分子量、InChIKey、CID，
结构图直接用 PubChem 的 PNG 接口渲染，无需本地安装 RDKit。
少数 TCMSP 特有成分（如 乙酰茯苓酸、泽泻萜醇 E）在 PubChem 无收录，会显示“暂无结构图”。

## 中药用途

`herb_info.py` 里按药名维护：性味、归经、功效主治、方中角色（君/臣/佐/使）、分类。
新增药材时按同样格式加一条即可，应用会自动显示。

## 核心靶点的蛋白信息与三维结构

```bash
python fetch_target_info.py     # 从 UniProt 抓取靶点蛋白信息 → target_info.json
```

数据来源：UniProt REST（人源、Swiss-Prot 已审阅）+ AlphaFold DB 预测结构。
应用内新增「靶点详情」页，展示：

- **蛋白名 / UniProt 登录号 / 物种 / 氨基酸长度**，并给出 UniProt 页面链接；
- **亚细胞定位**（细胞膜、细胞核、细胞质、线粒体、分泌蛋白等，来自 UniProt 注释）；
- **蛋白三维结构**：用 3Dmol.js 页面内渲染 AlphaFold 预测模型（PDB 内容内联，可旋转缩放），
  并提供 PDB 文件下载；
- **核心靶点的亚细胞定位分布**：柱状图统计核心靶点在各细胞区室的分布。

结构文件按需下载并缓存，不会大量占用仓库体积。

## 录入更多中药（含 TCMSP 全部 499 味）

TCMSP 旧版数据库已下线，无法直接抓取；请用下面任一方式拿到成分表后导入：

1. 从 TCMSP / HERB / SymMap 等数据库导出 CSV 或 Excel；
2. 执行导入（表头中英文均可，自动识别）：

```bash
python import_tcmsp.py 你的文件.csv            # 与内置数据合并
python import_tcmsp.py 你的文件.xlsx --replace  # 只用导入的数据
python data_prep.py
```

需要的列：中药名、成分名，可选：成分编号、英文名、OB、DL、靶点（多个用 `;` 或 `,` 分隔）、SMILES。

## 当前示例数据规模

9 味中药（黄芪、熟地黄、山茱萸、山药、茯苓、泽泻、牡丹皮、丹参、甘草）、
53 个活性成分、158 个靶点、740 条关系；53 个成分中 49 个已补全 PubChem 结构。
第一轮 CytoNCA 筛选得到 32 个核心靶点（含 TNF、IL6、AKT1、TP53）。

## 部署到 Streamlit Community Cloud

1. 推送到 GitHub（需包含 `app.py`、`data_prep.py`、`topology.py`、`herb_info.py`、
   `network.json`、`compound_structures.json`、`requirements.txt`）；
2. 打开 <https://share.streamlit.io> 连接仓库，入口文件选 `app.py`；
3. 部署完成后获得公开链接。
