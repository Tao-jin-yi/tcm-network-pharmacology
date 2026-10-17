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
| `formulas.py` | 复方库（14 首方剂）：组成剂量、功用主治、君臣佐使方解、配伍特点、药理与使用注意 | 复方说明页 |
| `resolve_links.py` | 用 Europe PMC 为文献补全 DOI / PMCID / 开放获取全文链接 | 文献直达 |
| `herb_info.py` | 中药用途：性味归经、功效主治、君臣佐使、分类 | 用途数据 |
| `import_tcmsp.py` | 导入 TCMSP 导出的 CSV/Excel（支持录入全部 499 味中药） | 批量入库 |
| `fetch_target_info.py` | 用 UniProt 补全靶点蛋白信息（蛋白名、登录号、亚细胞定位、功能） | 蛋白数据 |
| `target_info.json` | 上一步的产物，被 `app.py` 直接读取 | 蛋白数据缓存 |
| `fetch_literature.py` | 从 PubMed 抓取中药药理文献并抽取机制关键词与涉及靶点 | 文献机制 |
| `literature.json` | 上一步的产物（题目 / 期刊 / 年份 / 摘要 / 机制 / 靶点） | 文献数据缓存 |
| `pkpd.py` | 药动学与量-效-时模型（Hill 方程、一室一级吸收 PK、曲线拟合） | 曲线计算 |
| `animal_view.py` | 细胞受体点位图、器官-靶点映射（UniProt 组织特异性 → 器官） | 动物模型可视化 |
| `anatomy.py` | 解剖层级动物 SVG（心脏四腔/肺叶/肝分叶/豆形肾/骨骼/大血管，3 个物种） | 解剖可视化 |
| `cell_pathway.py` | 细胞水平通路图（细胞外配体 → 膜受体 → 胞质级联 → 入核 → 效应） | 细胞水平可视化 |
| `make_video.py` | 通路级联动画视频（`--cell` 生成细胞水平版本；MP4，失败回退 GIF） | 教学视频 |
| `pathways.py` | 12 条信号通路拓扑（配体→受体→激酶→转录因子→效应）与靶点归类 | 机制通路图 |
| `requirements.txt` | 依赖清单 | 部署 |
| `.streamlit/config.toml` | Streamlit 运行配置（headless、CORS、主题） | 部署 |
| `public_tunnel.py` | 一键把本机应用映射成公网地址（Cloudflare Quick Tunnel） | 公开访问 |
| `smoke_test.py` | 界面冒烟测试（AppTest 跑各方剂与极端阈值，检查整页崩溃） | 回归测试 |

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

## 公开访问

### 方式一：本机临时公网地址（免服务器、免注册，推荐先用这个）

用 Cloudflare Quick Tunnel 把本机 8501 端口映射成公网 HTTPS 地址：

```powershell
winget install --id Cloudflare.cloudflared -e   # 只需装一次
python public_tunnel.py                          # 自动拉起应用 + 隧道
```

脚本会打印形如 `https://xxxx.trycloudflare.com` 的地址（同时写入 `public_url.txt`），
任何人用浏览器打开即可访问。

- 不需要公网 IP、不需要路由器端口映射；
- **地址每次重启都会变化**，电脑关机 / 睡眠后失效；
- 按 `Ctrl+C` 停止隧道（会一并关闭隧道进程）。

### 方式二：部署到 Streamlit Community Cloud（固定网址，长期可访问）

仓库已在 GitHub（`Tao-jin-yi/tcm-network-pharmacology`），直接部署即可：

1. 打开 <https://share.streamlit.io/>，用 **GitHub 账号**登录；
2. 点 **New app** → 选仓库 `Tao-jin-yi/tcm-network-pharmacology`、
   分支 `main`、主文件 `app.py`；
3. URL 子域名填想要的短名（例如 `tcm-network`），点 **Deploy**；
4. 几十秒后得到固定地址：`https://tcm-network.streamlit.app`。

部署要点（已配好，无需改动）：

- 依赖走 `requirements.txt`（含 `matplotlib`、`imageio-ffmpeg`，视频区可正常播放 / 生成）；
- `.streamlit/config.toml`：headless、CORS、主题；
- 数据文件（`network.json` / `target_info.json` / `literature.json` /
  `compound_structures.json` / `videos/*.mp4`）均已入库，云端直接可读；
- 云端是 Linux，视频生成较慢，建议只点"生成单条"，不要点"生成全部"（会超时）。

> 注意：无论哪种方式，页面上的成分-靶点数据都来自数据库预测，
> 公开分享时请保留页面内的"预测结果、非临床结论"说明。

## 质量与稳定性

改动后建议先跑冒烟测试（用 Streamlit 官方 AppTest 真机跑一遍界面）：

```bash
python smoke_test.py
```

覆盖 8 首方剂切换（含无成分数据的方剂）+ 4 个极端阈值（OB/DL 拉到 0 与上限，
网络可能被清空），检查是否有未捕获异常导致整页崩溃。

已修过的典型问题：切换方剂后 `picked_gene` 未定义（空数据路径）、
配伍协同页缩进错误导致 `info` 未定义、空 `selectbox` options 崩溃、
`set & list` 类型错误、`use_container_width` 弃用（新版已改为 `width="stretch"`）。

已知技术债（尚未处理）：

- `st.components.v1.html` 官方公告 2026-06-01 后移除，后续需替换为 `st.iframe`
  （涉及解剖图、细胞通路图、3Dmol 结构，替换前需逐个验证渲染）；
- 14 首方剂中只有「六味地黄丸」6/6 味药材有成分-靶点数据，其余需先用
  `import_tcmsp.py` 导入数据才能做网络分析（侧栏已标注覆盖率）；
- 页面在窄屏（手机）上多列布局较拥挤，未做移动端专门适配。

## 界面功能

- **网络图**：渲染“中药-成分-靶点”网络，点击节点高亮邻居、可拖拽布局，
  支持切换 `cose / fcose / breadthfirst / circle / concentric / grid / cola / dagre`。
- **配伍协同分析**：选一味药查看靶点覆盖率、独有靶点，与其他药做交集与 Jaccard 相似度比较；
  同时显示该药性味归经、功效主治与方中角色（君/臣/佐/使）。
- **核心靶点**：CytoNCA 六指标迭代筛选，输出每轮节点数/靶点数与核心靶点的拓扑指标表（可下载 CSV）。
- **成分明细**：OB / DL、所属中药、靶点数、**分子式、分子量、SMILES**，并展示 PubChem 的 2D 结构式图片。

## 复方库与复方说明（14 首方剂）

应用内第一个标签页「复方说明」（`formulas.py`），侧栏「分析方剂」可切换方剂，
切换后网络图、配伍分析、核心靶点、通路与动物模型页都会只分析该方剂所含药材
相关的成分与靶点。

方剂库（共 14 首）：

| 方剂 | 组成（味） | 功用 |
| --- | --- | --- |
| 本项目示例方（六味地黄丸加黄芪·丹参·甘草） | 9 | 滋补肝肾、益气活血、利湿泄浊 |
| 六味地黄丸 | 6 | 滋补肝肾 |
| 知柏地黄丸 | 8 | 滋阴降火 |
| 杞菊地黄丸 | 8 | 滋肾养肝明目 |
| 金匮肾气丸（肾气丸） | 8 | 温补肾阳、化气行水 |
| 参芪地黄汤 | 8 | 益气养阴、补肾泄浊 |
| 四君子汤 | 4 | 益气健脾 |
| 补中益气汤 | 8 | 补中益气、升阳举陷 |
| 当归补血汤 | 2 | 补气生血 |
| 生脉散 | 3 | 益气生津、敛阴止汗 |
| 血府逐瘀汤 | 11 | 活血化瘀、行气止痛 |
| 逍遥散 | 8 | 疏肝解郁、养血健脾 |
| 二陈汤 | 6 | 燥湿化痰、理气和中 |
| 银翘散 | 10 | 辛凉透表、清热解毒 |

每首方剂给出：组成与剂量（g）、功用、主治、辨证要点、用法、
君臣佐使方解、配伍特点、现代药理研究方向、使用注意。

**数据覆盖率提示**：只有药材在数据库中有成分-靶点数据时才能做网络分析。
例如六味地黄丸 6/6 味有数据（可完整分析），而银翘散等多数方剂暂无成分数据，
页面会明确提示"仅展示方解"，可用 `import_tcmsp.py` 导入数据后再分析。

其中本项目示例方的完整说明如下：

- **方名 / 来源**：六味地黄丸加味（黄芪 · 丹参 · 甘草）——底方为《小儿药证直诀》
  六味地黄丸，加黄芪、丹参、甘草。
- **组成**：黄芪、熟地黄、山茱萸、山药、丹参、茯苓、泽泻、牡丹皮、甘草（9 味），
  给出每味药的常规剂量（g）、君臣佐使、性味归经、分类与功效主治，
  性味归经等字段自动取自 `herb_info.py`。
- **功用主治**：滋补肝肾、益气活血、利湿泄浊、兼清虚热；
  主治肝肾阴虚兼气虚血瘀证，并给出辨证要点与舌脉。
- **方解（君臣佐使）**：君药黄芪、熟地黄（气阴并补）；
  臣药山茱萸、山药（助滋肾填精）、丹参（活血通络）；
  佐药茯苓、泽泻、牡丹皮（「三泻」，渗湿浊、清虚热、散瘀滞）；
  使药甘草（调和诸药）。
- **配伍特点**：三补三泻补而不滞、气阴并补标本兼顾、寒温相制、
  益气活血与滋补肝肾相须。
- **现代药理研究方向**：抗炎免疫（NF-κB / TLR4-NLRP3）、抗氧化（Nrf2/HO-1）、
  抗凋亡、糖脂代谢（AMPK）、抗纤维化（TGF-β/Smad）、血管保护（HIF-1/VEGF）、
  自噬与铁死亡——均标注为文献报道方向。
- **使用注意**：熟地黄滋腻碍胃、甘草大剂量久服的水钠潴留与低钾风险、
  丹参与抗凝药的出血风险、黄芪升阳的禁忌、特殊人群与外感期间的用法。
- **说明与数据来源**：剂量与功效按《中国药典》2020 版与方剂学通行表述整理；
  平台结果属网络药理学预测，不能替代实验与临床证据。

> 用药须由医师辨证加减，本页内容用于学习与研究，不构成诊疗建议。

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

## 文献机制（PubMed）

```bash
python fetch_literature.py              # 按方中药材抓取（默认每味 20 篇）
python fetch_literature.py --herb 丹参   # 只抓一味药
python fetch_literature.py --reparse    # 不重新抓取，仅重算机制与靶点
```

用 NCBI E-utilities 检索并下载摘要，再用词典从摘要中抽取 18 类机制 / 通路
（NF-κB、PI3K/AKT、MAPK、JAK/STAT、TLR4/NLRP3、Nrf2/HO-1、AMPK/mTOR 自噬、铁死亡 …），
同时把 `TNF-α`、`IL-6`、`caspase-3` 等写法归一化为标准基因符号后与本方靶点取交集。

界面中给出：文献数、机制分布柱状图、靶点命中柱状图、文献列表与摘要，
以及**一键直达原文**按钮区。

### 文献一键直达（解决"打不开全文"）

PubMed 只给摘要页，找全文很麻烦。新增 `resolve_links.py` 用 Europe PMC REST
按 PMID 解析出 DOI、PMCID、是否开放获取与全文地址：

```bash
python resolve_links.py            # 补全缺失链接（已解析过的跳过）
python resolve_links.py --force    # 全部重查
```

界面中选中某篇文献后，直接给出按钮（有几项显示几项）：
**PubMed 摘要 / Europe PMC / DOI 原文页 / PMC 全文页 / 开放获取 PDF / 学术检索**，
并标注该文是否为开放获取（Open Access）。

当前解析结果：180 篇中 179 篇解析成功，**178 篇有 DOI，73 篇为开放获取（可直接打开全文）**；
其余为订阅制期刊，可通过"学术检索"找免费版本。

## 药动学（PK）与量-效-时（PD）曲线

模型见 `pkpd.py`：

- 量-效：`E(D) = E0 + Emax·D^h / (EC50^h + D^h)`（Hill 方程）
- 药-时：一级吸收一室模型 `C(t) = F·D·ka/(V(ka-ke))·(e^(-ke·t) - e^(-ka·t))`
- 时-效：把 `C(t)` 代入 Hill 方程（PK 驱动）

两种用法：

1. **参数模拟**：调 D / F / V / ka / ke 与 E0 / Emax / EC50 / h，实时得到
   血药浓度-时间曲线、剂量-效应曲线、效应-时间曲线、多剂量 PK 比较，
   并输出 Cmax、Tmax、AUC₀₋t、消除半衰期。
2. **上传实验数据拟合**：CSV 上传后选自变量 / 因变量列，用 `scipy.curve_fit`
   反推参数并给出 R² 与散点+拟合线图。

> 模拟曲线基于输入参数，属推演结果；拟合曲线才代表真实实验数据，请把实测 CSV 上传使用。

## 实验动物模型可视化

在「药动学 / 量-效-时」页底部，可编辑：动物品系、造模方法、给药途径、终点指标，
以及分组方案（组别 / 动物数 / 造模 / 剂量 / 给药），界面给出：

- 模型要素卡片（品系 → 造模 → 给药 → 终点指标）
- 实验时间轴甘特图（适应期 → 造模期 → 给药干预 → 取材检测）
- 分组剂量柱状图，分组方案可导出 CSV

默认内容是可编辑模板，请按真实实验方案修改。

## 动物模型模拟实验（器官 / 受体 / 通路）

应用内「动物模型模拟实验」页分三层可视化：

1. **动物器官点位图**：小鼠（C57BL/6）、大鼠（SD）、兔（新西兰白兔）三套 SVG 轮廓，
   器官圆点大小与颜色深浅 = 该器官相关靶点数量，选中器官后右侧列出靶点、
   蛋白名、功能类别、亚细胞定位与组织特异性原文。
   - 器官-靶点映射依据 **UniProt 组织特异性注释**（17 个器官类别，164 条记录中 117 条有注释），
     未标注的蛋白不会出现在图里。
2. **受体点位（细胞分区图）**：把靶点按 UniProt 亚细胞定位放到
   细胞外（分泌/配体）、细胞膜（受体）、细胞质（激酶/信号）、细胞核（转录因子）、
   线粒体、内质网/高尔基体 六个区室，可切换显示单个靶点或整个靶点集合。
   - 靶点类别（膜受体 / 核内转录因子 / 蛋白激酶 / 分泌型细胞因子 / 蛋白酶 …）
     由蛋白名与定位自动判定。
3. **作用机制通路图**：`pathways.py` 里整理 12 条通路
   （NF-κB、PI3K/AKT、MAPK、JAK/STAT、TLR4-NLRP3、Nrf2/HO-1、AMPK/mTOR 自噬、
   TGF-β/Smad、Wnt/β-catenin、HIF-1/VEGF、线粒体凋亡、铁死亡），
   用 st-link-analysis 渲染「配体 → 受体 → 激酶 → 转录因子 → 效应」拓扑，
   选中靶点会自动跳到它参与的通路。

> 器官位置是示意性布局，用于表达"靶点-器官"关系，不是解剖学精确作图；
> 通路图是人工整理的通用通路骨架，用于说明机制走向。

## 解剖层级动物图与通路教学视频

### 解剖层级动物图

「动物模型模拟实验」页的动物图已从圆点示意升级为解剖层级 SVG
（`anatomy.py`）：体壁轮廓、头部、四肢、尾按物种形态绘制，内脏按真实形态
（心脏四腔与主动脉弓、左右肺叶与支气管、肝分叶、胃、脾、胰腺、小肠盘曲、
盲肠（兔特别发达）、结肠框、豆形肾含肾门与输尿管、膀胱、生殖腺、胸腺、淋巴结），
并叠加骨骼层（脊柱、肋骨、骨盆、股骨/胫骨与关节）与主动脉/腔静脉。
器官按靶点数着色，选中器官高亮并引出标注线。

> 依据解剖学比例绘制的矢量示意图，用于表达器官层级关系，不等同于影像学图像。

### 细胞水平通路视图与动画

「动物模型模拟实验」页的"作用机制通路图"支持两种视图：

- **细胞水平**（默认，`cell_pathway.py`）：把通路节点放进细胞结构——
  细胞外（配体/细胞因子）→ 细胞膜（受体，跨膜蛋白）→ 细胞质（接头/激酶级联，
  含线粒体）→ 细胞核（核膜内侧转录因子 → 核内基因转录 → 生物学效应），
  箭头表现跨膜与入核过程；红框高亮当前选中靶点，★ 标注本方成分作用靶点。
- **网络拓扑**：原来的 st-link-analysis 网络布局（breadthfirst/dagre/cose/grid）。

细胞水平动画（与静态图同一套坐标）：

```bash
python make_video.py --pathway "NF-κB 通路" --cell   # 单条细胞水平动画
python make_video.py --all --cell                    # 全部 12 条细胞水平动画
python make_video.py --all --overview --cell         # 全部 + 总览
```

动画中有一个"信号小点"沿箭头移动，直观演示信号跨膜、胞质级联与入核转导。

### 成分 → 靶点 → 通路 三层联动动画（通路教学视频的默认形态）

```bash
python make_video.py --pathway "NF-κB 通路"   # 单条（默认就是三层联动）
python make_video.py --all                    # 全部 12 条
python make_video.py --all --cell             # 细胞水平版本（不含成分层）
python make_video.py --all --plain            # 纯级联版本（无细胞结构）
```

- 顶部紫色条带为本方成分（从 network.json 的 成分→靶点 边反查，按命中该通路
  靶点数排序取前 6 个，标注所属中药）；
- 紫色虚线箭头指向各成分在细胞内作用的靶点（膜受体 / 胞质激酶 / 核内因子），
  直观呈现"成分如何进入细胞、作用于哪一层"；
- 随后信号沿级联逐级推进至核内效应，★ 标注本方成分作用靶点。

### 器官 → 细胞 → 通路 两级联动

「动物模型模拟实验」页选中器官后，会自动列出该器官表达的本方靶点参与的
信号通路（按命中靶点数排序），并直接给出该通路的细胞水平视图，
★ 只标注在该器官表达的靶点——把"药物作用于哪个器官 → 器官内哪些细胞
分子被调控 → 触发什么通路"串成一条线。

### 视频输出与依赖

- 输出 MP4 到 `videos/`（用 imageio-ffmpeg 自带的 ffmpeg 编码 H.264；
  编码失败自动回退 GIF），另有 `00_通路总览.mp4` 逐条列出 12 条通路及其关键节点。
- 应用内「动物模型模拟实验」页底部可直接播放 / 下载，并可一键重新生成。
- 依赖：`matplotlib`、`imageio-ffmpeg`（已写入 requirements.txt）。
- 三种形态对照：

  | 参数 | 画面内容 |
  | --- | --- |
  | 默认 / `--triad` | 成分条带 + 细胞结构 + 级联（三层联动） |
  | `--cell` | 细胞结构 + 级联（不含成分层） |
  | `--plain` | 纯级联拓扑（无细胞结构） |

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
