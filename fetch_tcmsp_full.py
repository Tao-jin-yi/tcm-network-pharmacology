"""抓取 TCMSP 全库 502 味药的「成分-靶点-疾病」数据（四层数据源）。

数据结构（每味药）：
    ingredients : 成分表（MOL_ID、名称、OB、MW、AlogP、Caco-2、BBB、Hdon 等完整 ADME）
    targets     : 成分-靶点表（MOL_ID -> 靶点名/ID/DrugBank ID/SVM 预测分/RF 预测分）
    diseases    : 靶点-疾病表（靶点 -> 疾病名/疾病 ID）

特性：
    - 增量保存 tcmsp_full.json，中断后重跑自动续传；
    - token 失效自动刷新重试；
    - 限速 1.2 秒/请求，总耗时约 25-30 分钟。

用法：
    python fetch_tcmsp_full.py            # 抓全部 502 味
    python fetch_tcmsp_full.py --limit 3  # 只抓前 3 味（测试）
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import lxml.html
import requests

HERE = Path(__file__).parent
CATALOG_FILE = HERE / "tcmsp_herb_catalog.json"
OUTPUT_FILE = HERE / "tcmsp_full.json"

BASE = "https://old.tcmsp-e.com"
ROOT = f"{BASE}/tcmspsearch.php"
DELAY = 1.2          # 每次请求间隔（秒），避免触发封禁
MAX_RETRY = 3

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                  "(KHTML, like Gecko) Chrome/120.0 Safari/537.36",
    "Accept-Encoding": "gzip, deflate",
}


class TcmspClient:
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update(HEADERS)
        self.token = ""

    def get(self, url: str) -> str:
        r = self.session.get(url, timeout=60)
        r.raise_for_status()
        return r.content.decode("utf-8", errors="replace")

    def refresh_token(self) -> None:
        self.get(f"{BASE}/tcmsp.php")
        html = self.get(ROOT)
        tokens = lxml.html.fromstring(html).xpath(
            '//form[@id="SearchForm"]//input[@name="token"]/@value')
        if not tokens:
            raise RuntimeError("无法获取 TCMSP token（页面结构可能已变更）")
        self.token = tokens[0]

    def query(self, params: dict) -> str:
        for attempt in range(MAX_RETRY):
            if not self.token:
                self.refresh_token()
            try:
                return self.get(f"{ROOT}?{'&'.join(f'{k}={v}' for k, v in params.items())}")
            except Exception:
                self.token = ""   # 失效则刷新重试
                time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"查询失败：{params}")


def extract_json_array(text: str, start: int) -> str:
    """从 start（指向 '['）开始做括号配平扫描，返回完整 JSON 数组文本。

    名称里可能含 ']'（如化合物名），正则截断会得到残缺 JSON，
    故用逐字符扫描：跳过字符串字面量、按深度配平括号。
    """
    depth = 0
    in_str = False
    esc = False
    for i in range(start, len(text)):
        ch = text[i]
        if in_str:
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
            continue
        if ch == '"':
            in_str = True
        elif ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
            if depth == 0:
                return text[start:i + 1]
    return ""


def extract_grids(html: str) -> dict[str, list]:
    """从详情页 HTML 中提取 grid（成分）/grid2（成分-靶点）/grid3（靶点-疾病）。"""
    result = {"ingredients": [], "targets": [], "diseases": []}
    for pattern, key in (("grid", "ingredients"), ("grid2", "targets"),
                         ("grid3", "diseases")):
        decl = html.find(f'$("#{pattern}").kendoGrid')
        if decl < 0:
            continue
        dm = re.search(r"data:\s*\[", html[decl:decl + 2000])
        if not dm:
            continue
        arr_start = decl + dm.end() - 1
        arr_text = extract_json_array(html, arr_start)
        if arr_text:
            try:
                result[key] = json.loads(arr_text)
            except json.JSONDecodeError:
                pass
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0, help="只抓前 N 味（测试用）")
    parser.add_argument("--retry-zeros", action="store_true",
                        help="补抓此前存为空结果的药材（偶发首请求页面不完整）")
    args = parser.parse_args()

    catalog = json.loads(CATALOG_FILE.read_text(encoding="utf-8"))
    if args.limit:
        catalog = catalog[:args.limit]

    store: dict = {}
    if OUTPUT_FILE.exists():
        store = json.loads(OUTPUT_FILE.read_text(encoding="utf-8"))
        print(f"已有 {len(store)} 味，续传")

    if args.retry_zeros:
        zero = [cn for cn, h in store.items()
                if not (h.get("ingredients") or h.get("targets") or h.get("diseases"))]
        for cn in zero:
            del store[cn]
        if zero:
            print(f"补抓 {len(zero)} 味空结果：{zero[:10]}{'...' if len(zero) > 10 else ''}")

    client = TcmspClient()
    t0 = time.time()
    done = 0
    for i, herb in enumerate(catalog):
        cn = herb["herb_cn_name"]
        en = herb["herb_en_name"]
        if cn in store and store[cn].get("fetched"):
            continue
        try:
            html = client.query({"qr": requests.utils.quote(en),
                                 "qsr": "herb_en_name",
                                 "token": client.token})
            data = extract_grids(html)
            # 偶发：token 签发后首个详情请求返回不完整页面，空结果自动重试一次
            if not (data["ingredients"] or data["targets"] or data["diseases"]):
                time.sleep(2)
                client.token = ""
                html = client.query({"qr": requests.utils.quote(en),
                                     "qsr": "herb_en_name",
                                     "token": client.token})
                data = extract_grids(html)
        except Exception as e:
            print(f"[{i + 1}/{len(catalog)}] {cn}: 失败 {e}", flush=True)
            time.sleep(3)
            client.token = ""
            continue
        store[cn] = {
            "cn": cn,
            "pinyin": herb.get("herb_pinyin", ""),
            "en": en,
            "category": herb.get("child_cn_name", ""),
            "ingredients": data["ingredients"],
            "targets": data["targets"],
            "diseases": data["diseases"],
            "fetched": True,
        }
        done += 1
        # 每 10 味写盘一次（断点续传）
        if done % 10 == 0 or i == len(catalog) - 1:
            OUTPUT_FILE.write_text(json.dumps(store, ensure_ascii=False),
                                   encoding="utf-8")
        n_ing = len(data["ingredients"])
        n_tgt = len(data["targets"])
        n_dis = len(data["diseases"])
        print(f"[{i + 1}/{len(catalog)}] {cn}: 成分{n_ing} 靶点对{n_tgt} 疾病对{n_dis}"
              f"  ({time.time() - t0:.0f}s)", flush=True)
        time.sleep(DELAY)

    OUTPUT_FILE.write_text(json.dumps(store, ensure_ascii=False), encoding="utf-8")

    # 汇总
    n_ing = sum(len(h["ingredients"]) for h in store.values())
    n_tgt = sum(len(h["targets"]) for h in store.values())
    n_dis = sum(len(h["diseases"]) for h in store.values())
    ok = sum(1 for h in store.values() if h.get("fetched"))
    print()
    print(f"完成 {ok}/{len(catalog)} 呗；成分记录 {n_ing}，成分-靶点 {n_tgt}，靶点-疾病 {n_dis}")
    print(f"已保存 {OUTPUT_FILE}")
    print("下一步：python data_prep.py 重建 network.json，然后 streamlit run app.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
