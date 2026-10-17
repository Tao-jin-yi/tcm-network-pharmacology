"""界面冒烟测试：用 Streamlit 官方 AppTest 跑关键场景，检查是否有未捕获异常。

覆盖：切换各首方剂（含无成分数据的方剂）、OB / DL 拉到极值（网络为空）、
默认全方。用于在改动后快速发现"整页崩溃"类问题。

用法：
    python smoke_test.py          # 全部通过时退出码 0
"""

from __future__ import annotations

import logging
import sys

logging.disable(logging.WARNING)

from streamlit.testing.v1 import AppTest

FORMULAS = ["银翘散", "生脉散", "血府逐瘀汤", "六味地黄丸", "知柏地黄丸",
            "当归补血汤", "二陈汤", "全方（全部药材）"]
EXTREMES = [("ob", 100.0), ("ob", 0.0), ("dl", 1.0), ("dl", 0.0)]


def find_option(widget, prefix: str):
    """在下拉选项中按前缀找（现选项带覆盖率后缀，如「六味地黄丸　可分析 6/6 味」）。"""
    for option in widget.options:
        if str(option).startswith(prefix):
            return str(option)
    return None


def run_case(formula: str, slider=None) -> tuple[int, str]:
    app = AppTest.from_file("app.py", default_timeout=300)
    app.run()
    box = next((s for s in app.sidebar.selectbox if s.label == "选择复方"), None)
    if box is None:
        return 1, "未找到「选择复方」下拉框"
    option = find_option(box, formula)
    if option is None:
        return 1, f"下拉框中没有 {formula}"
    box.select(option).run()
    if slider:
        key, value = slider
        idx = 0 if key == "ob" else 1
        app.sidebar.slider[idx].set_value(value).run()
    if app.exception:
        return len(app.exception), str(app.exception[0].value)[:200]
    return 0, ""


def main() -> int:
    failures = 0
    for formula in FORMULAS:
        count, detail = run_case(formula)
        status = "OK" if count == 0 else f"失败({count})"
        print(f"{formula:<12} {status} {detail}")
        failures += count
    for key, value in EXTREMES:
        count, detail = run_case("全方（全部药材）", (key, value))
        status = "OK" if count == 0 else f"失败({count})"
        print(f"{key}={value:<6} {status} {detail}")
        failures += count
    print("\n全部通过 ✅" if failures == 0 else f"\n发现 {failures} 处异常 ❌")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
