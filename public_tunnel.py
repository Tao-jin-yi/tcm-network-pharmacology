"""一键把本地 Streamlit 应用变成公网可访问地址（Cloudflare Quick Tunnel）。

原理：cloudflared 在本地与 Cloudflare 边缘之间建立加密隧道，
不需要路由器端口映射，也不需要公网 IP。
免费、免注册；但网址是随机域名，**每次重启都会变**，
且电脑关机 / 睡眠后地址失效。

前置：
    winget install --id Cloudflare.cloudflared -e

用法：
    python public_tunnel.py            # 默认转发 http://localhost:8501
    python public_tunnel.py --port 8600

启动后脚本会打印公网地址，并写入 public_url.txt。
按 Ctrl+C 停止（会同时关闭隧道）。
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent
URL_FILE = HERE / "public_url.txt"
CLOUDFLARED = [
    r"C:\Program Files (x86)\cloudflared\cloudflared.exe",
    r"C:\Program Files\Cloudflare\cloudflared\cloudflared.exe",
    "cloudflared",
]


def find_cloudflared() -> str:
    for candidate in CLOUDFLARED:
        if candidate == "cloudflared":
            return candidate
        if Path(candidate).exists():
            return candidate
    return ""


def streamlit_alive(port: int) -> bool:
    try:
        import requests
        resp = requests.get(f"http://localhost:{port}/_stcore/health", timeout=5)
        return resp.status_code == 200
    except Exception:
        return False


def start_streamlit(port: int) -> subprocess.Popen:
    return subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", "app.py",
         "--server.headless=true", f"--server.port={port}"],
        cwd=str(HERE), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8501)
    args = parser.parse_args()
    port = args.port

    exe = find_cloudflared()
    if not exe:
        print("未找到 cloudflared，请先安装：")
        print("    winget install --id Cloudflare.cloudflared -e")
        return

    streamlit_proc = None
    if not streamlit_alive(port):
        print(f"正在启动 Streamlit（端口 {port}）…")
        streamlit_proc = start_streamlit(port)
        for _ in range(60):
            time.sleep(1)
            if streamlit_alive(port):
                break
        else:
            print("Streamlit 启动超时，请检查 app.py")
            return
        print("Streamlit 已就绪")
    else:
        print(f"检测到已在运行的 Streamlit（端口 {port}）")

    print("正在建立公网隧道…")
    tunnel = subprocess.Popen(
        [exe, "tunnel", "--url", f"http://localhost:{port}"],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="ignore", bufsize=1)

    url = ""
    deadline = time.time() + 60
    while time.time() < deadline and tunnel.poll() is None:
        line = tunnel.stdout.readline() if tunnel.stdout else ""
        if line:
            match = re.search(r"https://[a-zA-Z0-9\-]+\.trycloudflare\.com", line)
            if match:
                url = match.group(0)
                break
        time.sleep(0.2)

    if not url:
        print("未能获取公网地址，隧道输出如下：")
        if tunnel.stdout:
            print(tunnel.stdout.read()[:2000])
        tunnel.terminate()
        return

    URL_FILE.write_text(url, encoding="utf-8")
    print()
    print("=" * 62)
    print(f"公网访问地址：{url}")
    print("（已写入 public_url.txt；此地址每次重启都会变化）")
    print("按 Ctrl+C 停止")
    print("=" * 62)

    try:
        while tunnel.poll() is None:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\n正在关闭隧道…")
        tunnel.terminate()
        if streamlit_proc:
            streamlit_proc.terminate()
        print("已停止")


if __name__ == "__main__":
    main()
