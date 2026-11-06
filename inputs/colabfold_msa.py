#!/usr/bin/env python3
import hashlib, io, sys, tarfile, time
from pathlib import Path
import requests

HOST = "https://api.colabfold.com"
HDR = {"User-Agent": "mica-bench2026 academic benchmark (colabfold API client)"}
# 服务器结果里的两个 a3m，改成 CryoZeta（Protenix 的目录结构）要读的文件名。
NAMES = {"uniref.a3m": "uniref100_hits.a3m", "bfd.mgnify30.metaeuk30.smag30.a3m": "mmseqs_other_hits.a3m"}


def submit(seq):
    # mode "env" = colabfold_search 的默认：UniRef30 加 ColabFold 的 environmental database。
    # 服务器忙（RATELIMIT 等）时等 120 秒再交。
    while True:
        r = requests.post(f"{HOST}/ticket/msa", data={"q": f">101\n{seq}\n", "mode": "env"}, headers=HDR, timeout=60).json()
        if r.get("status") not in ("RATELIMIT", "UNKNOWN", "MAINTENANCE"):
            return r
        time.sleep(120)


if len(sys.argv) < 3:
    sys.exit("usage: colabfold_msa.py <out_root> <sequence> [<sequence> ...]")
out = Path(sys.argv[1])
for seq in sys.argv[2:]:
    d = out / hashlib.sha256(seq.encode()).hexdigest()  # 每条序列一个文件夹，名字是序列的 sha256；已有就跳过
    if d.exists():
        print("exists", d.name)
        continue
    t0 = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    r = submit(seq)
    while r.get("status") in ("PENDING", "RUNNING"):  # 每 30 秒问一次好了没有
        time.sleep(30)
        r = requests.get(f"{HOST}/ticket/{r['id']}", headers=HDR, timeout=60).json()
    if r.get("status") != "COMPLETE":
        print("FAILED", d.name, r)
        continue
    raw = requests.get(f"{HOST}/result/download/{r['id']}", headers=HDR, timeout=600).content
    # 先写进 .part 文件夹，全部写完再改名：中途中断不会留下不完整的结果
    tmp = d.with_name(d.name + ".part")
    tmp.mkdir(parents=True)
    (tmp / "colabfold_result.tar.gz").write_bytes(raw)
    with tarfile.open(fileobj=io.BytesIO(raw)) as tf:
        for src, dst in NAMES.items():
            # 服务器用 null 字符分隔多条结果，去掉以后其余内容不变
            (tmp / dst).write_text(tf.extractfile(src).read().decode().replace("\x00", ""))
    (tmp / "request.txt").write_text(f"server {HOST}\nmode env\nsubmitted {t0}\nticket {r['id']}\nsequence {seq}\n")
    tmp.rename(d)
    print("done", d.name, len(seq))
