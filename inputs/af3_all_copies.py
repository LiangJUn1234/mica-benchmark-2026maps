#!/usr/bin/env python3
import copy, json, sys
from pathlib import Path

if len(sys.argv) != 4:
    sys.exit("usage: af3_all_copies.py <mica_json> <data_json> <out_json>")
mica_json, data_json, out = map(Path, sys.argv[1:])
assert not out.exists(), out

# mica_json 是 MICA 自带的 utils/fasta_to_AF3_json.py 写出的 AlphaFold Server 格式：每份 copy 一个 proteinChain。
# 本地 AF3 3.0.1 不接受这种格式，所以改用等价输入：单 copy 那次的 data pipeline 输出（同一序列的 MSA 和 template），
# 只把 protein 的 id 从 "A" 换成 MICA 列出的全部 chain id。AF3 看到 id 是一个列表，就会预测这么多份 copy。
mica = json.loads(mica_json.read_text())
chains = [s["proteinChain"] for s in mica[0]["sequences"]]
ids = [c["id"] for c in chains]
assert len(mica) == 1 and {c["count"] for c in chains} == {1} and len({c["sequence"] for c in chains}) == 1

d = json.loads(data_json.read_text())
assert len(d["sequences"]) == 1 and d["sequences"][0]["protein"]["id"] == "A"
assert d["sequences"][0]["protein"]["sequence"] == chains[0]["sequence"]
new = copy.deepcopy(d)
new["sequences"][0]["protein"]["id"] = ids

# 确认除了 id 以外，新旧两份完全相同。
a, b = copy.deepcopy(d), copy.deepcopy(new)
a["sequences"][0]["protein"].pop("id"); b["sequences"][0]["protein"].pop("id")
assert a == b

out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(new, indent=2))
print(f"{out}: {len(ids)} copies {ids}")
