#!/usr/bin/env python3
import csv, sys
import requests

if len(sys.argv) != 3:
    sys.exit("usage: extract_sequences.py <id_list.tsv> <out.fasta>")
ids = sorted({r["pdb"].upper() for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")})
# RCSB 的 GraphQL 查询：每个 entry 的每个 polymer entity 取 canonical 一字母序列。
GQL = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id polymer_entities{rcsb_id
 entity_poly{type pdbx_seq_one_letter_code_can}}}}"""
out, found = [], set()
for i in range(0, len(ids), 200):  # 每次请求 200 个 entry
    r = requests.post("https://data.rcsb.org/graphql", timeout=120,
                      json={"query": GQL, "variables": {"ids": ids[i:i + 200]}})
    r.raise_for_status()
    for e in r.json()["data"]["entries"]:
        found.add(e["rcsb_id"])
        # 只要 protein（polypeptide(L)），FASTA 头是 entity id，例如 >7ABC_1
        out += [f">{p['rcsb_id']}\n{p['entity_poly']['pdbx_seq_one_letter_code_can']}"
                for p in e["polymer_entities"] if p["entity_poly"]["type"] == "polypeptide(L)"]
open(sys.argv[2], "w").write("\n".join(out) + "\n")
print(f"{len(found)}/{len(ids)} entries found; {len(out)} protein entities; missing: {sorted(set(ids) - found)}")
