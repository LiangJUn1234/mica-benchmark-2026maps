#!/usr/bin/env python3
import json, sys, time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import requests

if len(sys.argv) != 3:
    sys.exit("usage: query_candidates.py <out_dir> <freeze_date YYYY-MM-DD>")
OUT, FREEZE = Path(sys.argv[1]), sys.argv[2]
RAW = OUT / "raw"; RAW.mkdir(parents=True, exist_ok=True)
S = requests.Session()

def get_json(url, cache, **kw):
    # 每次 API 的回应都存进 out_dir/raw/；再跑时直接读存档，所以同一份存档一定得到同一张表。
    f = RAW / cache
    if f.exists():
        return json.loads(f.read_text())
    for i in range(5):  # 失败时等 1、2、4、8、16 秒再试
        r = S.post(url, timeout=120, **kw) if "json" in kw else S.get(url, timeout=60)
        if r.status_code == 200:
            f.write_text(r.text)
            return r.json()
        if r.status_code == 404:
            return None
        time.sleep(2 ** i)
    raise RuntimeError(f"{url} -> {r.status_code}")

# 第 1 步，RCSB 搜索：cryo-EM、2026-01-01 到 freeze_date 之间发布、含 protein、resolution ≤ 6 Å。
def node(attr, op, val):
    return {"type": "terminal", "service": "text", "parameters": {"attribute": attr, "operator": op, "value": val}}
query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
    node("exptl.method", "exact_match", "ELECTRON MICROSCOPY"),
    node("rcsb_accession_info.initial_release_date", "range",
         {"from": "2026-01-01", "to": FREEZE, "include_lower": True, "include_upper": True}),
    node("rcsb_entry_info.polymer_entity_count_protein", "greater", 0),
    node("rcsb_entry_info.resolution_combined", "less_or_equal", 6.0)]},
    "return_type": "entry", "request_options": {"return_all_hits": True}}
hits = get_json("https://search.rcsb.org/rcsbsearch/v2/query", f"search_{FREEZE}.json", json=query)
ids = sorted(h["identifier"] for h in hits["result_set"])
print(f"step1 RCSB hits: {len(ids)}", flush=True)

# 第 2a 步，RCSB GraphQL 取每个 entry 的信息，每次 200 个。
GQL = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
 rcsb_accession_info{initial_release_date} rcsb_entry_info{resolution_combined}
 rcsb_entry_container_identifiers{emdb_ids} em_experiment{reconstruction_method}
 polymer_entities{rcsb_id entity_poly{type pdbx_seq_one_letter_code_can}
  rcsb_polymer_entity_container_identifiers{auth_asym_ids}}}}"""
entries = []
for i in range(0, len(ids), 200):
    d = get_json("https://data.rcsb.org/graphql", f"gql_{FREEZE}_{i}.json",
                 json={"query": GQL, "variables": {"ids": ids[i:i + 200]}})
    entries += d["data"]["entries"]

rows, fasta = [], []
for e in entries:
    prot = [p for p in e["polymer_entities"] if p["entity_poly"]["type"] == "polypeptide(L)"]
    emdb = e["rcsb_entry_container_identifiers"]["emdb_ids"] or []
    row = {"pdb": e["rcsb_id"], "pdb_release": e["rcsb_accession_info"]["initial_release_date"][:10],
           "res": (e["rcsb_entry_info"]["resolution_combined"] or [None])[0],
           "n_emdb": len(emdb), "emdb": emdb[0] if len(emdb) == 1 else ";".join(emdb),
           "recon": (e["em_experiment"] or {}).get("reconstruction_method"),
           "n_prot_entity": len(prot),
           # protein residue 数 = 每种 protein 的序列长度 × 沉积的 chain 数（份数），再加总
           "n_prot_res": sum(len(p["entity_poly"]["pdbx_seq_one_letter_code_can"])
                             * len(p["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]) for p in prot),
           "n_other_entity": len(e["polymer_entities"]) - len(prot)}
    # reason 记第一条不符合的条件；空字符串表示还在候选里
    row["reason"] = ("emdb_count!=1" if row["n_emdb"] != 1 else
                     "not_single_particle" if row["recon"] != "SINGLE PARTICLE" else "")
    rows.append(row)
    for p in prot:
        chains = ",".join(p["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"])
        fasta.append(f">{p['rcsb_id']}|Chains {chains}\n{p['entity_poly']['pdbx_seq_one_letter_code_can']}")

# 第 2b 步，对还在候选里的 entry 查 EMDB：map 发布日期、方法、contour level、有没有反向链接到这个 PDB。
def emdb_info(row):
    d = get_json(f"https://www.ebi.ac.uk/emdb/api/entry/{row['emdb'].split('-')[1]}", f"{row['emdb']}.json")
    if d is None:
        return {"reason": "emdb_missing"}
    c = d["map"].get("contour_list", {}).get("contour", [{}])
    pdbs = [x.get("pdb_id", "").upper() for x in d.get("crossreferences", {}).get("pdb_list", {}).get("pdb_reference", [])]
    return {"map_release": d["admin"]["key_dates"]["map_release"][:10],
            "emdb_method": d["structure_determination_list"]["structure_determination"][0].get("method"),
            "contour0": c[0].get("level"),  # MICA 的 emdb_extractor.py 用第一个 contour；没有的话 MICA 无法 dock
            "contour_primary": next((x.get("level") for x in c if x.get("primary")), None),
            "pdb_in_emdb": row["pdb"] in pdbs}

live = [r for r in rows if not r["reason"]]
with ThreadPoolExecutor(8) as ex:  # 8 个请求同时发
    for r, info in zip(live, ex.map(emdb_info, live)):
        r.update(info)
for r in live:
    r["reason"] = r.get("reason") or (
        "map_release_not_2026" if not ("2026-01-01" <= r["map_release"] <= FREEZE) else
        "emdb_not_singleParticle" if r["emdb_method"] != "singleParticle" else
        "pdb_not_linked_in_emdb" if not r["pdb_in_emdb"] else "")
for r in rows:
    # 没有被排除的 entry 按 resolution 分层：≤ 4 Å 是 main，其余是 4–6 Å
    r["tier"] = "" if r["reason"] else ("main_le4A" if r["res"] <= 4.0 else "tier2_4to6A")

cols = ["pdb", "emdb", "pdb_release", "map_release", "res", "tier", "reason", "n_emdb", "recon", "emdb_method",
        "contour0", "contour_primary", "pdb_in_emdb", "n_prot_entity", "n_prot_res", "n_other_entity"]
with open(OUT / "candidates_2026.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
(OUT / "proteins_2026.fasta").write_text("\n".join(fasta) + "\n")

print("reasons:", dict(Counter(r["reason"] or "eligible" for r in rows)))
print("tiers:", dict(Counter(r["tier"] for r in rows if r["tier"])))
