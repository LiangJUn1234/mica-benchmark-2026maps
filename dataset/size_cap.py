#!/usr/bin/env python3
import csv, hashlib, json, re, sys, time, urllib.request
from pathlib import Path

if len(sys.argv) not in (4, 5):
    sys.exit("usage: size_cap.py <eligible_flags.tsv> <manifest.tsv> <out_dir> [cache_dir]")
EL, MAN, OUT = (Path(a) for a in sys.argv[1:4])
RAW = Path(sys.argv[4]) if len(sys.argv) == 5 else OUT / "raw"
# 上限 11,109 就是 7ZQB 组装体 1 的 modeled residue 数（MICA 论文 test set 里最大的蛋白），这里一起报出来对照。
CAP, REF = 11109, "7ZQB"
OUT.mkdir(parents=True, exist_ok=True); RAW.mkdir(parents=True, exist_ok=True)
for f in ("asm1_all.tsv", "sampled90.tsv", "summary.txt", "SHA256SUMS"):
    assert not (OUT / f).exists(), f"refusing to overwrite {OUT / f}"
el = list(csv.DictReader(open(EL), delimiter="\t"))
man = {r["pdb"].upper(): r for r in csv.DictReader(open(MAN), delimiter="\t")}
Q = """{entries(entry_ids:%s){rcsb_id
 polymer_entities{entity_poly{type pdbx_seq_one_letter_code_can} rcsb_polymer_entity_container_identifiers{asym_ids auth_asym_ids}}
 assemblies{rcsb_assembly_info{assembly_id polymer_entity_instance_count_protein polymer_monomer_count modeled_polymer_monomer_count}
  pdbx_struct_assembly_gen{asym_id_list oper_expression}}}}"""

def gql(ids, cache):
    f = RAW / cache
    if f.exists(): return json.loads(f.read_text())
    for k in range(5):
        try:
            req = urllib.request.Request("https://data.rcsb.org/graphql", data=json.dumps({"query": Q % json.dumps(ids)}).encode(),
                                         headers={"Content-Type": "application/json"})
            d = json.load(urllib.request.urlopen(req, timeout=180))
            if "errors" not in d: f.write_text(json.dumps(d)); return d
            print("errors", cache, d["errors"][:2])
        except Exception as x:
            print("retry", cache, k, x)
        time.sleep(5 * (k + 1))
    raise RuntimeError(cache)

# 组装体 1 里，一条 chain 被复制几次 = oper_expression 里 operator 的个数。
# 例如 '(1-60)' 是 60 次；'(1-5)(6-10)' 是两组相乘，5 × 5 = 25 次。
def n_ops(expr):
    n = 1
    for g in re.findall(r"\(([^)]*)\)", expr) or [expr]:
        n *= sum((int(m[2]) - int(m[1]) + 1) if (m := re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", p)) else 1 for p in g.split(","))
    return n

# 和 query_candidates.py 的算法相同（序列长度 × 份数），只是份数换成组装体 1 里的份数，而不是沉积的 chain 数。
def count(e):
    a1 = next((a for a in e["assemblies"] or [] if a["rcsb_assembly_info"]["assembly_id"] == "1"), None)
    inst = {}
    for g in (a1 or {}).get("pdbx_struct_assembly_gen") or []:
        for s in g["asym_id_list"]: inst[s] = inst.get(s, 0) + n_ops(g["oper_expression"])
    prot = [p for p in e["polymer_entities"] or [] if p["entity_poly"]["type"] == "polypeptide(L)"]
    ids = lambda p: p["rcsb_polymer_entity_container_identifiers"]
    L = lambda p: len(p["entity_poly"]["pdbx_seq_one_letter_code_can"])
    info = (a1 or {}).get("rcsb_assembly_info") or {}
    return dict(dep_recomputed=sum(L(p) * len(ids(p)["auth_asym_ids"]) for p in prot),
                asm1_prot_res=sum(L(p) * sum(inst.get(s, 0) for s in ids(p)["asym_ids"]) for p in prot) if a1 else "",
                asm1_prot_copies=sum(inst.get(s, 0) for p in prot for s in ids(p)["asym_ids"]) if a1 else "",
                rcsb_asm1_prot_instances=info.get("polymer_entity_instance_count_protein", ""),
                rcsb_asm1_polymer_monomers=info.get("polymer_monomer_count", ""),
                rcsb_asm1_modeled_polymer_monomers=info.get("modeled_polymer_monomer_count", ""))

ids = sorted({r["pdb"].upper() for r in el} | {REF}); E = {}
for i in range(0, len(ids), 100):
    for e in gql(ids[i:i + 100], f"gql_{i // 100:03d}.json")["data"]["entries"] or []:
        E[e["rcsb_id"].upper()] = e
ref = count(E[REF]); CNT = list(ref)
# dep_match、copies_match 是自查：沉积的 residue 数要和 sample.py 的 n_prot_res 一样，组装体的份数要和 RCSB 自己给的一样。
cols = ["pdb", "emdb", "tier", "passes_identity", "stratum", "n_prot_res", "dep_recomputed", "dep_match", "asm1_prot_res",
        "asm1_prot_copies", "rcsb_asm1_prot_instances", "copies_match", "within_cap_dep", "within_cap_asm1", "n_other_entity",
        "rcsb_asm1_polymer_monomers", "rcsb_asm1_modeled_polymer_monomers"]
rows = []
for r in el:
    p = r["pdb"].upper(); c = count(E[p]) if p in E else dict.fromkeys(CNT, "")
    o = {k: r.get(k, "") for k in ("pdb", "emdb", "tier", "passes_identity", "n_prot_res", "n_other_entity")}
    o.update(c, stratum=man[p]["stratum"] if p in man else "", missing=int(p not in E))
    o["dep_match"] = int(str(c["dep_recomputed"]) == str(int(float(r["n_prot_res"])))) if p in E else ""
    o["copies_match"] = int(c["asm1_prot_copies"] == c["rcsb_asm1_prot_instances"]) if p in E and c["asm1_prot_copies"] != "" else ""
    o["within_cap_dep"] = int(float(r["n_prot_res"]) <= CAP)
    o["within_cap_asm1"] = int(c["asm1_prot_res"] <= CAP) if c.get("asm1_prot_res", "") != "" else ""
    rows.append(o)
w = lambda path, rs: (path.write_text("\t".join(cols) + "\n" + "".join("\t".join(str(x.get(k, "")) for k in cols) + "\n" for x in rs)))
w(OUT / "asm1_all.tsv", rows)
s90 = sorted((x for x in rows if x["stratum"]), key=lambda x: (x["stratum"], x["emdb"])); w(OUT / "sampled90.tsv", s90)

pool = lambda tier, ident: [x for x in rows if x["tier"] == tier and x["passes_identity"] == ident]
L = [f"cap {CAP}; counting = FASTA length (pdbx_seq_one_letter_code_can) of polypeptide(L) entities x copies; deposited = auth chains, asm1 = assembly 1 operators",
     f"reference {REF}: deposited {ref['dep_recomputed']}, asm1 {ref['asm1_prot_res']} (copies {ref['asm1_prot_copies']}), RCSB asm1 polymer_monomer_count {ref['rcsb_asm1_polymer_monomers']}, modeled_polymer_monomer_count {ref['rcsb_asm1_modeled_polymer_monomers']}",
     f"eligible entries {len(rows)}; missing from RCSB now {sum(x['missing'] for x in rows)}; deposited count != frozen n_prot_res {sum(x['dep_match'] == 0 for x in rows)}; "
     f"asm1 copies != RCSB protein instances {sum(x['copies_match'] == 0 for x in rows)}; no assembly 1 {sum(x['asm1_prot_res'] == '' for x in rows if not x['missing'])}"]
for name, tier, ident in (("main", "main_le4A", "True"), ("train_similar", "main_le4A", "False"), ("tier2_4to6A", "tier2_4to6A", "True")):
    P = pool(tier, ident); d = [x for x in P if x["within_cap_dep"] == 1]; a = [x for x in P if x["within_cap_asm1"] == 1]
    L.append(f"pool {name}: within cap by deposited {len(d)}, by assembly 1 {len(a)}; in deposited pool but over by asm1 {sum(x['within_cap_asm1'] == 0 for x in d)}; "
             f"over by deposited but within by asm1 {sum(x['within_cap_dep'] == 0 and x['within_cap_asm1'] == 1 for x in P)}")
L.append("sampled entries over the cap by assembly 1: " + ", ".join(f"{x['emdb']} ({x['pdb']}, {x['stratum']}) {x['asm1_prot_res']}" for x in s90 if x["within_cap_asm1"] == 0))
L.append("eligible entries missing from RCSB now (assembly 1 unknown, not counted as within): " + (", ".join(f"{x['pdb']} ({x['tier']}, identity {x['passes_identity']}, deposited {x['n_prot_res']})" for x in rows if x["missing"]) or "none"))
L.append("asm1 copies != RCSB protein instances (other protein types not counted, as in query_candidates.py): " + (", ".join(f"{x['pdb']} {x['asm1_prot_copies']}/{x['rcsb_asm1_prot_instances']} asm1 {x['asm1_prot_res']}" for x in rows if x["copies_match"] == 0) or "none"))
L.append("sampled entries with deposited != frozen n_prot_res or copies mismatch: " + (", ".join(x["emdb"] for x in s90 if x["dep_match"] == 0 or x["copies_match"] == 0) or "none"))
L.append("largest asm1 counts among sampled within cap: " + ", ".join(f"{x['emdb']} {x['asm1_prot_res']}" for x in sorted((x for x in s90 if x["within_cap_asm1"] == 1), key=lambda x: -x["asm1_prot_res"])[:5]))
(OUT / "summary.txt").write_text("\n".join(L) + "\n"); print("\n".join(L))
sha = lambda f: hashlib.sha256((OUT / f).read_bytes()).hexdigest()
(OUT / "SHA256SUMS").write_text("".join(f"{sha(f)}  {f}\n" for f in ("asm1_all.tsv", "sampled90.tsv", "summary.txt")))
