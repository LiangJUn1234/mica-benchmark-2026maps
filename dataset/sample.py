#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd

if len(sys.argv) != 5:
    sys.exit("usage: sample.py <candidates_2026.tsv> <hits.m8> <training_union_pdb_ids.tsv> <out_dir>")
CAND, HITS, UNION, OUT = (Path(a) for a in sys.argv[1:5])
# 11,109 = MICA 论文 test set 里最大蛋白的 residue 数；20261020 是抽样的 seed；0.25 = 25% sequence identity。
CAP, SEED, ID_CUT = 11109, 20261020, 0.25
N = {"main": 50, "tier2_4to6A": 20, "train_similar": 20}

c = pd.read_csv(CAND, sep="\t")
e = c[c.tier.notna()].copy()  # 有 tier 的才是通过前面筛选的 entry
m = pd.read_csv(HITS, sep="\t", header=None,
                names="q t fident alnlen qcov tcov evalue".split())
m = m[m.fident > ID_CUT]
# hit 表里的 id 是 entity id（例如 9ABC_1），取 "_" 前面的部分就是 PDB id
m["pdb"], m["tpdb"] = m.q.str.split("_").str[0], m.t.str.split("_").str[0]
src = pd.read_csv(UNION, sep="\t").set_index("pdb").sources
# 每个候选 entry：最高的 identity，以及命中的是哪些方法的 training set。
# 命中的 target 不在 union 名单里，说明它是某个过时 entry 的替代条目，记为 replaces_obsolete。
hits = m.groupby("pdb").apply(lambda g: (g.fident.max(), ",".join(sorted(
    {s for t in g.tpdb.unique() for s in str(src.get(t, "replaces_obsolete")).split(",")}))))
e["max_id_union"] = e.pdb.map(lambda p: hits[p][0] if p in hits else 0.0)
e["hit_training_sets"] = e.pdb.map(lambda p: hits[p][1] if p in hits else "")
e["passes_identity"] = e.max_id_union <= ID_CUT
e["within_size_cap"] = e.n_prot_res <= CAP

# 三个抽样池。train_similar = 2–4 Å 里和 training set 相近（identity 超过 25%）的 entry。
pools = {"main": (e.tier == "main_le4A") & e.passes_identity & e.within_size_cap,
         "tier2_4to6A": (e.tier == "tier2_4to6A") & e.passes_identity & e.within_size_cap,
         "train_similar": (e.tier == "main_le4A") & ~e.passes_identity & e.within_size_cap}
rng = np.random.default_rng(SEED)
e["stratum"] = ""
for name, mask in pools.items():  # 三个池共用同一个随机数生成器，顺序固定，所以结果可以重现
    ids = sorted(e.loc[mask, "pdb"])
    pick = rng.choice(ids, size=min(N[name], len(ids)), replace=False)
    e.loc[e.pdb.isin(pick), "stratum"] = name

OUT.mkdir(parents=True, exist_ok=True)
e.to_csv(OUT / "eligible_flags.tsv", sep="\t", index=False)
cols = ["stratum", "emdb", "pdb", "pdb_release", "map_release", "res", "contour0", "n_prot_entity", "n_prot_res",
        "n_other_entity", "max_id_union", "hit_training_sets"]
e[e.stratum != ""].sort_values(["stratum", "res"])[cols].to_csv(OUT / "manifest.tsv", sep="\t", index=False)

lines = [f"freeze rules: identity>{ID_CUT} vs union => flagged; n_prot_res<={CAP}; seed={SEED}"]
for tier in ["main_le4A", "tier2_4to6A"]:
    t = e[e.tier == tier]
    lines.append(f"{tier}: eligible {len(t)} | pass identity {t.passes_identity.sum()} | "
                 f"pass identity & size {(t.passes_identity & t.within_size_cap).sum()} | "
                 f"excluded by size (among identity-pass) {(t.passes_identity & ~t.within_size_cap).sum()}")
for name, mask in pools.items():
    s = e[e.stratum == name]
    lines.append(f"stratum {name}: pool {mask.sum()} -> sampled {len(s)} | res median {s.res.median():.2f} "
                 f"[{s.res.min()}-{s.res.max()}] | residues median {s.n_prot_res.median():.0f} "
                 f"[{s.n_prot_res.min()}-{s.n_prot_res.max()}]")
(OUT / "freeze_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
