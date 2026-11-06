import argparse, json, os, re
from collections import Counter
from pathlib import Path
import gemmi, numpy as np, pandas as pd
from scipy.spatial import cKDTree

PROTEIN = (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)


def copies(cif):
    """reference 里每种 protein chain（entity）有几份 copy；返回 (最多的份数, protein entity 数)。"""
    st = gemmi.read_structure(str(cif)); st.setup_entities()
    n = Counter(st.get_entity_of(ch.get_polymer()).name for ch in st[0]
                if st.get_entity_of(ch.get_polymer()) and ch.get_polymer().check_polymer_type() in PROTEIN)
    return max(n.values()), len(n)


def ca_xyz(path):
    st = gemmi.read_structure(str(path)); st.setup_entities()
    return np.array([r.find_atom("CA", "*").pos.tolist() for ch in st[0] for r in ch.get_polymer() if r.find_atom("CA", "*")])


def n_chains(path):
    """数 model 里有 protein Cα 的 chain 条数（先去掉 ligand 和水，避免把钙离子 CA 算进去）。"""
    st = gemmi.read_structure(str(path)); st.setup_entities(); st.remove_ligands_and_waters()
    return sum(1 for ch in st[0] if any(r.find_atom("CA", "*") for r in ch))


def max_copies_fasta(fa):
    """test_2025 的 FASTA 头写着 'Chains A, B[auth C], ...'，数每条记录列了几条 chain，取最大。"""
    heads = [h for h in Path(fa).read_text().split("\n") if h.startswith(">")]
    return max(len(re.findall(r"[A-Za-z0-9]+(?:\[auth [^\]]+\])?", m.group(1))) if (m := re.search(r"Chains? ([^|]+)", h)) else 1
               for h in heads)


def main():
    env = lambda var, default: os.environ.get(var) or default
    bench = env("BENCH_ROOT", os.path.join(os.getcwd(), "work"))
    scores = env("SCORES_DIR", os.path.join(bench, "scores"))
    ap = argparse.ArgumentParser(description="Numbers in the report text that the other analysis scripts do not print.")
    ap.add_argument("--input", default=os.path.join(scores, "tables", "v2_main46", "input.tsv"), help="collect_scores.py output")
    ap.add_argument("--camatch", default=os.path.join(scores, "camatch", "cells.tsv"), help="whole_target_camatch.py output")
    ap.add_argument("--manifest", default=env("MANIFEST", os.path.join(bench, "dataset", "manifest.tsv")))
    ap.add_argument("--refs", default=env("REFS_DIR", os.path.join(bench, "refs")))
    ap.add_argument("--sets", default=str(Path(__file__).resolve().parent / "map_sets.json"))
    ap.add_argument("--truth", default=os.path.join(env("CALIB_DIR", os.path.join(bench, "calib")), "truth_2025", "Test_Dataset", "2025_Test"),
                    help="test_2025 references of the MICA paper (as in calibrate_test2025.py)")
    ap.add_argument("--allcopies", default=os.path.join(bench, "outputs_allcopies"),
                    help="RESULTS_DIR of the EMD-56389 run with every copy given to AF3")
    a = ap.parse_args()

    d = pd.read_csv(a.input, sep="\t"); d = d[(d.variant == "main") & (d.aligner == "mmalign")]
    man = pd.read_csv(a.manifest, sep="\t").set_index("emdb")
    main46 = json.loads(Path(a.sets).read_text())["main46"]
    res = man.res
    tier2 = sorted(e for e in d.emdb.unique() if man.stratum[e] == "tier2_4to6A")
    train = sorted(e for e in d.emdb.unique() if man.stratum[e] == "train_similar")
    tm = d.pivot(index="emdb", columns="method", values="tm")          # 行 = map，列 = 方法；没出 model 的格子在 input.tsv 里已经是 0
    ratio = d.assign(r=d.model_len / d.ref_len).pivot(index="emdb", columns="method", values="r")   # model 长度 / reference 长度

    # 0) Fig. 1 图注：三组 map 的 resolution 平均值（范围），以及 2–4 Å 那 46 张的 protein residue 数。
    for name, es in (("2-4 A", main46), ("4-6 A", tier2), ("training-similar", train)):
        print(f"{name}: {len(es)} maps, resolution mean {res[es].mean():.2f} A ({res[es].min()}-{res[es].max()})")
    n46 = man.n_prot_res[main46]
    print(f"2-4 A: {n46.min()}-{n46.max()} protein residues (mean {n46.mean():.0f})")

    # 1) 每组 map 上各方法的 mean、median、TM-score > 0.9 的张数。只列在这组每张 map 上都跑了的方法。
    for name, es in (("2-4 A (46)", main46), ("4-6 A (19)", tier2), ("training-similar (20)", train)):
        x = tm.loc[es].dropna(axis=1)
        t = pd.DataFrame(dict(mean=x.mean(), median=x.median(), n_above_0_9=(x > 0.9).sum())).sort_values("mean", ascending=False)
        print(f"\n== {name}: TM-score by method\n{t.round(3).to_string()}")
    print(f"\nMICA > AF3-docked on {(tm.loc[main46, 'mica'] > tm.loc[main46, 'af3d']).sum()} of 46 maps")
    r = ratio.loc[main46].dropna(axis=1).mean().sort_values()
    print("Model length / reference length, mean over the 46 maps:", r.round(2).to_dict())

    # 2) copy 数：取 reference 里份数最多的那种 chain。
    cp = {e: copies(Path(a.refs) / f"{man.pdb[e]}.cif") for e in main46 + tier2}
    ncp = lambda e: cp[e][0]
    one = [e for e in main46 if ncp(e) == 1]; two3 = [e for e in main46 if 2 <= ncp(e) <= 3]; multi = [e for e in main46 if ncp(e) > 1]
    lo = [e for e in main46 if ncp(e) >= 4 and res[e] < 3.3]; hi = [e for e in main46 if ncp(e) >= 4 and res[e] >= 3.3]
    diff = lambda m1, m2, es: (tm.loc[es, m1] - tm.loc[es, m2])
    print(f"\n== Copies (2-4 A): 1 copy {len(one)}, 2-3 copies {len(two3)}, >= 4 copies {len(lo) + len(hi)}"
          f" (better than 3.3 A {len(lo)}, 3.3 A or worse {len(hi)})")
    print(f"1 copy: MICA mean {tm.loc[one, 'mica'].mean():.3f}; MICA - CryoZeta {diff('mica', 'cz', one).mean():+.3f}; "
          f"MICA - EModelX(+AF) {diff('mica', 'emx', one).mean():+.3f}, MICA higher on {(diff('mica', 'emx', one) > 0).sum()} of {len(one)}")
    print(f"Multi-copy ({len(multi)}): MICA - CryoZeta {diff('mica', 'cz', multi).mean():+.3f}")
    g = diff("mica", "cz", two3).sort_values()
    print(f"2-3 copies: the three maps with the largest MICA - CryoZeta gap {g.index[:3].tolist()} carry {g.iloc[:3].sum() / g.sum():.0%} of it")
    print(f">= 4 copies, MICA - CryoZeta summed: better than 3.3 A {diff('mica', 'cz', lo).sum():+.2f}, 3.3 A or worse {diff('mica', 'cz', hi).sum():+.2f}")
    print(f">= 4 copies, better than 3.3 A: MICA > 0.9 on {(tm.loc[lo, 'mica'] > 0.9).sum()} of {len(lo)}; "
          f"AF3-docked model covers {ratio.loc[lo, 'af3d'].mean():.1%} of the reference on average")
    # 从好于 3.3 Å 到 3.3 Å 或更差，各方法平均 TM-score 降了多少（Discussion 里 MICA、ModelAngelo、E3-CryoFold、EModelX(+AF) 的比较）。
    x = tm.loc[lo + hi].dropna(axis=1)
    t = pd.DataFrame({"better than 3.3 A": x.loc[lo].mean(), "3.3 A or worse": x.loc[hi].mean()}); t["drop"] = t.iloc[:, 0] - t.iloc[:, 1]
    print(f">= 4 copies, mean TM-score by resolution:\n{t.round(3).to_string()}")
    e = "EMD-66113"
    print(f"{e}: {ncp(e)} copies, {man.n_prot_res[e]} residues (largest of the 46: {man.n_prot_res[main46].idxmax()}), "
          f"AF3 covers {ratio.loc[e, 'af3d']:.1%}, MICA {tm.loc[e, 'mica']:.3f}")
    t2hi = [e for e in tier2 if ncp(e) >= 4]; gap = diff("mica", "dm", tier2)
    print(f"4-6 A: {len(t2hi)} maps with >= 4 copies carry {gap[t2hi].sum() / gap.sum():.0%} of the MICA - DiffModeler gap; "
          f"MICA builds {ratio.loc[t2hi, 'mica'].mean():.0%} of the reference residues on them")
    tc = {Path(f).stem: max_copies_fasta(f) for f in sorted(Path(a.truth).glob("*/*.fasta"))}
    print(f"test_2025: {len(tc)} maps, >= 4 copies: {[k for k, v in tc.items() if v >= 4]}")

    # 3) 只建出一部分 copy：default Cα match 比 whole-target 高 5 分以上。
    c = pd.read_csv(a.camatch, sep="\t"); c = c[c.emdb.isin(main46)]
    part = lambda m: set(c[(c.method == m) & (c.ca_match_default - c.ca_match_whole >= 5)].emdb)
    # EMD-65039 不算：它的 reference 只沉积了组装体 4 条 chain 里的 2 条，MICA 那条没配上的 chain 建在没沉积的 copy 上，
    # 差值来自 reference 不完整，不是 MICA 少建。
    p = sorted(part("mica") - {"EMD-65039"}); p4 = [e for e in p if ncp(e) >= 4]
    print(f"\n== Partial copies (default - whole-target Ca match >= 5): MICA {len(p)} maps {p}; CryoZeta on these: {sorted(part('cz') & set(p))}")
    print(f"of them >= 4 copies: {len(p4)}, all 3.3 A or worse: {set(p4) <= set(hi)}, all one protein entity: {all(cp[e][1] == 1 for e in p4)}")

    # 4) Fig. 4c：EMD-64462 中间一圈（离环轴 < 60 Å）里，model 的 Cα 有多少离 reference 的 Cα 超过 3 Å。
    e = "EMD-64462"; R = ca_xyz(Path(a.refs) / f"{man.pdb[e]}.cif"); c0 = R.mean(0)
    w, v = np.linalg.eigh(np.cov((R - c0).T))
    ax = v[:, 2] if abs(w[0] - w[1]) < abs(w[1] - w[2]) else v[:, 0]   # 环轴 = 三个主轴里特征值和另外两个不一样的那个
    rad = lambda X: np.linalg.norm((X - c0) - np.outer((X - c0) @ ax, ax), axis=1)
    for m in ("mica", "cz"):
        M = ca_xyz(d[(d.method == m) & (d.emdb == e)].model.iloc[0]); central = rad(M) < 60
        far = cKDTree(R).query(M[central], distance_upper_bound=3.0)[0] > 3.0
        print(f"{e} {m}: central ring {far.sum()} of {central.sum()} Ca more than 3 A from the reference")

    # 5) EMD-56389：把 15 份 copy 都给 AF3 以后重跑的 MICA 和 AF3-docked。
    e = "EMD-56389"; A = Path(a.allcopies)
    full = {m: pd.read_csv(A / m / e / "metrics.tsv", sep="\t").iloc[0] for m in ("mica", "af3d")}
    one_model = lambda m: d[(d.method == m) & (d.emdb == e)].model.iloc[0]
    w56 = c[(c.method == "mica") & (c.emdb == e)].iloc[0]
    print(f"\n== {e}: {ncp(e)} copies; one copy given to AF3: MICA {tm.loc[e, 'mica']:.3f} ({n_chains(one_model('mica'))} chains, "
          f"default Ca match {w56.ca_match_default}, whole-target {w56.ca_match_whole}), EModelX(+AF) {tm.loc[e, 'emx']:.3f} "
          f"({n_chains(one_model('emx'))} chains)")
    print(f"every copy given to AF3: MICA {full['mica'].tm:.3f} ({n_chains(full['mica'].model)} chains), AF3-docked {full['af3d'].tm:.3f}, "
          f"difference {full['mica'].tm - full['af3d'].tm:.3f}")


if __name__ == "__main__":
    main()
