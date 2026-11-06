import argparse, json, os, sys, tempfile
from functools import partial
from multiprocessing import Pool
from pathlib import Path
import gemmi, numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scoring"))  # 让这里能 import scoring/ 里的打分函数
from evaluate import check_tools, evaluate, protein
from whole_target_camatch import cc

# CryoACE 只跑了 2–4 Å 的 46 张里、不超过它声明的大小上限（2,000 个 protein residue）的 22 张，所以单独在这 22 张上和 MICA、CryoZeta、EModelX(+AF) 比。
SEED, BOOT = 20261020, 10000  # bootstrap 的 seed 和次数，和 paired_tests.py 一样
MAX_RES = 2000  # CryoACE 声明支持的大小上限：protein residue 数
OTHERS = ("mica", "cz", "emx")  # 和 CryoACE 比较的方法：MICA、CryoZeta、EModelX(+AF)


def chains(path):
    """数一个结构里有几条 protein chain，以及 chain 断了几处（相邻 Cα 距离超过 4.2 Å 算断开）。"""
    with tempfile.TemporaryDirectory() as t:
        protein(path, Path(t) / "p.cif"); st = gemmi.read_structure(str(Path(t) / "p.cif"))
    out = []
    for ch in st[0]:
        X = np.array([[a.pos.x, a.pos.y, a.pos.z] for r in ch for a in r if a.name == "CA" and a.altloc in ("\0", "A")]).reshape(-1, 3)
        out.append(int((np.linalg.norm(np.diff(X, axis=0), axis=1) > 4.2).sum()) if len(X) > 1 else 0)
    return len(out), sum(out)


def label(model, results):
    """把 model 的路径写成相对路径，输出表里就不会出现本机的绝对路径。"""
    p = Path(model)
    for base in (Path(results).absolute().parent, Path(results).resolve().parent):
        if p.is_relative_to(base):
            return str(p.relative_to(base))
    return model


def one(ctx, cell):
    """给一张 map 上一个方法的 model 打分。出错时不中断，把错误写进 model 这一列。"""
    try:
        e, m = cell; d = Path(ctx["results"]) / m / e; ref = Path(ctx["refs"]) / f"{ctx['pdb'][e]}.cif"
        if not (d / "DONE").exists(): return dict(emdb=e, method=m, model="no model")
        mp = next(f.removeprefix("model=") for f in (d / "DONE").read_text().split() if "/" in f); r = evaluate(mp, ref, ("mmalign",))["mmalign"]
        with tempfile.TemporaryDirectory() as t:
            a, b = Path(t) / "ref.cif", Path(t) / "model.cif"; protein(ref, a); protein(mp, b)
            wca, wsq, *_ = cc(a, b, t, "test_unique_part_of_target_only=False")
        mc, mb = chains(mp); rc, _ = chains(ref)
        return dict(emdb=e, method=m, model=label(mp, ctx["results"]), tm=r["tm"], tm_by_model=r["tm_by_model"],
                    ca_match_default=r["ca_match"], seq_match_default=r["seq_match"], ca_match_whole=wca, seq_match_whole=wsq,
                    model_len=r["model_len"], ref_len=r["ref_len"], model_chains=mc, model_breaks=mb, ref_chains=rc)
    except BaseException as x:
        return dict(emdb=cell[0], method=cell[1], model=f"ERROR {type(x).__name__}: {x}")


def main():
    env = lambda var, default: os.environ.get(var) or default  # 优先用 config/paths.sh 里设的路径
    bench = env("BENCH_ROOT", os.path.join(os.getcwd(), "work"))
    scores = env("SCORES_DIR", os.path.join(bench, "scores"))
    manifest = env("MANIFEST", os.path.join(env("DATASET_DIR", os.path.join(bench, "dataset")), "manifest.tsv"))
    ap = argparse.ArgumentParser(description="CryoACE within its size range compared with MICA, CryoZeta and EModelX(+AF).")
    ap.add_argument("out_dir"); ap.add_argument("n_proc", nargs="?", type=int, default=6)
    ap.add_argument("--sets", default=str(Path(__file__).resolve().parent / "map_sets.json"))
    ap.add_argument("--manifest", default=manifest, help="the sampled entries (default $MANIFEST)")
    ap.add_argument("--results", default=env("RESULTS_DIR", os.path.join(bench, "outputs")),
                    help="method results (default $RESULTS_DIR)")
    ap.add_argument("--refs", default=env("REFS_DIR", os.path.join(bench, "refs")), help="reference structures (default $REFS_DIR)")
    ap.add_argument("--input", default=os.path.join(scores, "tables", "v2_main46", "input.tsv"),
                    help="per-map table of collect_scores.py (default $SCORES_DIR/tables/v2_main46/input.tsv)")
    ap.add_argument("--camatch", default=os.path.join(scores, "camatch", "cells.tsv"),
                    help="output of whole_target_camatch.py (default $SCORES_DIR/camatch/cells.tsv)")
    a = ap.parse_args()
    check_tools(("mmalign",))
    S = json.loads(Path(a.sets).read_text()); man = pd.read_csv(a.manifest, sep="\t").set_index("emdb")
    ENTS = sorted(e for e in S["main46"] if man.n_prot_res[e] <= MAX_RES)
    hn = sorted(set(ENTS) & set(S["header_normalized"]))
    assert not hn, f"maps with normalized headers in scope; their results are not under --results: {hn}"
    ctx, O, A = dict(results=a.results, refs=a.refs, pdb=man.pdb.to_dict()), Path(a.out_dir), Path(a.results) / "ace"

    O.mkdir(parents=True, exist_ok=True)
    scope = pd.DataFrame([dict(emdb=e, n_prot_res=int(man.n_prot_res[e]), ace_done=(A / e / "DONE").exists(), ace_dir=(A / e).is_dir()) for e in ENTS])
    scope.to_csv(O / "scope.tsv", sep="\t", index=False)
    print(f"scope {len(scope)} entries; ace DONE {int(scope.ace_done.sum())}; without DONE: {scope[~scope.ace_done].emdb.tolist()}", flush=True)
    with Pool(a.n_proc) as p:
        rows = list(p.imap_unordered(partial(one, ctx), [(e, m) for e in ENTS for m in ("ace", "mica")]))
    t = pd.DataFrame(rows).sort_values(["method", "emdb"]); t.to_csv(O / "cells_scored.tsv", sep="\t", index=False, float_format="%.5g")
    assert not t.model.str.startswith("ERROR").any(), t[t.model.str.startswith("ERROR")]

    # 检查：重新打的 MICA 分数要和存档的 metrics.tsv 完全一样。
    mi = t[(t.method == "mica") & (t.model != "no model")]
    fr = pd.DataFrame([pd.read_csv(Path(a.results) / "mica" / e / "metrics.tsv", sep="\t").iloc[0].to_dict() | {"emdb": e} for e in mi.emdb])
    g = mi.merge(fr[["emdb", "tm", "tm_by_model", "ca_match", "seq_match"]], on="emdb", suffixes=("", "_stored"))
    chk = pd.DataFrame(dict(emdb=g.emdb, tm=g.tm == g.tm_stored, tm_by_model=g.tm_by_model == g.tm_by_model_stored,
                            ca_match=g.ca_match_default == g.ca_match, seq_match=g.seq_match_default == g.seq_match))
    chk.to_csv(O / "mica_rescore_check.tsv", sep="\t", index=False)
    print("MICA rescored = stored:", {c: f"{int(chk[c].sum())}/{len(chk)}" for c in ("tm", "tm_by_model", "ca_match", "seq_match")},
          f"(cells with a MICA model {len(mi)})", flush=True)

    # 其他方法的分数直接读逐图分数表；交不出 model 的记 0 分。
    d = pd.read_csv(a.input, sep="\t"); d = d[(d.variant == "main") & (d.aligner == "mmalign") & d.emdb.isin(ENTS)].copy()
    d.loc[d.failed == 1, ["tm", "tm_by_model", "ca_match"]] = 0
    w5 = pd.read_csv(a.camatch, sep="\t")[["emdb", "method", "ca_match_whole"]]
    d = d.merge(w5, on=["emdb", "method"], how="left"); d["ca_match_whole"] = d.ca_match_whole.fillna(0.0)
    ace = t[t.method == "ace"].set_index("emdb"); rows = []
    for treat in ("zero", "drop"):
        v = ace[["tm", "ca_match_default", "ca_match_whole"]].astype(float)
        v = v.reindex(ENTS).fillna(0.0) if treat == "zero" else v.dropna()
        rows.append(dict(method="ace", treatment=treat, n=len(v), tm=v.tm.mean(), ca_match_default=v.ca_match_default.mean(),
                         ca_match_whole=v.ca_match_whole.mean()))
    for m, x in d.groupby("method"):
        rows.append(dict(method=m, treatment="frozen", n=x.emdb.nunique(), tm=x.tm.mean(), ca_match_default=x.ca_match.mean(),
                         ca_match_whole=x.ca_match_whole.mean()))
    pd.DataFrame(rows).to_csv(O / "means_22.tsv", sep="\t", index=False, float_format="%.4f")

    # CryoACE 减每个方法的 TM-score 配对差，以及高、低、相同的张数。
    W = d.pivot(index="emdb", columns="method", values="tm"); pr, wn = [], []
    for treat in ("zero", "drop"):
        at = ace.tm.astype(float).reindex(ENTS); at = at.fillna(0.0) if treat == "zero" else at
        for m in OTHERS:
            x = (at - W[m]).dropna(); diff = x.to_numpy(); bs = np.random.default_rng(SEED).choice(diff, (BOOT, len(diff))).mean(1)
            pr.append(dict(treatment=treat, ref="ace", method=m, n=len(diff), mean_diff=diff.mean(), ci_low=np.percentile(bs, 2.5),
                           ci_high=np.percentile(bs, 97.5)))
            wn.append(dict(treatment=treat, vs=m, n=len(x), ace_win=int((x > 0).sum()), ace_loss=int((x < 0).sum()), tie=int((x == 0).sum()),
                           ace_win_tol01=int((x >= 0.01).sum()), ace_loss_tol01=int((x <= -0.01).sum()), tie_tol01=int((x.abs() < 0.01).sum())))
    pd.DataFrame(pr).to_csv(O / "ace_paired_diff.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(wn).to_csv(O / "ace_wins.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False)); print(pd.DataFrame(pr).round(3).to_string(index=False))
    print(pd.DataFrame(wn).to_string(index=False))
    print(t[t.method == "ace"][["emdb", "tm", "ca_match_default", "ca_match_whole", "model_chains", "ref_chains", "model_breaks"]].to_string(index=False))


if __name__ == "__main__":
    main()
