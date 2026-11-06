import argparse, json, os, re, subprocess, sys, tempfile
from functools import partial
from multiprocessing import Pool
from pathlib import Path
import gemmi, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I 不会把脚本所在目录加进 import 路径，手动加上
from evaluate import PHENIX, check_tools, protein

# phenix.chain_comparison 默认（test_unique_part_of_target_only=True）报的 Cα match，是「只算一份 copy」和「算整个 target」两个值里较大的那个，
# 所以只建出一部分 copy 的 model 也可能拿高分。这里对每格跑两次：一次默认，一次 =False（whole target），两个值并排存。
# 方法：AF3-docked baseline、MICA、ModelAngelo、EModelX(+AF)、CryoAtom、EMProt、E3-CryoFold、EM3DFold、CryoZeta
METHODS = ["af3d", "mica", "ma", "emx", "cat", "emp", "e3cf", "em3d", "cz"]
COLS = ["emdb", "method", "root", "ca_match_default", "seq_match_default", "ca_match_whole", "seq_match_whole",
        "pct_unique_true", "pct_unique_false", "target_res_used_default", "target_res_original", "frozen_ca_match",
        "frozen_seq_match", "ref_protein_chains", "ref_protein_entities", "ref_max_chains_per_entity"]


def cc(r, q, t, extra):
    txt = subprocess.run(["bash", "-c", f"source {PHENIX} && phenix.chain_comparison {r} {q} {extra}"],
                         capture_output=True, text=True, cwd=t).stdout.splitlines()
    g = lambda pat: next((re.search(pat, l)[1] for l in txt if re.search(pat, l)), "")
    if any("MATCH(%)" in l for l in txt):  # 默认那次会打印表格，读法和 evaluate.py 一样
        h = next(i for i, l in enumerate(txt) if "MATCH(%)" in l)
        row = next(l.split() for l in txt[h + 1:] if l.strip())
        ca, sq = float(row[7]), float(row[9])
    else:  # =False 那次不打印表格，只有一行 'All residues near target: ... Seq match (%): S  % Found: F'
        x = next((re.search(r"Seq match \(%\):\s*([\d.]+)\s+% Found:\s*([\d.]+)", l) for l in txt if "All residues near target" in l), None)
        if x is None and any("Residues far from target" in l for l in txt):  # model 没有任何 residue 靠近 target 时没有这一行，记 0, 0
            ca, sq = 0.0, 0.0
        else:
            ca, sq = float(x[2]), float(x[1])
    return ca, sq, g(r"unique_part_of_target_only=True:\s*([\d.]+)"), g(r"unique_part_of_target_only=False:\s*([\d.]+)"), \
        g(r"Used residues in target:\s*(\d+)"), g(r"Original residues in target:\s*(\d+)")


def one(ctx, cell):
    """出错时不抛异常，把错误写进这一行：子进程里抛出的 StopIteration 会让主进程的 imap_unordered 循环悄悄结束。"""
    try:
        return _one(ctx, cell)
    except BaseException as x:
        return [cell[0], cell[1], "ERROR", f"{type(x).__name__}: {x}"] + [""] * (len(COLS) - 4)


def _one(ctx, cell):
    e, m = cell
    h = e in ctx["header_normalized"]
    root = "results_hdrnorm" if h else "results"  # header 规范化过的 map 用规范化后的结果和平移后的 reference
    d = Path(ctx["hdrnorm"] if h else ctx["results"]) / m / e
    ref = Path(ctx["shifted_refs"] if h else ctx["refs"]) / f"{ctx['pdb'][e]}.cif"
    mp = next(f.removeprefix("model=") for f in (d / "DONE").read_text().split() if "/" in f)
    fr = pd.read_csv(d / "metrics.tsv", sep="\t").iloc[0]  # 存档的分数：默认那次要和它一样
    st = gemmi.read_structure(str(Path(ctx["refs"]) / f"{ctx['pdb'][e]}.cif"))
    st.setup_entities()
    ents = pd.Series([st.get_entity_of(c.get_polymer()).name for c in st[0] if st.get_entity_of(c.get_polymer()) and
                      c.get_polymer().check_polymer_type() in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD)]).value_counts()
    with tempfile.TemporaryDirectory() as t:
        r, q = Path(t) / "ref.cif", Path(t) / "model.cif"
        protein(ref, r)
        protein(mp, q)
        a1, s1, pt, pf, used, orig = cc(r, q, t, "")
        a2, s2, *_ = cc(r, q, t, "test_unique_part_of_target_only=False")
    return [e, m, root, a1, s1, a2, s2, pt, pf, used, orig, float(fr.ca_match), float(fr.seq_match),
            int(ents.sum()), len(ents), int(ents.max())]


def main():
    env = lambda var, default: os.environ.get(var) or default  # 优先用 config/paths.sh 设的路径
    bench = env("BENCH_ROOT", os.path.join(os.getcwd(), "work"))
    manifest = env("MANIFEST", os.path.join(env("DATASET_DIR", os.path.join(bench, "dataset")), "manifest.tsv"))
    ap = argparse.ArgumentParser(description="Whole-target versus default Cα match of phenix.chain_comparison.")
    ap.add_argument("out_dir"); ap.add_argument("n_proc", nargs="?", type=int, default=6)
    ap.add_argument("--sets", default=str(Path(__file__).resolve().parents[1] / "analysis/map_sets.json"))
    ap.add_argument("--manifest", default=manifest, help="the sampled entries (default $MANIFEST)")
    ap.add_argument("--results", default=env("RESULTS_DIR", os.path.join(bench, "outputs")),
                    help="results on the original maps (default $RESULTS_DIR)")
    ap.add_argument("--hdrnorm", default=env("RESULTS_HDRNORM_DIR", os.path.join(bench, "outputs_hdrnorm")),
                    help="results on the header-normalized maps (default $RESULTS_HDRNORM_DIR)")
    ap.add_argument("--refs", default=env("REFS_DIR", os.path.join(bench, "refs")), help="reference structures (default $REFS_DIR)")
    ap.add_argument("--shifted-refs", default=env("REFS_HDRNORM_DIR", os.path.join(bench, "refs_hdrnorm")),
                    help="translated references of the header-normalized maps (default $REFS_HDRNORM_DIR)")
    a = ap.parse_args()
    check_tools()
    S = json.loads(Path(a.sets).read_text())
    man = pd.read_csv(a.manifest, sep="\t").set_index("emdb")
    ctx = dict(header_normalized=list(S["header_normalized"]), pdb=man.pdb.to_dict(), results=a.results,
               hdrnorm=a.hdrnorm, refs=a.refs, shifted_refs=a.shifted_refs)
    O = Path(a.out_dir)
    O.mkdir(parents=True, exist_ok=True)
    P = O / "cells.part"
    # 每算完一格就追加一行到 cells.part；重跑时跳过已经算过的格子
    done = set(map(tuple, pd.read_csv(P, sep="\t", header=None).iloc[:, :2].values.tolist())) if P.exists() else set()
    cells = [(e, m) for e in S["main46"] for m in METHODS
             if (Path(a.hdrnorm if e in S["header_normalized"] else a.results) / m / e / "DONE").exists() and (e, m) not in done]
    print(f"{len(cells)} cells to do ({len(done)} done before)", flush=True)
    with Pool(a.n_proc) as p, open(P, "a") as f:
        for row in p.imap_unordered(partial(one, ctx), cells):
            f.write("\t".join(map(str, row)) + "\n")
            f.flush()
    t = pd.read_csv(P, sep="\t", header=None, names=COLS)
    t.to_csv(O / "cells.tsv", sep="\t", index=False)
    print(f"wrote {len(t)} rows")


if __name__ == "__main__":
    main()
