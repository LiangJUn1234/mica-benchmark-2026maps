import argparse, os, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

METRICS = ["tm", "aligned_len", "seq_id", "ca_match", "seq_match", "ca_quality"]  # MICA 论文的六个指标
KEY = ["variant", "aligner"]
TIER2 = "tier2_4to6A"
emd = lambda s: "EMD-" + re.sub(r"\D", "", str(s))  # 统一写成 EMD-<编号>


def load(path):
    df = pd.read_csv(path, sep="\t")
    df = df.assign(emdb=df.emdb.map(emd), variant=df.get("variant", "main"), aligner=df.get("aligner", "mmalign"))
    df["failed"] = df.failed.fillna(0).astype(int)
    df.loc[df.failed == 1, METRICS] = 0  # 没交出 model 的格子每项记 0
    return df


def main():
    env = lambda var, default: os.environ.get(var) or default  # 优先用 config/paths.sh 设的路径
    bench = env("BENCH_ROOT", os.path.join(os.getcwd(), "work"))
    manifest = env("MANIFEST", os.path.join(env("DATASET_DIR", os.path.join(bench, "dataset")), "manifest.tsv"))
    ap = argparse.ArgumentParser(description="Paired tests of a reference method against every other method.")
    ap.add_argument("--input", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--manifest", default=manifest, help="columns emdb, stratum (default $MANIFEST)")
    ap.add_argument("--ref", default="mica"); ap.add_argument("--main", default="main")
    ap.add_argument("--boot", type=int, default=10000); ap.add_argument("--seed", type=int, default=20261020)
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        ap.error(f"{out} exists")

    df = load(a.input)
    man = pd.read_csv(a.manifest, sep="\t").assign(emdb=lambda m: m.emdb.map(emd))[["emdb", "stratum"]]
    outside = sorted(set(df.emdb) - set(man.emdb))
    df = df.merge(man, on="emdb")
    dup = df.duplicated(["method", "emdb"] + KEY)
    assert not dup.any(), f"duplicate rows: {df[dup][['method', 'emdb'] + KEY].values.tolist()}"
    methods = [a.ref] + sorted(set(df.method) - {a.ref})
    rng0 = lambda: np.random.default_rng(a.seed)  # 每做一个比较都用同一个 seed 重新开始，结果和比较的先后顺序无关

    tests = []
    # 两组 map 按 manifest 的 stratum 分：main = 2–4 Å，tier2 = 4–6 Å。training-similar 那 20 张不做检验。
    for name, strata in [("main", [a.main]), ("tier2", [TIER2])]:
        d = df[df.stratum.isin(strata)]
        if d.empty:
            continue
        for (var, alg), g in d.groupby(KEY):
            ref = g[g.method == a.ref].set_index("emdb")
            for m in methods[1:]:
                # 只用两个方法都有分数的 map，逐图配对
                pair = ref.join(g[g.method == m].set_index("emdb")[METRICS], how="inner", rsuffix="_o")
                for k in METRICS:
                    diff = (pair[k] - pair[f"{k}_o"]).to_numpy(float)
                    r = dict(set=name, variant=var, aligner=alg, ref=a.ref, method=m, metric=k, n=len(diff),
                             mean_ref=pair[k].mean(), mean_other=pair[f"{k}_o"].mean(), mean_diff=np.nan,
                             ci_low=np.nan, ci_high=np.nan, t=np.nan, p_t=np.nan, W=np.nan, p_wilcoxon=np.nan)
                    if len(diff) == 0:
                        tests.append(r)
                        continue
                    rng = rng0()
                    r["mean_diff"] = diff.mean()
                    # bootstrap：从配对差里有放回地抽 a.boot 次（默认 10,000），每次算平均，取 2.5% 和 97.5% 分位数当 95% CI
                    bs = rng.choice(diff, (a.boot, len(diff))).mean(1)
                    # 差值全部相同时 t 检验算不出来，全为 0 时 Wilcoxon 算不出来，这两种情况留空
                    if len(diff) > 1 and np.any(diff != diff[0]):
                        r["t"], r["p_t"] = stats.ttest_rel(pair[k], pair[f"{k}_o"])
                    if np.any(diff != 0):
                        r["W"], r["p_wilcoxon"] = stats.wilcoxon(diff)
                    r["ci_low"], r["ci_high"] = np.percentile(bs, [2.5, 97.5])
                    tests.append(r)

    tests = pd.DataFrame(tests)
    if not tests.empty:
        for c in ["p_t", "p_wilcoxon"]:
            tests[f"{c}_holm"] = np.nan
            # Holm 校正：同一组 map、同一个指标下，MICA 对其他各方法的 p 值放在一起校正
            for _, g in tests.groupby(["set", "variant", "aligner", "metric"]):
                ok = g[c].notna()
                if ok.any():
                    tests.loc[g.index[ok], f"{c}_holm"] = multipletests(g.loc[ok, c], method="holm")[1]
    tests.to_csv(out, sep="\t", index=False, float_format="%.6g")
    print(f"{len(df)} rows, methods {methods}; maps outside the manifest ignored: {len(outside)}; "
          f"ref {a.ref}, boot {a.boot}, seed {a.seed}; wrote {len(tests)} tests to {out}")


if __name__ == "__main__":
    main()
