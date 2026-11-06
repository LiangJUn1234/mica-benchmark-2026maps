import argparse, os, sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python -I 不会把脚本所在目录加进 import 路径，手动加上
from evaluate import COLS, check_tools, evaluate

# 用本 benchmark 的打分流程，重打 MICA 论文公布的 test_2025 model（Zenodo record 15756654），和论文 Table S3 比。
# MICA 给了两个 model 文件（refine 前、refine 后），论文没说 Table S3 用的是哪个，所以两个都打，MM-align 和 US-align 各打一次。
VARIANTS = {"all_atom": "MICA_all_atom_model.pdb", "refined": "MICA_all_atom_model_real_space_refined_000.pdb"}


def run(models, truth, task):
    emdb, pdb, v = task
    res = evaluate(Path(models) / emdb / f"{emdb}_{pdb}_{VARIANTS[v]}", Path(truth) / emdb / f"{pdb}.pdb", ("mmalign", "usalign"))
    return [dict(emdb=emdb, pdb=pdb, model_file=v, aligner=a, **m) for a, m in res.items()]


def main():
    env = lambda var, default: os.environ.get(var) or default  # 优先用 config/paths.sh 设的路径
    calib = env("CALIB_DIR", os.path.join(env("BENCH_ROOT", os.path.join(os.getcwd(), "work")), "calib"))
    ap = argparse.ArgumentParser(description="Rescore MICA's test_2025 models and compare with Table S3.")
    ap.add_argument("out_dir"); ap.add_argument("n_proc", type=int)
    ap.add_argument("--models", default=os.path.join(calib, "Results/Results/2025_Test"),
                    help="default $CALIB_DIR/Results/Results/2025_Test")
    ap.add_argument("--truth", default=os.path.join(calib, "truth_2025/Test_Dataset/2025_Test"),
                    help="default $CALIB_DIR/truth_2025/Test_Dataset/2025_Test")
    ap.add_argument("--table", default=str(Path(__file__).resolve().parent / "table_s3_mica_test2025.tsv"))
    a = ap.parse_args()
    check_tools(("mmalign", "usalign"))
    s3 = pd.read_csv(a.table, sep="\t", dtype=str)  # 先按论文印的原样读成文字，下面再转成数字
    out = Path(a.out_dir)
    out.mkdir(parents=True)
    tasks = [(r.emdb, r.pdb, v) for r in s3.itertuples() for v in VARIANTS]
    with ProcessPoolExecutor(a.n_proc) as ex:
        ours = pd.DataFrame([row for rows in ex.map(partial(run, a.models, a.truth), tasks) for row in rows])
    paper = s3.drop(columns=["pdb"]).set_index("emdb").astype(float).add_prefix("paper_")
    df = ours.join(paper, on="emdb")
    for c in [c for c in COLS if f"paper_{c}" in df]:  # diff = 我们的分数 − 论文的分数
        df[f"diff_{c}"] = df[c] - df[f"paper_{c}"]
    df.to_csv(out / "calib_2025.tsv", sep="\t", index=False)
    lines = []
    for (v, al), g in df.groupby(["model_file", "aligner"]):
        d = g.diff_tm.abs()
        lines.append(f"{v:9s} {al:8s} |dTM|<=0.02: {(d <= 0.02).sum()}/{len(d)}  mean|dTM|={d.mean():.4f}  max={d.max():.4f} "
                     f"| mean|d| Ca_match={g.diff_ca_match.abs().mean():.2f} "
                     f"seq_match={g.diff_seq_match.abs().mean():.2f} model_len={g.diff_model_len.abs().mean():.1f}")
    (out / "calib_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
