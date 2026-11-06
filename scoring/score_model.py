import argparse, hashlib, os, sys, urllib.request
from pathlib import Path
import gemmi, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))  # 用 python -I 运行时，脚本所在目录不在 import 路径里；手动加上才能 import 同目录的 evaluate.py
from evaluate import COLS, evaluate


def sha(f):
    return hashlib.sha256(Path(f).read_bytes()).hexdigest()


def main():
    env = lambda var, default: os.environ.get(var) or default  # 优先用 config/paths.sh 设的路径
    bench = env("BENCH_ROOT", os.path.join(os.getcwd(), "work"))
    ap = argparse.ArgumentParser(description="Score one (method, map) result with evaluate.py.")
    ap.add_argument("results_dir"); ap.add_argument("method"); ap.add_argument("emdb"); ap.add_argument("pdb")
    ap.add_argument("final", choices=["0", "1"], help="1: write a failed row if there is no DONE")
    ap.add_argument("offsets", nargs="?", help="offsets.tsv of the header-normalized maps")
    ap.add_argument("--refs", default=env("REFS_DIR", os.path.join(bench, "refs")), help="reference directory (default $REFS_DIR)")
    ap.add_argument("--shifted-refs", default=env("REFS_HDRNORM_DIR", os.path.join(bench, "refs_hdrnorm")),
                    help="directory for translated references (default $REFS_HDRNORM_DIR)")
    a = ap.parse_args()
    m, e, p, final = a.method, a.emdb, a.pdb, a.final
    R = Path(a.results_dir) / m / e
    # 已经打过分就跳过；没有 DONE（还没跑完或跑失败）也跳过。
    # final=1 表示不会再重跑：这时没有 DONE 就写一行 failed=1，分数记 0。
    if (R / "metrics.tsv").exists() or not ((R / "DONE").exists() or final == "1"):
        sys.exit(print(f"{m} {e}: skipped"))  # print 返回 None，所以这一行是打印后正常退出（退出码 0）
    ref = Path(a.refs) / f"{p}.cif"
    if not ref.exists():
        # 本地没有 reference 就从 RCSB 下载，并把 sha256 记下来
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{p}.cif", ref)
        open(ref.parent / "refs_added.sha256", "a").write(f"{sha(ref)}  {ref.name}\n")
    # header 规范化过的 map（offsets.tsv 里有这个 PDB）：reference 要平移 ref_shift，才能和规范化后的 map 对齐。
    # 平移后的 reference 另存一份，原来的不改。
    if a.offsets and (o := pd.read_csv(a.offsets, sep="\t").set_index("pdb")).index.isin([p]).any():
        dst = Path(a.shifted_refs) / f"{p}.cif"
        if not dst.exists():
            st = gemmi.read_structure(str(ref))
            v = gemmi.Position(*o.loc[p, ["ref_shift_x", "ref_shift_y", "ref_shift_z"]].astype(float))
            for md in st:
                for ch in md:
                    for res in ch:
                        for at in res:
                            at.pos = at.pos + v
            dst.parent.mkdir(parents=True, exist_ok=True)
            st.make_mmcif_document().write_file(str(dst))
            open(dst.parent / "refs.sha256", "a").write(f"{sha(dst)}  {dst.name}  from {sha(ref)}\n")
        ref = dst
    # model 路径 = DONE 里第一个带 "/" 的字段（去掉开头的 'model='）；没有 DONE 时用一个不存在的路径，evaluate 会记 0 分
    done = (R / "DONE").read_text().split() if (R / "DONE").exists() else []
    model = next((f.removeprefix("model=") for f in done if "/" in f), "/nonexistent")
    row = evaluate(model, ref, ("mmalign",))["mmalign"] | dict(variant="main", aligner="mmalign", failed=int(not done),
                                                                model=model, ref=f"{ref.name} {sha(ref)[:16]}")
    R.mkdir(parents=True, exist_ok=True)
    # 先写 .part 再改名：中途中断也不会留下写了一半的 metrics.tsv
    pd.DataFrame([row])[["variant", "aligner", "failed"] + COLS + ["model", "ref"]].to_csv(
        R / "metrics.tsv.part", sep="\t", index=False)
    os.replace(R / "metrics.tsv.part", R / "metrics.tsv")
    print(f"{m} {e}: metrics.tsv written (failed={int(not done)})")


if __name__ == "__main__":
    main()
