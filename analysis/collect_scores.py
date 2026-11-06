import argparse, hashlib, json, os
from pathlib import Path
import pandas as pd

PADFIX_EMDB = ("EMD-47284", "EMD-70558")  # 这两张 map 上 ModelAngelo 的 model 有 padding 造成的偏移


def main():
    env = lambda var, default: os.environ.get(var) or default  # 优先用 config/paths.sh 设的路径
    bench = env("BENCH_ROOT", os.path.join(os.getcwd(), "work"))
    ap = argparse.ArgumentParser(description="Collect metrics.tsv files into per-map tables, one per analysis set.")
    ap.add_argument("out")
    ap.add_argument("--sets", default=str(Path(__file__).resolve().parent / "map_sets.json"))
    ap.add_argument("--results", default=env("RESULTS_DIR", os.path.join(bench, "outputs")),
                    help="results on the original maps (default $RESULTS_DIR)")
    ap.add_argument("--hdrnorm", default=env("RESULTS_HDRNORM_DIR", os.path.join(bench, "outputs_hdrnorm")),
                    help="results on the header-normalized maps (default $RESULTS_HDRNORM_DIR)")
    ap.add_argument("--padfix", default=env("PADFIX_DIR", os.path.join(bench, "padfix")),
                    help="pad-offset-corrected ModelAngelo scores (default $PADFIX_DIR), or none")
    a = ap.parse_args()
    D = Path(env("DATASET_DIR", os.path.join(bench, "dataset")))  # map_sets.json 里排除表的路径是相对于这个目录的

    cfg, out = json.load(open(a.sets)), Path(a.out)
    out.mkdir(parents=True)
    # 读一个结果目录下所有 <方法>/<EMDB>/metrics.tsv，从路径里取出方法和 EMDB 两列
    read = lambda root: pd.concat([pd.read_csv(p, sep="\t").assign(method=p.parts[-3], emdb=p.parts[-2])
                                   for p in sorted(Path(root).glob("*/EMD-*/metrics.tsv"))]
                                  or [pd.DataFrame(columns=["method", "emdb", "failed", "ref_len"])], ignore_index=True)
    orig, hdr = read(a.results), read(a.hdrnorm)

    # ModelAngelo 在这两张 map 上的那一行，整行换成坐标移回去以后的分数（位置不变）；换下来的旧行和新行一起另存一个表
    if a.padfix != "none":
        sup = []
        for e in PADFIX_EMDB:
            f = Path(a.padfix) / "ma" / e / "metrics.tsv"
            new = pd.read_csv(f, sep="\t").assign(method="ma", emdb=e)
            i = hdr.index[(hdr.method == "ma") & (hdr.emdb == e)]
            assert len(i) == 1 and len(new) == 1 and list(new.columns) == list(hdr.columns), (e, len(i))
            old = hdr.loc[i].copy()
            for c in hdr.columns:
                hdr.at[i[0], c] = new.at[0, c]
            for ver, d, src in (("uncorrected (ModelAngelo output as is)", old, Path(a.hdrnorm) / "ma" / e / "metrics.tsv"),
                                ("pad-offset corrected (main table)", new, f)):
                sup.append(d.assign(version=ver, source=str(src), source_sha256=hashlib.sha256(src.read_bytes()).hexdigest()))
        pd.concat(sup, ignore_index=True).to_csv(out / "S_ma_padfix_uncorrected.tsv", sep="\t", index=False)

    hn = cfg["header_normalized"]
    # 每个分析集（报告里只有 v2_main46）：header 规范化过的 map 取规范化后的结果，其余取原来的结果；
    # 再去掉组装体 1 超过大小上限的 map（size_cap/sampled90.tsv 里 within_cap_asm1 = 0）和 exclude_emdb 里的 map。
    # 失败的格子保留 failed = 1，paired_tests.py 把它的每项分数记 0。
    for s in cfg["sets"]:
        d = pd.concat([orig[~orig.emdb.isin(hn)], hdr[hdr.emdb.isin(hn)]], ignore_index=True)
        if "exclude" in s:
            x = pd.read_csv(D / s["exclude"]["file"], sep="\t")
            d = d[~d.emdb.isin(x[x[s["exclude"]["column"]] == s["exclude"]["value"]].emdb)]
        if "exclude_emdb" in cfg:
            d = d[~d.emdb.isin(cfg["exclude_emdb"])]
        (out / s["name"]).mkdir()
        d.to_csv(out / s["name"] / "input.tsv", sep="\t", index=False)
        print(f"{s['name']}: {len(d)} rows, {d.emdb.nunique()} maps -> {out / s['name'] / 'input.tsv'}")


if __name__ == "__main__":
    main()
