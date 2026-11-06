import os, re, subprocess, sys, tempfile
from pathlib import Path
import gemmi

# 工具路径来自 config/paths.sh；没设的话用 paths.example.sh 里的占位路径，check_tools 会报出来。
PHENIX = os.environ.get("PHENIX_ENV") or "/path/to/phenix-1.20.1-4487/phenix_env.sh"
# "mmalign" = MM-align（benchmark 的分数全用它）；"usalign" = US-align -mm 1（只在校准时用）
ALIGN = {"mmalign": [os.environ.get("MMALIGN") or "/path/to/MMalign"],
         "usalign": [os.environ.get("USALIGN") or "/path/to/USalign", "-mm", "1"]}
TOOL_VAR = {"mmalign": "MMALIGN", "usalign": "USALIGN"}
COLS = ["tm", "tm_by_model", "aligned_len", "seq_id", "ca_match", "seq_match", "ca_quality", "model_len", "ref_len"]


def check_tools(aligners=()):
    """PHENIX_ENV 或要用的 aligner 没设、找不到时直接报错停下。"""
    need = [("PHENIX_ENV", PHENIX)] + [(TOOL_VAR[a], ALIGN[a][0]) for a in aligners]
    bad = [name for name, path in need if not path or not Path(path).is_file()]
    if bad:
        raise RuntimeError(f"not set or not found: {', '.join(bad)} (run `source config/paths.sh`)")


def protein(src, dst):
    """只留第一个 model 里的 protein（去掉氢、ligand、水和核酸），写成 mmCIF；返回有 Cα 的 residue 数。
    写成 mmCIF 是因为 PDB 格式最多只能放 62 条 chain，大的组装体放不下。"""
    st = gemmi.read_structure(str(src))
    while len(st) > 1:
        del st[len(st) - 1]
    st.setup_entities()  # 没有 entity 记录的 PDB 文件（例如 MICA 的输出）要先做这一步，才能正确去掉 ligand
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    for ch in st[0]:
        for i in reversed(range(len(ch))):
            info = gemmi.find_tabulated_residue(ch[i].name)
            if info is None or not info.is_amino_acid():
                del ch[i]
    st.remove_empty_chains()
    st.setup_entities()
    st.make_mmcif_document().write_file(str(dst))
    return sum(1 for ch in st[0] for r in ch if r.find_atom("CA", "*"))


def evaluate(model, ref, aligners=("mmalign",)):
    """一个 model 对一个 reference 打分，返回 {aligner: {指标: 值}}。"""
    check_tools(aligners)
    with tempfile.TemporaryDirectory() as t:
        m, r = Path(t) / "model.cif", Path(t) / "ref.cif"
        ref_len = protein(ref, r)
        if not Path(model).exists() or protein(model, m) == 0:
            return {a: dict.fromkeys(COLS, 0) | {"ref_len": ref_len} for a in aligners}  # 没有 model 或是空的：每项都记 0
        model_len = protein(model, m)
        # Cα match 和 sequence match：phenix.chain_comparison 默认参数，reference 当 target；
        # 读表头含 "MATCH(%)" 那一行下面的第一行，第 8 列是 FOUND(%)，第 10 列是 SEQ MATCH(%)。
        cc = subprocess.run(["bash", "-c", f"source {PHENIX} && phenix.chain_comparison {r} {m}"],
                            capture_output=True, text=True, cwd=t).stdout.splitlines()
        h = next(i for i, l in enumerate(cc) if "MATCH(%)" in l)
        row = next(l.split() for l in cc[h + 1:] if l.strip())
        ca_match, seq_match = float(row[7]), float(row[9])
        out = {}
        for a in aligners:
            txt = subprocess.run(ALIGN[a] + [str(m), str(r)], capture_output=True, text=True).stdout
            al = re.search(r"Aligned length=\s*(\d+).*?Seq_ID=n_identical/n_aligned=\s*([\d.]+)", txt)
            # MM-align 输出两个 TM-score：Structure_2（reference）归一化的是 tm，Structure_1（model）归一化的是 tm_by_model
            tm = {s: v for v, s in re.findall(r"TM-score=\s*([\d.]+) \(normalized by length of (Structure_\d)", txt)}
            # ca_quality = Cα match × model 长度 / reference 长度，这样算能重现 MICA 论文 Table S3
            out[a] = dict(tm=float(tm["Structure_2"]), tm_by_model=float(tm["Structure_1"]),
                          aligned_len=int(al[1]), seq_id=float(al[2]), ca_match=ca_match, seq_match=seq_match,
                          ca_quality=ca_match * model_len / ref_len, model_len=model_len, ref_len=ref_len)
        return out


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4) or (len(sys.argv) == 4 and sys.argv[3] not in ALIGN):
        sys.exit("usage: python evaluate.py <model> <reference> [mmalign|usalign]")
    a = sys.argv[3] if len(sys.argv) > 3 else "mmalign"
    res = evaluate(sys.argv[1], sys.argv[2], (a,))[a]
    print("\t".join(COLS))
    print("\t".join(str(round(res[c], 4)) for c in COLS))
