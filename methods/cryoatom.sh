#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
# 先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径，用来检查输出
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 2 ] || { echo "usage: $0 <EMDB> <PDB>" >&2; exit 2; }
# KEEP=1 时加 -k，CryFold 就不删 see_alpha_output/ 和 CryNet_round_x/；默认 0
emdb=$1; pdb=$2; n=${emdb#EMD-}; KEEP=${KEEP:-0}
# CryoAtom 用论文版本：作者网站上的 CryFold v1.0（CryoAtom 论文的 Reporting Summary 写的是 1.0）
P=$ENVS_DIR/cryoatom; W=$WEIGHTS_DIR/cryoatom; I=$INPUTS_DIR

OUT=$RESULTS_DIR/cat/$emdb
mkdir -p "$OUT" || exit 1

SCR=$(mktemp -d "$SCRATCH_DIR/cat_$emdb.XXXXXX") || exit 1; WD=$SCR/$emdb
mkdir -p $WD
echo "$emdb $pdb; work dir $WD (kept after the run)"
export TORCH_HOME=$W/torch HF_HOME=$W/hf XDG_CACHE_HOME=$W/xdg PYTHONNOUSERSITE=1
gunzip -c $I/$pdb/emd_$n.map.gz > $WD/emd_$n.map; FA=$I/$pdb/$pdb.fasta

# README 的命令，参数全是默认
OPT=(); [ "$KEEP" = 1 ] && OPT=(-k)
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
cd $WD && /usr/bin/time -v $P/bin/build -s $FA -v emd_$n.map -o cat_$n "${OPT[@]}" > $OUT/run.log 2>&1; rc=$?
echo "build exit $rc"
F=$WD/cat_$n/cat_$n.cif   # 最终 model 是 README 说的 <out>/<out>.cif，不是 _raw 那个
# 一行 Python：数氨基酸残基上的 Cα 个数；文件不存在或读不了时记 0
ca=$($PY -I -c "import gemmi,sys;s=gemmi.read_structure(sys.argv[1]);print(sum(1 for c in s[0] for r in c if gemmi.find_tabulated_residue(r.name) and gemmi.find_tabulated_residue(r.name).is_amino_acid() and r.find_atom('CA','*')))" $F 2>/dev/null || echo 0)
echo "final model $F: protein CA $ca"
[ -f $F ] && cp $F $OUT/ && cp $F $OUT/model.cif

# 成功 = 最终 model 至少有一个 protein Cα
# ca 不是数字时比较会报错，2>/dev/null 不显示错误，结果按失败处理
if [ "$ca" -ge 1 ] 2>/dev/null; then echo "model=$OUT/model.cif exit_code=$rc" > $OUT/DONE; echo "$emdb DONE"
else { echo "no final model with protein CA (exit_code=$rc, protein_CA=$ca)"; tail -50 $OUT/run.log; } > $OUT/FAILED; echo "$emdb FAILED"; exit 1; fi
