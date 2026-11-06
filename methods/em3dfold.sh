#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
# 先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径，用来检查输出
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 2 ] || { echo "usage: $0 <EMDB> <PDB>" >&2; exit 2; }
EMDB=$1; PDB=$2; N=${EMDB#EMD-}; KEEP=${KEEP:-0}   # KEEP=1 时加 --keep-temp-files
P=$ENVS_DIR/em3dfold_v2
# weights 在 $W/dl/weights：EM3DFold v1.2、RiNALMo、ESM-2 650M，用作者的 scripts/download_with_mamba.sh 下载
W=$WEIGHTS_DIR/em3dfold

OUT=$RESULTS_DIR/em3d/$EMDB
mkdir -p "$OUT" || exit 1
fail() {  # 写 FAILED（原因加 run.log 最后 50 行）然后退出
  { echo "FAILED: $1"; echo "--- last 50 lines of run.log ---"; tail -50 $OUT/run.log 2>/dev/null; } > $OUT/FAILED
  echo "FAILED: $1"; exit 1
}

SCR=$(mktemp -d "$SCRATCH_DIR/em3d_$EMDB.XXXXXX") || exit 1
X=$SCR/$EMDB
case $X in *results*) fail "work dir path contains 'results': $X";; esac   # 路径里有 results 就停
mkdir -p $X/out || fail "cannot create work dir $X"
echo "$EMDB $PDB; work dir $X (kept after the run)"

MAPGZ=$INPUTS_DIR/$PDB/emd_$N.map.gz
FASTA=$INPUTS_DIR/$PDB/$PDB.fasta
[ -s $MAPGZ ] && [ -s $FASTA ] || fail "missing input $MAPGZ or $FASTA"
gzip -dc $MAPGZ > $X/emd_$N.map || fail "cannot decompress $MAPGZ"
cp $FASTA $X/$PDB.fasta

# README 的命令 em3dfold build，参数全是默认
export EM_WEIGHTS_DIR=$W/dl/weights TORCH_HOME=$W/torch HF_HOME=$W/hf XDG_CACHE_HOME=$W/xdg
export PYTHONNOUSERSITE=1
export LD_LIBRARY_PATH=$P/lib:${LD_LIBRARY_PATH:-}   # 不加的话在 GPU 节点上 llvmlite 载入不了自己的共享库
OPT=(); [ "$KEEP" = 1 ] && OPT=(--keep-temp-files)
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
cd $X && /usr/bin/time -v $P/bin/em3dfold build --map emd_$N.map --protein $PDB.fasta --output out --device 0 "${OPT[@]}" \
  > $OUT/run.log 2>&1
RC=$?
echo "em3dfold build exit $RC"

# 最终 model 是 README 说的 out/output.cif（原名和 model.cif 各存一份），EM3DFold 自己的 log 也保存
[ -e $X/out/output.cif ] && cp $X/out/output.cif $OUT/output.cif && cp $X/out/output.cif $OUT/model.cif
[ -e $X/out/run.log ] && cp $X/out/run.log $OUT/em3dfold_run.log

# 成功 = 退出码 0，而且最终 model 至少有一个 protein Cα
NCA=$($PY -I -c "
import sys, gemmi
st = gemmi.read_structure(sys.argv[1])
aa = lambda n: (gemmi.find_tabulated_residue(n) is not None and gemmi.find_tabulated_residue(n).is_amino_acid())
print(sum(1 for ch in st[0] for r in ch if r.find_atom('CA', '*') and aa(r.name)))  # same protein rule as in scoring
" $OUT/output.cif 2>/dev/null || echo 0)
echo "protein CA $NCA"
# ${NCA:-0} = NCA 没有值时当作 0；2>/dev/null 让 NCA 不是数字时不显示比较的错误，按失败处理
[ "$RC" = 0 ] && [ "${NCA:-0}" -ge 1 ] 2>/dev/null || fail "exit code $RC, protein CA in output.cif: ${NCA:-0}"
echo "$OUT/output.cif exit_code $RC protein_CA $NCA" > $OUT/DONE
echo "DONE: $(cat $OUT/DONE)"
