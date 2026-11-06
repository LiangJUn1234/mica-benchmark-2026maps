#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
REPO_DIR=$(cd "$(dirname "$0")/.." && pwd)
# 先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径，用来放 AF3 model 和检查输出
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 2 ] || { echo "usage: $0 <EMDB> <PDB>" >&2; exit 2; }
EMDB=$1; PDB=$2; N=${EMDB#EMD-}; KEEP=${KEEP:-0}   # KEEP=1 时加 --keep_temp_files
P=$ENVS_DIR/emprot_flex
WT=$WEIGHTS_DIR/emprot
export TORCH_HOME=$WT/torch HF_HOME=$WT/hf XDG_CACHE_HOME=$WT/xdg
# 把环境的 lib/ 加进 LD_LIBRARY_PATH：README 的 troubleshooting 写了，不加的话在 GPU 节点上 llvmlite 载入失败
export PYTHONNOUSERSITE=1 LD_LIBRARY_PATH=$P/lib:${LD_LIBRARY_PATH:-}

OUT=$RESULTS_DIR/emp/$EMDB
mkdir -p "$OUT" || exit 1

SCR=$(mktemp -d "$SCRATCH_DIR/emp_$EMDB.XXXXXX") || exit 1
W=$SCR/$EMDB
case $W in *results*) echo "work dir contains 'results': $W"; exit 2;; esac   # 路径里有 results 就停
mkdir -p $W/out   # EMProt 可以接受已经存在的输出目录
echo "$EMDB $PDB; work dir $W (kept after the run)"

# 输入：map、FASTA，以及每种 chain 一个 AF3 单链 model
gunzip -c $INPUTS_DIR/$PDB/emd_$N.map.gz > $W/emd_$N.map || { echo "map decompress failed" | tee $OUT/FAILED; exit 1; }
cp $INPUTS_DIR/$PDB/$PDB.fasta $W/
$PY -I $REPO_DIR/inputs/af3_json.py layout $INPUTS_DIR $AF3_DIR $W $PDB || { echo "af3_json.py layout failed" | tee $OUT/FAILED; exit 1; }
CHAINS=$(cd $W && ls templates/*/ranked_0.pdb | sort -V)   # sort -V 按 entity 编号排序（_2 在 _10 前面），按这个顺序交给 --chain

# README 例 3 的命令（--seq 加 --chain），参数全是默认
cd $W
OPT=""; [ "$KEEP" = 1 ] && OPT="--keep_temp_files"
# 栈上限调成不限：EMProt 自带的 dockx 在大 map 上用默认 8 MB 栈会 segmentation fault；
# 在一张两种设置都能跑的 map 上核对过，输出逐字节相同
ulimit -s unlimited; echo "ulimit -s $(ulimit -s)"
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
/usr/bin/time -v $P/bin/emprot build --map emd_$N.map --seq $PDB.fasta --chain $CHAINS --device 0 --output out $OPT > $OUT/run.log 2>&1
RC=$?
echo "emprot build exit $RC"

# 最终 model 是 README 说的 out/output.cif，原名保存一份，再存一份叫 model.cif
FINAL=$W/out/output.cif
[ -s $FINAL ] && cp $FINAL $OUT/output.cif && cp $FINAL $OUT/model.cif

# 成功 = 有最终 model 而且至少一个 protein Cα（不只看退出码）
# 下面 <<'PY' 到 PY 之间是一段直接交给 Python 执行的代码：数氨基酸残基上的 Cα
NCA=$($PY -I - $OUT/model.cif 2>/dev/null <<'PY' || echo 0
import gemmi, sys
def aa(name):
    t = gemmi.find_tabulated_residue(name)
    return t is not None and t.is_amino_acid()
m = gemmi.read_structure(sys.argv[1])[0]
print(sum(1 for c in m for r in c if aa(r.name) and r.find_atom("CA", "*")))
PY
)
echo "protein CA atoms in final model: ${NCA:-0}"
if [ "${NCA:-0}" -ge 1 ]; then
  echo "$OUT/output.cif exit $RC" > $OUT/DONE; echo "DONE $EMDB"
else
  { echo "no final model with protein CA (exit $RC)"; tail -50 $OUT/run.log; } > $OUT/FAILED; echo "FAILED $EMDB"; exit 1
fi
