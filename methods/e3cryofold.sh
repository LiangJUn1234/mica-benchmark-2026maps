#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
# 先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径，用来检查输出
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 2 ] || { echo "usage: $0 <EMDB> <PDB>" >&2; exit 2; }
EMDB=$1; PDB=$2
# 用的是作者更新后的仓库，weights 是它 release 里的 models.zip（论文当时的 weights 链接已经失效）
P=$ENVS_DIR/e3cryofold_v2; R=$TOOLS_DIR/E3-CryoFold
# save_name 里不能有 '_' 或 '.'：它的 run_pulchra 从文件名里解析 chain id。所以写成 <编号>-<pdb 小写>
N=${EMDB#EMD-}; NAME=$N-$(echo $PDB | tr 'A-Z' 'a-z')
export TORCH_HOME=$WEIGHTS_DIR/e3cryofold/torch HF_HOME=$WEIGHTS_DIR/e3cryofold/hf \
       XDG_CACHE_HOME=$WEIGHTS_DIR/e3cryofold/xdg PYTHONNOUSERSITE=1

OUT=$RESULTS_DIR/e3cf/$EMDB
mkdir -p "$OUT" || exit 1

# 工作目录里放一份代码（不复制最上层的 .git 和 models/；chroma/models 要留着），models/ 用链接指向 weights
SCR=$(mktemp -d "$SCRATCH_DIR/e3cf_$EMDB.XXXXXX") || exit 1
WD=$SCR/$EMDB
mkdir -p $WD/inputs $WD/outputs
echo "$EMDB $PDB; work dir $WD (kept after the run)"
rsync -a --exclude /.git --exclude /models $R/ $WD/E3-CryoFold/ && ln -s $R/models $WD/E3-CryoFold/models
gzip -dc $INPUTS_DIR/$PDB/emd_$N.map.gz > $WD/inputs/emd_$N.map && cp -a $INPUTS_DIR/$PDB/$PDB.fasta $WD/inputs/

# README 的 pre_align 命令，参数默认（--t 0.1；默认会跑 PULCHRA）。
# --save_dir 结尾要有 '/'：它的代码用字符串直接拼路径。
cd $WD/E3-CryoFold
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
/usr/bin/time -v $P/bin/python inference.py --map_path $WD/inputs/emd_$N.map --fasta_path $WD/inputs/$PDB.fasta \
  --save_dir $WD/outputs/ --save_name $NAME --protocol pre_align --t 0.1 > $OUT/run.log 2>&1
RC=$?
echo "inference.py exit $RC"
F=$WD/outputs/${NAME}_all_atom_model.pdb
# 一段 Python：数氨基酸残基上的 Cα 个数；读不了时记 0
NCA=$($PY -I -c "
import gemmi, sys
try: m = gemmi.read_structure(sys.argv[1])[0]
except Exception: print(0); sys.exit()
aa = lambda r: (gemmi.find_tabulated_residue(r.name) is not None and gemmi.find_tabulated_residue(r.name).is_amino_acid())
print(sum(1 for c in m for r in c if r.find_atom('CA', '*') and aa(r)))" $F 2>/dev/null || echo 0)
[ -s $F ] && cp -an $F $OUT/ && cp -an $F $OUT/model.pdb

# 成功 = 最终 model 至少有一个 protein Cα
if [ "$NCA" -gt 0 ] 2>/dev/null; then   # NCA 不是数字时比较会报错，2>/dev/null 不显示错误，结果按失败处理
  printf "final_model %s\nexit_code %s\n" "$OUT/$(basename $F)" "$RC" > $OUT/DONE; echo "DONE $EMDB CA=$NCA rc=$RC"
else
  { echo "reason: final model missing or without protein CA (exit_code $RC, CA $NCA)"; echo "--- last 50 lines of run.log"; tail -50 $OUT/run.log; } > $OUT/FAILED
  echo "FAILED $EMDB rc=$RC"; exit 1
fi
