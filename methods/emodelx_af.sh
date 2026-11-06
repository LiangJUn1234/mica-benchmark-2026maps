#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
REPO_DIR=$(cd "$(dirname "$0")/.." && pwd)
# 先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径，用来放 AF3 template 和检查输出
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 3 ] || { echo "usage: $0 <EMDB> <PDB> <resolution>" >&2; exit 2; }
EMDB=$1; PDB=$2; RES=$3; N=${EMDB#EMD-}
P=$ENVS_DIR/emodelx_b; R=$TOOLS_DIR/EModelX_new; PX=$PHENIX_ENV
export TORCH_HOME=$WEIGHTS_DIR/emodelx/torch HF_HOME=$WEIGHTS_DIR/emodelx/hf XDG_CACHE_HOME=$WEIGHTS_DIR/emodelx/xdg \
       PYTHONNOUSERSITE=1

OUT=$RESULTS_DIR/emx/$EMDB
mkdir -p "$OUT" || exit 1
fail() { { echo "reason: $1"; tail -50 $OUT/run.log 2>/dev/null; } > $OUT/FAILED; echo "FAILED $EMDB: $1"; exit 1; }

# 工作目录（路径里不能有 'results'）；outputs 要先建好，EModelX 自己不建
SCR=$(mktemp -d "$SCRATCH_DIR/emx_$EMDB.XXXXXX") || exit 1
W=$SCR/$EMDB
mkdir -p $W/inputs $W/outputs
echo "$EMDB $PDB res $RES; work dir $W (kept after the run)"

# 输入：map、FASTA，以及当 template 用的 AF3 model（templates/<PDB>_<n>/ranked_0.pdb）
gunzip -c $INPUTS_DIR/$PDB/emd_$N.map.gz > $W/inputs/emd_$N.map || fail "map decompress"
cp $INPUTS_DIR/$PDB/$PDB.fasta $W/inputs/ || fail "fasta copy"
$PY -I $REPO_DIR/inputs/af3_json.py layout $INPUTS_DIR $AF3_DIR $W $PDB > $W/layout.log 2>&1 || fail "af3_json.py layout: $(tail -1 $W/layout.log)"
# EModelX 要在自己的仓库目录里跑（用相对路径找 phenix.sh、phenix.eff、./models），所以复制一份；weights 用链接，不复制
rsync -a --exclude=.git --exclude=models $R/ $W/EModelX_new/ && ln -s $R/models $W/EModelX_new/models || fail "repo copy"

# README 的写法，--protocol=temp_flex 加 --template_dir 指向我们的 AF3 model（不用 --download_afdb，
# 这样每个方法拿到同一批 AF3 model）；其余参数默认，包括 PULCHRA 和 Phenix refinement。
cd $W/EModelX_new
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
/usr/bin/time -v $P/bin/python run.py --protocol=temp_flex --EM_map=$W/inputs/emd_$N.map --fasta=$W/inputs/$PDB.fasta \
  --template_dir=$W/templates --output_dir=$W/outputs --run_pulchra --pulchra_path=modules/pulchra304/src/pulchra \
  --run_phenix --phenix_act=$PX --resolution=$RES > $OUT/run.log 2>&1
RC=$?
# run.py 把 phenix.real_space_refine 放到后台（phenix.sh 里是 `$3 &`）就马上返回，所以要等它跑完
sleep 20
while pgrep -u ${USER:-$(id -un)} -f -- "real_space_refine.*$W/outputs" > /dev/null; do sleep 30; done
echo "run.py exit $RC"

# 最终 model = README 流程最后的 refine 后 model；refine 前的也保存
F=$(ls $W/outputs/*'_EModelX(+AF)_all_atom_model_real_space_refined_000.pdb' 2>/dev/null | head -1)
NCA=0
if [ -n "$F" ]; then
  # 一段 Python：数氨基酸残基上的 Cα 个数
  NCA=$($PY -I -c "
import gemmi, sys
n = 0
for c in gemmi.read_structure(sys.argv[1])[0]:
    for r in c:
        t = gemmi.find_tabulated_residue(r.name)
        if t and t.is_amino_acid(): n += sum(a.name == 'CA' for a in r)
print(n)" "$F" 2>/dev/null || echo 0)
  cp -n -- "$F" $OUT/ && cp -n -- "$F" $OUT/model.pdb
  cp -n -- $W/outputs/*'_EModelX(+AF)_all_atom_model.pdb' $OUT/ 2>/dev/null
fi
echo "final model ${F:-none}  protein CA $NCA"

# 成功 = 有最终 model 而且至少一个 protein Cα（不只看退出码）
if [ "$NCA" -ge 1 ]; then printf 'model %s\nrun.py_exit %s\n' "$OUT/model.pdb" $RC > $OUT/DONE; echo "DONE $EMDB"
else fail "no final model with protein CA (run.py exit $RC)"; fi
