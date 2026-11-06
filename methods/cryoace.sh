#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
REPO_DIR=$(cd "$(dirname "$0")/.." && pwd)
# 先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径，用来准备输入和检查输出
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 4 ] || { echo "usage: $0 <EMDB> <PDB> <resolution> <contour>" >&2; exit 2; }
EMDB=$1; PDB=$2; RES=$3; CONTOUR=$4; N=${EMDB#EMD-}
# CryoACE（Boltz-1 的分支）要一个已经对齐进 map 的初始 model，这里用 MICA 第 1 步 dock 好的 AF3 model，由 cryoace_inputs.py 整理：
#   VARIANT=B（默认）：每种 chain 的全长 AF3 chain，叠到每份 copy 里 dock 进去最大的那一段上
#   VARIANT=A：dock 进去的片段原样使用
VARIANT=${VARIANT:-B}
[[ $VARIANT =~ ^[AB]$ ]] || { echo "VARIANT must be A or B"; exit 2; }
P=$ENVS_DIR/cryoace; Wt=$WEIGHTS_DIR/cryoace
CKPT=$Wt/dl/cryoace_boltz1_mapfit_alpha500_with_conditioners.ckpt
AF3D=$RESULTS_DIR/af3d

OUT=$RESULTS_DIR/ace/$EMDB
# 要先跑 mica_dock.sh：af3d/<EMDB>/DONE 里第一个像路径的字段（去掉 'model='）就是 dock 好的 model。
# 没有的话什么都不写就退出（这不算 CryoACE 失败）。
DOCK=$(awk 'NR==1{for(i=1;i<=NF;i++){f=$i; sub(/^model=/,"",f); if(f ~ /\//){print f; exit}}}' $AF3D/$EMDB/DONE 2>/dev/null)
[ -s "$DOCK" ] || { echo "no docked model in $AF3D/$EMDB/DONE (run mica_dock.sh first)"; exit 3; }
mkdir -p "$OUT" || exit 1

J=$(mktemp -d "$SCRATCH_DIR/ace_$EMDB.XXXXXX") || exit 1; WK=$J/$EMDB
mkdir -p $WK/in $WK/out $J/keops
echo "$EMDB $PDB res $RES contour $CONTOUR variant $VARIANT; docked model $DOCK; work dir $WK (kept after the run)"
export TORCH_HOME=$Wt/torch HF_HOME=$Wt/hf XDG_CACHE_HOME=$Wt/xdg MPLCONFIGDIR=$Wt/mpl PYTHONNOUSERSITE=1
# KeOps 编译到工作目录里（KeOps 2.3 不看 XDG，不设的话会编译到 ~/.cache）
export KEOPS_CACHE_FOLDER=$J/keops
# cryoace_libcompat 里放一个不带版本号的 libnvrtc.so，给 KeOps 的 -lnvrtc 用（PyTorch 只带了 libnvrtc.so.13）
export LIBRARY_PATH=$ENVS_DIR/cryoace_libcompat:${LIBRARY_PATH:-}
export LD_LIBRARY_PATH=$ENVS_DIR/cryoace_libcompat:$P/lib/python3.10/site-packages/nvidia/cu13/lib:${LD_LIBRARY_PATH:-}

# 输入：map、FASTA、dock 好的 model；variant B 还要放好 AF3 model；然后生成对齐的初始 model 和 Boltz 要的 FASTA
gunzip -c $INPUTS_DIR/$PDB/emd_$N.map.gz > $WK/in/emd_$N.map
cp -a $INPUTS_DIR/$PDB/$PDB.fasta $DOCK $WK/in/
[ $VARIANT = B ] && $PY -I $REPO_DIR/inputs/af3_json.py layout $INPUTS_DIR $AF3_DIR $WK/af3 $PDB > $WK/in/af3_layout.log 2>&1
$PY -I $REPO_DIR/methods/cryoace_inputs.py $WK/in/$PDB.fasta $WK/in/$(basename $DOCK) $VARIANT \
  $WK/af3/templates $WK/in $EMDB > $WK/in/prep.log 2>&1; PREP=$?
cp $WK/in/prep.log $OUT/

# README "Using The Checkpoint" 的命令，参数默认，加 --use_msa_server（运行时向 ColabFold 服务器要 MSA；没设 seed）
RC=99
if [ $PREP = 0 ]; then
  cd $WK/out
  # /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
  /usr/bin/time -v $P/bin/boltz predict $WK/in/$EMDB.fasta --density_map $WK/in/emd_$N.map \
    --aligned_model $WK/in/aligned_$VARIANT.cif --res $RES --thresh $CONTOUR --checkpoint $CKPT \
    --cryoace_density_conditioning --cryoace_atom_profile_conditioning --use_msa_server \
    --out_dir $WK/out --cache $Wt/boltz_cache > $OUT/run.log 2>&1
  RC=$?
fi
echo "cryoace_inputs.py exit $PREP; boltz exit $RC"

# 最终 model（Boltz 的 predictions/<name>/<name>_model_0.cif）要存在且至少有一个 protein Cα。
# GPU 内存不够时 Boltz 只在 log 里写 "ran out of memory, skipping batch"，不出 model；
# benchmark 把这种情况算资源失败，换内存更大的 GPU（H100/H200）重跑。
F=$WK/out/boltz_results_$EMDB/predictions/$EMDB/${EMDB}_model_0.cif
NCA=$($PY -I -c "import gemmi,sys
try:
    m=gemmi.read_structure(sys.argv[1])[0]; print(sum(1 for c in m for r in c if r.find_atom('CA','*') and gemmi.find_tabulated_residue(r.name) and gemmi.find_tabulated_residue(r.name).is_amino_acid()))
except Exception: print(0)" $F)
[ "$NCA" -ge 1 ] && cp -a $F $OUT/ && cp -a $F $OUT/model.cif
if [ "$NCA" -ge 1 ]; then
  echo "$OUT/model.cif boltz_exit=$RC protein_CA=$NCA" > $OUT/DONE; echo "DONE: $(cat $OUT/DONE)"
else
  { echo "no final model with protein CA (boltz exit $RC, prep exit $PREP)"; tail -50 $OUT/run.log 2>/dev/null || tail -50 $WK/in/prep.log; } > $OUT/FAILED
  echo "FAILED: see $OUT/FAILED"; exit 1
fi
