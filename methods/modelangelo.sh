#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
[ $# -eq 2 ] || { echo "usage: $0 <EMDB> <PDB>" >&2; exit 2; }
# KEEP=1 时加 --keep-intermediate-results，只是让 ModelAngelo 最后不删中间文件；默认 0，和 README 一样
EMDB=$1; PDB=$2; KEEP=${KEEP:-0}
P=$ENVS_DIR/modelangelo
W=$WEIGHTS_DIR/modelangelo

OUT=$RESULTS_DIR/ma/$EMDB
mkdir -p "$OUT" || exit 1
# fail：写 FAILED（原因加 run.log 最后 50 行）然后退出
fail() { { echo "FAILED: $1"; echo "--- last 50 lines of run.log ---"; tail -50 $OUT/run.log 2>/dev/null; } > $OUT/FAILED; echo "FAILED $EMDB: $1"; exit 1; }

ROOT=$(mktemp -d "$SCRATCH_DIR/ma_$EMDB.XXXXXX") || exit 1
WD=$ROOT/$EMDB
mkdir -p $WD
echo "$EMDB $PDB; work dir $WD (kept after the run)"
# ModelAngelo 从 TORCH_HOME 读 weights；事先确认 weights 在，否则它会自己下载到别的地方
export PYTHONNOUSERSITE=1 TORCH_HOME=$W HF_HOME=$W/hf XDG_CACHE_HOME=$ROOT/xdg_cache
[ -f $W/hub/checkpoints/model_angelo_v1.0/nucleotides/success.txt ] && [ -f $W/hub/checkpoints/esm1b_t33_650M_UR50S.pt ] \
  || fail "weights missing under TORCH_HOME=$W (refusing to let ModelAngelo download into another location)"

N=${EMDB#EMD-}; MAPGZ=$INPUTS_DIR/$PDB/emd_$N.map.gz; FASTA=$INPUTS_DIR/$PDB/$PDB.fasta
[ -s $MAPGZ ] && [ -s $FASTA ] || fail "input missing: $MAPGZ or $FASTA"
gunzip -c $MAPGZ > $WD/emd_$N.map || fail "gunzip $MAPGZ"

# README 的命令 model_angelo build，参数全是默认
KEEPFLAG=""; [ "$KEEP" = 1 ] && KEEPFLAG="--keep-intermediate-results"
CMD="$P/bin/model_angelo build -v emd_$N.map -pf $FASTA -o out $KEEPFLAG"
echo "command: $CMD"
cd $WD
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
/usr/bin/time -v $CMD > $OUT/run.log 2>&1
RC=$?
echo "model_angelo exit $RC"

# 最终 model 是 README 说的 out/out.cif，另存一份叫 model.cif
[ -f out/out.cif ] && cp -n out/out.cif $OUT/out.cif && cp -n out/out.cif $OUT/model.cif
[ -f out/model_angelo.log ] && cp -n out/model_angelo.log $OUT/model_angelo.log
# 一段 Python（ModelAngelo 环境里的 Biopython）：数 protein Cα 个数和 chain 数
read NCA NCH < <($P/bin/python -c "
import sys
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
m = MMCIFParser(QUIET=True).get_structure('m', sys.argv[1])[0]
print(sum(1 for r in m.get_residues() if is_aa(r, standard=False) and 'CA' in r), len(m))
" out/out.cif 2>/dev/null)
echo "protein CA ${NCA:-0}; chains ${NCH:-0}"   # ${NCA:-0} = NCA 没有值时当作 0

# 成功 = 退出码 0、有最终 model、至少一个 protein Cα，而且 log 里有 "refinement round 3"（最后一轮跑完了）
[ $RC -eq 0 ] || fail "model_angelo exit code $RC"
[ -s out/out.cif ] || fail "final model out/out.cif missing"
[ "${NCA:-0}" -ge 1 ] || fail "final model has no protein CA"
grep -q 'refinement round 3' out/model_angelo.log 2>/dev/null || fail "'refinement round 3' not in model_angelo.log"
printf 'model %s\nexit_code %s\n' $OUT/model.cif $RC > $OUT/DONE
echo "DONE $EMDB -> $OUT/model.cif"
# 需要 padding 的 map 上，坐标要再用 modelangelo_pad_offset.py 和 modelangelo_apply_offset.py 移回去（benchmark 里是 EMD-47284、EMD-70558）。
