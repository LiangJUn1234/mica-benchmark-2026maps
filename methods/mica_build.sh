#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
# 后面会把 MICA 环境放到 PATH 最前面，所以先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 3 ] || { echo "usage: $0 <EMDB> <PDB> <resolution>" >&2; exit 2; }
emdb=$1; pdb=$2; res=$3; n=${emdb#EMD-}   # 去掉 "EMD-" 前缀，例如 EMD-56389 → 56389
P=$ENVS_DIR/mica_v5; PX=$PHENIX_ENV

# MICA 第 2 步：在 mica_dock.sh dock 好的 entry 目录上建模（README 的 2.6，run.py --run_phenix）。
# 分两步跑和直接跑 MICA_pipeline.sh 结果相同（benchmark 前核对过，TM-score 差 0.0000）。
OUT=$RESULTS_DIR/mica/$emdb
mkdir -p "$OUT" || exit 1
# 从 af3d 的 DONE 里取出 handoff= 后面的路径，也就是第 1 步打包的 entry 目录
H=$(sed -n 's/.*handoff=\([^ ]*\).*/\1/p' $RESULTS_DIR/af3d/$emdb/DONE 2>/dev/null)
[ -n "$H" ] && [ -f "$H" ] || { echo "no docking output ($RESULTS_DIR/af3d/$emdb/DONE missing or without handoff)" | tee $OUT/FAILED; exit 1; }

# 工作目录里放一份 MICA 和解开的 entry 目录（路径里不能有 'results'，原因见 mica_dock.sh）。
W=$(mktemp -d "$SCRATCH_DIR/mica_build_$emdb.XXXXXX") || exit 1; E=$W/e/$n
echo "$emdb $pdb res $res; work dir $W (kept after the run)"
mkdir -p $W/e && cp -a $TOOLS_DIR/MICA $W/MICA && mkdir $E && tar -xf $H -C $E && gunzip -c $INPUTS_DIR/$pdb/emd_$n.map.gz > $E/emd_$n.map

cd $W/MICA && export PATH=$P/bin:$PATH PYTHONNOUSERSITE=1
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
/usr/bin/time -v python run.py -m $E/emd_$n.map -f $E/$pdb.fasta -a $E/AF3_results -p modules/pulchra304/src/pulchra \
  --run_phenix -x $PX -r $res > $W/run.log 2>&1; RUN=$?
# run.py 退出时 Phenix 的 real_space_refine 可能还在后台跑：每 30 秒查一次，最多等 12 小时。
for i in $(seq 1 1440); do pgrep -u ${USER:-$(id -un)} -f "real_space_refine.*$W" > /dev/null || break; sleep 30; done
echo "run.py exit $RUN; Phenix wait ~$((i * 30)) s"
# 最终 model 用 Phenix refine 之前的那个（和 MICA 论文 Table S3 一样）；refine 后的 FR 也一起保存。
F=$W/MICA/output/${n}_${pdb}_MICA_all_atom_model.pdb; FR=${F%.pdb}_real_space_refined_000.pdb
# 一行 Python：去掉 ligand 和水之后数 protein Cα 的个数；文件不存在或读不了时记 0
NCA=$($PY -I -c "import gemmi,sys; s=gemmi.read_structure(sys.argv[1]); s.setup_entities(); s.remove_ligands_and_waters(); print(sum(1 for c in s[0] for r in c if r.find_atom('CA','*')))" $F 2>/dev/null || echo 0)
echo "final model $F: protein CA $NCA; refined present: $([ -f $FR ] && echo yes || echo no)"
cp $W/run.log $OUT/

# 成功 = 最终 model 至少有一个 protein Cα。
if [ "$NCA" -gt 0 ] 2>/dev/null; then   # NCA 不是数字时比较会报错，2>/dev/null 不显示错误，结果按失败处理
  cp $F $OUT/ && cp $F $OUT/model.pdb && { [ -f $FR ] && cp $FR $OUT/ || true; }
  echo "model=$OUT/model.pdb exit=$RUN" > $OUT/DONE; echo "DONE mica $n"
else
  { echo "no final model with protein CA (run.py exit $RUN)"; tail -50 $W/run.log; } > $OUT/FAILED; echo "FAILED mica $n"; exit 1
fi
