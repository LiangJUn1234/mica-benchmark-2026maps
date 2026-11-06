#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
REPO_DIR=$(cd "$(dirname "$0")/.." && pwd)
# 后面会把 MICA 环境放到 PATH 最前面，所以先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 4 ] || { echo "usage: $0 <EMDB> <PDB> <resolution> <contour>" >&2; exit 2; }
emdb=$1; pdb=$2; res=$3; contour=$4; n=${emdb#EMD-}   # ${emdb#EMD-} 去掉 "EMD-" 前缀，例如 EMD-56389 → 56389
P=$ENVS_DIR/mica_v5; PX=$PHENIX_ENV

# MICA 第 1 步：把每种 chain 的 AF3 model dock 进 map（README 的 2.4 和 2.5 两步，只用 CPU 和 Phenix）。
# dock 好的 model 同时就是 benchmark 里的 AF3-docked baseline；整个 entry 目录交给 mica_build.sh 做第 2 步。
OUT=$RESULTS_DIR/af3d/$emdb
mkdir -p "$OUT" || exit 1

# MICA 会把 AF3_results 路径里的 'results' 改成 'PDBs'，所以这个工作目录的路径里不能有 'results'。
W=$(mktemp -d "$SCRATCH_DIR/mica_dock_$emdb.XXXXXX") || exit 1; E=$W/e/$n   # mktemp -d 新建一个不重名的工作目录，跑完保留
echo "$emdb $pdb res $res contour $contour; work dir $W (kept after the run)"
mkdir -p $W/e && cp -a $TOOLS_DIR/MICA $W/MICA && mkdir $E
gunzip -c $INPUTS_DIR/$pdb/emd_$n.map.gz > $E/emd_$n.map && cp $INPUTS_DIR/$pdb/$pdb.fasta $E/
$PY -I $REPO_DIR/inputs/af3_json.py layout $INPUTS_DIR $AF3_DIR $E $pdb; LAYOUT=$?   # $? = 上一条命令的退出码

cd $W/MICA && export PATH=$P/bin:$PATH PYTHONNOUSERSITE=1
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
/usr/bin/time -v python utils/process_AF3_results.py -f $E/$pdb.fasta -a $E/AF3_results > $W/run.log 2>&1; PROC=$?
# dock_in_map.py 会同时开 min(0.6 × CPU 核数, 48) 个 Phenix 进程；benchmark 里每次分配 48 核。
/usr/bin/time -v python utils/dock_in_map.py -m $E/emd_$n.map -c $contour -r $res -f $E/$pdb.fasta -a $E/AF3_results \
  -x $PX >> $W/run.log 2>&1; DOCK=$?
echo "layout exit $LAYOUT; process_AF3_results exit $PROC; dock_in_map exit $DOCK"
F=$E/${n}_af3_docked.pdb
# 一行 Python：去掉 ligand 和水之后数 protein Cα 的个数；文件不存在或读不了时记 0
NCA=$($PY -I -c "import gemmi,sys; s=gemmi.read_structure(sys.argv[1]); s.setup_entities(); s.remove_ligands_and_waters(); print(sum(1 for c in s[0] for r in c if r.find_atom('CA','*')))" $F 2>/dev/null || echo 0)
echo "docked model $F: protein CA $NCA"
cp $W/run.log $OUT/

# 成功 = dock 好的 model 至少有一个 protein Cα（退出码只记录，不作判断）。
if [ "$NCA" -gt 0 ] 2>/dev/null; then   # NCA 不是数字时比较会报错，2>/dev/null 不显示错误，结果按失败处理
  cp $F $OUT/ && cp $F $OUT/model.pdb
  # 把整个 entry 目录（FASTA、AF3_results、docked model）打包交给 mica_build.sh；
  # map 不打包，因为 mica_build.sh 会从 $INPUTS_DIR 重新解压
  H=$OUT/dock_handoff.tar
  tar -C $E -cf $H --exclude=emd_$n.map .
  echo "model=$OUT/model.pdb exit=0 handoff=$H" > $OUT/DONE; echo "DONE af3d $n"
else
  { echo "no docked model with protein CA (layout $LAYOUT process $PROC dock $DOCK)"; tail -50 $W/run.log; } > $OUT/FAILED
  echo "FAILED af3d $n"; exit 1
fi
