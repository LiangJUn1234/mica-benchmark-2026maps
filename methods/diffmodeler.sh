#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
REPO_DIR=$(cd "$(dirname "$0")/.." && pwd)
# 先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径，用来准备输入和检查输出
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 4 ] || { echo "usage: $0 <EMDB> <PDB> <resolution> <contour>" >&2; exit 2; }
EMDB=$1; PDB=$2; RES=$3; CONTOUR=$4; N=${EMDB#EMD-}
P=$ENVS_DIR/diffmodeler_b
RP=$TOOLS_DIR/DiffModeler

OUT=$RESULTS_DIR/dm/$EMDB
mkdir -p "$OUT" || exit 1
fail() {  # 写 FAILED（原因加 run.log 最后 50 行）然后退出
  { echo "FAILED: $1"; echo "--- last 50 lines of run.log ---"; tail -50 $OUT/run.log 2>/dev/null; } > $OUT/FAILED
  echo "FAILED: $1"; exit 1
}
# --contour 照 README 的建议取 0.5 × 作者 contour，作者 contour 用 EMDB 推荐的 contour level
CH=$(awk -v c="$CONTOUR" 'BEGIN{printf "%.6g", 0.5 * c}')
echo "$EMDB $PDB res $RES contour $CONTOUR; --contour $CH (README: 0.5 x author contour)"

SCR=$(mktemp -d "$SCRATCH_DIR/dm_$EMDB.XXXXXX") || exit 1
X=$SCR/$EMDB
case $X in *results*) fail "work dir path contains 'results': $X";; esac   # af3_json.py layout 不接受含 results 的路径
mkdir -p $X/out $X/dm_in || fail "cannot create work dir $X"
echo "work dir $X (kept after the run)"

# 单链结构用每种 chain 排第一的 AF3 model，代替 README 里从 AlphaFold DB / EBI 搜的结构；mode 0 不用数据库，也不连外部服务。
MAPGZ=$INPUTS_DIR/$PDB/emd_$N.map.gz
FASTA=$INPUTS_DIR/$PDB/$PDB.fasta
[ -s $MAPGZ ] && [ -s $FASTA ] || fail "missing input $MAPGZ or $FASTA"
gzip -dc $MAPGZ > $X/emd_$N.map || fail "cannot decompress $MAPGZ"
cp $FASTA $X/$PDB.fasta
$PY -I $REPO_DIR/inputs/af3_json.py layout $INPUTS_DIR $AF3_DIR $X/af3 $PDB > $X/layout.log 2>&1 || fail "af3_json.py layout: $(tail -1 $X/layout.log)"
# 下面 <<'EOF' 到 EOF 之间是一段 Python：写 -M 要的对照表，每行 "<单链文件> <这种 chain 的全部 chain 名>"。
# chain 名用 RCSB FASTA 头里的 label id（"Chains A, B[auth C]" → A B），所以 DiffModeler 知道每种 chain 有几份 copy。
$PY -I - $X/$PDB.fasta $X/af3/templates $X/dm_in $X/dm_info.txt <<'EOF' || fail "cannot build DiffModeler single-chain inputs"
import re, shutil, sys
from pathlib import Path
fasta, tpl, out, info = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
lines = []
for b in Path(fasta).read_text().split(">")[1:]:
    head = b.strip().split("\n")[0]
    eid = head.split("|")[0]
    m = re.match(r"Chains? (.*)", head.split("|")[1].strip())
    chains = [c.strip().split("[")[0].strip() for c in m.group(1).split(",")]
    shutil.copy(tpl / eid / "ranked_0.pdb", out / f"{eid}.pdb")
    lines.append(f"{eid}.pdb " + " ".join(chains))
Path(info).write_text("\n".join(lines) + "\n")
print("\n".join(lines))
EOF
cp $X/dm_info.txt $OUT/dm_info.txt
# main.py 每次启动都会在仓库里 compileall VESPER_CUDA，所以复制一份仓库到工作目录再跑；weights（best_model/）用链接
mkdir -p $X/DiffModeler && ( cd $RP && tar --exclude=./best_model -cf - . ) | ( cd $X/DiffModeler && tar -xf - ) \
  && ln -s $RP/best_model $X/DiffModeler/best_model || fail "cannot copy the DiffModeler repo to $X"

# README 的 template mode（--mode=0），作者的 config 不改；diffusion 的临时文件写在 --output 下，合并后自己删掉
export PATH=$P/bin:$PATH LD_LIBRARY_PATH=$P/lib:${LD_LIBRARY_PATH:-} PYTHONNOUSERSITE=1
export TORCH_HOME=$WEIGHTS_DIR/diffmodeler/torch XDG_CACHE_HOME=$WEIGHTS_DIR/diffmodeler/xdg
# /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
cd $X/DiffModeler && /usr/bin/time -v $P/bin/python main.py --mode=0 -F=$X/emd_$N.map -P=$X/dm_in -M=$X/dm_info.txt \
  --config=config/diffmodeler.json --contour=$CH --gpu=0 --resolution=$RES --output=$X/out > $OUT/run.log 2>&1
RC=$?
echo "main.py exit $RC"

# 最终 model 是 README 说的 <output>/DiffModeler.cif（原名和 model.cif 各存一份）；带 fit score 的那份也保存
[ -e $X/out/DiffModeler.cif ] && cp $X/out/DiffModeler.cif $OUT/DiffModeler.cif && cp $X/out/DiffModeler.cif $OUT/model.cif
[ -e $X/out/DiffModeler_fitscore.cif ] && cp $X/out/DiffModeler_fitscore.cif $OUT/

# 成功 = 退出码 0，而且最终 model 至少有一个 protein Cα（和打分时判断 protein 的规则相同）
NCA=$($PY -I -c "
import sys, gemmi
st = gemmi.read_structure(sys.argv[1])
aa = lambda n: (gemmi.find_tabulated_residue(n) is not None and gemmi.find_tabulated_residue(n).is_amino_acid())
print(sum(1 for ch in st[0] for r in ch if r.find_atom('CA', '*') and aa(r.name)))
" $OUT/model.cif 2>/dev/null || echo 0)
echo "protein CA $NCA"
# ${NCA:-0} = NCA 没有值时当作 0；2>/dev/null 让 NCA 不是数字时不显示比较的错误，按失败处理
[ "$RC" = 0 ] && [ "${NCA:-0}" -ge 1 ] 2>/dev/null || fail "exit code $RC, protein CA in DiffModeler.cif: ${NCA:-0}"
echo "$OUT/model.cif exit_code $RC protein_CA $NCA" > $OUT/DONE
echo "DONE: $(cat $OUT/DONE)"
