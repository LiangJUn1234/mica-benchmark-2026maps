#!/bin/bash
source "$(dirname "$0")/../config/paths.sh" || exit 1
set -uo pipefail
# 先记下 BENCH_PYTHON（装了 gemmi 的 Python）的完整路径，用来准备输入和检查输出
PY=$(command -v "$BENCH_PYTHON") || { echo "BENCH_PYTHON not found: $BENCH_PYTHON" >&2; exit 1; }
[ $# -eq 4 ] || { echo "usage: $0 <EMDB> <PDB> <resolution> <contour>" >&2; exit 2; }
EMDB=$1; PDB=$2; RES=$3; CONTOUR=$4
# 用作者的 Apptainer image 和作者没改过的 wrapper 脚本（run_cryozeta.sh、run_cryozeta_large.sh），README 默认参数
CZ=$WEIGHTS_DIR/cryozeta; SIF=$CZ/sif/CryoZeta.sif; ASSETS=$CZ/assets
N=${EMDB#EMD-}; NAME=$(echo $PDB | tr A-Z a-z)

OUT=$RESULTS_DIR/cz/$EMDB
mkdir -p "$OUT" || exit 1

SCR=$(mktemp -d "$SCRATCH_DIR/cz_$EMDB.XXXXXX") || exit 1
W=$SCR/$EMDB
[[ $W == *results* ]] && { echo "scratch path contains 'results'" > $OUT/FAILED; exit 6; }   # 路径里有 results 就停
mkdir -p $W/inputs $W/msa $W/out $W/cache/xdg $W/tmp && cd $W
echo "$EMDB $PDB res $RES contour $CONTOUR; work dir $W (kept after the run)"

# 输入：map、FASTA，以及写好 MSA 路径的 input.json。
# 下面 <<'PY' 到 PY 之间是一段直接交给 Python 执行的代码，最后打印 "模式 residue数"。
gunzip -c $INPUTS_DIR/$PDB/emd_$N.map.gz > $W/inputs/emd_$N.map && cp $INPUTS_DIR/$PDB/$PDB.fasta $W/inputs/
MODE=$($PY -I - "$W" "$PDB" "$N" "$RES" "$CONTOUR" "$MSA_DIR" <<'PY'
import hashlib, json, os, shutil, sys
w, pdb, n, res, contour, msa_root = sys.argv[1:]
seqs, cur = [], None
for l in open(f"{w}/inputs/{pdb}.fasta"):
    l = l.strip()
    if l.startswith(">"):
        cur = [l, ""]; seqs.append(cur)
    elif cur:
        cur[1] += l
ents, tokens = [], 0
for head, s in seqs:
    # RCSB 的 FASTA 头是 ">9ABC_1|Chains A, B, C|..."：数出几条 chain，就是这种 chain 的份数（count）。
    # CryoZeta 因此知道每种 chain 有几份 copy，这一点和只拿一份 AF3 model 的方法不同。
    chains = head.split("|")[1].replace("Chains ", "").replace("Chain ", "").split(", ")
    h = hashlib.sha256(s.encode()).hexdigest()
    src, dst = f"{msa_root}/{h}", f"{w}/msa/{h}"   # colabfold_msa.py 事先算好的 MSA，文件夹名是序列的 sha256
    if not os.path.isdir(dst):
        os.mkdir(dst)
        for f in ("uniref100_hits.a3m", "mmseqs_other_hits.a3m"):
            shutil.copy2(f"{src}/{f}", dst)  # 缺 MSA 会报错，整个脚本记 FAILED
    ents.append({"proteinChain": {"sequence": s, "count": len(chains),
                 "msa": {"precomputed_msa_dir": dst, "pairing_db": "uniref100"}}})
    tokens += len(s) * len(chains)
json.dump([{"name": pdb.lower(), "modelSeeds": [], "map_path": f"{w}/inputs/emd_{n}.map", "resolution": float(res),
            "contour_level": float(contour), "sequences": ents}], open(f"{w}/inputs/input.json", "w"), indent=1)
print(("large" if tokens > 2800 else "standard"), tokens)
PY
) || { echo "input preparation failed (missing MSA?)" > $OUT/FAILED; exit 7; }
set -- $MODE; MODE=$1; TOKENS=$2   # 把 "模式 residue数" 拆成两个变量
echo "mode $MODE (protein residues incl. copies = $TOKENS; README: > ~2,800 -> large-complex mode)"

# 所有 cache（包括现场编译的 extension）都放在这次运行的工作目录里，所以每次运行都为自己的 GPU 重新编译。
export PYTHONNOUSERSITE=1 CRYOZETA_CACHE_DIR=$W/cache XDG_CACHE_HOME=$W/cache/xdg TORCH_HOME=$W/cache/torch \
       HF_HOME=$W/cache/hf TMPDIR=$W/tmp APPTAINER_CACHEDIR=$W/cache/apptainer
# image 的 %environment 会在读入外面的环境变量之后再设这四个变量，所以要用 APPTAINERENV_* 传进去（它在 %environment 之后生效）。
# WC 用真实路径，在 wrapper 挂载的当前目录里面。
WC=$(readlink -f $W)/cache/cryozeta
mkdir -p $WC/torch_extensions $WC/inductor $WC/triton
export APPTAINERENV_CRYOZETA_CACHE_DIR=$WC APPTAINERENV_TORCH_EXTENSIONS_DIR=$WC/torch_extensions \
       APPTAINERENV_TORCHINDUCTOR_CACHE_DIR=$WC/inductor APPTAINERENV_TRITON_CACHE_DIR=$WC/triton
G=${CUDA_VISIBLE_DEVICES:-0}; G=${G%%,*}   # 用 CUDA_VISIBLE_DEVICES 里的第一张卡（没设就是 0）；${G%%,*} 取第一个逗号前面的部分
if [ $MODE = standard ]; then
  # 标准模式（≤ 2,800 residue，copy 都算上）：run_cryozeta.sh 跑两次，就是 README 默认的 combined 流程。
  # 第 1 次：detection 加 standard model（--skip-combine）；第 2 次：interpolation model 加 combine（--skip-detection）。
  echo "=== run_cryozeta.sh: detection + standard model" >> $OUT/run.log
  # /usr/bin/time -v 记下运行时间和内存峰值，和输出一起写进 run.log
  /usr/bin/time -v bash $CZ/sif/run_cryozeta.sh --sif $SIF --gpu $G --skip-combine inputs/input.json out $ASSETS >> $OUT/run.log 2>&1
  E1=$?; E2=skipped
  if [ $E1 = 0 ]; then
    echo "=== run_cryozeta.sh: interpolation model + combine" >> $OUT/run.log
    /usr/bin/time -v bash $CZ/sif/run_cryozeta.sh --sif $SIF --gpu $G --interpolation --skip-detection inputs/input.json out $ASSETS >> $OUT/run.log 2>&1
    E2=$?
  fi
  FINAL=$W/out/$NAME/CryoZeta-Final/${NAME}_sample_0.cif
else
  # 大结构模式（> 2,800 residue，README 的 Large Structure Inference）：run_cryozeta_large.sh，参数默认。
  # REGISTRATION=svd 加作者提供的 --registration svd；benchmark 里只用在默认跑法因 TEASER++ 每个 sample 都没给出结果
  # （"TEASER++ NPZ file not found"）而停下的 entry 上；这几张用 svd 也没有出 model，所以冻结结果里没有一个分数来自 svd。
  REGOPT=""; [ -n "${REGISTRATION:-}" ] && { REGOPT="--registration $REGISTRATION"; echo "registration=$REGISTRATION"; }
  echo "=== run_cryozeta_large.sh (defaults${REGOPT:+, $REGOPT})" >> $OUT/run.log
  /usr/bin/time -v bash $CZ/sif/run_cryozeta_large.sh --sif $SIF --gpu $G $REGOPT inputs/input.json out $ASSETS >> $OUT/run.log 2>&1
  E1=$?; E2=n/a
  FINAL=$W/out/combined.cif
fi
echo "exit codes: $E1 $E2"

# 保存最终 model；成功 = 最终 model 至少有一个 protein Cα
NCA=0
if [ -s $FINAL ]; then
  cp -n $FINAL $OUT/ && cp -n $FINAL $OUT/model.cif
  # 和别的脚本不同：查不到的残基名按 UNK 算，UNK 也算氨基酸，所以 CryoZeta 写的未知残基也计入
  NCA=$($PY -I -c "
import gemmi, sys
st = gemmi.read_structure(sys.argv[1]); st.setup_entities()
print(sum(1 for ch in st[0] for r in ch if r.find_atom('CA', '*') and (gemmi.find_tabulated_residue(r.name) or gemmi.find_tabulated_residue('UNK')).is_amino_acid()))
" $FINAL 2>/dev/null || echo 0)
fi
echo "final model $FINAL | protein CA $NCA"
if [ "$NCA" -gt 0 ] 2>/dev/null; then   # NCA 不是数字时比较会报错，2>/dev/null 不显示错误，结果按失败处理
  printf "final_model %s\nmode %s\nexit_codes %s %s\nprotein_CA %s\n" "$OUT/model.cif" $MODE $E1 $E2 $NCA > $OUT/DONE
  echo "DONE $EMDB"
else
  { echo "no final model with protein CA ($FINAL); exit codes $E1 $E2"; tail -50 $OUT/run.log; } > $OUT/FAILED
  echo "FAILED $EMDB"; exit 1
fi
