# benchmark 写出的所有东西都放在这里（默认是当前目录下的 work/）。路径用绝对路径（脚本会 cd 到各自的工作目录），
# 而且不能含 "results" 这个词（MICA 会把它改掉）。
# `: "${VAR:=默认值}"` 的意思是：VAR 已经有值就保留，没有才用默认值。所以可以只对一次运行改路径，
# 例如 header 规范化过的 map，输入从第二组文件夹读，结果也写进第二组：
#   INPUTS_DIR=$INPUTS_HDRNORM_DIR RESULTS_DIR=$RESULTS_HDRNORM_DIR bash methods/modelangelo.sh EMD-12345 9ABC
: "${BENCH_ROOT:=$PWD/work}"; export BENCH_ROOT
: "${DATASET_DIR:=$BENCH_ROOT/dataset}"; export DATASET_DIR
: "${MANIFEST:=$DATASET_DIR/manifest.tsv}"; export MANIFEST
# 每个 entry 一个文件夹：$INPUTS_DIR/<PDB>/emd_<编号>.map.gz 和 <PDB>.fasta
: "${INPUTS_DIR:=$BENCH_ROOT/inputs}"; export INPUTS_DIR
# header 规范化过的 map 和它们的 offsets.tsv（inputs/normalize_header.py 写的）
: "${INPUTS_HDRNORM_DIR:=$BENCH_ROOT/inputs_hdrnorm}"; export INPUTS_HDRNORM_DIR
: "${REFS_DIR:=$BENCH_ROOT/refs}"; export REFS_DIR
# 平移到规范化 map 坐标系里的 reference（scoring/score_model.py 写的）
: "${REFS_HDRNORM_DIR:=$BENCH_ROOT/refs_hdrnorm}"; export REFS_HDRNORM_DIR
: "${AF3_DIR:=$BENCH_ROOT/af3}"; export AF3_DIR
# 事先算好的 ColabFold MSA，每条序列一个文件夹（inputs/colabfold_msa.py 写，methods/cryozeta.sh 读）
: "${MSA_DIR:=$BENCH_ROOT/msa_colabfold}"; export MSA_DIR
# 方法输出和分数，每个方法每张 map 一个文件夹：$RESULTS_DIR/<方法简称>/<EMDB>/
: "${RESULTS_DIR:=$BENCH_ROOT/outputs}"; export RESULTS_DIR
# 在 header 规范化过的 map 上跑的结果
: "${RESULTS_HDRNORM_DIR:=$BENCH_ROOT/outputs_hdrnorm}"; export RESULTS_HDRNORM_DIR
# ModelAngelo 按 padding 偏移移回去以后的 model 和分数：$PADFIX_DIR/ma/<EMDB>/
: "${PADFIX_DIR:=$BENCH_ROOT/padfix}"; export PADFIX_DIR
: "${SCORES_DIR:=$BENCH_ROOT/scores}"; export SCORES_DIR
# MICA 论文公布的 test_2025 数据（Zenodo record 15756654），scoring/calibrate_test2025.py 重新打分用
: "${CALIB_DIR:=$BENCH_ROOT/calib}"; export CALIB_DIR
# 每次方法运行的工作目录建在这下面，跑完保留
: "${SCRATCH_DIR:=${TMPDIR:-/tmp}}"; export SCRATCH_DIR
# 装了 gemmi 的 Python 3，方法脚本用它准备输入、检查输出
: "${BENCH_PYTHON:=python3}"; export BENCH_PYTHON
# 下载好的方法仓库：MICA、EModelX_new、E3-CryoFold、DiffModeler（版本见 methods/versions.tsv）
: "${TOOLS_DIR:=/path/to/tools}"; export TOOLS_DIR
# 各方法的 conda 环境，每个环境一个文件夹（见 methods/envs/）
: "${ENVS_DIR:=/path/to/conda/envs}"; export ENVS_DIR
# 各方法的 weights 和 cache，每个方法一个文件夹
: "${WEIGHTS_DIR:=/path/to/weights}"; export WEIGHTS_DIR
: "${AF3_IMAGE:=/path/to/alphafold3.sif}"; export AF3_IMAGE
# AlphaFold3 的模型参数要直接向 Google DeepMind 申请
: "${AF3_WEIGHTS_DIR:=/path/to/af3/params}"; export AF3_WEIGHTS_DIR
: "${AF3_DB_DIR:=/path/to/af3/databases}"; export AF3_DB_DIR
# Phenix 1.20.1-4487 的 phenix_env.sh：打分时用 chain_comparison，MICA 和 EModelX 用它做 real-space refinement
: "${PHENIX_ENV:=/path/to/phenix-1.20.1-4487/phenix_env.sh}"; export PHENIX_ENV
# US-align 只在校准时用；benchmark 的 TM-score 全用 MM-align
: "${USALIGN:=/path/to/USalign}"; export USALIGN
: "${MMALIGN:=/path/to/MMalign}"; export MMALIGN
