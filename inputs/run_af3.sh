#!/bin/bash
source "$(dirname "$0")/../config/paths.sh"
set -uo pipefail
[ $# -eq 1 ] || { echo "usage: $0 <json_list.txt>" >&2; exit 2; }
mapfile -t jobs < "$1"   # 每行一个 AF3 输入 JSON 的路径（af3_json.py json 写的），要在 $AF3_DIR 下面

# AF3 3.0.1 官方容器，除了路径和分两步跑的开关以外，参数全是默认：
# max_template_date = 2021-09-30，5 个 diffusion sample，10 次 recycle；seed 用 JSON 里的 1。
for J in "${jobs[@]}"; do
  [ -n "$J" ] || continue
  N=$(basename "$J" .json); echo "json $J"   # N = 作业名，也就是文件名去掉 .json

  # 第 1 步：只跑 data pipeline（MSA 和 template 搜索，只用 CPU）。输出已存在就跳过。
  if [ -e "$AF3_DIR/msa/$N/${N}_data.json" ]; then
    echo "already done: $AF3_DIR/msa/$N"
  else
    # -B 把这两个目录挂进容器，容器里才看得到
    singularity exec -B "$AF3_DIR,$AF3_DB_DIR" "$AF3_IMAGE" \
      python /app/alphafold3/run_alphafold.py --json_path="$J" --output_dir="$AF3_DIR/msa" \
      --db_dir="$AF3_DB_DIR" --run_inference=false
    echo "af3 msa exit $?"
  fi

  # 第 2 步：只跑 inference（GPU），输入是第 1 步写的 data.json。输出已存在就跳过。
  D=$AF3_DIR/msa/$N/${N}_data.json; echo "data $D"
  if [ -e "$AF3_DIR/out/$N/${N}_model.cif" ]; then
    echo "already done: $AF3_DIR/out/$N"
  else
    singularity exec --nv -B "$AF3_DIR,$AF3_WEIGHTS_DIR" "$AF3_IMAGE" \
      python /app/alphafold3/run_alphafold.py --json_path="$D" --output_dir="$AF3_DIR/out" \
      --model_dir="$AF3_WEIGHTS_DIR" --run_data_pipeline=false
    echo "af3 inference exit $?"
  fi
done
