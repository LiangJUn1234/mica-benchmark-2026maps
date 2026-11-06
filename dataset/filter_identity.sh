#!/bin/bash
source "$(dirname "$0")/../config/paths.sh"
set -euo pipefail
[ $# -eq 3 ] || { echo "usage: $0 <candidates.fasta> <training.fasta> <out.m8>" >&2; exit 2; }
query=$1; target=$2; out=$3; tmp=${out%.m8}_tmp   # ${out%.m8} = 去掉结尾的 .m8；MMseqs2 的临时文件放在 <out>_tmp/
mkdir -p "$(dirname "$out")"

# 这里只出 hit 表，不做判断；「identity 超过 25% 就算和 training set 相近」这条规则在 sample.py 里用。
# -s 7.5 是高灵敏度，-e 1e-3 是 E-value 上限；其余都是 MMseqs2 默认（不设 coverage 或 identity 门槛）。
mmseqs easy-search "$query" "$target" "$out" "$tmp" -s 7.5 --threads 8 \
  --format-output query,target,fident,alnlen,qcov,tcov,evalue -e 1e-3
