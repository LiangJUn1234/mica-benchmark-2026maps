#!/usr/bin/env python3
import json, sys
from pathlib import Path
import gemmi


def shifted(src, dst, v):
    """每个 model 的每个原子都加上 v；写文件用 gemmi，和 score_model.py 写平移后的 reference 用的是同一个 writer。"""
    st = gemmi.read_structure(str(src)); g = gemmi.Position(*v)
    [setattr(at, "pos", at.pos + g) for md in st for ch in md for r in ch for at in r]
    dst.parent.mkdir(parents=True, exist_ok=True); assert not dst.exists(), dst
    st.make_mmcif_document().write_file(str(dst)); return dst


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: modelangelo_apply_offset.py <offset.json> <model> <out.cif>")
    # offset.json = modelangelo_pad_offset.py 对 ModelAngelo 实际用的那张 map 输出的一行
    off = [json.loads(l) for l in Path(sys.argv[1]).read_text().splitlines() if l.startswith("{")]
    if len(off) != 1 or "shift" not in off[0]:
        sys.exit(f"{sys.argv[1]}: expected one line of modelangelo_pad_offset.py output with a 'shift' field")
    print(shifted(sys.argv[2], Path(sys.argv[3]), off[0]["shift"]))
    # 给移回去的 model 打分之前，要在它旁边手写一个 DONE（score_model.py 从 DONE 读 model 路径），例如
    #   printf 'model %s\n' "$PADFIX_DIR/ma/EMD-47284/model_padfix.cif" > "$PADFIX_DIR/ma/EMD-47284/DONE"
    #   python3 scoring/score_model.py "$PADFIX_DIR" ma EMD-47284 9DXD 0 "$INPUTS_HDRNORM_DIR/offsets.tsv"
