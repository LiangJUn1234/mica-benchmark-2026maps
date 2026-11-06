#!/usr/bin/env python3
import csv, gzip, hashlib, shutil, sys, tempfile
from pathlib import Path
import mrcfile, numpy as np

if len(sys.argv) != 4:
    sys.exit("usage: normalize_header.py <list.tsv> <inputs_dir> <out_dir>")
lst, inp, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
assert not (out / "offsets.tsv").exists(), "refusing to overwrite"
rows = []
for r in csv.DictReader(open(lst), delimiter="\t"):
    n, p = r["emdb"].split("-")[1], r["pdb"]; src, dst = inp / p / f"emd_{n}.map.gz", out / p / f"emd_{n}.map.gz"
    assert not dst.exists(), dst; dst.parent.mkdir(parents=True)
    with tempfile.NamedTemporaryFile(suffix=".map") as t, mrcfile.open(src, permissive=True) as m:
        h, d, v = m.header, np.asarray(m.data), m.voxel_size
        # 非标准 header：MAPC/MAPR/MAPS（数组的哪一维对应 x、y、z）不是 1,2,3，或者 NSTART、ORIGIN 不是 0。
        # 有的方法读这种 header 会把原点读错，所以这里改成标准写法；voxel 的数值一个都不改，只是换轴的顺序。
        crs = [int(h.mapc), int(h.mapr), int(h.maps)]
        start = dict(zip(crs, [int(h.nxstart), int(h.nystart), int(h.nzstart)]))
        vox = np.array([float(v.x), float(v.y), float(v.z)])
        # 原来的 map 在空间里的起点 = ORIGIN + NSTART × voxel size（Å）
        off = np.array([float(h.origin.x), float(h.origin.y), float(h.origin.z)]) + np.array([start[1], start[2], start[3]]) * vox
        zyx = np.transpose(d, [[crs[2], crs[1], crs[0]].index(a) for a in (3, 2, 1)])  # 把数组的轴排成 z, y, x
        with mrcfile.new(t.name, overwrite=True) as o:
            o.set_data(np.ascontiguousarray(zyx).astype(d.dtype)); o.voxel_size = tuple(vox); o.header.origin = (0.0, 0.0, 0.0)
        with open(t.name, "rb") as f, gzip.open(dst, "wb") as g:
            shutil.copyfileobj(f, g)
    shutil.copy2(inp / p / f"{p}.fasta", dst.parent)
    # 新 map 的起点移到了 0，所以 reference 要平移 ref_shift = −off 才能和新 map 对齐；这是直接算出来的，没有拟合。
    rows.append([r["emdb"], p, ",".join(map(str, crs)), ",".join(str(start[a]) for a in (1, 2, 3)), *np.round(-off, 4),
                 hashlib.sha256(src.read_bytes()).hexdigest(), hashlib.sha256(dst.read_bytes()).hexdigest()])
(out / "offsets.tsv").write_text("emdb\tpdb\taxes_in\tnstart_xyz\tref_shift_x\tref_shift_y\tref_shift_z\tsrc_sha256\tout_sha256\n"
                                 + "".join("\t".join(map(str, x)) + "\n" for x in rows))
print(*("\t".join(map(str, x[:7])) for x in rows), sep="\n")
