#!/usr/bin/env python3
import json, sys, types
import mrcfile, numpy as np
import model_angelo.utils.grid as G

# 要在 ModelAngelo 自己的环境里跑（用到 model_angelo 这个包）。
# ModelAngelo v1.0.18 建模前先把 map 补成立方体（make_cubic：边长 = max(最长边, 128)，取偶数），
# 坐标写在补过的盒子里，所以需要补边的 map 上 model 会整体偏移。这里直接调用 ModelAngelo 自己的
# load_mrc 和 make_cubic 算出每个轴前面补了几格 c，修正量 shift = −c × voxel（Å），由 modelangelo_apply_offset.py 加回去。
if len(sys.argv) < 3:
    sys.exit("usage: modelangelo_pad_offset.py {header|full} <map> [<map> ...]")
# 确认 ModelAngelo 输出坐标的 (x, y, z) 对应数组的第 (2, 1, 0) 维：这行就是它 c_alpha/inference.py 第 40 行的写法
lat = np.flip(G.get_lattice_meshgrid_np(4, no_shift=True), -1)
assert lat[1, 2, 3].tolist() == [3, 2, 1], lat[1, 2, 3]
mode, maps = sys.argv[1], sys.argv[2:]
assert mode in ("header", "full"), mode


def header_only(fn, *a, **k):
    """header 模式：只读 header，数据换成一个不占内存的全 0 数组（形状和真数据一样），大 map 也能很快算完。"""
    with mrcfile.open(fn, header_only=True) as m:
        h, vs = m.header.copy(), m.voxel_size.copy()
    data = np.lib.stride_tricks.as_strided(np.zeros(1, np.float32), shape=(int(h.nz), int(h.ny), int(h.nx)), strides=(0, 0, 0))
    return types.SimpleNamespace(header=h, voxel_size=vs, data=data)


if mode == "header":
    G.mrc = types.SimpleNamespace(open=header_only)  # 只替换 load_mrc 读文件那一步，load_mrc 和 make_cubic 本身不改
for fn in maps:
    with mrcfile.open(fn, header_only=True) as m:
        h = m.header
        hd = dict(nxyz=[int(h.nx), int(h.ny), int(h.nz)], axes=[int(h.mapc), int(h.mapr), int(h.maps)],
                  nstart=[int(h.nxstart), int(h.nystart), int(h.nzstart)], origin=[float(h.origin.x), float(h.origin.y), float(h.origin.z)])
    try:
        o = G.load_mrc(fn)
    except RuntimeError as x:
        print(json.dumps(dict(map=fn, mode=mode, **hd, error=str(x)))); continue
    nbox, c, _ = G.make_cubic(o.grid)   # c = 每个数组轴前面补的格数
    s = int(nbox.shape[0]); del nbox
    c_xyz = np.flip(np.asarray(c), -1)  # 数组轴顺序 (z, y, x) 换成 (x, y, z)
    print(json.dumps(dict(map=fn, mode=mode, **hd, grid_shape=list(o.grid.shape), s=s, c_axes=np.asarray(c).tolist(),
                          c_xyz=c_xyz.tolist(), voxel=o.voxel_size, shift=(-c_xyz * o.voxel_size).tolist())), flush=True)
