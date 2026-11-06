#!/usr/bin/env python3
import sys
import gemmi

if len(sys.argv) < 7:
    sys.exit("usage: cryoace_inputs.py <rcsb.fasta> <docked.pdb> <A|B> <templates_dir> <out_dir> <name>")
fa, dock, v, tdir, out, name = sys.argv[1:7]
# 每种 chain（entity）：(entity id, 序列, 份数)。份数 = FASTA 头 "Chains A, B, C" 里列了几条 chain。
ents = []
for b in open(fa).read().split(">")[1:]:
    head, *seq = b.strip().split("\n")
    chains = head.split("|")[1].replace("Chains ", "").replace("Chain ", "").split(", ")
    ents.append((head.split("|")[0], "".join(seq), len(chains)))
# CryoACE 的 load_cif 只比较 chain 名的第一个字符，所以 chain id 只能用单个字符，最多 62 个
ids = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
if sum(n for _, _, n in ents) > 62:
    sys.exit("ERROR: more than 62 chain copies")
copies = {e: [[] for _ in range(n)] for e, _, n in ents}
# 把 dock 进来的每一段分给某种 chain 的某一份 copy：从最长的一段开始；
# 按序列判断属于哪种 chain，再放进第一份 residue 编号和已有片段不重叠的 copy。
for ch in sorted(gemmi.read_structure(dock)[0], key=len, reverse=True):
    code = gemmi.one_letter_code([r.name for r in ch]).upper()
    e = next((e for e, s, _ in ents if all(0 < r.seqid.num <= len(s) and s[r.seqid.num - 1] == c
                                           for r, c in zip(ch, code))), None)
    nums = {r.seqid.num for r in ch}
    slot = next((c for c in copies.get(e, []) if not nums & {r.seqid.num for p in c for r in p}), None)
    print(f"piece {ch.name} {min(nums)}-{max(nums)} -> {e}, copy {copies[e].index(slot) + 1 if slot is not None else 'NONE'}")
    if slot is not None:
        slot.append(ch)
st = gemmi.Structure()
st.add_model(gemmi.Model("1"))
recs = []
for e, s, _ in ents:
    for i, pieces in enumerate(copies[e]):
        cid = next(ids)
        recs.append(f">{cid}|protein\n{s}\n")  # Boltz 的 FASTA：每份 copy 一条记录，顺序和 aligned model 的 chain 一样
        if not pieces:
            print(f"copy {e}#{i + 1} ({cid}): no docked piece, absent from aligned model")
            continue
        ch = gemmi.Chain(cid)
        if v == "A":
            # A：这份 copy 的几段原样拼起来（可能不完整）
            for r in sorted((r for p in pieces for r in p), key=lambda r: r.seqid.num):
                ch.add_residue(r)
        else:
            # B：拿这种 chain 的全长 AF3 model，用 Cα 叠到这份 copy 里最长的那一段上，所以放进去的 copy 都是完整的
            big = max(pieces, key=len)
            full = gemmi.read_structure(f"{tdir}/{e}/ranked_0.pdb")[0][0]
            fix = {r.seqid.num: r["CA"][0].pos for r in big if r.find_atom("CA", "*")}
            sup = gemmi.superpose_positions([fix[r.seqid.num] for r in full if r.seqid.num in fix],
                                            [r["CA"][0].pos for r in full if r.seqid.num in fix])
            for r in full:
                ch.add_residue(r)
            for r in ch:
                for a in r:
                    a.pos = gemmi.Position(sup.transform.apply(a.pos))
            print(f"copy {e}#{i + 1} ({cid}): AF3 fitted on piece {big.name}, CA rmsd {sup.rmsd:.2f}")
        for r in ch:
            r.label_seq = r.seqid.num  # CryoACE 读文件时按 label_seq（序列位置）找 residue
        st[0].add_chain(ch)
        print(f"copy {e}#{i + 1} ({cid}): {len(st[0][cid])}/{len(s)} residues in aligned model")
st.setup_entities()
st.make_mmcif_document().write_file(f"{out}/aligned_{v}.cif")
open(f"{out}/{name}.fasta", "w").write("".join(recs))
