#!/usr/bin/env python3
import hashlib, json, re, sys, time
from pathlib import Path
import gemmi


def entities(inputs, pdb):
    """从 RCSB 的 FASTA 读出每种 chain（entity）：返回 (entity id，例如 9ABC_1；序列)。"""
    for b in (Path(inputs) / pdb / f"{pdb}.fasta").read_text().split(">")[1:]:
        head, *seq = b.strip().split("\n")
        yield head.split("|")[0], "".join(seq)


def job(seq):
    # AF3 的作业名用序列的 sha1 来起：同一条序列在不同 map 里只跑一次 AF3。AF3 会把作业名转成小写，所以这里本来就用小写。
    return "s" + hashlib.sha1(seq.encode()).hexdigest()[:16]


if len(sys.argv) < 4 or sys.argv[1] not in ("json", "leak", "layout"):
    sys.exit("usage: af3_json.py json <inputs_dir> <af3_dir> <PDB>...\n"
             "       af3_json.py leak <inputs_dir> <af3_dir>\n"
             "       af3_json.py layout <inputs_dir> <af3_dir> <run_dir> <PDB>")
cmd, inputs, af3 = sys.argv[1], sys.argv[2], Path(sys.argv[3])
if cmd == "json":
    # 每条不重复的序列写一个 AF3 输入，只有一份 copy（id "A"），一个 seed。
    (af3 / "json").mkdir(parents=True, exist_ok=True)
    for pdb in sys.argv[4:]:
        for eid, seq in entities(inputs, pdb):
            p = af3 / "json" / f"{job(seq)}.json"
            if not p.exists():
                p.write_text(json.dumps({"name": job(seq), "modelSeeds": [1], "dialect": "alphafold3", "version": 1,
                                         "sequences": [{"protein": {"id": "A", "sequence": seq}}]}, indent=1))
            print(eid, len(seq), p.name)
elif cmd == "leak":
    # 检查 AF3 用到的每个 template 都是 2021-09-30（AF3 默认的 template 截止日）以前发布的，避免用到测试结构本身。
    # 发布日期取 template mmCIF 里最早的 revision_date；没有日期的也算不通过。
    rows = []
    for f in sorted((af3 / "msa").glob("*/*_data.json")):
        for s in json.load(open(f))["sequences"]:
            for tp in s.get("protein", {}).get("templates") or []:
                d = re.findall(r"_pdbx_audit_revision_history.revision_date\s+(\d{4}-\d\d-\d\d)", tp["mmcif"])
                e = re.search(r"_entry.id\s+(\S+)", tp["mmcif"])
                rows.append((f.parent.name, e[1] if e else "?", min(d) if d else "?"))
    late = [r for r in rows if r[2] == "?" or r[2] > "2021-09-30"]
    out = af3 / f"gate_leak_{time.strftime('%Y%m%d_%H%M%S')}.tsv"
    out.write_text("job\ttemplate\trelease\n" + "".join("\t".join(r) + "\n" for r in rows))
    jobs = len(list((af3 / "msa").glob("*/*_data.json")))
    print(f"{jobs} AF3 jobs, {len(rows)} templates; latest release {max((r[2] for r in rows), default='-')}; "
          f"after 2021-09-30 or undated: {len(late)} -> template date check {'PASS' if not late else 'FAIL'}  ({out.name})", late[:5])
elif cmd == "layout":
    # 把每种 chain 排第一的 AF3 model 按两个方法要的文件名放好：
    #   AF3_results/<PDB>_<n>/fold_<pdb>_<n>_model_0.cif  给 MICA 的 process_AF3_results.py（它从文件名第 3 段读 entity 编号）
    #   templates/<PDB>_<n>/ranked_0.pdb                  给 EModelX 的 --template_dir
    run, pdb = Path(sys.argv[4]), sys.argv[5]
    assert "results" not in str(run), "MICA would rewrite 'results' inside run_dir"  # MICA 会把路径里的 results 改成 PDBs
    for eid, seq in entities(inputs, pdb):
        top = af3 / "out" / job(seq) / f"{job(seq)}_model.cif"
        a, t = run / "AF3_results" / eid, run / "templates" / eid
        a.mkdir(parents=True), t.mkdir(parents=True)  # 文件夹已存在就报错，不会覆盖
        (a / f"fold_{pdb.lower()}_{eid.split('_')[1]}_model_0.cif").write_bytes(top.read_bytes())
        gemmi.read_structure(str(top)).write_pdb(str(t / "ranked_0.pdb"))
        print(eid, job(seq), "->", a.name)
