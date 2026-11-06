# MICA benchmark on 2026 cryo-EM maps

Code for the benchmark in the accompanying report. The test dataset is on Zenodo: [DOI to be added].

## Repository structure

<pre>
<a href="environment.yml">environment.yml</a>
<a href="config">config/</a>
└── <a href="config/paths.example.sh">paths.example.sh</a>
<a href="dataset">dataset/</a>
├── <a href="dataset/query_candidates.py">query_candidates.py</a>
├── <a href="dataset/extract_sequences.py">extract_sequences.py</a>
├── <a href="dataset/training_union_pdb_ids.tsv">training_union_pdb_ids.tsv</a>
├── <a href="dataset/training_obsolete_replacements.tsv">training_obsolete_replacements.tsv</a>
├── <a href="dataset/filter_identity.sh">filter_identity.sh</a>
├── <a href="dataset/sample.py">sample.py</a>
└── <a href="dataset/size_cap.py">size_cap.py</a>
<a href="inputs">inputs/</a>
├── <a href="inputs/af3_json.py">af3_json.py</a>
├── <a href="inputs/run_af3.sh">run_af3.sh</a>
├── <a href="inputs/colabfold_msa.py">colabfold_msa.py</a>
├── <a href="inputs/normalize_header.py">normalize_header.py</a>
└── <a href="inputs/af3_all_copies.py">af3_all_copies.py</a>
<a href="methods">methods/</a>
├── <a href="methods/mica_dock.sh">mica_dock.sh</a>
├── <a href="methods/mica_build.sh">mica_build.sh</a>
├── <a href="methods/modelangelo.sh">modelangelo.sh</a>
├── <a href="methods/modelangelo_pad_offset.py">modelangelo_pad_offset.py</a>
├── <a href="methods/modelangelo_apply_offset.py">modelangelo_apply_offset.py</a>
├── <a href="methods/emodelx_af.sh">emodelx_af.sh</a>
├── <a href="methods/cryoatom.sh">cryoatom.sh</a>
├── <a href="methods/emprot.sh">emprot.sh</a>
├── <a href="methods/e3cryofold.sh">e3cryofold.sh</a>
├── <a href="methods/em3dfold.sh">em3dfold.sh</a>
├── <a href="methods/cryozeta.sh">cryozeta.sh</a>
├── <a href="methods/cryoace.sh">cryoace.sh</a>
├── <a href="methods/cryoace_inputs.py">cryoace_inputs.py</a>
├── <a href="methods/diffmodeler.sh">diffmodeler.sh</a>
├── <a href="methods/envs">envs/</a>
└── <a href="methods/versions.tsv">versions.tsv</a>
<a href="scoring">scoring/</a>
├── <a href="scoring/evaluate.py">evaluate.py</a>
├── <a href="scoring/score_model.py">score_model.py</a>
├── <a href="scoring/whole_target_camatch.py">whole_target_camatch.py</a>
├── <a href="scoring/calibrate_test2025.py">calibrate_test2025.py</a>
└── <a href="scoring/table_s3_mica_test2025.tsv">table_s3_mica_test2025.tsv</a>
<a href="analysis">analysis/</a>
├── <a href="analysis/map_sets.json">map_sets.json</a>
├── <a href="analysis/collect_scores.py">collect_scores.py</a>
├── <a href="analysis/paired_tests.py">paired_tests.py</a>
├── <a href="analysis/cryoace_compare.py">cryoace_compare.py</a>
└── <a href="analysis/report_numbers.py">report_numbers.py</a>
</pre>

## Files

**environment.yml**: Python environment for the dataset, input, scoring and analysis scripts.

### config/

- **paths.example.sh**: every path the scripts use. Copy it to `config/paths.sh`, fill in the local tools, environments and weights, and `source` it before running any script. Keep the word "results" out of all paths: MICA rewrites it in its working paths.

### dataset/

- **query_candidates.py**: 2026 cryo-EM PDB entries at ≤ 6 Å from the RCSB and EMDB APIs, with the reason each excluded entry fails.
- **extract_sequences.py**: protein sequences of a list of PDB entries, used to build the training-set FASTA.
- **training_union_pdb_ids.tsv**: union of the training sets of the benchmarked methods (9,425 PDB IDs), with the source of each.
- **training_obsolete_replacements.tsv**: entries that replace the 25 obsolete IDs in the union; their sequences are used in place of the obsolete ones.
- **filter_identity.sh**: MMseqs2 search of the candidate sequences against the training sequences.
- **sample.py**: applies the identity rule (> 25%) and the size cap, draws the seeded sample and writes `manifest.tsv`.
- **size_cap.py**: protein residues in biological assembly 1 of every eligible entry, checked against the size cap.

### inputs/

- **af3_json.py**: `json` writes one AlphaFold3 input per unique sequence; `leak` checks that every template is dated on or before 2021-09-30; `layout` places the AF3 models where MICA and EModelX(+AF) read them.
- **run_af3.sh**: AlphaFold3 3.0.1 on a list of input JSONs, data pipeline then inference.
- **colabfold_msa.py**: MSAs from the ColabFold server, given to CryoZeta.
- **normalize_header.py**: rewrites non-standard map headers without changing voxel values and records the shift for the reference in `offsets.tsv`.
- **af3_all_copies.py**: AF3 input with every copy of a chain, for the EMD-56389 rerun. MICA's `utils/fasta_to_AF3_json.py` writes AlphaFold Server JSON, which the local AF3 3.0.1 does not read; this script gives the single-copy MSA and templates all the chain IDs instead. The rerun:

```bash
J=s904744fa97cdceb0; A=$BENCH_ROOT/af3_allcopies; O=$BENCH_ROOT/outputs_allcopies
mkdir -p "$A" && (cd "$A" && python3 "$TOOLS_DIR/MICA/utils/fasta_to_AF3_json.py" -f "$INPUTS_DIR/9TX6/9TX6.fasta" -n 56389)
python3 inputs/af3_all_copies.py "$A/input/56389/AF3_JSON/9TX6_1.json" \
    "$AF3_DIR/msa/$J/${J}_data.json" "$A/msa/$J/${J}_data.json"
echo "$A/json/$J.json" > allcopies_jobs.txt && AF3_DIR=$A bash inputs/run_af3.sh allcopies_jobs.txt
AF3_DIR=$A RESULTS_DIR=$O bash methods/mica_dock.sh EMD-56389 9TX6 3.9 0.1
AF3_DIR=$A RESULTS_DIR=$O bash methods/mica_build.sh EMD-56389 9TX6 3.9
for m in af3d mica; do python3 scoring/score_model.py "$O" $m EMD-56389 9TX6 1; done
```

### methods/

Each runner takes the map's EMDB and PDB IDs, plus resolution and contour where the method needs them, runs the method in a work directory under `$SCRATCH_DIR`, and writes `$RESULTS_DIR/<code>/<EMDB>/` with `DONE` (model path) or `FAILED` (reason and log tail).

| Script | Method | Code | Arguments |
| --- | --- | --- | --- |
| mica_dock.sh | AF3-docked (MICA step 1) | af3d | EMDB PDB resolution contour |
| mica_build.sh | MICA (needs mica_dock.sh) | mica | EMDB PDB resolution |
| modelangelo.sh | ModelAngelo | ma | EMDB PDB |
| emodelx_af.sh | EModelX(+AF) | emx | EMDB PDB resolution |
| cryoatom.sh | CryoAtom | cat | EMDB PDB |
| emprot.sh | EMProt | emp | EMDB PDB |
| e3cryofold.sh | E3-CryoFold | e3cf | EMDB PDB |
| em3dfold.sh | EM3DFold | em3d | EMDB PDB |
| cryozeta.sh | CryoZeta | cz | EMDB PDB resolution contour |
| cryoace.sh | CryoACE (needs mica_dock.sh) | ace | EMDB PDB resolution contour |
| diffmodeler.sh | DiffModeler | dm | EMDB PDB resolution contour |

- **modelangelo_pad_offset.py**: the padding offset of a map, computed with ModelAngelo's own code (run in ModelAngelo's environment).
- **modelangelo_apply_offset.py**: moves a ModelAngelo model back by that offset.
- **cryoace_inputs.py**: CryoACE's FASTA and initial model from MICA's docked AF3 model.
- **envs/**: conda exports and pip freezes of the method environments.
- **versions.tsv**: repository, commit, weights and environment file of each method.

### scoring/

- **evaluate.py**: scores one model against its reference: TM-score with MM-align, Cα match and sequence match with phenix.chain_comparison, and the other metrics of the MICA paper.
- **score_model.py**: scores one method on one map and writes `metrics.tsv` next to its `DONE`.
- **whole_target_camatch.py**: Cα match with `test_unique_part_of_target_only=False`, next to the default.
- **calibrate_test2025.py**: rescores MICA's published test_2025 models and compares with Table S3 of the MICA paper.
- **table_s3_mica_test2025.tsv**: the Table S3 values used by the calibration.

### analysis/

- **map_sets.json**: maps with normalized headers, the 46 main maps, the excluded maps, and the analysis set.
- **collect_scores.py**: collects every `metrics.tsv` into one per-map table.
- **paired_tests.py**: MICA against each method: mean difference with bootstrap CI, paired t-test and Wilcoxon test, Holm-corrected.
- **cryoace_compare.py**: CryoACE against MICA, CryoZeta and EModelX(+AF) on the 22 main maps within CryoACE's size limit.
- **report_numbers.py**: the numbers in the report that the scripts above do not print: medians and maps above 0.9 per method, results by number of copies, maps where MICA built only some copies, the central-ring count of Fig. 4c, and the EMD-56389 rerun with every copy.

## Requirements

Python 3.11 ([environment.yml](environment.yml)); Phenix 1.20.1-4487, MM-align and US-align 20261019, MMseqs2 18.8cc5c, AlphaFold3 3.0.1. Method environments are in [methods/envs/](methods/envs) and versions in [methods/versions.tsv](methods/versions.tsv). Most method runs used one NVIDIA A100 GPU, 8–32 CPU cores, 64–128 GB memory and a 4 h time limit (more for large entries); MICA's docking step used 48 CPU cores, 24 GB and no GPU, and AF3 used one A100 for inference and no GPU for the data pipeline.
