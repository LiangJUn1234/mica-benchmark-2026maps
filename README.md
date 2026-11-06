# MICA benchmark on 2026 cryo-EM maps

Code for the benchmark in the accompanying report. Test dataset on Zenodo: [DOI to be added].

| In the report | Code |
| --- | --- |
| Test-set selection (Methods, Fig. 1a) | [dataset/](dataset) |
| AF3 models given to the methods (Table 1, note a) | [af3_json.py](inputs/af3_json.py), [run_af3.sh](inputs/run_af3.sh) |
| Header normalization; ModelAngelo padding (Procedure) | [normalize_header.py](inputs/normalize_header.py); [modelangelo_pad_offset.py](methods/modelangelo_pad_offset.py), [modelangelo_apply_offset.py](methods/modelangelo_apply_offset.py) |
| Running each method and the AF3-docked baseline (Procedure) | [methods/](methods) |
| Six metrics; whole-target Cα match (Evaluation Metrics) | [evaluate.py](scoring/evaluate.py); [whole_target_camatch.py](scoring/whole_target_camatch.py) |
| Calibration against Table S3 (Fig. 1b) | [calibrate_test2025.py](scoring/calibrate_test2025.py) |
| Mean scores and paired tests at 2–4 Å and 4–6 Å (Fig. 2, Fig. 3) | [collect_scores.py](analysis/collect_scores.py), [paired_tests.py](analysis/paired_tests.py) |
| CryoACE on its 22 maps (Fig. 3a) | [cryoace_compare.py](analysis/cryoace_compare.py) |
| Medians, training-similar means, results by number of copies, partial copies, numbers in the figure captions (Results, Fig. 1–4, Discussion) | [report_numbers.py](analysis/report_numbers.py) |
| EMD-56389 with every copy given to AF3 (Results, Fig. 5) | [af3_all_copies.py](inputs/af3_all_copies.py), then [mica_dock.sh](methods/mica_dock.sh), [mica_build.sh](methods/mica_build.sh); numbers from [report_numbers.py](analysis/report_numbers.py) |

All files, by folder (click a folder to open it):


<details><summary><b>config/</b> paths</summary>

| File | What it does |
| --- | --- |
| [paths.example.sh](config/paths.example.sh) | Every path the scripts use. |

</details>

<details><summary><b>dataset/</b> test-set selection</summary>

| File | What it does |
| --- | --- |
| [query_candidates.py](dataset/query_candidates.py) | 2026 cryo-EM PDB entries at ≤ 6 Å from RCSB and EMDB, with the reason each excluded entry fails. |
| [extract_sequences.py](dataset/extract_sequences.py) | Protein sequences of a list of PDB entries, used for the training-set FASTA. |
| [training_union_pdb_ids.tsv](dataset/training_union_pdb_ids.tsv) | PDB IDs in the training sets of the benchmarked methods (9,425), with their source. |
| [training_obsolete_replacements.tsv](dataset/training_obsolete_replacements.tsv) | Entries that replace the 25 obsolete IDs in that list. |
| [filter_identity.sh](dataset/filter_identity.sh) | MMseqs2 search of the candidate sequences against the training sequences. |
| [sample.py](dataset/sample.py) | Identity rule (> 25%), size cap and the seeded sample. |
| [size_cap.py](dataset/size_cap.py) | Protein residues in biological assembly 1, checked against the size cap. |

</details>

<details><summary><b>inputs/</b> map headers, AlphaFold3, MSAs</summary>

| File | What it does |
| --- | --- |
| [af3_json.py](inputs/af3_json.py) | AF3 inputs (one per unique sequence), the template-date check, and placement of AF3 models for the methods. |
| [run_af3.sh](inputs/run_af3.sh) | AlphaFold3 3.0.1. |
| [colabfold_msa.py](inputs/colabfold_msa.py) | ColabFold MSAs given to CryoZeta. |
| [normalize_header.py](inputs/normalize_header.py) | Rewrites non-standard map headers and records the shift for the reference. |
| [af3_all_copies.py](inputs/af3_all_copies.py) | AF3 input with every copy of a chain, for the EMD-56389 rerun (MICA's own script writes AlphaFold Server JSON, which the local AF3 3.0.1 does not read). |

</details>

<details><summary><b>methods/</b> one runner per method; results in &lt;code&gt;/&lt;EMDB&gt;/</summary>

| File | What it does |
| --- | --- |
| [mica_dock.sh](methods/mica_dock.sh) | MICA step 1: docks the AF3 models into the map; the docked model is the AF3-docked baseline (code af3d). |
| [mica_build.sh](methods/mica_build.sh) | MICA step 2: model building on the docked entry (mica). |
| [modelangelo.sh](methods/modelangelo.sh) | ModelAngelo (ma). |
| [modelangelo_pad_offset.py](methods/modelangelo_pad_offset.py) | ModelAngelo's padding offset of a map. |
| [modelangelo_apply_offset.py](methods/modelangelo_apply_offset.py) | Moves a ModelAngelo model back by that offset. |
| [emodelx_af.sh](methods/emodelx_af.sh) | EModelX(+AF) with our AF3 models as templates (emx). |
| [cryoatom.sh](methods/cryoatom.sh) | CryoAtom (cat). |
| [emprot.sh](methods/emprot.sh) | EMProt (emp). |
| [e3cryofold.sh](methods/e3cryofold.sh) | E3-CryoFold (e3cf). |
| [em3dfold.sh](methods/em3dfold.sh) | EM3DFold (em3d). |
| [cryozeta.sh](methods/cryozeta.sh) | CryoZeta (cz). |
| [cryoace.sh](methods/cryoace.sh) | CryoACE, starting from MICA's docked model (ace). |
| [cryoace_inputs.py](methods/cryoace_inputs.py) | CryoACE's FASTA and initial model from MICA's docked model. |
| [diffmodeler.sh](methods/diffmodeler.sh) | DiffModeler (dm). |
| [envs/](methods/envs) | Environment files (conda exports, pip freezes, the authors' yml) of the methods. |
| [versions.tsv](methods/versions.tsv) | Repository, commit, weights and environment of each method; versions of the tools and Python packages. |

</details>

<details><summary><b>scoring/</b> scores</summary>

| File | What it does |
| --- | --- |
| [evaluate.py](scoring/evaluate.py) | Scores one model against its reference: TM-score (MM-align), Cα match and sequence match (phenix.chain_comparison), and the other metrics of the MICA paper. |
| [score_model.py](scoring/score_model.py) | Scores one method on one map. |
| [whole_target_camatch.py](scoring/whole_target_camatch.py) | Whole-target Cα match next to the default one. |
| [calibrate_test2025.py](scoring/calibrate_test2025.py) | Rescores MICA's test_2025 models and compares with Table S3 of the MICA paper. |
| [table_s3_mica_test2025.tsv](scoring/table_s3_mica_test2025.tsv) | The Table S3 values used by the calibration. |

</details>

<details><summary><b>analysis/</b> tables and tests</summary>

| File | What it does |
| --- | --- |
| [map_sets.json](analysis/map_sets.json) | Map lists (normalized headers: five maps, of which EMD-70884 is not analysed; the 46 maps at 2–4 Å; excluded maps) and the analysis set. |
| [collect_scores.py](analysis/collect_scores.py) | All scores in one per-map table. |
| [paired_tests.py](analysis/paired_tests.py) | MICA against each method: paired difference, bootstrap CI, t-test, Wilcoxon, Holm. |
| [cryoace_compare.py](analysis/cryoace_compare.py) | CryoACE against MICA, CryoZeta and EModelX(+AF) on its 22 maps. |
| [report_numbers.py](analysis/report_numbers.py) | The other numbers in the report: medians, results by number of copies, partial copies, Fig. 4c, the EMD-56389 rerun. |

</details>
