# MICA benchmark on 2026 cryo-EM maps

Code for the benchmark in the accompanying report. Test dataset on Zenodo: [DOI to be added].

<table>
<tr><th>In the report</th><th>

Code $\color{red}{\text{(click to see what each file does)}}$

</th></tr>
<tr><td>Test-set selection (Methods, Fig. 1a)</td><td><details><summary>dataset/ (7 files)</summary>
<a href="dataset/query_candidates.py">query_candidates.py</a>: 2026 cryo-EM PDB entries at ≤ 6 Å from RCSB and EMDB, with the reason each excluded entry fails.<br>
<a href="dataset/extract_sequences.py">extract_sequences.py</a>: Protein sequences of a list of PDB entries, used for the training-set FASTA.<br>
<a href="dataset/training_union_pdb_ids.tsv">training_union_pdb_ids.tsv</a>: PDB IDs in the training sets of the benchmarked methods (9,425), with their source.<br>
<a href="dataset/training_obsolete_replacements.tsv">training_obsolete_replacements.tsv</a>: Entries that replace the 25 obsolete IDs in that list.<br>
<a href="dataset/filter_identity.sh">filter_identity.sh</a>: MMseqs2 search of the candidate sequences against the training sequences.<br>
<a href="dataset/sample.py">sample.py</a>: Identity rule (> 25%), size cap and the seeded sample.<br>
<a href="dataset/size_cap.py">size_cap.py</a>: Protein residues in biological assembly 1, checked against the size cap.
</details></td></tr>
<tr><td>AF3 models given to the methods (Table 1, note a)</td><td><details><summary>af3_json.py, run_af3.sh</summary>
<a href="inputs/af3_json.py">af3_json.py</a>: AF3 inputs (one per unique sequence), the template-date check, and placement of AF3 models for the methods.<br>
<a href="inputs/run_af3.sh">run_af3.sh</a>: AlphaFold3 3.0.1: data pipeline, then inference.
</details></td></tr>
<tr><td>Header normalization; ModelAngelo padding (Procedure)</td><td><details><summary>normalize_header.py, modelangelo_pad_offset.py, modelangelo_apply_offset.py</summary>
<a href="inputs/normalize_header.py">normalize_header.py</a>: Rewrites non-standard map headers and records the shift for the reference.<br>
<a href="methods/modelangelo_pad_offset.py">modelangelo_pad_offset.py</a>: ModelAngelo's padding offset of a map.<br>
<a href="methods/modelangelo_apply_offset.py">modelangelo_apply_offset.py</a>: Moves a ModelAngelo model back by that offset.
</details></td></tr>
<tr><td>Running each method and the AF3-docked baseline (Procedure, Table 1)</td><td><details><summary>method runners, environments, versions</summary>
In brackets: the folder that holds each method's results.<br>
<a href="methods/mica_dock.sh">mica_dock.sh</a>: MICA step 1: docks the AF3 models into the map; the docked model is the AF3-docked baseline (af3d).<br>
<a href="methods/mica_build.sh">mica_build.sh</a>: MICA step 2: model building on the docked entry (mica).<br>
<a href="methods/modelangelo.sh">modelangelo.sh</a>: ModelAngelo (ma).<br>
<a href="methods/emodelx_af.sh">emodelx_af.sh</a>: EModelX(+AF) with our AF3 models as templates (emx).<br>
<a href="methods/cryoatom.sh">cryoatom.sh</a>: CryoAtom (cat).<br>
<a href="methods/emprot.sh">emprot.sh</a>: EMProt (emp).<br>
<a href="methods/e3cryofold.sh">e3cryofold.sh</a>: E3-CryoFold (e3cf).<br>
<a href="methods/em3dfold.sh">em3dfold.sh</a>: EM3DFold (em3d).<br>
<a href="methods/cryozeta.sh">cryozeta.sh</a>: CryoZeta (cz).<br>
<a href="inputs/colabfold_msa.py">colabfold_msa.py</a>: ColabFold MSAs given to CryoZeta.<br>
<a href="methods/cryoace.sh">cryoace.sh</a>: CryoACE, starting from MICA's docked model (ace).<br>
<a href="methods/cryoace_inputs.py">cryoace_inputs.py</a>: CryoACE's FASTA and initial model from MICA's docked model.<br>
<a href="methods/diffmodeler.sh">diffmodeler.sh</a>: DiffModeler (dm).<br>
<a href="methods/envs">envs/</a>: Environment files (conda exports, pip freezes, the authors' yml) of the methods.<br>
<a href="methods/versions.tsv">versions.tsv</a>: Repository, commit, weights and environment of each method; versions of the tools and Python packages.
</details></td></tr>
<tr><td>Six metrics; whole-target Cα match (Evaluation Metrics)</td><td><details><summary>evaluate.py, score_model.py, whole_target_camatch.py</summary>
<a href="scoring/evaluate.py">evaluate.py</a>: Scores one model against its reference: TM-score (MM-align), Cα match and sequence match (phenix.chain_comparison), and the other metrics of the MICA paper.<br>
<a href="scoring/score_model.py">score_model.py</a>: Scores one method on one map.<br>
<a href="scoring/whole_target_camatch.py">whole_target_camatch.py</a>: Whole-target Cα match next to the default one.
</details></td></tr>
<tr><td>Calibration against Table S3 (Fig. 1b)</td><td><details><summary>calibrate_test2025.py, table_s3_mica_test2025.tsv</summary>
<a href="scoring/calibrate_test2025.py">calibrate_test2025.py</a>: Rescores MICA's test_2025 models and compares with Table S3 of the MICA paper.<br>
<a href="scoring/table_s3_mica_test2025.tsv">table_s3_mica_test2025.tsv</a>: The Table S3 values used by the calibration.
</details></td></tr>
<tr><td>Mean scores and paired tests at 2–4 Å and 4–6 Å (Fig. 2, Fig. 3)</td><td><details><summary>map_sets.json, collect_scores.py, paired_tests.py</summary>
<a href="analysis/map_sets.json">map_sets.json</a>: Map lists (normalized headers: five maps, of which EMD-70884 is not analysed; the 46 maps at 2–4 Å; excluded maps) and the analysis set.<br>
<a href="analysis/collect_scores.py">collect_scores.py</a>: All scores in one per-map table.<br>
<a href="analysis/paired_tests.py">paired_tests.py</a>: MICA against each method: paired difference, bootstrap CI, t-test, Wilcoxon, Holm.
</details></td></tr>
<tr><td>CryoACE on its 22 maps (Fig. 3a)</td><td><details><summary>cryoace_compare.py</summary>
<a href="analysis/cryoace_compare.py">cryoace_compare.py</a>: CryoACE against MICA, CryoZeta and EModelX(+AF) on its 22 maps.
</details></td></tr>
<tr><td>Medians, training-similar means, results by number of copies, partial copies, numbers in the figure captions (Results, Fig. 1–4, Discussion)</td><td><details><summary>report_numbers.py</summary>
<a href="analysis/report_numbers.py">report_numbers.py</a>: The other numbers in the report: medians, training-similar means, results by number of copies, partial copies, numbers in the figure captions, the EMD-56389 rerun.
</details></td></tr>
<tr><td>EMD-56389 with every copy given to AF3 (Results, Fig. 5)</td><td><details><summary>af3_all_copies.py</summary>
<a href="inputs/af3_all_copies.py">af3_all_copies.py</a>: AF3 input with every copy of a chain, for the EMD-56389 rerun (MICA's own script writes AlphaFold Server JSON, which the local AF3 3.0.1 does not read).<br>
Then <a href="methods/mica_dock.sh">mica_dock.sh</a> and <a href="methods/mica_build.sh">mica_build.sh</a>; numbers from <a href="analysis/report_numbers.py">report_numbers.py</a>.
</details></td></tr>
<tr><td>Paths used by every script</td><td><details><summary>paths.example.sh</summary>
<a href="config/paths.example.sh">paths.example.sh</a>: Every path the scripts use.
</details></td></tr>
</table>
