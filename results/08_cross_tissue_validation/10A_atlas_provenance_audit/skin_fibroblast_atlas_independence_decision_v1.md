# Step 10A atlas provenance and label audit decision

> This is a metadata/provenance gate. It is not a module validation result and does not change the Step 08C2 PIEZO2 decision.

## Official source inventory

- Official source inventory URL: https://raw.githubusercontent.com/haniffalab/skin_fibroblast_atlas/refs/heads/main/data/scrna_seq_data.csv
- Accession-resolved source rows: 30
- Rows overlapping GSE173252: 1
- Rows overlapping PRJNA607098: 1
- Independent candidate accession rows after exclusion: 28

The official inventory includes both GSE173252 and PRJNA607098. Therefore, the integrated atlas cannot be treated as a fully independent validation cohort unless these source datasets are explicitly excluded from the validation layer.

## Label finding

The Nature article uses F8 for fascia-like myofibroblasts, whereas the current Step 08B3 source-specific output uses F7: Fascia-like myofibroblast in the atlas celltype field. The label relationship is recorded as unresolved until the official object definition and source-specific mapping are audited.

## Sample-level provenance finding

The official scrna_seq_data.csv inventory is source-level and does not provide a sample/donor column. The existing obs/_index sample reconstruction applies to PRJNA607098 only and cannot be generalized to every atlas source.

## Gate A

- Gate A status: **PARTIAL**
- PASS requires an independent source after exclusion plus sample/donor provenance and a resolved target-label crosswalk.
- PARTIAL permits source-level exploratory mapping only; it does not upgrade the current evidence grade.

## Next action

1. Restrict any cross-tissue validation to source accessions that are not GSE173252 or PRJNA607098.
2. Audit source-specific sample/donor fields before computing any donor-level result.
3. Resolve the F7/F8 nomenclature before comparing fascia-like states.
4. Do not download the full atlas or run a new module test in this step.

## Output files

- Source map: ./results/08_cross_tissue_validation/10A_atlas_provenance_audit/skin_fibroblast_atlas_source_accession_map_v1.csv
- Label crosswalk: ./results/08_cross_tissue_validation/10A_atlas_provenance_audit/skin_fibroblast_atlas_F6_F7_F8_label_crosswalk_v1.csv
- Sample-field audit: ./results/08_cross_tissue_validation/10A_atlas_provenance_audit/skin_fibroblast_atlas_source_sample_field_audit_v1.csv
