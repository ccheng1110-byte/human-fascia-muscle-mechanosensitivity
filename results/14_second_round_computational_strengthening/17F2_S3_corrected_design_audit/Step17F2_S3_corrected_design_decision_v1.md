# Step 17F2 S3 corrected design audit

## GSE123100

- Corrected eligible primary samples: 11/19.
- Excluded clinical tissue samples: 8.
- Corrected stiffness levels: 1.1, 2.5, 4.2, 11.9, 34.3, 34.4, 50 kPa.
- Corrected gate: **PARTIAL_PROCEED_TO_PRIMARY_MATRIX_AUDIT_WITH_DESCRIPTIVE_DOSE_LIMIT**.
- The HTM series can support a cautious dose-response-form analysis, but most stiffness-by-culture-group cells contain one sample; formal inferential claims must remain limited.
- The eight clinical trabecular-meshwork tissue samples are excluded from the primary stiffness series and must not be mixed with cultured HTM samples.

## GSE276045

- Corrected eligible samples: 178/178.
- Cell-model conditions: hTERT, WT.
- Corrected stiffness levels retained: 0.5, 2, 16, 32, 2300000 kPa-equivalent.
- Corrected gate: **PARTIAL_PROCEED_TO_WI38_CELL_MODEL_MATRIX_AUDIT**.
- WT/hTERT is retained as a cell-model/genotype factor, not relabeled as proliferation status.
- The planned stiffness-by-proliferation confounding test is therefore NOT_ESTIMABLE from current metadata; the next analysis may only test stiffness direction by cell-model condition and timepoint, with this limitation explicit.

## S3 boundary

- Neither source is fascia-direct validation.
- No GWAS branch is activated by this audit.
- A passing targeted matrix audit will authorize descriptive or model-based cross-tissue checks only; it will not change the overall project evidence grade from CAUTION.

## Material Passport

- Input: Step 17F official GEO metadata and supplementary inventories.
- Transformation: title/characteristic-supported correction of cultured versus clinical samples, stiffness values, cell-model conditions, timepoints and replicate labels.
- New data downloaded: none.
- Next step: targeted audit of the selected official processed matrices using these corrected sample maps.
