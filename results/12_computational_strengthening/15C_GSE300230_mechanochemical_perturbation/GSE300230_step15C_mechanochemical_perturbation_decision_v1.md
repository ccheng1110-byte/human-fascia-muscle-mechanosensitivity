# Step 15C GSE300230 mechanochemical perturbation decision v1

- Perturbation gate: **PARTIAL_DIRECTIONAL_PLAUSIBILITY**.
- Overall manuscript evidence grade remains: **CAUTION**.
- Frozen candidates positive in both cell lines: 8/15.
- Frozen candidates with candidate-family q <= 0.05 in both lines: 0/15.

## Core-module primary contrast

- actomyosin_rho: ROBUST_BOTH_LINES (GM08401 Up, q=5.071e-09; GM09503 Up, q=1.197e-05).
- integrin_focal_adhesion: DIRECTIONALLY_CONCORDANT_ONE_LINE (GM08401 Up, q=0.0002666; GM09503 Up, q=0.1151).

## Interpretation boundary

This analysis tests mechanochemical response within two primary dermal fibroblast cell lines. It does not establish fascia specificity, population-level donor replication, an independent age effect, or causal relevance to human pain. Age is completely confounded with cell line. Null and discordant module/candidate results remain part of the evidence record.

## Material Passport

- Input: official GSE300230 16-sample processed raw-count matrix.
- Contract: step15_GSE300230_analysis_contract_v1.md.
- Model: line-stratified limma-voom with four mechanochemical groups plus batch.
- Multiplicity: BH within nine modules and within the frozen 15-gene candidate family.
- Reproducibility: all design, module, candidate, transcriptome and decision outputs are stored in this result directory.
