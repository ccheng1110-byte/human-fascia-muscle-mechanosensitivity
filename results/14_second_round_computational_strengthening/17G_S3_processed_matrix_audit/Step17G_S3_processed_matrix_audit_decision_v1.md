# Step 17G S3 processed-matrix audit

- This step downloaded only the selected official processed matrices through R and performed no biological hypothesis test.
- GSE123100 clinical tissue samples were excluded from the primary HTM stiffness-series mapping.
- GSE276045 WT/hTERT labels remain cell-model conditions; they are not treated as proliferation measurements.

## Source decisions

### GSE123100

- Matrix dimensions: 20298 x 18.
- Eligible metadata samples: 11; mapped eligible samples: 0.
- Resolved gene column: Symbol.
- Frozen candidate coverage: 14/15.
- Core module gate: TRUE.
- Decision: **HOLD_SAMPLE_MAPPING_OR_MATRIX_REVIEW**.

### GSE276045

- Matrix dimensions: 58096 x 179.
- Eligible metadata samples: 178; mapped eligible samples: 0.
- Resolved gene column: "Genes".
- Frozen candidate coverage: 0/15.
- Core module gate: FALSE.
- Decision: **HOLD_SAMPLE_MAPPING_OR_MATRIX_REVIEW**.

## S3 interpretation boundary

- A passing matrix audit authorizes the corresponding targeted expression analysis only.
- GSE123100 may support a cautious cross-tissue stiffness dose-response-form analysis, not fascia-direct validation or causal proof.
- GSE276045 may support a stiffness direction by WT/hTERT cell-model condition and timepoint; the proliferation-confounding question remains NOT_ESTIMABLE.
- Overall project evidence grade remains CAUTION.
