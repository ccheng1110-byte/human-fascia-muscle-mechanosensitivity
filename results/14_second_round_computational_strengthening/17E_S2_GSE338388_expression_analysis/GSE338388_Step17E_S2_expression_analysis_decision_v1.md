# Step 17E GSE338388 S2 frozen-contract expression analysis

- Analysis gate: **PASS_TO_S2_EVIDENCE_SYNTHESIS**.
- Interpretation class: **REGULATORY_AXIS_CONCORDANT**.
- Overall project evidence grade remains: **CAUTION**.
- This analysis is regulatory-axis cross-validation, not mechanical causality.

## Prespecified contrasts

- `tgfb_main_effect`: marginal TGF-beta exposure effect across TEAD levels.
- `tead_main_effect`: marginal TEAD-inhibition effect across TGF-beta levels.
- `tgfb_x_tead_interaction`: difference in the TGF-beta effect between TEAD states.

## Axis-focused interpretation

- Hippo/YAP/TAZ module under TEAD inhibition: Down; FDR=0.001461.
- TGF/fibrosis module under TGF-beta exposure: Up; FDR=0.02584.
- Expected-direction concordance: TRUE.
- Expected-direction FDR support: TRUE.

## Interpretation boundary

- A positive result supports TGF-beta/SMAD and TEAD-related regulatory-axis cross-validation within this vocal-fold fibroblast dataset.
- It does not establish mechanical loading or stiffness causality, YAP/TAZ-independent mechanical driving, fascia specificity, or pain causality.
- Candidate-level results are exploratory within the frozen 15-gene contract; module-level results are the primary S2 evidence.
- Null or discordant results are retained as boundary evidence and do not justify changing the frozen panel or relaxing thresholds.

## Material Passport

- Input: official GSE338388 normalized expression matrix already downloaded during Step 17D3/17D4.
- Provenance/design: corrected official GEO sample design and validated Yap164 suffix mapping from Step 17D4.
- Preprocessing: log2(x+1) for non-negative normalized expression; duplicate gene identifiers collapsed by mean where needed.
- Model: limma linear model with four factorial condition means; no batch term was introduced because the locked 2x2 design has three replicates per cell and no validated independent batch factor.
- Multiplicity: genome-wide BH for all genes, frozen-family BH for candidate genes, and module-level FDR from camera.
- No FASTQ, RAW archive, H5AD, or additional large file was downloaded in Step 17E.
