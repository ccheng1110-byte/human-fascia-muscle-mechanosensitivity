## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-24
- Verification Status: ANALYZED
- Version Label: PRJNA607098_F7_sample_level_v1

## Step 08C2 paired sample-level F7 validation

### Frozen analysis

- Statistical units: 12 reconstructed SRS samples; cells are aggregated within sample.
- Primary contrast: within-sample mean expression in F7 minus pooled non-F7 cells.
- Directional test: one-sided one-sample Wilcoxon signed-rank test across 12 paired differences.
- Multiplicity: Benjamini-Hochberg correction across the frozen 15-candidate family.
- Strict gene support requires median difference > 0, at least 10/12 positive samples, and BH FDR <= 0.05.
- PIEZO1, TRPV4 and PKD2 are secondary context genes and do not enter the frozen gate.

### Pre-specified gate

- At least 8/15 frozen candidates must pass strict sample-level support.
- At least 4/6 integrin/focal-adhesion candidates must pass.
- PIEZO2 must pass.

### Results

- Frozen candidates supported: 10/15.
- Integrin/focal-adhesion candidates supported: 5/6.
- PIEZO2 supported: FALSE.
- Step 08C2 gate passed: FALSE.
- Evidence grade: CAUTION.

### Evidence boundary

- Passing upgrades the result from cell-descriptive CAUTION to sample-replicated CAUTION+, not to causal evidence.
- SRS accessions are sample units; donor independence still requires explicit accession-to-donor verification.
- The contrast uses normalized atlas expression, not raw-count pseudobulk, and pooled non-F7 cells are composition-dependent.
- The data are observational and disease-source specific; protein abundance, channel activity, force response and causality remain untested.

### Statistical fallacy scan (11/11 checked)

- Simpson's paradox: reduced by within-sample pairing, but non-F7 state composition remains a sensitivity issue.
- Ecological fallacy: sample-level expression is not generalized to individual clinical outcomes.
- Berkson's paradox: surgical and atlas inclusion may induce selection bias.
- Collider bias: no adjusted causal model is fitted.
- Base-rate neglect: sample and cell counts are retained in the outputs.
- Regression to the mean: no repeated extreme-value selection is used.
- Survivorship bias: atlas QC may exclude low-quality cells or samples.
- Look-elsewhere effect: controlled within the frozen 15-gene family by BH FDR.
- Garden of forking paths: direction, family, support rule and gate were frozen before 08C2 extraction.
- Correlation is not causation: causal claims are prohibited.
- Reverse causality: expression-state association cannot establish mechanosensor activation direction.
