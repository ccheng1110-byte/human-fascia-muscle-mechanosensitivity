## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-23
- Verification Status: ANALYZED
- Version Label: PRJNA607098_F7_celltype_column_v3

## Step 08B3 final PRJNA607098 F7 source-specific validation

### Label-column correction audit

- Step 08B v1 searched for F7 in `celltype_skinspecific_nomenclature`, where PRJNA607098 has no F7 label; v1 is invalid as a biological negative result.
- Step 08B2 tested the auxiliary F6/Myofibroblast population (1,714 cells) in that column and is retained as a sensitivity analysis.
- Step 08B3 uses the intended atlas `celltype` column and the frozen target `F7: Fascia-like myofibroblast`.
- The frozen 15 genes, contrast rule and gate are unchanged.

### Pre-specified gate

Proceed to donor-level Step 08C only if all conditions are met:

- At least 8 of 15 frozen candidates show descriptive support.
- At least 4 of 6 integrin/focal-adhesion candidates show support.
- PIEZO2 shows support.

### Results

- Frozen candidates supported: 13/15.
- Integrin/focal-adhesion candidates supported: 6/6.
- PIEZO2 supported: TRUE.
- Step 08C gate passed: TRUE.
- Known expression-chunk volume: 0 MiB; 18 chunk size(s) unavailable after non-blocking network audit.
- Required files were streamed by the user-run R workflow; no 50 MiB per-file limit was applied.

### Evidence decision

- Evidence grade remains CAUTION even if the gate passes.
- This step isolates PRJNA607098 from GSE173252, improving dataset independence.
- The public atlas lacks sample/donor labels, so cells are not treated as replicates and no p values are calculated.
- A passing gate supports investment in Step 08C donor/sample-level reconstruction.
- Mechanosensor activity, protein abundance and causality remain untested.

### Statistical fallacy scan (11/11 checked)

- Simpson's paradox: not testable without donor strata; retained as a caution.
- Ecological fallacy: avoided; no cell-state result is generalized to individual patients.
- Berkson's paradox: possible surgical/atlas selection bias; caution.
- Collider bias: no covariate-adjusted causal model was fitted.
- Base-rate neglect: cell counts and expression prevalence are reported.
- Regression to the mean: not applicable to this cross-sectional screen.
- Survivorship bias: atlas QC selection may exclude low-quality cells; caution.
- Look-elsewhere effect: reduced by using a frozen 15-gene candidate panel.
- Garden of forking paths: gate and contrasts were fixed before 08B results.
- Correlation is not causation: causal claims are prohibited.
- Reverse causality: state-expression direction cannot establish mechanism.
