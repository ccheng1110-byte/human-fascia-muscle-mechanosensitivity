# Step 10C frozen-contract validation decision

## Material Passport

- Origin: local Step 10B analysis contract and the validated PRJNA607098 Zarr extraction helper.
- Statistical unit: 12 reconstructed PRJNA607098 sample IDs.
- Target: exact atlas label `F7: Fascia-like myofibroblast`.
- Comparator: within-sample pooled non-F7 cells.
- Gene sets: config/mechanotransduction_module_registry_v2.csv.
- Candidate panel: config/frozen_candidate_panel_v2.csv.
- No raw H5AD download: only targeted Zarr chunks were streamed.

## Results

- Frozen candidate support: 10/15.
- Integrin/focal-adhesion support: 5/6.
- PIEZO2 supported: FALSE.
- Core module gate: TRUE.
- Leave-one-sample-out gate: TRUE.
- Competition robustness gate: FALSE.
- Bounded support gate: FALSE.
- Evidence grade: CAUTION.

## Interpretation boundary

- This is a sample-level observational validation, not a causal test.
- Step 10A Gate A remains PARTIAL because the integrated atlas contains GSE173252 and PRJNA607098 and lacks a universal donor field.
- A positive result cannot be called independent replication or direct F8 replication while the F7/F8 crosswalk remains unresolved.
- PIEZO2 is evaluated as a channel-specific branch. A negative result narrows the conclusion but does not erase the integrin–actomyosin program.
- Competitor residual checks are robustness audits, not causal adjustment models.
