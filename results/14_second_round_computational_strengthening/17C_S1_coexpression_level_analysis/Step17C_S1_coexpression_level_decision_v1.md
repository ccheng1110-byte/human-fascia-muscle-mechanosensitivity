# Step 17C S1 co-expression-level analysis

## Decision

- Step status: **COMPLETED_WITH_EVIDENCE_BOUNDARY**.
- The analysis estimates within-sample co-detection patterns; it does not prove a cell-intrinsic mechanism.
- PRJNA607098 and GSE130973 were analyzed separately and were not pooled into one statistical sample.

## Analysis contract

- Primary unit: sample/subject-level summary.
- Cell-level co-detection: descriptive only.
- Null construction: within-sample permutation preserving gene-level detection margins; library-size stratification was used when a valid field was available.
- Primary pair families: integrin × actomyosin, integrin × mechanosensor, and integrin × hippo candidate genes.
- Descriptive support flag: positive median excess co-detection, at least 75% positive sample summaries, and BH FDR ≤ 0.05 within source/pair-family. This is not an evidence-grade upgrade rule.

## Evidence boundary

- `WITHIN_SAMPLE_COEXPRESSION_SIGNAL` must be written as co-expression-level evidence only.
- It cannot be described as cell-intrinsic causality, direct mechanosensitivity, or fascia specificity.
- Detection-rate, library-depth, and cell-state composition limitations remain part of the interpretation.
- The project evidence grade remains CAUTION.

- PRJNA607098 source cells: 34569; target-state cells: 30803.
- PRJNA607098 samples analyzed: 12.
- GSE130973 subjects analyzed: 5.
- GSE130973 missing targeted genes: .

## Material Passport

- New large source file downloaded: none; only targeted Zarr chunks were streamed.
- PRJNA607098 output includes the targeted chunk manifest and gene-coverage inventory.
- Next step: Step 17D GSE338388 provenance/design audit for TGFβ exposure × TEAD inhibition regulatory-axis cross-validation.
