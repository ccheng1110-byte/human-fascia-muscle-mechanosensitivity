# Step 17D4 GSE338388 processed-matrix sample-column re-audit

## Decision

- Gate 17D4: **PASS_TO_S2_EXPRESSION_ANALYSIS**.
- This step audited the official processed matrix and performed no expression hypothesis test.
- No FASTQ, RAW archive, H5AD, or unprocessed sequencing file was downloaded.

## Matrix and design checks

- Matrix file: GSE338388_Cnt-K975-TGF-TGFandK975-12h_normalizedExpression.csv.gz (existing).
- Matrix dimensions: 19791 rows x 16 columns.
- Resolved gene identifier column: Name.
- Resolved Yap164 sample columns: 12.
- Sample-column mapping gate: TRUE.
- Complete TGF-beta x TEAD 2x2 matrix design: TRUE.
- Three matrix replicates in each combination: TRUE.
- Frozen candidate coverage: 15/15.
- Registry coverage: 178/181.
- Core-module coverage gate (>=80% per core module): TRUE.

## S2 interpretation boundary

- A passing audit authorizes the next expression analysis only.
- The next analysis must report TGF-beta main effect, TEAD-inhibition main effect, and interaction separately.
- Because GSE338388 contains no mechanical loading or stiffness factor, it remains regulatory-axis cross-validation rather than mechanical causality.
- Missing genes are not replaced and the overall project evidence grade remains CAUTION.

## Material Passport

- Source: official NCBI GEO GSE338388 processed expression matrix.
- Transformation: matrix parsing, sample-column mapping, frozen candidate/module coverage audit.
- New source type downloaded: one processed matrix candidate only.
- Next step: S2 frozen-contract two-factor expression analysis, if the gate is PASS.
