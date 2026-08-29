# Step 17D GSE338388 provenance and design audit

## Decision

- Gate 17D: **HOLD_NEEDS_MANUAL_DESIGN_REVIEW**.
- This step audited official GEO metadata and supplementary-file inventory only.
- No processed expression matrix, FASTQ, RAW archive, or H5AD was downloaded.
- No biological hypothesis test was performed.

## Locked design requirement

- Required design: TGFβ exposure (not exposed/exposed) × TEAD inhibition (not inhibited/inhibited).
- This dataset contains no mechanical loading or stiffness factor.
- A passing audit therefore supports regulatory-axis cross-validation only, not mechanical causality or fascia-direct replication.

## Reconstructed design

- GEO sample records: 12.
- TGFβ levels detected: ambiguous.
- TEAD levels detected: ambiguous, tead_inhibited.
- Complete 2×2 combination table: FALSE.
- Supplementary files listed: 2.
- Listing error, if any: none.

## S2 boundaries

- The TGFβ main effect, TEAD-inhibition main effect, and interaction must be reported separately in the next step.
- Frozen candidate and module coverage must be audited before expression analysis; missing genes are not replaced.
- A TEAD-associated module response cannot be described as proof that YAP/TAZ independently drives a mechanical program.
- The project evidence grade remains CAUTION regardless of this gate.

## Material Passport

- Source: official NCBI GEO record GSE338388.
- Transformation: sample metadata parsing, factor reconstruction, 2×2 design gate, and supplementary inventory.
- Frozen candidate panel recorded: 15 genes; expression coverage remains pending the next matrix audit.
- Next step: S2 processed-matrix audit and frozen module/candidate coverage check.
