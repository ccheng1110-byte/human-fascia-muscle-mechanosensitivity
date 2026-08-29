# Step 17D2 corrected GSE338388 design audit

## Decision

- Gate 17D2: **PASS_TO_S2_MATRIX_AUDIT**.
- The preliminary automatic parser was corrected after manual review.
- No expression matrix was downloaded and no hypothesis test was performed.

## Why the preliminary parser failed

- GEO's generic treatment protocol field states that the experiment used both TGF-beta1 and K-975 for every sample.
- That field describes the protocol, not the sample-specific treatment assignment.
- The valid sample-specific fields are the title and the characteristic `treatment`.

## Corrected design

- Samples: 12.
- TGF-beta levels: exposed, not_exposed.
- TEAD levels: tead_inhibited, tead_not_inhibited.
- Complete 2×2 design: TRUE.
- Three replicates in each combination: TRUE.
- Control = no TGF-beta exposure and no TEAD inhibition.
- K-975 = TEAD inhibition without TGF-beta exposure.
- TGF-beta1 = TGF-beta exposure without TEAD inhibition.
- TGF-beta1 + K-975 = both factors present.

## S2 interpretation boundary

- GSE338388 has no mechanical loading or stiffness factor.
- A passing design gate supports TGF-beta/SMAD versus TEAD-related regulatory-axis cross-validation only.
- It cannot prove mechanical causality, YAP/TAZ-independent mechanical driving, or fascia-direct replication.
- Frozen candidate and module coverage must be checked in the processed matrix before any expression analysis.
- Overall project evidence grade remains CAUTION.

## Material Passport

- Input: official GSE338388 sample metadata snapshot from Step 17D.
- Transformation: manual-review-supported correction using sample-specific title and treatment characteristic.
- New data downloaded: none.
- Next step: Step 17E S2 processed-matrix audit and frozen panel/module coverage check.
