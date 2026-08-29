## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-23
- Verification Status: ANALYZED
- Version Label: PRJNA607098_sample_metadata_audit_v1

## Step 08C1 PRJNA607098 sample-metadata audit

### Frozen eligibility rule

- At least 6 samples must each contain at least 20 F7 cells and 20 non-F7 cells.
- Missing sample IDs must account for no more than 1% of source cells.
- No single sample may account for more than 40% of source cells.

### Results

- Resolved sample field: `obs/_index`.
- Sample-ID provenance: reconstructed from obs/_index final suffix; format defined by the official atlas integration code (barcode + '_' + sample_id).
- Sample-index format valid for all PRJNA607098 cells: True.
- PRJNA607098 cells: 34569.
- Distinct sample IDs: 12.
- Public-record expected sample count: 12.
- Exact expected sample count: True.
- Eligible paired sample-state units: 12/12.
- Missing sample-ID fraction: 0.000000.
- Largest sample fraction: 0.315919.
- Proceed to Step 08C2: True.

### Evidence boundary

- This audit checks sample-level feasibility only; it performs no hypothesis test.
- Recovered SRS accessions are sample-level units, not automatically independent donors.
- Step 08C2 must retain the accession provenance and verify donor independence before donor-level claims.
- Step 08C2 must aggregate within sample before inference; individual cells are not replicates.
- A successful audit does not test protein abundance, channel activity or causality.
