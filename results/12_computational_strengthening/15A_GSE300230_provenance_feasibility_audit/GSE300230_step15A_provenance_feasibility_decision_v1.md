# Step 15A GSE300230 provenance and feasibility audit

## Scope

- Official GEO metadata and supplementary-file inventory were audited.
- Processed text matrices were downloaded when available.
- No FASTQ, BAM, H5AD, or GEO RAW.tar archive was downloaded.
- No condition-specific expression comparison or biological hypothesis test was performed.

## Gate

- Step 15A gate: **HOLD**.
- Proceed to Step 15B contract freeze: **FALSE**.
- GSM samples: 24.
- Resolved biological units: 2.
- Mechanical levels: tension.
- TGF-beta levels: absent | present.
- Age levels: older | young.
- Best frozen candidate coverage: 15/15.

## Interpretation boundary

A PASS or PARTIAL gate establishes feasibility only. It does not support mechanotransduction, specificity, or causality. Step 15B must freeze the factor coding, biological-unit model, selected matrix, aliases, endpoints, and multiplicity rules before Step 15C examines expression contrasts.

## Material Passport

- Origin: official NCBI GEO GSE300230 metadata and supplementary directory.
- Local frozen inputs: frozen_candidate_panel_v2.csv and mechanotransduction_module_registry_v2.csv.
- Transformation: metadata flattening, conservative factor reconstruction, file classification, and gene-identifier coverage audit.
- Integrity boundary: unresolved metadata fields require explicit Step 15B adjudication; they must not be inferred from expression results.
- Supplementary listing warning: none.
