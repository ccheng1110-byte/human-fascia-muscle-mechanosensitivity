## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-24
- Verification Status: ANALYZED
- Version Label: GSE273293_PRJNA1206333_sample_recovery_09A_v1

## Step 09A PRJNA1206333 sample-recovery metadata audit

### Scope

- This step downloads NCBI RunInfo metadata only.
- No SRA, FASTQ, BAM or expression matrix is downloaded.
- No expression hypothesis is tested.

### Metadata result

- Metadata gate passed: TRUE.
- Runs / experiments / SRA samples / BioSamples: 14 / 14 / 14 / 14.
- Clinical labels: 10 GMC and 4 nonfibrotic controls.
- Estimated total raw archive size: 617.25 GiB.
- Every raw run is larger than 50 MB; raw download is not part of Step 09A.

### Local GSE273293 result

- Local GEO series-metadata rows: 1.
- Local supplementary TAR sample prefixes: GSM8425828_GMC01.
- The current local GEO files do not support 14-sample donor-level inference.

### Gate decision

- Gate A metadata component: PASS.
- Gate A processed-expression component: BLOCKED/UNKNOWN.
- Full raw-data download: DO NOT PROCEED.
- Proceed to Step 09B: search for a sample-resolved processed object,
  barcode-to-sample map, or author-supplied processed matrices.
- If Step 09B fails, retain GSE273293 as pooled descriptive deep-fascia
  evidence and do not represent it as a 14-donor expression analysis.

### Evidence boundary

- Fourteen accession records establish sample availability, not processed
  cell-level donor labels in the local matrix.
- Library names support the 10 GMC + 4 control design but do not by
  themselves validate every clinical variable reported in the paper.
- No disease effect, mechanotransduction effect or causal claim is tested.
