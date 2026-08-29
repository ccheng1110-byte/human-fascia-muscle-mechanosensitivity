# Step 11A independent external source screening

## Scope

- Only accession-resolved sources not overlapping GSE173252 or PRJNA607098 were screened.
- This step downloaded GEO metadata only; no expression matrix, H5AD, or raw sequencing file was downloaded.
- Priority scores are triage rules, not evidence of replication.

## Screening result

- Candidate GSE sources audited: 27.
- High-priority manual object audits: 1.
- Metadata audit failures requiring follow-up: 0.

## Selection rule for Step 11B

Prioritize sources with sample-level donor/subject fields, explicit fibroblast or myofibroblast labels, and single-cell evidence. Before expression analysis, manually verify the object-level cell-state field, source accession, sample/donor mapping, and availability of the frozen candidate genes.

## Material Passport

- Origin: corrected Step 10A source map.
- Transformation: GEO series metadata were flattened into a source-level priority table and a sample-level metadata snapshot.
- No expression values were read or tested.
- Reproducibility: GEO metadata cache is under `.runtime/geo_step11a/geo_metadata_cache`; all exported tables are in the result directory.
- Integrity boundary: keyword-based priority is only a screening aid and requires manual verification before validation.
