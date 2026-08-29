## Material Passport

- Origin Skill: academic-research-suite / experiment-agent + nature-downloader
- Origin Mode: validate / lawful OA source audit
- Origin Date: 2026-08-24
- Verification Status: ANALYZED
- Version Label: GSE273293_processed_data_recovery_09B_v1

## Step 09B processed-data recovery audit

### Objective

Determine whether the 14 PRJNA1206333 clinical samples can be recovered in a
sample-resolved processed expression object without downloading and
reprocessing approximately 617.25 GiB of raw sequence data.

### Sources audited

- Local GSE273293 GEO archive and series metadata.
- NCBI GSE273293 / PRJNA1141379.
- NCBI PRJNA1206333 and the Step 09A sample-accession map.
- The article's Data Availability and Code Availability statements.
- Targeted web/GitHub discovery using the accession, BioProject, DOI, exact
  title and author-name combinations.
- Publisher/PMC records for the two electronic supplementary DOCX files.
- PMC open-access package inventory.

### Findings

1. The local GEO archive contains one matrix triplet only:
   `GSM8425828_GMC01_{barcodes,genes,matrix}`.
2. The local GEO metadata contains one row and therefore cannot support 14
   sample-level expression inference.
3. PRJNA1206333 resolves the intended 10 GMC + 4 control design into 14
   one-to-one raw Run/Experiment/SRA-sample/BioSample records, but it does not
   expose a sample-resolved processed expression object in the audited public
   NCBI records.
4. The paper's Data Availability statement points only to GSE273293.
5. The paper says that original code will be deposited at GitHub, but gives no
   repository URL. No matching study-specific repository was found in the
   targeted search as of 2026-08-24.
6. Two public supplementary DOCX files exist (approximately 8.5 MB and 185.4
   KB), but they have not been downloaded because Supporting Information
   requires explicit user confirmation. Their ability to recover a
   barcode-to-sample map therefore remains unknown.

### Gate decision

- Step 09A metadata component: **PASS**.
- Direct public processed-object recovery: **NOT FOUND**.
- Supplementary-material inspection: **PENDING USER CONFIRMATION**.
- Full raw-data reprocessing: **NOT AUTHORIZED / NOT RECOMMENDED**.
- Current Step 09B state: **PARTIAL BLOCK — SI review and/or author request**.

### Preferred recovery sequence

1. With explicit approval, download and inspect both supplementary DOCX files
   for sample metadata, per-sample cell counts, barcode mapping, processed
   object links, repository links or code references.
2. If no cell-level sample mapping is present, send the prepared corresponding-
   author request for one of the following:
   - a sample-resolved Seurat/SingleCellExperiment/AnnData object;
   - per-sample 10X count matrices;
   - a `cell_barcode -> GMC/Con sample_id` table;
   - cell-type/subcluster annotations and preprocessing scripts.
3. If neither route succeeds, close Gate A as `processed-expression failed` and
   retain GSE273293 only as pooled descriptive deep-fascia evidence.

### Claim boundary

- Fourteen raw accession records do not imply that the one-GSM processed matrix
  contains recoverable donor labels.
- Per-sample composition tables in a supplement would support abundance
  summaries but would not reconstruct sample-level gene expression unless cell
  barcodes or per-sample matrices are also supplied.
- No biological hypothesis is tested in Step 09B.

## Source anchors

- GSE273293: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273293
- PRJNA1141379: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1141379
- PRJNA1206333: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1206333
- Article: https://doi.org/10.1186/s12967-024-05889-y
- PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC11834283/
- PMC OA inventory: https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC11834283

