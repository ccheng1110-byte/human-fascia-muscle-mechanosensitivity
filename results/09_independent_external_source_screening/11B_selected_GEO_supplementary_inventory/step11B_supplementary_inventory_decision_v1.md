# Step 11B selected-source supplementary inventory

## Scope

- Selected sources: GSE175817, GSE130973, GSE202352, GSE228421, and GSE138669.
- Only FTP directory listings were read.
- No supplementary file was downloaded.

## Next selection rule

Prioritize a processed single-cell object with explicit cell-state labels and a sample/donor field. A raw archive, a matrix without cell-state labels, or a repeated-measures source without donor mapping is not sufficient for independent donor-level replication.

GSE228421 must be analysed as a repeated-patient/longitudinal design if selected; its GSM count must not be treated as the number of independent biological samples.

## Material Passport

- Origin: Step 11A metadata-only screening and official NCBI GEO FTP directory listings.
- Transformation: file names and candidate URLs were inventoried; no expression values were read.
- Reproducibility: all URLs and classification flags are exported in the result directory.
- Integrity boundary: file-name classification is preliminary and requires inspection of the selected object before download.
