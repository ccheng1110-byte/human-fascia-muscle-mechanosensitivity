## Material Passport

- Origin Skill: academic-research-suite / experiment-agent; nature-downloader
- Origin Mode: validate / lawful public-data audit
- Origin Date: 2026-08-23
- Verification Status: SOURCE_AUDITED
- Version Label: step08_external_validation_design_v1

## Step 08 validation target

The immediate target is the 2025 human skin fibroblast atlas from Steele et al. The atlas contains both the discovery source used in Steps 03-07 (`GSE173252`) and a second Dupuytren source (`PRJNA607098`). Only the latter can contribute independent-dataset evidence.

The claim under test is deliberately narrow:

> The frozen Step-07 shortlist represents a fascia/myofibroblast-enriched transcriptional state, with selected mechanotransduction candidates, rather than a universal or proven causal force-sensing program.

## Local search outcome

A separate local archive directory was searched recursively on 2026-08-23. No `.h5ad`, `.h5seurat`, `.loom`, `.zarr`, `adata_webportal`, `cellatlas`, `PRJNA607098`, or matching atlas record was found. Files larger than 1 GB in that directory were unrelated GWAS/LAVA resources.

## Official source audit

- Article: https://www.nature.com/articles/s41590-025-02267-8
- Author repository: https://github.com/haniffalab/skin_fibroblast_atlas
- Atlas portal: https://cellatlas.io/studies/skin-fibroblast
- Full H5AD: https://storage.googleapis.com/haniffalab/skin-fibroblast/adata_webportal.h5ad
- Cloud Zarr: https://storage.googleapis.com/haniffalab/skin-fibroblast/zarr/adata_webportal.zarr/
- Independent Dupuytren source: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA607098

The full H5AD server header reports 27,231,223,312 bytes (approximately 25.36 GiB), ETag `e2cf64a7d04d6afe87f9278d046c3a46`, last modified 2024-12-04. It was not downloaded.

The consolidated Zarr metadata reports 357,276 cells by 36,601 genes. Its expression matrix is chunked as all 357,276 cells by 10 genes. This supports future gene-panel extraction without obtaining the whole H5AD. The public webportal observation table contains only:

- `GSE`
- `Patient_status`
- `disease_category_orig`
- `celltype`
- `lesional_vs_nonlesional`
- `celltype_skinspecific_nomenclature`

It does not expose a sample or donor identifier. Therefore, source-specific atlas analysis can be descriptive but cannot provide donor-level inferential statistics.

## Files already saved (all below 50 MB)

Author-provided dataset inventory, marker rankings, supplementary marker tables, and consolidated Zarr metadata are stored under:

`./data/external_validation/skin_fibroblast_atlas_2025`

The largest saved file is `DEGs_results.csv` at 14,643,361 bytes. SHA-256 values and source URLs are recorded in `metadata/step08_small_file_download_manifest.csv`.

## Validation ladder

1. **Step 08A — current:** analyze the authors' marker-rank table against the frozen Step-07 shortlist. This is a low-volume exploratory screen.
2. **Step 08B — only if 08A is positive:** extract the frozen genes from cloud Zarr and stratify by `GSE`, explicitly isolating `PRJNA607098` from `GSE173252`. This remains descriptive because donor labels are absent.
3. **Step 08C — evidence-upgrading route:** reconstruct sample-level matrices from `PRJNA607098` and perform donor/sample pseudobulk comparisons. Raw sequence volume is approximately 20 GB and would require user-managed downloads.
4. **Biological validation:** protein localization or perturbation is still required before any causal mechanosensitivity claim.

## 50 MB download rule

No project script may automatically download a file larger than 50 MB. If the full atlas is later chosen, the user-managed destination is:

`./data/external_validation/skin_fibroblast_atlas_2025/manual_downloads/adata_webportal.h5ad`

The recommended route is currently cloud-Zarr gene-panel extraction, not the full H5AD.
