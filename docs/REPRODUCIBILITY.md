# Reproducibility guide

## 1. Obtain the repository

Download or clone the repository and start R from the repository root. All copied scripts use relative project paths, so the working directory must be the repository root.

```r
setwd("path/to/human-fascia-muscle-mechanosensitivity")
normalizePath(getwd(), winslash = "/")
```

## 2. Install dependencies

Install CRAN and Bioconductor packages listed in `environment/requirements_R.txt`. Remote-Zarr scripts use `reticulate::py_require()` and provision their bounded Python dependencies automatically. Restart R before each remote-Zarr workflow if Python has already been initialised.

## 3. Retrieve external datasets

The public accession inventory is stored in `data/public_accession_inventory.csv`. Begin with the relevant acquisition or audit script rather than copying source data into Git.

Large files may require several gigabytes of free disk space. The original GSE173252 and GSE130973 processed objects are not included in this release.

## 4. Main analysis sequence

The numbered filenames preserve the historical execution order. The principal evidence chain is:

1. Steps 01–07: GSE173252 acquisition, audit, cluster reproduction, scoring, and specificity analysis.
2. Steps 08–10: PRJNA607098 atlas validation, sample reconstruction, frozen-contract analysis, and evidence synthesis.
3. Steps 11–12: independent-source screening, GSE130973 analysis, and first evidence closeout.
4. Steps 15C–15D: GSE300230 mechanical-tension strengthening and score-sensitivity audit.
5. Steps 17B–17H: second-round co-expression, TGF-beta/TEAD, and stiffness analyses.
6. Step 17F synthesis: final computational evidence resynthesis.
7. Steps 20B and 21B: final Figures 1–6 and figure source-data generation.

Run one script at a time and inspect its decision record before proceeding. For example:

```r
source("scripts/17E_analyze_GSE338388_S2_frozen_contract.R", encoding = "UTF-8")
```

## 5. Frozen analysis contracts

Do not alter candidate genes, module definitions, statistical thresholds, or gate rules when reproducing the confirmatory analyses. These are recorded in `config/` and in the associated decision files under `results/`.

## 6. Expected boundaries

Reproduction should retain the CAUTION evidence grade. The computational results support an actomyosin-centred, mechanotransduction-compatible programme with computational triangulation. They do not establish causal mechanotransduction, fascia specificity, cell-intrinsic regulation, or donor-replicated functional effects.

## 7. File verification

After downloading an archived release, verify files against `metadata/SHA256SUMS.txt`. The complete inventory, including byte size and SHA-256 hash, is stored in `metadata/file_manifest.csv`.
