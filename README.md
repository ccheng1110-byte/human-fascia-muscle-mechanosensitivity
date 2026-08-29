# Human fascia–muscle mechanosensitivity transcriptomic analysis

This repository contains the analysis code, frozen analysis contracts, machine-readable result tables, decision records, and figure source data supporting the manuscript:

> Cross-tissue transcriptomic evidence supports a mechanotransduction-compatible ECM–integrin–cytoskeletal programme in human fibroblast states

## Scientific scope

The project reanalyses publicly available human transcriptomic datasets to test whether an extracellular-matrix, integrin/focal-adhesion, and actomyosin/Rho programme is reproducible across fibroblast states. The strongest defensible interpretation is an evidence-bounded mechanistic hypothesis supported by computational triangulation. The analyses do not establish fascia specificity, cell-intrinsic causality, a universal PIEZO2 mechanism, or a donor-replicated functional mechanism.

## Repository contents

- `scripts/`: R, Python, and PowerShell scripts used for data acquisition, audit, analysis, evidence synthesis, figure generation, manuscript support, and reference checking.
- `config/`: frozen candidate panels, module registries, analysis contracts, gates, and sample maps.
- `results/`: machine-readable analysis outputs and decision records. Manuscript drafts and rendered review pages are intentionally excluded.
- `figures/`: final Figures 1–6 in PDF, PNG, SVG, and TIFF formats, with figure source-data tables and legends.
- `data/`: public accession inventory and explanation of why third-party raw data are not redistributed.
- `environment/`: R and Python dependency information.
- `metadata/`: release manifests and SHA-256 checksums.
- `docs/`: reproduction instructions, release notes, and Zenodo metadata guidance.

## Data availability

No new primary sequencing data were generated. The analyses reuse public NCBI GEO and BioProject resources listed in `data/public_accession_inventory.csv`. Large third-party raw and processed objects are not redistributed in this repository. They should be obtained from the original repositories using the listed accessions and the acquisition scripts in `scripts/`.

## Reproduction

Run the scripts from the repository root. On Windows with RStudio, open a new R session, set the working directory to this folder, and source the required scripts in the order described in `docs/REPRODUCIBILITY.md`.

Several workflows use `reticulate` to stream selected chunks from a public remote Zarr resource. These scripts declare their Python requirements before Python is initialised. Restart R before running such a script if `reticulate` has already initialised Python.

## Evidence and provenance

The frozen candidate panel and module definitions are stored under `config/`. Analysis decisions and boundary statements are retained under `results/`. Source data underlying the final figures are stored alongside the figure files in `figures/`.

## Licences

- Analysis code in `scripts/` is released under the MIT License; see `LICENSE`.
- Derived result tables, figures, configuration records, and documentation are released under the Creative Commons Attribution 4.0 International License; see `LICENSES/CC-BY-4.0.txt`.
- Third-party source datasets remain subject to the terms of their original repositories and are not relicensed here.

## Citation

Use the metadata in `CITATION.cff` when citing this software release. After Zenodo deposition, cite the archived release DOI shown on the Zenodo record.

## Contact

Correspondence concerning the associated manuscript may be directed to Jian Sun, Guangzhou University of Chinese Medicine.
