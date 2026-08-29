# FAIR and repository-readiness audit

## Findable

- Versioned title, creators, keywords, release date, and citation metadata are provided.
- Public source datasets are mapped to persistent GEO or BioProject accessions.
- A Zenodo DOI remains to be assigned after deposition.

## Accessible

- All release files use open, commonly readable formats.
- Large third-party datasets are accessed through their original public repositories.
- No sensitive or directly identifiable participant data are included.

## Interoperable

- Machine-readable outputs use CSV, CSV.GZ, JSON, or plain text.
- Figure source data are mapped to individual figures.
- Script and file manifests use repository-relative paths and SHA-256 hashes.

## Reusable

- Code and derived outputs have explicit licences.
- Frozen candidate panels, module definitions, analysis contracts, gates, and decision records are included.
- R and Python requirements and an execution guide are included.
- Third-party datasets are not relicensed.

## Remaining release actions

1. Create the public GitHub repository.
2. Publish a tagged GitHub release named `v1.0.0`.
3. Archive that release in Zenodo and obtain the version DOI.
4. Add the GitHub URL and Zenodo DOI to the manuscript Data Availability statement.
5. Optionally add the DOI and repository URL to `CITATION.cff` and the README in a follow-up commit.
