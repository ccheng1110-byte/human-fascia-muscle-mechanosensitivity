## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-24
- Verification Status: ANALYZED
- Version Label: PRJNA607098_08C2_interpretation_accession_audit_v1

## Step 08C2 interpretation

The frozen overall gate failed and must remain failed. Ten of fifteen candidates
and five of six integrin/focal-adhesion genes passed, but the pre-specified gate
also required PIEZO2 to pass.

PIEZO2 did not show a borderline failure. It was lower in F7 than pooled non-F7
cells in 9/12 samples, with median paired difference -0.1109, paired rank-biserial
effect -0.6923, one-sided Wilcoxon p = 0.9829 and BH FDR = 0.9829. Therefore the
external sample-level evidence does not support F7-specific PIEZO2 enrichment.

Module-level interpretation:

- actomyosin/Rho: 4/4 strict support;
- integrin/focal adhesion: 5/6 strict support;
- Hippo/YAP/TAZ: 1/2 strict support;
- mechanosensor channels: 0/3 strict support.

The primary biological direction should therefore be narrowed to an
integrin-focal-adhesion/actomyosin transcriptional state. PIEZO2 and other
mechanosensor channels remain heterogeneous exploratory candidates and must not
be presented as externally replicated F7 markers.

## Accession and sample-unit audit

NCBI PRJNA607098 metadata contains 12 SRA experiments and 12 BioSamples. The 12
atlas SRS identifiers map one-to-one to 12 distinct BioSample accessions, 12
distinct `isolate` identifiers, 12 SRX experiments and 12 SRR runs. The records
cover two sequencing batches (six samples each), ages 42-68 years, nine male and
three female hand-tissue samples.

This one-to-one mapping supports treating the 12 atlas SRS units as distinct
specimen/donor records for the paired analysis. It does not turn the
observational transcriptomic comparison into functional or causal evidence.

Primary sources:

- NCBI BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA607098
- NCBI SRA RunInfo: https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA607098
- NCBI E-utilities SRA metadata: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/

## Evidence decision

- Frozen combined claim: CAUTION; gate failed.
- Integrin/focal-adhesion and actomyosin subclaim: independently sample-replicated transcriptomic support.
- PIEZO2/F7-specific channel subclaim: not externally replicated and directionally contradicted in this source.
- Protein abundance, channel activity, force response and causality remain untested.

## Statistical fallacy scan (11/11 checked)

- Simpson's paradox: the earlier median-across-states summary masked the actual
  sample-weighted non-F7 comparison; Step 08C2 is preferred.
- Ecological fallacy: no individual clinical outcome is inferred from aggregated expression.
- Berkson's paradox: surgical and atlas inclusion remain potential selection biases.
- Collider bias: no adjusted causal model was fitted.
- Base-rate neglect: sample and cell counts were retained.
- Regression to the mean: no extreme-value repeat-selection design was used.
- Survivorship bias: atlas QC may exclude low-quality cells or samples.
- Look-elsewhere effect: BH correction was applied to the frozen 15-gene family.
- Garden of forking paths: the failed PIEZO2 gate was not changed after observing results.
- Correlation is not causation: causal claims are prohibited.
- Reverse causality: cross-sectional expression cannot establish mechanotransduction direction.
