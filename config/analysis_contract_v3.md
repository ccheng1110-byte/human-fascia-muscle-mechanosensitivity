# Human Fascia–Muscle Mechanosensitivity

## Step 10B analysis contract v3

**Freeze date:** 2026-08-25  
**Purpose:** freeze the confirmatory analysis before any new expression-level result is inspected.  
**Applies to:** GSE173252 discovery, PRJNA607098 atlas validation, and any later independent source.

## 1. Current evidence boundary

Step 10A correctly resolved 30 official atlas source rows. The official inventory contains both GSE173252 and PRJNA607098, so only 28 source rows remain after source-level exclusion. The atlas source inventory has no donor/sample field. Therefore:

- the integrated atlas is not a fully independent donor-level replication of the current study;
- source-level mapping may be used for provenance and exploratory validation;
- donor-level confirmation must rely on a source with sample metadata and no overlap with either discovery or current validation data;
- the paper label F8 “fascia-like myofibroblast” and the atlas label F7 “Fascia-like myofibroblast” remain a crosswalk issue. The exact atlas label may be analysed, but it must not be described as a proven one-to-one F8 replication.

## 2. Frozen mechanotransduction modules

The exact gene sets are stored in `config/mechanotransduction_module_registry_v2.csv` and are inherited from `data/gene_sets/mechanosensitivity_program_v1.csv`.

### 2.1 Core confirmatory modules

- `integrin_focal_adhesion`: extracellular-matrix force coupling and focal-adhesion transmission.
- `actomyosin_rho`: Rho/actomyosin contractility and cytoskeletal force generation.

These two modules are the primary mechanotransduction claim. They are evaluated at the sample level and are not interpreted as causal proof.

### 2.2 Mechanism-gate modules

- `mechanosensor_channels`: retained as a specificity gate, not as a mandatory prerequisite for the broad claim.
- `hippo_yap_taz`: retained as a downstream mechanotransduction module; a positive NF2 result is not equivalent to direct YAP/TAZ activation.

The negative PIEZO2 result in Step 08C2 is therefore a boundary on channel-specific interpretation, not a reason to discard the integrin–actomyosin hypothesis.

### 2.3 Explicit competing explanations

The following modules must be shown beside every confirmatory result:

- `ecm_remodeling`
- `tgf_fibrosis`
- `inflammation_ap1_nfkb`
- `hypoxia`
- `cell_cycle`

They are not optional covariates to be selected after seeing the data. They are fixed alternative explanations. The competition analysis is exploratory and must not be presented as causal adjustment in a 12-sample dataset.

## 3. Frozen candidate panel

The 15-gene panel in `config/frozen_candidate_panel_v2.csv` is inherited from Step 08C2 and cannot be expanded, substituted, or re-ranked after Step 10C results are seen.

- Integrin/focal adhesion: ITGAV, ITGB1, ITGA2, PARVA, ITGA5, FERMT2.
- Actomyosin/Rho: CFL1, LIMK1, CNN2, CDC42.
- Mechanosensor channels: PIEZO2, TMEM63B, PANX1.
- Hippo/YAP/TAZ downstream: NF2, TEAD1.

Step 08C2 supported 10/15 candidates, 5/6 integrin/focal-adhesion candidates, and did not support PIEZO2. These are frozen prior results, not a permission to claim independent replication.

## 4. Frozen target and comparator states

### Primary contrast

`PRJNA607098`: exact atlas-defined F7/Fascia-like myofibroblast label versus pooled non-F7 cells within the same sample.

The target label must be read from the validated atlas field; no relabelling from paper prose is allowed. Each eligible sample contributes one target mean, one comparator mean, and their difference.

### Sensitivity contrasts

1. F7 versus `Reticular: universal`, only when both states meet the cell-count rule.
2. F7 versus the pooled non-F7 fibroblast-state comparator, only when state labels and sample coverage are sufficient.

These contrasts are sensitivity analyses. They do not replace the primary pooled non-F7 contrast.

### Cell and sample exclusions

- Exclude cells failing the dataset’s pre-existing quality/doublet rules.
- Do not treat cells as biological replicates.
- Keep the 12 eligible paired samples for the primary contrast if the validated Step 08C2 construction is reproduced.
- For sensitivity contrasts, require at least 50 cells in each group within a sample; report excluded samples explicitly.
- Never mix GSE173252 cells with PRJNA607098 cells in a single donor-level test.

## 5. Frozen statistical contract

1. Calculate gene-level and module-level values within each sample and state.
2. Compute `F7 - comparator` for every sample.
3. Report median difference, IQR, number positive/negative/zero, paired rank-biserial effect, one-sided paired Wilcoxon p value, and BH FDR within the frozen candidate family or module family.
4. Use the 15 candidate genes as a confirmatory family; do not add genes because they are significant in Step 10. Candidate-level strict support retains the Step 08C2 rule: positive median, at least 10/12 positive samples, and BH-adjusted one-sided q <= 0.05.
5. For core modules, require a positive median and at least 9/12 positive samples for a module-level pass. A p value cannot rescue a directionally inconsistent module.
6. Perform leave-one-sample-out reanalysis for the two core modules and the primary composite. A sign reversal downgrades robustness.
7. Cell-level p values are descriptive only and cannot be used to upgrade evidence.

## 6. Competition analysis

For each core module, report the unadjusted paired difference and then a pre-specified robustness audit against each competitor module. The audit should include:

- correlation of core and competitor sample-level differences;
- high-versus-low competitor descriptive stratification using the frozen median split;
- a simple leave-one-sample-out residual or rank-based sensitivity check.

The signal is considered residual only if the core direction remains positive and at least 8/12 samples remain positive in the robustness audit. If the signal attenuates to near zero or reverses after an explicit competitor audit, the result is reclassified as a shared ECM/fibrosis/state signal.

Because n = 12 is small, no high-dimensional multivariable model, machine-learning classifier, or post-hoc covariate selection is allowed in the confirmatory layer.

## 7. Evidence-grade rules

- `CAUTION`: any source overlap, unresolved label crosswalk, failed candidate specificity, or failed leave-one-out robustness remains.
- `SUPPORTED_WITH_BOUNDARY`: both core modules pass sample-level and robustness gates, but source independence or F8/F7 mapping remains partial.
- `INDEPENDENT_REPLICATION_SUPPORTED`: requires a non-overlapping source with sample/donor metadata, a valid label crosswalk, and the same frozen core/competition contract passing.
- `MECHANOSENSOR_SPECIFICITY_SUPPORTED`: reserved for a separately positive channel-specific result; it is not implied by a positive integrin/actomyosin result.

The current Step 08C2 result remains `CAUTION` because the integrated atlas is not independent and PIEZO2 is not supported. The intended direction is therefore a narrower, stronger claim: fascia-like fibroblast states show a reproducible integrin–actomyosin-associated mechanotransduction program, while channel-specific PIEZO2 evidence is unresolved or negative.

## 8. Next data strategy

1. Complete Step 10C with the frozen PRJNA607098 contract.
2. Complete Step 10D as a provenance-aware evidence synthesis; do not count the integrated atlas as a second cohort.
3. Prioritize one genuinely non-overlapping source with donor/sample metadata over downloading the 617 GiB raw PRJNA1206333 collection.
4. Use processed matrices or targeted gene-panel extraction whenever possible; reserve full raw download for a pre-specified necessity that cannot be met by processed data.
5. If no independent donor-level source is recoverable, publish the work as a carefully bounded reanalysis/triangulation study and state the missing independent replication as the main limitation.

## 9. Timeline and risk response

| Phase | Deliverable | Target duration | Triggered fallback |
|---|---|---:|---|
| Step 10C | Frozen PRJNA607098 sample-level reanalysis | 1 day | Reuse validated Step 08C2 objects; no new labels |
| Step 10D | Provenance and evidence-grade synthesis | 1 day | Keep Gate A/B partial and retain CAUTION |
| Step 11 | One non-overlapping donor-level source | 3–5 days | Search source inventory; avoid >50 GiB raw download |
| Step 12 | Final figures, limitations, and manuscript evidence table | 2–3 days | Report descriptive atlas validation if independent source fails |

Risk responses:

- **Source overlap:** remove the source from confirmatory replication and retain it only for exploratory mapping.
- **F8/F7 label mismatch:** report the exact atlas label and downgrade direct anatomical equivalence.
- **PIEZO2 negative:** focus the hypothesis on integrin–focal adhesion and actomyosin/Rho; do not force a channel mechanism.
- **Low sample coverage:** retain primary pooled comparator, mark sensitivity contrasts as unavailable, and report coverage.
- **SSL/reticulate/download failure:** use cached metadata and small targeted chunks; do not redownload the full atlas.
- **No independent cohort:** stop escalation of evidence grade and convert the final output to a transparent, hypothesis-generating reanalysis.

## Material Passport

- **Material type:** configuration and analysis-contract Markdown/CSV files.
- **Source:** local project evidence, especially Step 08C2 paired validation and corrected Step 10A provenance audit.
- **Transformations:** frozen modules and candidate panel were transcribed from the existing v1 gene-set registry and Step 08C2 candidate statistics; no new expression values were computed.
- **Computational environment:** Windows/RStudio project at `.`; UTF-8 text.
- **Key assumptions:** PRJNA607098 remains the current sample-level validation source; the integrated atlas is partially overlapping; the exact F7 label is retained.
- **Reproducibility:** all later scripts must read these files rather than recreate gene lists or thresholds inline.
- **Integrity status:** Step 10B is a pre-analysis freeze; it does not itself upgrade the evidence grade.
