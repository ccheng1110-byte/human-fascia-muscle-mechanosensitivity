# Step 15A2 GSE300230 gate adjudication v1

## Adjudicated outcome

- Original automated gate: `HOLD`.
- Adjudicated gate: `PARTIAL_PROCEED_WITH_STRICT_BOUNDARIES`.
- Proceed to Step 15B: `TRUE`.

## Why the automated gate failed

1. The parser treated `NoTension` as ambiguous because the string also contains `Tension`. Official GEO characteristics explicitly distinguish `tension` from `no tension`; this is a deterministic parsing defect rather than unresolved metadata.
2. GSE300230 contains two primary dermal fibroblast cell lines: GM08401 (75-year-old donor) and GM09503 (10-year-old donor). The original feasibility rule expected at least four biological units and therefore failed.

## Corrected design interpretation

- The 16-sample `GSE300230_raw_counts_tensionTGFb.csv.gz` matrix contains a complete 2 x 2 mechanical-tension-by-TGF-beta design within each cell line, with two recorded batches/replicates per condition.
- Mechanical condition has two levels: `relaxed` and `tension`.
- TGF-beta has two levels: `absent` and `present`.
- The frozen 15-gene candidate panel is covered 15/15.
- Donor age is perfectly confounded with cell line. No population-level age effect, donor-general age interaction, or age-specific mechanism may be inferred.
- With only two cell lines, perturbation results establish at most external mechanochemical plausibility. They do not constitute independent donor-level replication and cannot by themselves upgrade the overall evidence above `CAUTION`.

## Authorized next step

Step 15B may freeze a stratified analysis contract that:

1. analyzes GM08401 and GM09503 separately;
2. treats batch 1/2 as a blocking factor;
3. uses the TGF-beta-absent tension-versus-relaxed contrast as the primary mechanical contrast;
4. retains the TGF-beta-present mechanical contrast and mechanical-by-TGF-beta interaction as prespecified secondary contrasts;
5. reports cross-line directional concordance without treating two lines as a population sample;
6. preserves all frozen modules, candidates, competitor programs, and null results.

## Material Passport

- Source: official GSE300230 GEO sample characteristics, downloaded processed matrices, and Step 15A audit artifacts.
- Transformations: deterministic correction of `no tension` parsing and feasibility-rule adjudication.
- Expression values inspected: matrix headers and gene identifiers only; no condition-specific expression result was examined.
- Integrity status: prospective analysis remains untested at this adjudication point.

