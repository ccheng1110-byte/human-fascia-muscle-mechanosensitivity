# Step 15E computational-strengthening evidence resynthesis v1

## Final decision

- Strengthening package status: **COMPLETE**.
- Overall evidence grade: **CAUTION**.
- Strengthening interpretation class: **ACTOMYOSIN_PERTURBATION_PLAUSIBILITY_WITH_PROLIFERATION_COMPETITION**.
- Manuscript-level upgrade: **not authorized**.

## Evidence domains

### 1. Observational fascia-like state association

The discovery and atlas analyses support an association between fascia-like fibroblast states and a broader ECM/integrin–cytoskeletal transcriptional program. The atlas source has provenance overlap and incomplete F7/F8 equivalence; candidate specificity and competition robustness are not uniformly positive. This domain remains `CAUTION`.

### 2. Integrin/focal-adhesion support

Across the GSE300230 primary mechanical contrast, the integrin/focal-adhesion module was up in both cell lines, but module-family q values were 0.00027 for GM08401 and 0.115 for GM09503. The score sensitivity analysis retained positive direction in both lines, while the candidate-level result was only 8/15 positive in both lines and 0/15 q<=0.05 in both. This is directional, context-dependent support, not a universal integrin mechanism.

### 3. Actomyosin perturbation plausibility

The actomyosin/Rho module was up in both cell lines for tension versus relaxed conditions without exogenous TGF-beta, with module-family q values of 5.07e-09 and 1.20e-05. The positive direction was stable across z-mean, z-median, rank-mean and centered-mean scoring and remained positive in leave-one-gene-out analyses. This is the strongest new domain, but it is restricted to cultured dermal fibroblasts and two cell lines.

### 4. Competition and specificity

The cell-cycle module was also up in both cell lines for the primary mechanical contrast, with q values of 0.00225 and 5.12e-12. Expression-matched negative controls did not eliminate the actomyosin or integrin score effects, but they do not distinguish actomyosin from a biologically coupled proliferation response. Cell-intrinsic versus composition-associated explanations were `NOT_ESTIMABLE` because the input was a gene-by-sample cultured-fibroblast matrix without cell-level composition data.

### 5. Mechanosensor-specific evidence

The mechanosensor channel result was not sufficient for a channel-specific claim. PIEZO2 remained inconsistent across sources and must not be presented as validated activation. Hippo/YAP/TAZ evidence is also insufficient for a direct YAP/TAZ activation claim.

### 6. Age, fascia specificity and causality

GSE300230 contains one older and one young cell line, so age is perfectly confounded with cell line. No independent age effect may be inferred. The cultured dermal-fibroblast perturbation is not fascia tissue and does not contain a pain phenotype. Therefore fascia specificity, donor-population generalization and causal relevance remain unestablished.

## Revised primary hypothesis

Human fascia-like fibroblast states are associated with a context-dependent ECM/integrin–actomyosin transcriptional program, and mechanical tension can induce an actomyosin-associated transcriptional response in cultured human dermal fibroblasts. The integrin response is less uniform, mechanosensor-specific PIEZO2 evidence is unresolved, and the observed response may be partly coupled to proliferation-associated transcription. These findings are mechanistic plausibility evidence rather than proof of fascia-specific or pain-relevant causality.

## What can be claimed

- A reproducible actomyosin-associated response to mechanical tension is biologically plausible in the analyzed cultured human dermal fibroblast models.
- The original fascia-like cell-state association is compatible with, but does not prove, a mechanotransduction-related program.
- Integrin/focal-adhesion signals are directionally supported but not uniformly specific across candidates, sources or contexts.

## What cannot be claimed

- Direct activation of PIEZO2 or another specific mechanosensor.
- A fascia-specific mechanosensitivity mechanism.
- An independent human age effect.
- Cell-intrinsic rather than composition-associated causation from GSE300230.
- Causal relevance to pain, fibromyalgia, or clinical mechanosensitivity.

## Manuscript action

The manuscript should be revised by adding GSE300230 as a bounded perturbation-plausibility analysis, explicitly displaying cell-cycle competition, and replacing broad mechanistic language with the revised hypothesis above. The overall evidence label remains `CAUTION`.

## Material Passport

- Inputs: Step 12C frozen claims; Step 15C mechanochemical perturbation; Step 15D score-sensitivity and negative-control audits; prior source-provenance decisions.
- New expression computation in this synthesis: none.
- Transformation: domain-wise evidence separation, claim-boundary adjudication and manuscript action mapping.
- Integrity status: complete with explicit unresolved debts for proliferation competition, cell composition, age confounding, fascia specificity and causality.

