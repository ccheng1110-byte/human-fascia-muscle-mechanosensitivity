# Step 16A functional-validation protocol lock

Date: 2026-08-26  
Status: **Protocol planning gate; no experimental results generated**  
Current computational evidence grade: **CAUTION**

## Decision and stage boundary

The computational phase is closed for the current claim set. The next phase is a focused, donor-replicated functional validation rather than another unrestricted public-data screen. This document locks the scientific logic, endpoint hierarchy, analysis unit, and decision rules. Exact reagent catalogue numbers, instrument settings, and the final mechanical dose must be filled in and approved before confirmatory data collection.

Deferred from this stage: direct clinical pain testing, a population-level age analysis, a claim of F7/F8 equivalence, and any universal PIEZO2 mechanism. PIEZO2 remains a secondary context branch.

## One-sentence argument

The existing data justify testing whether a controlled mechanical challenge produces a reproducible fibroblast functional phenotype that depends on integrin/focal-adhesion or Rho/actomyosin signalling, but they do not justify assuming that such a phenotype is fascia-specific or independent of proliferation and generic fibroblast activation.

## Revised hypothesis

In donor-matched human fascia-relevant fibroblast preparations, a prespecified mechanical challenge will alter a functional mechanical phenotype. The effect should be reproducible across independent donors and attenuated by a prespecified integrin/focal-adhesion or Rho/actomyosin perturbation after viability, proliferation, ECM/TGF and inflammatory state are measured. A positive result will support mechanistic plausibility; it will not by itself prove pain relevance or fascia specificity.

## Experimental design

### Biological unit and replication

- Minimum: **3 independent human donors**; preferred: 4–6 donors if recruitment and material allow.
- Donor, not cell, image field or well, is the biological replicate.
- Use one predeclared passage window and record tissue source, donor code, age band, sex if available, passage, culture duration and batch.
- Use at least three technical wells per donor × condition. Technical wells improve precision but do not increase the inferential sample size.
- Randomize donor material across plates and mechanical-condition positions. Use plate/batch as a blocking factor.
- Blind image or morphology scoring to donor and condition labels until QC is complete.

### Minimal confirmatory matrix

Use a 2 × 3 design:

| Mechanical condition | Vehicle | Integrin/focal-adhesion perturbation | Rho/ROCK–actomyosin perturbation |
|---|---:|---:|---:|
| Relaxed/reference | required | required | required |
| Prespecified mechanical challenge | required | required | required |

The mechanical challenge should be one controlled stiffness or loading condition selected before the confirmatory run. If the platform requires a pilot to establish a non-toxic range, the pilot is feasibility work and is not pooled with the confirmatory analysis. Perturbation concentrations, exposure time, solvent, vehicle and target-engagement checks must be fixed before unblinding.

Recommended perturbation logic:

1. One intervention directed at the ITGB1/ITGAV–focal-adhesion axis, with matched vehicle or non-targeting control.
2. One intervention directed at Rho/ROCK or actomyosin contractility, with matched vehicle or non-targeting control.
3. A perturbation is interpretable only if it does not produce unacceptable toxicity or a large nonspecific loss of cell attachment.

Include a matched state-control arm for ECM/TGF or fibroblast activation if material permits. This arm is a confounder control, not a second primary hypothesis.

## Endpoint hierarchy

### Primary endpoint

Select one functional mechanical readout before the confirmatory experiment:

- Preferred: traction-force or force-response measurement.
- Feasible fallback: standardized collagen-gel contraction or mechanically induced morphology index.

The primary quantity is the donor-level mechanical response difference between challenge and relaxed conditions under vehicle. The primary mechanistic contrast is the change in that response under each prespecified perturbation. Do not switch the primary readout after seeing the data.

### Orthogonal and control endpoints

- YAP/TAZ nuclear-to-cytoplasmic localization or an equivalent orthogonal mechanotransduction readout.
- F-actin stress-fibre and focal-adhesion features, using predefined segmentation and scoring rules.
- Viability and cell attachment.
- Proliferation, such as EdU or Ki-67, measured in the same donor and condition structure.
- ECM/TGF/fibrosis and inflammatory state markers, selected before data inspection.
- PIEZO2 expression or functional response as a secondary context-specific observation; positive and null results must both be reported.

The actomyosin-associated transcriptional panel may be measured as a secondary molecular readout, but it cannot replace the primary functional endpoint. Candidate genes should remain separated from module-level interpretation.

## Analysis contract

1. Calculate a donor-level summary for each condition from passing technical wells; do not use cells as independent replicates.
2. Primary contrast: `challenge − relaxed` under vehicle for the selected functional endpoint.
3. Mechanism contrast: compare the primary mechanical response under each perturbation with the vehicle response, preferably as a mechanical-condition × perturbation interaction.
4. Report donor-level effect sizes, 95% confidence intervals, donor-wise directions, technical variation and missingness. Report p-values only alongside these quantities.
5. Use donor as a blocking or random effect when the number of donors supports a mixed model; otherwise use paired donor-level contrasts and label inference as preliminary.
6. Apply multiplicity control to prespecified secondary endpoints. The primary endpoint and the two mechanistic contrasts must remain identifiable and must not be selected after results are viewed.
7. Repeat the primary analysis after excluding predefined QC failures and after a sensitivity adjustment for proliferation or viability. These sensitivity analyses are not permission to remove inconvenient donors.

### Proposed QC and exclusion rules to lock before data collection

- Exclude only documented contamination, instrument failure, failed image acquisition, or a well-level viability/attachment failure defined before unblinding.
- A provisional operational threshold is viability below 70% of the matched reference well; the laboratory must confirm or amend this threshold before the run.
- Retain a donor only when the minimum number of technical wells passes QC in every primary condition; otherwise report the donor as incomplete and provide the prespecified sensitivity analysis.
- Never exclude a donor solely because its direction differs from the group direction.

## Evidence-upgrade gates

The project should move beyond the current CAUTION framing only if all four gates are met:

**Gate A — donor replication.** The primary functional response is directionally concordant in at least 3 independent donors and is not driven by a single donor.

**Gate B — alternative-state control.** The response persists after viability and proliferation checks and is not adequately explained by a generic ECM/TGF or inflammatory-state shift.

**Gate C — perturbation dependence.** At least one prespecified integrin/focal-adhesion or Rho/actomyosin perturbation attenuates the mechanical response with acceptable viability and attachment.

**Gate D — orthogonal convergence.** At least one independent readout, such as YAP/TAZ localization, focal-adhesion remodelling or contractility, changes in a direction consistent with the primary phenotype.

If Gate A passes but B–D do not, retain **CAUTION** and describe a reproducible association only. If A–B pass but C fails, retain the actomyosin/integrin result as biological plausibility without pathway-dependence. A successful functional validation still does not establish pain causality or universal PIEZO2 involvement without additional tissue- and clinical-level evidence.

## Timeline and go/no-go points

- **Week 1:** confirm donor/material availability, platform, perturbation reagents, blinding, randomization and preregistered endpoints.
- **Week 2:** non-confirmatory feasibility pilot to verify mechanical range, imaging quality, viability and target engagement.
- **Weeks 3–6:** confirmatory donor-replicated experiment with locked 2 × 3 matrix and matched controls.
- **Weeks 7–8:** repeat or expand donors only if the predefined QC or replication rule requires it; do not change hypotheses because of one discordant donor.
- **Weeks 9–10:** donor-level analysis, sensitivity audit, evidence decision and manuscript update.

Stop or redesign before confirmatory collection if the mechanical platform cannot produce a stable, non-toxic challenge, if a perturbation has no interpretable target engagement, or if donor material cannot support the minimum replication. Do not replace the primary endpoint with a more favourable molecular marker.

## Risk-response plan

| Risk | Prespecified response |
|---|---|
| No functional mechanical response | Verify platform and non-toxic range in the pilot; if confirmed, report a null result and downgrade the mechanistic candidate. |
| Response tracks proliferation or viability | Strengthen quiescence/state controls, perform the locked sensitivity analysis, and do not claim pathway specificity. |
| Strong donor heterogeneity | Expand toward 6 donors if feasible, report donor-wise effects, and avoid pooling away discordance. |
| Perturbation causes nonspecific toxicity | Run a dose/target-engagement feasibility check and retain only interpretable doses; do not treat cytotoxicity as pathway dependence. |
| PIEZO2 is absent or inconsistent | Report it as a null or context-specific secondary result; the primary program is not failed by this alone. |
| Equipment or reagent delay | Preserve the computational CAUTION claim set and revise the timeline; do not substitute an unvalidated endpoint without a documented amendment. |

## Expected outputs

- Preregistered or internally signed protocol and randomization record.
- Donor- and well-level raw measurement table with QC flags.
- Primary functional endpoint analysis and mechanistic interaction analysis.
- Orthogonal readout and competing-state control summaries.
- Gate A–D decision record and updated evidence grade.

## Material Passport

- Inputs: Step 12B research-direction lock; Step 12C functional-validation specification; Step 15E evidence resynthesis; Step 15D sensitivity and negative-control audit.
- Transformation: computational claim boundaries converted into a feasible, donor-level functional-validation protocol.
- New data generated: none.
- Current status: **ready for laboratory/platform feasibility confirmation; not yet ready for unamended confirmatory data collection**.
