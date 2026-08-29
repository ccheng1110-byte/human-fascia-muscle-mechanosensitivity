## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: run
- Origin Date: 2026-08-23
- Verification Status: ANALYZED
- Version Label: mechanosensitivity_scoring_v1

## GSE173252 mechanosensitivity scoring result

### Execution status

- Status: completed
- Singlets analyzed: 16662
- Predicted doublets excluded: 326
- Frozen gene-set rows: 182
- Gene-set rows absent from this matrix: 1
- Random seed: 20260823

### Paired DD Myofib versus FB module evidence

- ecm_remodeling: median paired Myofib-FB score difference = 0.446; positive in 100% of DD samples
- actomyosin_rho: median paired Myofib-FB score difference = 0.293; positive in 100% of DD samples
- cell_cycle: median paired Myofib-FB score difference = 0.189; positive in 100% of DD samples
- integrin_focal_adhesion: median paired Myofib-FB score difference = 0.100; positive in 100% of DD samples
- primary_mechanotransduction_composite: median paired Myofib-FB score difference = 0.094; positive in 100% of DD samples
- hypoxia: median paired Myofib-FB score difference = 0.047; positive in 67% of DD samples
- tgf_fibrosis: median paired Myofib-FB score difference = 0.007; positive in 67% of DD samples
- hippo_yap_taz: median paired Myofib-FB score difference = -0.089; positive in 33% of DD samples
- mechanosensor_channels: median paired Myofib-FB score difference = -0.119; positive in 0% of DD samples
- inflammation_ap1_nfkb: median paired Myofib-FB score difference = -0.943; positive in 0% of DD samples

### Top tier-A primary mechanotransduction candidates

- CNN1 [actomyosin_rho]: median Myofib-FB log2CPM difference = 3.924; expressing Myofib cells = 23.7%
- ITGA8 [integrin_focal_adhesion]: median Myofib-FB log2CPM difference = 3.144; expressing Myofib cells = 41.7%
- TAGLN [actomyosin_rho]: median Myofib-FB log2CPM difference = 2.917; expressing Myofib cells = 87.7%
- ACTA2 [actomyosin_rho]: median Myofib-FB log2CPM difference = 2.396; expressing Myofib cells = 62.5%
- TPM1 [actomyosin_rho]: median Myofib-FB log2CPM difference = 2.064; expressing Myofib cells = 94.5%
- TPM2 [actomyosin_rho]: median Myofib-FB log2CPM difference = 2.038; expressing Myofib cells = 87.7%
- ITGB1 [integrin_focal_adhesion]: median Myofib-FB log2CPM difference = 1.973; expressing Myofib cells = 98.0%
- ITGAV [integrin_focal_adhesion]: median Myofib-FB log2CPM difference = 1.336; expressing Myofib cells = 75.0%
- BCAR1 [integrin_focal_adhesion]: median Myofib-FB log2CPM difference = 1.253; expressing Myofib cells = 44.3%
- PIEZO2 [mechanosensor_channels]: median Myofib-FB log2CPM difference = 1.236; expressing Myofib cells = 27.0%
- TEAD2 [hippo_yap_taz]: median Myofib-FB log2CPM difference = 1.219; expressing Myofib cells = 29.9%
- LIMK1 [actomyosin_rho]: median Myofib-FB log2CPM difference = 1.127; expressing Myofib cells = 22.7%
- CNN2 [actomyosin_rho]: median Myofib-FB log2CPM difference = 1.025; expressing Myofib cells = 47.1%
- TEAD3 [hippo_yap_taz]: median Myofib-FB log2CPM difference = 1.012; expressing Myofib cells = 14.0%
- PFN1 [actomyosin_rho]: median Myofib-FB log2CPM difference = 0.872; expressing Myofib cells = 97.1%

### Interpretation boundary

- Scores represent transcriptional programs, not direct measurements of cellular mechanical sensitivity.
- DD1-DD3 are the only samples with a sufficiently represented Myofib population; therefore the paired analysis has n = 3 biological samples.
- Effect size and direction consistency are reported. No inferential p-value is claimed because n = 3 cannot support a stable formal test.
- ECM, TGF/fibrosis, inflammation, hypoxia and cell-cycle modules are explicit competing explanations and must be interpreted alongside the four primary modules.
- A Myofib-enriched transcript is a candidate marker, not evidence that the gene causally mediates mechanosensitivity.

### Frozen source anchors

- GO/MSigDB response to mechanical stimulus: GO:0009612
- Reactome Integrin cell surface interactions: R-HSA-216083
- Reactome shear-stress response: R-HSA-9860931
- Reactome YAP1/WWTR1 transcription: R-HSA-2032785
- Reactome RHO-ROCK signaling: R-HSA-5627117
- Reactome extracellular matrix organization: R-HSA-1474244
- Reactome TGF-beta receptor signaling: R-HSA-170834
