## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-23
- Verification Status: ANALYZED
- Version Label: mechanosensitivity_specificity_audit_v1

## GSE173252 mechanosensitivity specificity and donor-robustness audit

### Audit design

- Biological unit: DD donor/sample (DD1-DD3).
- Primary comparison: Myofib versus FB, VSMC and Pericyte separately within each DD sample.
- Robust specificity requires a positive Myofib-comparator log2CPM difference for all three identities in all three DD samples.
- Leave-one-donor-out analysis repeats direction checks after omitting DD1, DD2 or DD3.
- No cell-level p values are calculated.

### Module-level identity specificity

- ecm_remodeling: Myofib_dominant_all_DD; median Myofib rank = 1; minimum gap from best other identity = 0.256
- cell_cycle: Myofib_dominant_two_of_three_DD; median Myofib rank = 1; minimum gap from best other identity = -0.599
- integrin_focal_adhesion: Myofib_dominant_all_DD; median Myofib rank = 1; minimum gap from best other identity = 0.059
- hypoxia: Not_Myofib_dominant; median Myofib rank = 2; minimum gap from best other identity = -0.019
- primary_mechanotransduction_composite: Not_Myofib_dominant; median Myofib rank = 2; minimum gap from best other identity = -0.035
- tgf_fibrosis: Not_Myofib_dominant; median Myofib rank = 2; minimum gap from best other identity = -0.135
- mechanosensor_channels: Not_Myofib_dominant; median Myofib rank = 2; minimum gap from best other identity = -0.199
- actomyosin_rho: Not_Myofib_dominant; median Myofib rank = 3; minimum gap from best other identity = -0.663
- hippo_yap_taz: Not_Myofib_dominant; median Myofib rank = 4; minimum gap from best other identity = -0.349
- inflammation_ap1_nfkb: Not_Myofib_dominant; median Myofib rank = 4; minimum gap from best other identity = -1.283

### Step-06 candidates passing all-identity specificity

- ITGAV [integrin_focal_adhesion]: minimum pairwise log2CPM difference = 1.318; leave-one-donor-out all-pairs positive = TRUE
- ITGB1 [integrin_focal_adhesion]: minimum pairwise log2CPM difference = 1.002; leave-one-donor-out all-pairs positive = TRUE
- ITGA2 [integrin_focal_adhesion]: minimum pairwise log2CPM difference = 0.870; leave-one-donor-out all-pairs positive = TRUE
- PIEZO2 [mechanosensor_channels]: minimum pairwise log2CPM difference = 0.728; leave-one-donor-out all-pairs positive = TRUE
- PARVA [integrin_focal_adhesion]: minimum pairwise log2CPM difference = 0.655; leave-one-donor-out all-pairs positive = TRUE
- NF2 [hippo_yap_taz]: minimum pairwise log2CPM difference = 0.533; leave-one-donor-out all-pairs positive = TRUE
- ITGA5 [integrin_focal_adhesion]: minimum pairwise log2CPM difference = 0.514; leave-one-donor-out all-pairs positive = TRUE
- CFL1 [actomyosin_rho]: minimum pairwise log2CPM difference = 0.423; leave-one-donor-out all-pairs positive = TRUE
- LIMK1 [actomyosin_rho]: minimum pairwise log2CPM difference = 0.235; leave-one-donor-out all-pairs positive = TRUE
- TMEM63B [mechanosensor_channels]: minimum pairwise log2CPM difference = 0.225; leave-one-donor-out all-pairs positive = TRUE
- PANX1 [mechanosensor_channels]: minimum pairwise log2CPM difference = 0.129; leave-one-donor-out all-pairs positive = TRUE
- CNN2 [actomyosin_rho]: minimum pairwise log2CPM difference = 0.125; leave-one-donor-out all-pairs positive = TRUE
- CDC42 [actomyosin_rho]: minimum pairwise log2CPM difference = 0.109; leave-one-donor-out all-pairs positive = TRUE
- TEAD1 [hippo_yap_taz]: minimum pairwise log2CPM difference = 0.088; leave-one-donor-out all-pairs positive = TRUE
- FERMT2 [integrin_focal_adhesion]: minimum pairwise log2CPM difference = 0.007; leave-one-donor-out all-pairs positive = TRUE

### Focus mechanosensors

- PIEZO2: robust_Myofib_over_all_three_identities; expressing Myofib cells = 27.0%; median Myofib-FB difference = 1.236; minimum difference over all identities/donors = 0.728
- TRPV4: robust_Myofib_over_all_three_identities; expressing Myofib cells = 2.8%; median Myofib-FB difference = 0.630; minimum difference over all identities/donors = 0.619
- TMEM63B: robust_Myofib_over_all_three_identities; expressing Myofib cells = 17.5%; median Myofib-FB difference = 0.781; minimum difference over all identities/donors = 0.225
- PKD2: robust_over_FB_but_not_all_identities; expressing Myofib cells = 35.8%; median Myofib-FB difference = 0.629; minimum difference over all identities/donors = -0.380
- PIEZO1: donor_or_identity_dependent; expressing Myofib cells = 8.3%; median Myofib-FB difference = -0.834; minimum difference over all identities/donors = -1.102

### Mandatory interpretation limits

- Myofib specificity does not establish disease specificity because controls do not contain a comparable Myofib population.
- Transcriptional enrichment does not demonstrate force sensing, channel activity, protein abundance or causal function.
- Contractile genes such as ACTA2, TAGLN and CNN1 may describe myofibroblast differentiation rather than a distinct mechanosensitive state.
- PIEZO2 and other low-to-moderate prevalence genes require donor-level prevalence checks and independent dataset validation.
- With only three DD donors, leave-one-out analysis is a fragility diagnostic rather than formal statistical validation.
