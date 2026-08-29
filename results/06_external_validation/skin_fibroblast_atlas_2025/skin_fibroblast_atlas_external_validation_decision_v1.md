## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-23
- Verification Status: SCREENED
- Version Label: skin_fibroblast_atlas_marker_rank_screen_v1

## Step 08A external atlas marker-rank screen

### Frozen question

Do the 15 frozen Step-07 A/B candidates rank preferentially in the authors' F7 fascia-like myofibroblast marker list?

### Screening rule

Exploratory support requires F7 rank <= 3191 (top 10%) and a better F7 rank than the median rank across non-F7 states.
No cell-level p values are used and no new differential-expression test is run.

### Result

- Supported candidates: 11/15.
- Supported genes: ITGB1, ITGAV, CFL1, FERMT2, ITGA5, PIEZO2, TMEM63B, PARVA, CDC42, LIMK1, NF2.
- PIEZO2 F7 marker rank: 435 of 31908.
- Independent DD source PRJNA607098 listed by the atlas: TRUE.
- Discovery source GSE173252 also listed by the atlas: TRUE.

### Evidence decision

- Evidence grade remains CAUTION.
- This marker list was calculated for atlas cell states pooled across sources; it does not isolate PRJNA607098 from GSE173252.
- The public webportal object exposes dataset and cell-state labels but no sample/donor label.
- Therefore this step is an external atlas screen, not donor-level independent validation.
- A positive screen justifies Step 08B: source-specific PRJNA607098 expression extraction.
- Formal evidence upgrading requires Step 08C donor/sample-level pseudobulk or experimental validation.

### Mandatory interpretation limits

- Marker rank does not establish mechanosensor activity, protein abundance or causality.
- Atlas integration and cell-state annotation can induce source-dependent classification effects.
- The F7 state may reflect fascia/myofibroblast identity rather than disease-specific force sensing.
- The 50 MB automatic-download ceiling remains in force.
