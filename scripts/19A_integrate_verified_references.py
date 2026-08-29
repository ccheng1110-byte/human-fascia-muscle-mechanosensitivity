from pathlib import Path
import csv


PROJECT = Path(".")
INPUT_MD = PROJECT / "results/11_manuscript_preparation/18B_second_round_manuscript_QA/18B_second_round_QA_revised_manuscript_v3.md"
OUTPUT_DIR = PROJECT / "results/11_manuscript_preparation/19A_reference_integrated_manuscript"
OUTPUT_MD = OUTPUT_DIR / "19A_reference_integrated_manuscript_v4.md"
REGISTRY_CSV = OUTPUT_DIR / "19A_reference_registry_v1.csv"
MAP_CSV = OUTPUT_DIR / "19A_citation_insertion_map_v1.csv"
DECISION_MD = OUTPUT_DIR / "19A_reference_integration_decision_v1.md"


REFERENCES = [
    {
        "number": 1,
        "citation": "Nardone G, Oliver-De La Cruz J, Vrbsky J, et al. YAP regulates cell mechanics by controlling focal adhesion assembly. Nature Communications. 2017;8:15321.",
        "doi": "10.1038/ncomms15321",
        "url": "https://doi.org/10.1038/ncomms15321",
        "role": "YAP–focal-adhesion mechanobiology",
        "verification": "Official Nature article / DOI",
    },
    {
        "number": 2,
        "citation": "Ren B, et al. RAP2 mediates mechanoresponses of the Hippo pathway. Nature. 2018;560:655-660.",
        "doi": "10.1038/s41586-018-0444-0",
        "url": "https://doi.org/10.1038/s41586-018-0444-0",
        "role": "Hippo-pathway mechanotransduction",
        "verification": "DOI / published article record",
    },
    {
        "number": 3,
        "citation": "Oakes PW, Wagner E, Brand CA, et al. Optogenetic control of RhoA reveals zyxin-mediated elasticity of stress fibres. Nature Communications. 2017;8:15817.",
        "doi": "10.1038/ncomms15817",
        "url": "https://doi.org/10.1038/ncomms15817",
        "role": "RhoA–actomyosin force transmission",
        "verification": "Official Nature article / DOI",
    },
    {
        "number": 4,
        "citation": "Hong SP, Yang MJ, Cho H, et al. Distinct fibroblast subsets regulate lacteal integrity through YAP/TAZ-induced VEGF-C in intestinal villi. Nature Communications. 2020;11:4102.",
        "doi": "10.1038/s41467-020-17886-y",
        "url": "https://doi.org/10.1038/s41467-020-17886-y",
        "role": "Fibroblast state and YAP/TAZ context",
        "verification": "Official Nature article / DOI",
    },
    {
        "number": 5,
        "citation": "Yuan Y, Li M, Yao Q, Qi Y, Ke B. Piezo2 mechanically regulate scleral fibroblast differentiation by activating relA/RhoA pathway. Experimental Eye Research. 2026;267:110955.",
        "doi": "10.1016/j.exer.2026.110955",
        "url": "https://doi.org/10.1016/j.exer.2026.110955",
        "role": "Context-specific PIEZO2 evidence",
        "verification": "DOI / published article record",
    },
    {
        "number": 6,
        "citation": "Chen H, Qu J, Huang X, et al. Mechanosensing by the alpha6-integrin confers an invasive fibroblast phenotype and mediates lung fibrosis. Nature Communications. 2016;7:12564.",
        "doi": "10.1038/ncomms12564",
        "url": "https://doi.org/10.1038/ncomms12564",
        "role": "Matrix stiffness, integrin, fibroblast phenotype",
        "verification": "PubMed PMID 27535718 / PMC4992155 / DOI",
    },
    {
        "number": 7,
        "citation": "Zhou DW, Fernandez-Yague MA, Holland EN, et al. Force-FAK signaling coupling at individual focal adhesions coordinates mechanosensing and microtissue repair. Nature Communications. 2021;12:2359.",
        "doi": "10.1038/s41467-021-22602-5",
        "url": "https://doi.org/10.1038/s41467-021-22602-5",
        "role": "Force–FAK coupling at focal adhesions",
        "verification": "PubMed PMID 33883558 / DOI",
    },
    {
        "number": 8,
        "citation": "Shiu JY, Aires L, Lin Z, Vogel V. Nanopillar force measurements reveal actin-cap-mediated YAP mechanotransduction. Nature Cell Biology. 2018;20:262-271.",
        "doi": "10.1038/s41556-017-0030-y",
        "url": "https://doi.org/10.1038/s41556-017-0030-y",
        "role": "Actin-cap, integrin and YAP nuclear signaling",
        "verification": "PubMed PMID 29403039 / official Nature article / DOI",
    },
    {
        "number": 9,
        "citation": "Crowell HL, Soneson C, Germain PL, et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. Nature Communications. 2020;11:6077.",
        "doi": "10.1038/s41467-020-19894-4",
        "url": "https://doi.org/10.1038/s41467-020-19894-4",
        "role": "Sample-level inference in multi-sample single-cell studies",
        "verification": "PubMed PMID 33257685 / PMC7705760 / DOI",
    },
    {
        "number": 10,
        "citation": "Squair JW, Gautier M, Kathe C, et al. Confronting false discoveries in single-cell differential expression. Nature Communications. 2021;12:5692.",
        "doi": "10.1038/s41467-021-25960-2",
        "url": "https://doi.org/10.1038/s41467-021-25960-2",
        "role": "Biological replicate and pseudoreplication control",
        "verification": "PubMed PMID 34584091 / PMC8479118 / DOI",
    },
    {
        "number": 11,
        "citation": "Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Research. 2012;40(17):e133.",
        "doi": "10.1093/nar/gks461",
        "url": "https://doi.org/10.1093/nar/gks461",
        "role": "Competitive gene-set testing and inter-gene correlation",
        "verification": "Official Oxford Academic article / DOI",
    },
    {
        "number": 12,
        "citation": "Dobie R, et al. Deciphering mesenchymal drivers of human Dupuytren's disease at single-cell level. Journal of Investigative Dermatology. 2022;142(1):114-123.e8.",
        "doi": "10.1016/j.jid.2021.05.030",
        "url": "https://doi.org/10.1016/j.jid.2021.05.030",
        "role": "GSE173252 discovery resource",
        "verification": "PubMed PMID 34274346 / GEO GSE173252",
    },
    {
        "number": 13,
        "citation": "Steele L, et al. A single-cell and spatial genomics atlas of human skin fibroblasts reveals shared disease-related fibroblast subtypes across tissues. Nature Immunology. 2025;26(10):1807-1820.",
        "doi": "10.1038/s41590-025-02267-8",
        "url": "https://doi.org/10.1038/s41590-025-02267-8",
        "role": "PRJNA607098 skin fibroblast atlas",
        "verification": "PubMed PMID 40993240 / official Nature article / project GitHub",
    },
    {
        "number": 14,
        "citation": "Tie J, Chen D, Guo J, et al. Transcriptome-wide study of the response of human trabecular meshwork cells to the substrate stiffness increase. Journal of Cellular Biochemistry. 2020;121(5-6):3112-3123.",
        "doi": "10.1002/jcb.29578",
        "url": "https://doi.org/10.1002/jcb.29578",
        "role": "GSE123100 cross-tissue stiffness resource",
        "verification": "PubMed PMID 32115746 / GEO GSE123100",
    },
    {
        "number": 15,
        "citation": "Kaiser AM, Selahi A, Kong W, Ruby JG. Substrate stiffness dictates unique paths towards proliferative arrest in WI-38 cells. GeroScience. 2026;48(3):4173-4193.",
        "doi": "10.1007/s11357-025-01858-5",
        "url": "https://doi.org/10.1007/s11357-025-01858-5",
        "role": "GSE276045 WI-38 stiffness resource",
        "verification": "PMC13356154 / PMID 40974519 / DOI",
    },
]


def replace_once(text, old, new):
    if old not in text:
        raise RuntimeError(f"Expected manuscript text not found: {old[:120]}")
    return text.replace(old, new, 1)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    text = INPUT_MD.read_text(encoding="utf-8")

    text = replace_once(
        text,
        "Focal adhesions connect integrin–ECM interactions to the actin cytoskeleton, while mechanical cues can be transmitted to transcriptional regulators such as YAP and TAZ (Nardone et al., 2017; Ren et al., 2018). In stromal and fibroblast-related contexts, mechanical stretch, matrix stiffness, and cytoskeletal perturbation can alter YAP/TAZ localization and transcriptional outputs (Hong et al., 2020).",
        "Focal adhesions connect integrin–ECM interactions to the actin cytoskeleton, while mechanical cues can be transmitted to transcriptional regulators such as YAP and TAZ (Nardone et al., 2017; Ren et al., 2018; Chen et al., 2016; Zhou et al., 2021). In stromal and fibroblast-related contexts, mechanical stretch, matrix stiffness, and cytoskeletal perturbation can alter YAP/TAZ localization and transcriptional outputs (Hong et al., 2020; Shiu et al., 2018).",
    )
    text = replace_once(
        text,
        "The discovery resource was GSE173252, for which a processed mesenchyme object and supplementary files were available locally under the project data directory. The primary external validation used the PRJNA607098 representation of the 2025 human skin fibroblast atlas.",
        "The discovery resource was GSE173252, for which a processed mesenchyme object and supplementary files were available locally under the project data directory (Dobie et al., 2022). The primary external validation used the PRJNA607098 representation of the 2025 human skin fibroblast atlas (Steele et al., 2025).",
    )
    text = replace_once(
        text,
        "Module-level summaries and candidate-gene evidence were aggregated at the sample/cluster level rather than treating individual cells as independent biological replicates.",
        "Module-level summaries and candidate-gene evidence were aggregated at the sample/cluster level rather than treating individual cells as independent biological replicates, consistent with established guidance for replicated multi-sample single-cell inference (Crowell et al., 2020; Squair et al., 2021).",
    )
    text = replace_once(
        text,
        "This separation was maintained because coordinated module positivity does not imply that every individual candidate gene is conserved across tissues or datasets.",
        "This separation was maintained because coordinated module positivity does not imply that every individual candidate gene is conserved across tissues or datasets. Competitive module testing used CAMERA, which accounts for inter-gene correlation in competitive gene-set inference (Wu and Smyth, 2012).",
    )
    text = replace_once(
        text,
        "For S3, GSE123100 was restricted to a descriptive cultured HTM stiffness-form analysis, whereas GSE276045 was modeled as a WI-38 stiffness slope with WT/hTERT cell-model interaction and timepoint adjustment.",
        "For S3, GSE123100 was restricted to a descriptive cultured HTM stiffness-form analysis, whereas GSE276045 was modeled as a WI-38 stiffness slope with WT/hTERT cell-model interaction and timepoint adjustment; these source studies provide cross-tissue context rather than fascia-specific replication (Tie et al., 2020; Kaiser et al., 2026).",
    )
    text = replace_once(
        text,
        "This interpretation is consistent with the biological expectation that mechanical information is distributed across matrix attachment, force transmission, cytoskeletal state, and downstream transcriptional responses (Nardone et al., 2017; Ren et al., 2018).",
        "This interpretation is consistent with the biological expectation that mechanical information is distributed across matrix attachment, force transmission, cytoskeletal state, and downstream transcriptional responses (Nardone et al., 2017; Ren et al., 2018; Shiu et al., 2018).",
    )
    text = replace_once(
        text,
        "This concordance makes integrin-linked force transmission a rational priority for functional testing.",
        "This concordance makes integrin-linked force transmission a rational priority for functional testing, in line with experimental evidence that integrin signaling can couple matrix stiffness to fibroblast phenotypes and fibrosis-related behavior (Chen et al., 2016; Zhou et al., 2021).",
    )
    text = replace_once(
        text,
        "The WI-38 stiffness analysis is informative because actomyosin/Rho direction was supported in both WT and hTERT cell-model slopes after timepoint adjustment.",
        "The WI-38 stiffness analysis is informative because actomyosin/Rho direction was supported in both WT and hTERT cell-model slopes after timepoint adjustment, consistent with the source study's demonstration that the mechanical environment can redirect WI-38 proliferative and transcriptional trajectories (Kaiser et al., 2026).",
    )
    text = replace_once(
        text,
        "The HTM analysis provides a complementary dose-response form but is descriptive and has incomplete frozen-candidate coverage.",
        "The HTM analysis provides a complementary dose-response form but is descriptive and has incomplete frozen-candidate coverage; the original HTM study likewise reported stiffness-associated transcriptomic changes in an ocular connective-tissue model rather than fascia (Tie et al., 2020).",
    )

    ref_start = text.index("## References")
    references = ["## References", ""]
    for ref in REFERENCES:
        references.append(f"{ref['number']}. {ref['citation']} DOI: {ref['doi']}")
    references.append("")
    text = text[:ref_start] + "\n".join(references)
    OUTPUT_MD.write_text(text, encoding="utf-8")

    with REGISTRY_CSV.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["number", "citation", "doi", "url", "role", "verification"])
        writer.writeheader()
        writer.writerows(REFERENCES)

    insertion_rows = [
        ["Introduction", "ECM/integrin–cytoskeletal and YAP/TAZ mechanotransduction rationale", "1;2;6;7;8", "Mechanistic context", "Background evidence only; not causal evidence for the present reanalysis"],
        ["Introduction", "Fibroblast state and context-dependent PIEZO2 rationale", "4;5", "Fibroblast-state and channel context", "PIEZO2 is retained as context-specific"],
        ["Methods: Data sources and provenance audit", "Discovery and primary atlas source provenance", "12;13", "Dataset-source attribution", "Source citation does not establish label equivalence"],
        ["Methods: Candidate scoring", "Sample-level inference and biological replicate discipline", "9;10", "Statistical rationale", "Supports unit-of-inference choice; does not change project gates"],
        ["Methods: Candidate scoring", "CAMERA competitive module testing", "11", "Method citation", "Supports inter-gene-correlation-aware competitive testing"],
        ["Methods: Second-round strengthening", "Cross-tissue stiffness source context", "14;15", "Dataset-source attribution", "Cross-tissue plausibility only; no fascia direct replication"],
        ["Discussion: Multicomponent program", "Distributed force transmission and YAP readout", "1;2;8", "Mechanistic context", "Compatible with the hypothesis, not proof of causality"],
        ["Discussion: Integrin support", "Integrin-linked stiffness and focal-adhesion force transmission", "6;7", "Mechanistic context", "Rationalizes functional prioritization; specificity remains unresolved"],
        ["Discussion: Second-round evidence", "WI-38 and HTM cross-tissue interpretation", "14;15", "Dataset-source attribution", "The manuscript retains the cross-tissue boundary"],
    ]
    with MAP_CSV.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["manuscript_location", "claim_or_topic", "references", "support_role", "interpretation_boundary"])
        writer.writerows(insertion_rows)

    decision = """# Reference integration decision\n\n- The manuscript was expanded from 5 to 15 references.\n- Ten additions were selected from primary research or established statistical-method papers and checked against DOI, PubMed, PMC, official publisher, GEO, or project-source records.\n- Citations were inserted only where they support background mechanism, statistical design, or source-data provenance.\n- The reference list does not add a formal publication for GSE338388 because a source paper was not verified; the accession remains explicitly identified in the manuscript.\n- The evidence grade remains CAUTION. Additional references improve context and traceability but do not replace donor-replicated functional perturbation.\n\nOutputs:\n- `19A_reference_integrated_manuscript_v4.md`\n- `19A_reference_registry_v1.csv`\n- `19A_citation_insertion_map_v1.csv`\n"""
    DECISION_MD.write_text(decision, encoding="utf-8")
    print(f"Reference integration completed: {OUTPUT_MD}")
    print(f"Reference count: {len(REFERENCES)}")
    print(f"Registry: {REGISTRY_CSV}")
    print(f"Citation map: {MAP_CSV}")


if __name__ == "__main__":
    main()
