"""Prepare a BMC Genomics Research article manuscript from the style-polished draft.

This script performs format-level transformations only:
* converts author--year citations to Vancouver-style numbered citations;
* reorders the verified reference list by first citation;
* applies the BMC Genomics Research article section structure;
* adds title-page/declarations placeholders where author information is unavailable.

It does not change reported numerical results, evidence grades or scientific claims.
"""

from __future__ import annotations

import re
import unicodedata
from pathlib import Path


PROJECT = Path(r".")
INPUT_MD = PROJECT / "results/11_manuscript_preparation/22B_style_polished_manuscript_v7.md"
OUTPUT_DIR = PROJECT / "results/11_manuscript_preparation/23A_BMC_Genomics_format"
OUTPUT_MD = OUTPUT_DIR / "23A_BMC_Genomics_manuscript_v1.md"
AUDIT_MD = OUTPUT_DIR / "23A_BMC_Genomics_format_audit_v1.md"


FIGURE_INSERTIONS = {
    "### Provenance and primary validation structure": {
        "path": "../20B_reviewer_response_figures/20B_Figure1_evidence_workflow.png",
        "caption": (
            "Figure 1. Evidence workflow and inferential units. The staged design separates "
            "discovery, paired atlas validation, independent subject-level triangulation, "
            "direct mechanical perturbation, regulatory/stiffness context and bounded synthesis. "
            "Each source retains its own biological or inferential denominator."
        ),
    },
    "### Core modules were consistently positive in PRJNA607098": {
        "path": "../20B_reviewer_response_figures/20B_Figure2_sample_subject_module_contrasts.png",
        "caption": (
            "Figure 2. Primary-module contrasts at the correct inferential unit. PRJNA607098 "
            "points are 12 reconstructed SRS sample units comparing F7 with pooled non-F7 cells "
            "within sample. GSE130973 points are five subjects comparing marker-defined candidate "
            "fibroblast clusters with other clusters within subject. Black summaries are medians "
            "and IQRs; source-specific free y-axis scales are intentional."
        ),
    },
    "### Competition analysis limited integrin specificity": {
        "path": "../20B_reviewer_response_figures/20B_Figure4_competition_residuals.png",
        "caption": (
            "Figure 3. Competition residuals expose the specificity boundary. Within PRJNA607098, "
            "each core-module contrast was regressed on each competitor-module contrast. The "
            "prespecified primary gate requires at least 8 of 12 residuals above zero; a secondary "
            "sensitivity summary requires at least 9 of 12. Integrin versus ECM remodelling and "
            "three actomyosin comparisons meet only the minimum 8/12 boundary."
        ),
    },
    "### Independent skin data supported the modules but not candidate uniformity": {
        "path": "../20B_reviewer_response_figures/20B_Figure3_frozen_candidate_raw_contrasts.png",
        "caption": (
            "Figure 4. Frozen candidate effects are source dependent. All frozen candidates are "
            "shown without post hoc removal. Points are sample- or subject-level target-minus-" 
            "comparator contrasts; black diamonds and bars show median and IQR. Source separation "
            "preserves distinct expression scales and denominators."
        ),
    },
    "### Direct mechanical tension supported actomyosin/Rho but also cell cycle": {
        "path": "../20B_reviewer_response_figures/20B_Figure5_GSE300230_tension_response.png",
        "caption": (
            "Figure 5. Direct mechanical tension activates actomyosin and cell-cycle programmes. "
            "Model-estimated module-score contrasts compare tension with relaxed collagen without "
            "TGF-beta in each GSE300230 cell line. Labels give CAMERA BH q values. The concurrent "
            "cell-cycle response is an interpretation boundary; the two cell lines are not "
            "independent donors or age replicates."
        ),
    },
    "### Module-level support was more reproducible than candidate-gene support": {
        "path": "../21B_evidence_hierarchy_figure/21B_Figure6_evidence_hierarchy.png",
        "caption": (
            "Figure 6. Module-level support is more reproducible than individual-gene support "
            "across sources. The matrix summarises the frozen evidence hierarchy across primary "
            "paired validation, independent skin triangulation, direct tension, TEAD-linked "
            "regulatory cross-validation and cross-tissue stiffness analysis. Partial or discordant "
            "cells retain co-active competitor programmes, limited replication and source-dependent "
            "candidates. Colours are descriptive evidence classes, not a universal grading scale."
        ),
    },
}


# BMC Genomics follows Index Medicus/MEDLINE-style journal abbreviations in
# Vancouver references. The bibliographic fields were already verified before
# this format-only conversion; this map changes presentation, not identity.
JOURNAL_ABBREVIATIONS = {
    "British Journal of Sports Medicine": "Br J Sports Med",
    "Histochemistry and Cell Biology": "Histochem Cell Biol",
    "Journal of Cellular Physiology": "J Cell Physiol",
    "Nature Reviews Molecular Cell Biology": "Nat Rev Mol Cell Biol",
    "Journal of Cell Biology": "J Cell Biol",
    "Nature Communications": "Nat Commun",
    "Nature Cell Biology": "Nat Cell Biol",
    "Nature Medicine": "Nat Med",
    "Nature Immunology": "Nat Immunol",
    "Journal of Biomechanics": "J Biomech",
    "American Journal of Physiology-Lung Cellular and Molecular Physiology": "Am J Physiol Lung Cell Mol Physiol",
    "Experimental Eye Research": "Exp Eye Res",
    "Journal of Investigative Dermatology": "J Invest Dermatol",
    "Communications Biology": "Commun Biol",
    "Proceedings of the National Academy of Sciences of the United States of America": "Proc Natl Acad Sci U S A",
    "Nucleic Acids Research": "Nucleic Acids Res",
    "Journal of Cellular Biochemistry": "J Cell Biochem",
}


def normalise_key(value: str) -> str:
    value = unicodedata.normalize("NFKD", value)
    value = "".join(ch for ch in value if not unicodedata.combining(ch))
    return re.sub(r"[^A-Za-z0-9]", "", value).casefold()


def split_reference_lines(reference_block: str) -> list[str]:
    return [
        line.strip()
        for line in reference_block.splitlines()
        if re.match(r"^\d+\.\s+", line.strip())
    ]


def reference_key(reference: str) -> tuple[str, int]:
    numberless = re.sub(r"^\d+\.\s+", "", reference)
    author_part = numberless.split(". ", 1)[0]
    first_author = author_part.split()[0]
    year_match = re.search(r"\b(20\d{2})\b", numberless)
    if not year_match:
        raise ValueError(f"Could not identify publication year in reference: {reference}")
    return normalise_key(first_author), int(year_match.group(1))


def apply_journal_abbreviation(reference: str) -> str:
    formatted = reference
    for full_name, abbreviation in JOURNAL_ABBREVIATIONS.items():
        formatted = formatted.replace(f"*{full_name}*", f"*{abbreviation}*")
    return formatted


def parse_citation_group(group: str) -> list[tuple[str, int]]:
    """Parse a semicolon-separated author-year parenthetical citation."""
    keys: list[tuple[str, int]] = []
    for part in group.split(";"):
        years = re.findall(r"\b(20\d{2})\b", part)
        if not years:
            continue
        first_year_position = re.search(r"\b20\d{2}\b", part)
        assert first_year_position is not None
        author_text = part[: first_year_position.start()].strip(" ,")
        if not author_text:
            continue
        surname = author_text.split()[0]
        for year in years:
            keys.append((normalise_key(surname), int(year)))
    return keys


def convert_citations(body: str, ref_lookup: dict[tuple[str, int], str]) -> tuple[str, list[tuple[str, int]], list[str]]:
    citation_pattern = re.compile(r"\(([^()\n]*\b20\d{2}\b[^()\n]*)\)")
    citation_order: list[tuple[str, int]] = []
    unresolved: list[str] = []

    def replace(match: re.Match[str]) -> str:
        group = match.group(1)
        keys = parse_citation_group(group)
        if not keys:
            return match.group(0)
        numbers: list[int] = []
        for key in keys:
            if key not in ref_lookup:
                unresolved.append(group)
                continue
            if key not in citation_order:
                citation_order.append(key)
            numbers.append(citation_order.index(key) + 1)
        if not numbers:
            return match.group(0)
        unique_numbers = list(dict.fromkeys(numbers))
        return "[" + ",".join(str(number) for number in unique_numbers) + "]"

    return citation_pattern.sub(replace, body), citation_order, unresolved


def reorder_references(reference_lines: list[str], citation_order: list[tuple[str, int]]) -> tuple[list[str], list[tuple[str, int]]]:
    by_key: dict[tuple[str, int], str] = {}
    for reference in reference_lines:
        key = reference_key(reference)
        if key in by_key:
            raise ValueError(f"Duplicate reference key: {key}")
        by_key[key] = reference

    ordered_keys = list(citation_order)
    uncited = [reference_key(reference) for reference in reference_lines if reference_key(reference) not in ordered_keys]
    ordered_keys.extend(uncited)

    output: list[str] = []
    for number, key in enumerate(ordered_keys, start=1):
        reference = by_key[key]
        text = re.sub(r"^\d+\.\s+", "", reference)
        text = apply_journal_abbreviation(text)
        output.append(f"{number}. {text}")
    return output, uncited


def extract_abstract(text: str) -> tuple[dict[str, str], int, int]:
    abstract_start = text.index("## Abstract")
    introduction_start = text.index("## Introduction")
    block = text[abstract_start:introduction_start]
    matches = list(re.finditer(r"^###\s+(.+?)\s*$", block, flags=re.MULTILINE))
    sections: dict[str, str] = {}
    for index, match in enumerate(matches):
        start = match.end()
        end = matches[index + 1].start() if index + 1 < len(matches) else len(block)
        sections[match.group(1).strip()] = block[start:end].strip()
    return sections, abstract_start, introduction_start


def add_figure_markers(body: str) -> str:
    lines = body.splitlines()
    output: list[str] = []
    current_heading: str | None = None

    def flush_figure() -> None:
        if current_heading not in FIGURE_INSERTIONS:
            return
        info = FIGURE_INSERTIONS[current_heading]
        output.extend(["", f"![{current_heading[4:]}]({info['path']})", f"**{info['caption']}**", ""])

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("### ") or stripped.startswith("## "):
            flush_figure()
            current_heading = stripped
        output.append(line)
    flush_figure()
    return "\n".join(output).rstrip() + "\n"


def build_declarations(data_text: str, ethics_text: str) -> str:
    return f"""## Declarations

### Ethics approval and consent to participate

{ethics_text.strip()}

This secondary analysis used only de-identified public datasets and recruited no new participants or human specimens. Dataset-specific original ethics committee names, approval identifiers and consent wording should be added here if required by the source records or the journal.

### Consent for publication

Not applicable. No identifiable individual-level information is presented.

### Availability of data and materials

{data_text.strip()}

**Author action required before submission:** add the versioned public repository URL and archival identifier for the analysis code, frozen registries, source inventories, decision records and machine-readable result tables. Publicly available datasets should be cited with their persistent accessions or identifiers in accordance with BMC Genomics policy.

### Competing interests

[AUTHOR TO COMPLETE: state all financial and non-financial competing interests, or write “The authors declare that they have no competing interests.”]

### Funding

[AUTHOR TO COMPLETE: list all sources of funding and describe the role of each funder, if applicable.]

### Authors' contributions

[AUTHOR TO COMPLETE: describe each author's contribution using author initials.]

### Acknowledgements

[AUTHOR TO COMPLETE, or write “Not applicable.” Obtain permission from anyone named in this section.]
"""


def build_markdown() -> tuple[str, dict[str, object]]:
    text = INPUT_MD.read_text(encoding="utf-8")
    title = re.search(r"^#\s+(.+?)\s*$", text, flags=re.MULTILINE)
    if not title:
        raise ValueError("Title not found")

    abstract_sections, _, introduction_start = extract_abstract(text)
    required_abstract_sections = {"Background", "Methods", "Results", "Conclusions"}
    missing = required_abstract_sections.difference(abstract_sections)
    if missing:
        raise ValueError(f"Missing abstract section(s): {sorted(missing)}")

    references_start = text.index("## References")
    core_source = text[introduction_start:references_start]
    data_match = re.search(r"^## Data and code availability\s*$\n(.*?)(?=^## Ethics statement\s*$)", core_source, flags=re.MULTILINE | re.DOTALL)
    ethics_match = re.search(r"^## Ethics statement\s*$\n(.*)$", core_source, flags=re.MULTILINE | re.DOTALL)
    if not data_match or not ethics_match:
        raise ValueError("Data availability or ethics section not found")
    scientific_body = core_source[: data_match.start()].rstrip()
    data_text = data_match.group(1).strip()
    ethics_text = ethics_match.group(1).strip()

    reference_lines = split_reference_lines(text[references_start + len("## References"):])
    ref_lookup = {reference_key(reference): reference for reference in reference_lines}
    scientific_body, citation_order, unresolved = convert_citations(scientific_body, ref_lookup)
    if unresolved:
        raise ValueError(f"Unresolved author-year citation(s): {sorted(set(unresolved))}")
    ordered_references, uncited = reorder_references(reference_lines, citation_order)

    scientific_body = scientific_body.replace("## Introduction", "## Background", 1)
    scientific_body = scientific_body.replace("## Conclusion", "## Conclusions", 1)
    # The source draft mentions the competition panel before the independent-skin
    # panel. Renumber those two figures so first mention order is 1, 2, 3, 4, 5, 6.
    scientific_body = scientific_body.replace("(Fig. 4)", "(Fig. __COMPETITION_FIG__)")
    scientific_body = scientific_body.replace("(Fig. 3)", "(Fig. 4)")
    scientific_body = scientific_body.replace("(Fig. __COMPETITION_FIG__)", "(Fig. 3)")
    scientific_body = add_figure_markers(scientific_body)

    abstract_background = abstract_sections["Background"].strip() + "\n\n" + abstract_sections["Methods"].strip()
    abstract_results = abstract_sections["Results"].strip()
    # The source draft places the keyword line after the Conclusions block; it is
    # emitted once as the journal-level Keywords section below.
    abstract_conclusions = re.split(r"\n\*\*Keywords:\*\*", abstract_sections["Conclusions"], maxsplit=1)[0].strip()

    output = f"""# {title.group(1).strip()}

*Research Article*

**Authors:** [AUTHOR 1 FULL NAME]1, [AUTHOR 2 FULL NAME]2

1. [Department, Institution, City, Country]
2. [Department, Institution, City, Country]

**Corresponding author:** [FULL NAME, EMAIL ADDRESS]

## Abstract

### Background

{abstract_background}

### Results

{abstract_results}

### Conclusions

{abstract_conclusions}

## Keywords

fibroblast; fascia; mechanotransduction; integrin; focal adhesion; actomyosin; YAP/TAZ; single-cell RNA sequencing

{scientific_body.rstrip()}

## List of abbreviations

BH, Benjamini–Hochberg; CAMERA, correlation-adjusted mean-rank gene-set test; ECM, extracellular matrix; F7, fascia-like myofibroblast; H5AD, hierarchical data format for annotated data; IQR, interquartile range; RDS, R data serialization format; SRS, sequence read archive sample; TAZ, transcriptional co-activator with PDZ-binding motif; TEAD, TEA domain transcription factor; TGF-β, transforming growth factor beta; YAP, yes-associated protein.

{build_declarations(data_text, ethics_text).rstrip()}

## References

""" + "\n".join(ordered_references) + "\n"

    abstract_word_count = len(re.findall(r"\b[\w’/-]+\b", abstract_background + " " + abstract_results + " " + abstract_conclusions))
    audit = {
        "title": title.group(1).strip(),
        "abstract_word_count": abstract_word_count,
        "reference_count": len(ordered_references),
        "citation_count_in_order": len(citation_order),
        "uncited_reference_count": len(uncited),
        "unresolved_citations": sorted(set(unresolved)),
        "official_guideline": "https://link.springer.com/journal/12864/submission-guidelines/research-article",
    }
    return output, audit


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manuscript, audit = build_markdown()
    OUTPUT_MD.write_text(manuscript, encoding="utf-8", newline="\n")
    audit_lines = [
        "# BMC Genomics format audit",
        "",
        "This is a format/provenance audit of the BMC Genomics Research article package. Scientific results and evidence boundaries were not changed.",
        "",
        f"- Abstract word count: **{audit['abstract_word_count']}** (official limit: 350)",
        f"- Vancouver reference count: **{audit['reference_count']}**",
        f"- In-text citation keys assigned in first-citation order: **{audit['citation_count_in_order']}**",
        f"- Uncited reference count: **{audit['uncited_reference_count']}**",
        f"- Unresolved author-year citations: **{len(audit['unresolved_citations'])}**",
        "- Main section structure: Background; Methods; Results; Discussion; Conclusions",
        "- Declarations: all mandatory BMC Genomics Research article subheadings included",
        "- Submission placeholders retained: author details, correspondence, competing interests, funding, contributions, acknowledgements, public code/archive URL",
        "",
        "Official basis: [BMC Genomics Research article submission guidelines](https://link.springer.com/journal/12864/submission-guidelines/research-article).",
        "",
    ]
    AUDIT_MD.write_text("\n".join(audit_lines), encoding="utf-8", newline="\n")
    print(OUTPUT_MD)
    print(AUDIT_MD)


if __name__ == "__main__":
    main()
