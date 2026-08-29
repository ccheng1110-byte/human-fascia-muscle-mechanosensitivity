"""Build the reviewer-revised manuscript DOCX from Markdown and Figure 1-5 PNGs.

The document assembly reuses the project's established DOCX builder so the
Markdown remains the single source of manuscript text and all figure assets
are inserted at Results-level evidence boundaries.
"""

import os
import runpy
from pathlib import Path


PROJECT = Path(r".")
INPUT_MD = PROJECT / (
    "results/11_manuscript_preparation/20A_review_adoption_and_revision/"
    "20A_revised_manuscript_v5.md"
)
FIG_DIR = PROJECT / (
    "results/11_manuscript_preparation/20B_reviewer_response_figures"
)
OUTPUT_DIR = PROJECT / (
    "results/11_manuscript_preparation/20C_word_revised_manuscript_with_figures"
)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DOCX = OUTPUT_DIR / "20C_revised_manuscript_v5_with_figures.docx"

os.environ["MANUSCRIPT_INPUT_MD"] = str(INPUT_MD)
os.environ["MANUSCRIPT_OUTPUT_DIR"] = str(OUTPUT_DIR)
os.environ["MANUSCRIPT_OUTPUT_DOCX"] = str(OUTPUT_DOCX)

builder_path = PROJECT / "scripts/18C_build_word_manuscript_with_figures.py"
builder = runpy.run_path(str(builder_path))

figure_map = {
    "### Provenance and primary validation structure": {
        "path": FIG_DIR / "20B_Figure1_evidence_workflow.png",
        "caption": "Figure 1 | Evidence workflow and inferential units. The staged design separates discovery, paired atlas validation, independent subject-level triangulation, direct mechanical perturbation, regulatory/stiffness context and bounded synthesis. Each source retains its own biological or inferential denominator.",
    },
    "### Core modules were consistently positive in PRJNA607098": {
        "path": FIG_DIR / "20B_Figure2_sample_subject_module_contrasts.png",
        "caption": "Figure 2 | Primary-module contrasts at the correct inferential unit. PRJNA607098 points are 12 reconstructed SRS sample units comparing F7 with pooled non-F7 cells within sample. GSE130973 points are five subjects comparing marker-defined candidate fibroblast clusters with other clusters within subject. Black summaries are medians and IQRs; source-specific free y-axis scales are intentional.",
    },
    "### Competition analysis limited integrin specificity": {
        "path": FIG_DIR / "20B_Figure4_competition_residuals.png",
        "caption": "Figure 4 | Competition residuals expose the specificity boundary. Within PRJNA607098, each core-module contrast was regressed on each competitor-module contrast. Points are sample-level residuals and the frozen gate requires at least 8 of 12 residuals above zero. Integrin/focal adhesion fails against several competitors, whereas actomyosin/Rho is more robust.",
    },
    "### Independent skin-source triangulation was partial": {
        "path": FIG_DIR / "20B_Figure3_frozen_candidate_raw_contrasts.png",
        "caption": "Figure 3 | Frozen candidate effects are source dependent. All frozen candidates are shown without post hoc removal. Points are sample- or subject-level target-minus-comparator contrasts; black diamonds and bars show median and IQR. Source separation preserves distinct expression scales and denominators.",
    },
    "### Direct mechanical tension supported actomyosin/Rho but also cell cycle": {
        "path": FIG_DIR / "20B_Figure5_GSE300230_tension_response.png",
        "caption": "Figure 5 | Direct mechanical tension activates actomyosin and cell-cycle programs. Model-estimated module-score contrasts compare tension with relaxed collagen without TGF-beta in each GSE300230 cell line. Labels give CAMERA BH q values. The concurrent cell-cycle response is an interpretation boundary; the two cell lines are not independent donors or age replicates.",
    },
}
builder["FIGURES"] = figure_map
builder["build_doc"].__globals__["FIGURES"] = figure_map

path = builder["build_doc"]()
print(path)
