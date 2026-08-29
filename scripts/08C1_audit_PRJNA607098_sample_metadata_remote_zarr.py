"""Step 08C1: audit sample-level eligibility for PRJNA607098.

This helper reads only observation metadata from the public skin-fibroblast
atlas Zarr. It does not read the expression matrix and does not perform any
cell-level or sample-level hypothesis test.
"""

from __future__ import annotations

import csv
import os
import re
from collections import Counter

import fsspec
import numpy as np
import zarr


def _as_text(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values)
    if values.dtype.kind == "S":
        return np.asarray(
            [value.decode("utf-8", errors="replace") for value in values],
            dtype=str,
        )
    return values.astype(str)


def _decode_obs(root, path: str) -> np.ndarray:
    node = root[path]
    if hasattr(node, "keys") and "categories" in node and "codes" in node:
        categories = _as_text(np.asarray(node["categories"][:]))
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        values = np.full(codes.shape[0], "__MISSING__", dtype=object)
        valid = (codes >= 0) & (codes < categories.shape[0])
        values[valid] = categories[codes[valid]]
        return values.astype(str)
    return _as_text(np.asarray(node[:]))


def _first_existing_path(root, candidates: list[str]) -> str:
    for path in candidates:
        try:
            root[path]
            return path
        except (KeyError, AttributeError):
            continue
    raise RuntimeError(
        "None of the candidate observation fields were found: "
        + ", ".join(candidates)
    )


def _resolve_sample_id(
    root, source_mask: np.ndarray
) -> tuple[np.ndarray, str, bool, str]:
    """Resolve sample IDs without silently guessing from arbitrary barcodes.

    The web-portal Zarr omits the original ``sample_id`` column.  The atlas'
    official integration code constructed each observation name as
    ``barcode + '_' + sample_id`` before export, so the final underscore suffix
    can be used only after strict format checks.
    """
    direct_candidates = [
        "obs/sample_id",
        "obs/sampleID",
        "obs/Sample",
        "obs/sample",
    ]
    for path in direct_candidates:
        try:
            root[path]
        except (KeyError, AttributeError):
            continue
        return _decode_obs(root, path), path, True, "direct observation field"

    index_path = "obs/_index"
    try:
        observation_index = _decode_obs(root, index_path)
    except (KeyError, AttributeError) as exc:
        raise RuntimeError(
            "No direct sample field or obs/_index was found; sample-level "
            "analysis is not possible from this Zarr."
        ) from exc

    sample_id = np.full(observation_index.shape[0], "__MISSING__", dtype=object)
    format_valid = np.zeros(observation_index.shape[0], dtype=bool)
    for i, value in enumerate(observation_index.tolist()):
        if "_" not in value:
            continue
        barcode, candidate_sample = value.rsplit("_", 1)
        # PRJNA607098 uses 10x nucleotide barcodes followed by SRA study-run
        # sample accessions (SRS...).  Requiring both patterns prevents an
        # unrelated underscore suffix from being treated as a biological unit.
        valid_barcode = re.fullmatch(r"[ACGTN]{12,}(?:-\d+)?", barcode) is not None
        valid_sample = re.fullmatch(r"SRS\d+", candidate_sample) is not None
        if valid_barcode and valid_sample:
            sample_id[i] = candidate_sample
            format_valid[i] = True

    source_index_entries_valid = bool(np.all(format_valid[source_mask]))
    provenance = (
        "reconstructed from obs/_index final suffix; format defined by the "
        "official atlas integration code (barcode + '_' + sample_id)"
    )
    return sample_id.astype(str), index_path, source_index_entries_valid, provenance


def _write_rows(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_step08c1(
    output_dir: str,
    zarr_url: str = (
        "https://storage.googleapis.com/haniffalab/skin-fibroblast/"
        "zarr/adata_webportal.zarr"
    ),
    source_id: str = "PRJNA607098",
    target_state: str = "F7: Fascia-like myofibroblast",
    minimum_cells_per_state: int = 20,
    minimum_eligible_samples: int = 6,
    expected_public_samples: int = 12,
):
    os.makedirs(output_dir, exist_ok=True)

    filesystem = fsspec.filesystem("http")
    mapper = fsspec.mapping.FSMap(
        root=zarr_url,
        fs=filesystem,
        check=False,
        create=False,
    )
    root = zarr.open_consolidated(mapper, mode="r")

    gse = _decode_obs(root, "obs/GSE")
    celltype = _decode_obs(root, "obs/celltype")
    source_mask = gse == source_id
    if not np.any(source_mask):
        raise RuntimeError("No cells were found for source: " + source_id)

    sample_id, sample_path, sample_format_valid, sample_provenance = (
        _resolve_sample_id(root, source_mask)
    )
    if not (sample_id.shape[0] == gse.shape[0] == celltype.shape[0]):
        raise RuntimeError("Observation metadata arrays have inconsistent lengths.")

    source_samples = sample_id[source_mask]
    source_states = celltype[source_mask]
    total_source_cells = int(source_samples.shape[0])
    missing_sample_mask = (source_samples == "__MISSING__") | (source_samples == "")
    missing_sample_cells = int(np.sum(missing_sample_mask))

    valid_mask = ~missing_sample_mask
    valid_samples = source_samples[valid_mask]
    valid_states = source_states[valid_mask]
    sample_ids = sorted(np.unique(valid_samples).tolist())

    pair_counter = Counter(zip(valid_samples.tolist(), valid_states.tolist()))
    count_rows = []
    for (sample, state), cells in sorted(pair_counter.items()):
        count_rows.append(
            {
                "GSE": source_id,
                "sample_id": sample,
                "celltype": state,
                "cells": int(cells),
            }
        )

    count_output = os.path.join(
        output_dir, "PRJNA607098_sample_by_celltype_counts_v1.csv"
    )
    _write_rows(
        count_output,
        count_rows,
        ["GSE", "sample_id", "celltype", "cells"],
    )

    sample_rows = []
    for sample in sample_ids:
        mask = valid_samples == sample
        sample_states = valid_states[mask]
        total_cells = int(np.sum(mask))
        f7_cells = int(np.sum(sample_states == target_state))
        other_cells = total_cells - f7_cells
        eligible = bool(
            f7_cells >= int(minimum_cells_per_state)
            and other_cells >= int(minimum_cells_per_state)
        )
        sample_rows.append(
            {
                "GSE": source_id,
                "sample_id": sample,
                "total_cells": total_cells,
                "F7_cells": f7_cells,
                "other_state_cells": other_cells,
                "F7_fraction": f7_cells / total_cells if total_cells else float("nan"),
                "eligible_for_paired_sample_state_analysis": eligible,
            }
        )

    sample_output = os.path.join(
        output_dir, "PRJNA607098_sample_F7_eligibility_v1.csv"
    )
    _write_rows(
        sample_output,
        sample_rows,
        [
            "GSE",
            "sample_id",
            "total_cells",
            "F7_cells",
            "other_state_cells",
            "F7_fraction",
            "eligible_for_paired_sample_state_analysis",
        ],
    )

    n_samples = len(sample_rows)
    n_eligible = sum(
        row["eligible_for_paired_sample_state_analysis"] for row in sample_rows
    )
    missing_fraction = (
        missing_sample_cells / total_source_cells if total_source_cells else 1.0
    )
    largest_sample_fraction = (
        max(row["total_cells"] for row in sample_rows) / total_source_cells
        if sample_rows
        else 1.0
    )
    exact_expected_sample_count = n_samples == int(expected_public_samples)

    proceed_to_08c2 = bool(
        sample_format_valid
        and exact_expected_sample_count
        and n_eligible >= int(minimum_eligible_samples)
        and missing_fraction <= 0.01
        and largest_sample_fraction <= 0.40
    )

    decision_output = os.path.join(
        output_dir, "PRJNA607098_sample_metadata_audit_decision_v1.md"
    )
    report_lines = [
        "## Material Passport",
        "",
        "- Origin Skill: academic-research-suite / experiment-agent",
        "- Origin Mode: validate",
        "- Origin Date: 2026-08-23",
        "- Verification Status: ANALYZED",
        "- Version Label: PRJNA607098_sample_metadata_audit_v1",
        "",
        "## Step 08C1 PRJNA607098 sample-metadata audit",
        "",
        "### Frozen eligibility rule",
        "",
        f"- At least {minimum_eligible_samples} samples must each contain at least "
        f"{minimum_cells_per_state} F7 cells and {minimum_cells_per_state} non-F7 cells.",
        "- Missing sample IDs must account for no more than 1% of source cells.",
        "- No single sample may account for more than 40% of source cells.",
        "",
        "### Results",
        "",
        f"- Resolved sample field: `{sample_path}`.",
        f"- Sample-ID provenance: {sample_provenance}.",
        f"- Sample-index format valid for all PRJNA607098 cells: {sample_format_valid}.",
        f"- PRJNA607098 cells: {total_source_cells}.",
        f"- Distinct sample IDs: {n_samples}.",
        f"- Public-record expected sample count: {expected_public_samples}.",
        f"- Exact expected sample count: {exact_expected_sample_count}.",
        f"- Eligible paired sample-state units: {n_eligible}/{n_samples}.",
        f"- Missing sample-ID fraction: {missing_fraction:.6f}.",
        f"- Largest sample fraction: {largest_sample_fraction:.6f}.",
        f"- Proceed to Step 08C2: {proceed_to_08c2}.",
        "",
        "### Evidence boundary",
        "",
        "- This audit checks sample-level feasibility only; it performs no hypothesis test.",
        "- Recovered SRS accessions are sample-level units, not automatically independent donors.",
        "- Step 08C2 must retain the accession provenance and verify donor independence before donor-level claims.",
        "- Step 08C2 must aggregate within sample before inference; individual cells are not replicates.",
        "- A successful audit does not test protein abundance, channel activity or causality.",
    ]
    with open(decision_output, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(report_lines) + "\n")

    return {
        "sample_field": sample_path,
        "sample_provenance": sample_provenance,
        "sample_format_valid": sample_format_valid,
        "total_source_cells": total_source_cells,
        "n_samples": n_samples,
        "expected_public_samples": int(expected_public_samples),
        "exact_expected_sample_count": exact_expected_sample_count,
        "n_eligible_samples": n_eligible,
        "missing_sample_fraction": missing_fraction,
        "largest_sample_fraction": largest_sample_fraction,
        "proceed_to_08c2": proceed_to_08c2,
        "sample_celltype_counts": count_output,
        "sample_eligibility": sample_output,
        "decision": decision_output,
    }
