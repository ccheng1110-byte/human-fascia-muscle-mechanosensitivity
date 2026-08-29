"""Step 08B helper: source-specific gene-panel extraction from public Zarr.

The helper reads only the observation metadata, gene index, and Zarr expression
chunks containing the frozen gene panel. It records every required expression
chunk and its HTTP HEAD size before streaming it through the user-run R job.
"""

from __future__ import annotations

import csv
import os
import time
import urllib.error
import urllib.request
from collections import Counter

import fsspec
import numpy as np
import zarr


def _decode_categorical(root, path: str) -> np.ndarray:
    group = root[path]
    categories = np.asarray(group["categories"][:]).astype(str)
    codes = np.asarray(group["codes"][:], dtype=np.int64)
    values = np.full(codes.shape[0], "__MISSING__", dtype=object)
    valid = (codes >= 0) & (codes < categories.shape[0])
    values[valid] = categories[codes[valid]]
    return values.astype(str)


def _write_rows(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _safe_mean(x: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    return float(np.mean(x, dtype=np.float64))


def _safe_median(x: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    return float(np.median(x))


def _safe_prevalence(x: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    return float(np.mean(x > 0))


def _head_size(
    url: str,
    timeout_seconds: int = 120,
    max_attempts: int = 4,
) -> int | None:
    """Return remote size when available; network-size audit is non-blocking.

    Google Storage occasionally closes a TLS connection during a HEAD request.
    Retry transient failures, then use a one-byte range request. If both methods
    fail, return None because the user-run R workflow has no per-file size cap.
    """
    headers = {"User-Agent": "Step08B-source-specific-audit"}
    for attempt in range(max_attempts):
        try:
            request = urllib.request.Request(url, method="HEAD", headers=headers)
            with urllib.request.urlopen(
                request, timeout=timeout_seconds
            ) as response:
                value = response.headers.get("Content-Length")
            if value is not None:
                return int(value)
        except (urllib.error.URLError, TimeoutError, OSError):
            pass
        if attempt + 1 < max_attempts:
            time.sleep(min(2**attempt, 8))

    range_headers = dict(headers)
    range_headers["Range"] = "bytes=0-0"
    for attempt in range(max_attempts):
        try:
            request = urllib.request.Request(
                url, method="GET", headers=range_headers
            )
            with urllib.request.urlopen(
                request, timeout=timeout_seconds
            ) as response:
                content_range = response.headers.get("Content-Range")
                if content_range and "/" in content_range:
                    return int(content_range.rsplit("/", 1)[1])
                value = response.headers.get("Content-Length")
            if value is not None and int(value) > 1:
                return int(value)
        except (urllib.error.URLError, TimeoutError, OSError, ValueError):
            pass
        if attempt + 1 < max_attempts:
            time.sleep(min(2**attempt, 8))
    return None


def run_step08b(
    panel_genes,
    output_dir: str,
    manual_download_dir: str,
    zarr_url: str = (
        "https://storage.googleapis.com/haniffalab/skin-fibroblast/"
        "zarr/adata_webportal.zarr"
    ),
    state_label_field: str = "celltype_skinspecific_nomenclature",
    target_state: str = "F7: Fascia-like myofibroblast",
    output_version: str = "v1",
    audit_chunk_sizes: bool = True,
):
    panel_genes = [str(x) for x in panel_genes]
    panel_genes = list(dict.fromkeys(panel_genes))
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(manual_download_dir, exist_ok=True)

    filesystem = fsspec.filesystem("http")
    mapper = fsspec.mapping.FSMap(
        root=zarr_url,
        fs=filesystem,
        check=False,
        create=False,
    )
    root = zarr.open_consolidated(mapper, mode="r")

    gse = _decode_categorical(root, "obs/GSE")
    patient_status = _decode_categorical(root, "obs/Patient_status")
    disease_category = _decode_categorical(root, "obs/disease_category_orig")
    celltype = _decode_categorical(root, "obs/celltype")
    lesion_status = _decode_categorical(root, "obs/lesional_vs_nonlesional")
    skin_state = _decode_categorical(
        root, "obs/celltype_skinspecific_nomenclature"
    )
    state_vectors = {
        "celltype": celltype,
        "celltype_skinspecific_nomenclature": skin_state,
    }
    if state_label_field not in state_vectors:
        raise RuntimeError(
            "Unsupported state_label_field: "
            + state_label_field
            + ". Expected one of: "
            + ", ".join(sorted(state_vectors))
        )
    analysis_state = state_vectors[state_label_field]

    # Metadata counts are descriptive and do not treat cells as replicates.
    metadata_counter = Counter(
        zip(
            gse,
            patient_status,
            disease_category,
            celltype,
            lesion_status,
            skin_state,
        )
    )
    metadata_rows = []
    for key, count in sorted(metadata_counter.items()):
        metadata_rows.append(
            {
                "GSE": key[0],
                "Patient_status": key[1],
                "disease_category_orig": key[2],
                "celltype": key[3],
                "lesional_vs_nonlesional": key[4],
                "celltype_skinspecific_nomenclature": key[5],
                "cells": int(count),
            }
        )
    metadata_path = os.path.join(
        output_dir,
        f"skin_fibroblast_atlas_metadata_counts_{output_version}.csv",
    )
    _write_rows(
        metadata_path,
        metadata_rows,
        [
            "GSE",
            "Patient_status",
            "disease_category_orig",
            "celltype",
            "lesional_vs_nonlesional",
            "celltype_skinspecific_nomenclature",
            "cells",
        ],
    )

    gene_symbols = np.asarray(root["var/gene_symbol"][:]).astype(str)
    symbol_to_index = {}
    for index, symbol in enumerate(gene_symbols):
        if symbol not in symbol_to_index:
            symbol_to_index[symbol] = index

    found_genes = [gene for gene in panel_genes if gene in symbol_to_index]
    missing_genes = [gene for gene in panel_genes if gene not in symbol_to_index]
    if not found_genes:
        raise RuntimeError("None of the requested genes were found in the Zarr object.")

    x_array = root["X"]
    row_chunk, column_chunk = (int(x) for x in x_array.chunks)
    if row_chunk != int(x_array.shape[0]):
        raise RuntimeError(
            "Unexpected X row chunking; the pre-audited selective-read rule no longer applies."
        )

    gene_indices = [symbol_to_index[gene] for gene in found_genes]
    unique_chunk_indices = sorted({index // column_chunk for index in gene_indices})
    chunk_rows = []
    for chunk_index in unique_chunk_indices:
        relative_key = f"X/0.{chunk_index}"
        chunk_url = zarr_url.rstrip("/") + "/" + relative_key
        size_bytes = _head_size(chunk_url) if audit_chunk_sizes else None
        manual_path = os.path.join(
            manual_download_dir,
            "X",
            f"0.{chunk_index}",
        )
        if not audit_chunk_sizes:
            status = "REMOTE_STREAM_ALLOWED_SIZE_NOT_AUDITED"
        elif size_bytes is not None:
            status = "REMOTE_STREAM_ALLOWED"
        else:
            status = "REMOTE_STREAM_ALLOWED_SIZE_UNAVAILABLE"
        chunk_rows.append(
            {
                "zarr_key": relative_key,
                "url": chunk_url,
                "size_bytes": size_bytes,
                "size_mib": (
                    size_bytes / (1024**2) if size_bytes is not None else None
                ),
                "status": status,
                "manual_destination": manual_path,
            }
        )

    chunk_manifest_path = os.path.join(
        output_dir,
        f"skin_fibroblast_atlas_gene_panel_chunk_manifest_{output_version}.csv",
    )
    _write_rows(
        chunk_manifest_path,
        chunk_rows,
        [
            "zarr_key",
            "url",
            "size_bytes",
            "size_mib",
            "status",
            "manual_destination",
        ],
    )

    # Sort indices for predictable chunk access, then restore requested gene order.
    sorted_pairs = sorted(zip(gene_indices, found_genes))
    sorted_indices = [pair[0] for pair in sorted_pairs]
    sorted_genes = [pair[1] for pair in sorted_pairs]
    x_sorted = np.asarray(x_array.oindex[:, sorted_indices], dtype=np.float32)
    sorted_position = {gene: index for index, gene in enumerate(sorted_genes)}
    x_panel = np.column_stack(
        [x_sorted[:, sorted_position[gene]] for gene in found_genes]
    )
    del x_sorted

    target_sources = ["PRJNA607098", "GSE173252"]
    valid_states = sorted(
        state for state in np.unique(analysis_state) if state != "__MISSING__"
    )
    expression_rows = []
    for source in target_sources:
        source_mask = gse == source
        for state in valid_states:
            mask = source_mask & (analysis_state == state)
            n_cells = int(np.sum(mask))
            if n_cells == 0:
                continue
            values = x_panel[mask, :]
            for gene_index, gene in enumerate(found_genes):
                gene_values = values[:, gene_index]
                expression_rows.append(
                    {
                        "GSE": source,
                        "cell_state": state,
                        "gene": gene,
                        "cells": n_cells,
                        "mean_expression": _safe_mean(gene_values),
                        "median_expression": _safe_median(gene_values),
                        "percent_cells_expressing": 100.0
                        * _safe_prevalence(gene_values),
                    }
                )

    expression_path = os.path.join(
        output_dir,
        f"skin_fibroblast_atlas_source_state_gene_expression_summary_{output_version}.csv",
    )
    _write_rows(
        expression_path,
        expression_rows,
        [
            "GSE",
            "cell_state",
            "gene",
            "cells",
            "mean_expression",
            "median_expression",
            "percent_cells_expressing",
        ],
    )

    contrast_rows = []
    for source in target_sources:
        source_rows = [row for row in expression_rows if row["GSE"] == source]
        for gene in found_genes:
            gene_rows = [row for row in source_rows if row["gene"] == gene]
            target_rows = [row for row in gene_rows if row["cell_state"] == target_state]
            other_rows = [row for row in gene_rows if row["cell_state"] != target_state]
            if not target_rows:
                contrast_rows.append(
                    {
                        "GSE": source,
                        "gene": gene,
                        "F7_cells": 0,
                        "comparison_states": len(other_rows),
                        "F7_mean_expression": float("nan"),
                        "median_other_state_mean_expression": float("nan"),
                        "F7_minus_median_other_state_mean": float("nan"),
                        "F7_minus_best_other_state_mean": float("nan"),
                        "F7_percent_cells_expressing": float("nan"),
                        "median_other_state_percent_expressing": float("nan"),
                        "F7_minus_median_other_state_percent_expressing": float("nan"),
                        "F7_highest_mean_across_states": False,
                        "source_specific_descriptive_support": False,
                    }
                )
                continue

            target = target_rows[0]
            other_means = np.asarray(
                [row["mean_expression"] for row in other_rows], dtype=float
            )
            other_prevalence = np.asarray(
                [row["percent_cells_expressing"] for row in other_rows], dtype=float
            )
            median_other_mean = _safe_median(other_means)
            median_other_prevalence = _safe_median(other_prevalence)
            best_other_mean = float(np.max(other_means)) if other_means.size else float("nan")
            mean_difference = target["mean_expression"] - median_other_mean
            prevalence_difference = (
                target["percent_cells_expressing"] - median_other_prevalence
            )
            highest_mean = bool(
                other_means.size > 0 and target["mean_expression"] > best_other_mean
            )
            descriptive_support = bool(
                target["cells"] >= 20
                and other_means.size >= 2
                and mean_difference > 0
                and prevalence_difference > 0
            )
            contrast_rows.append(
                {
                    "GSE": source,
                    "gene": gene,
                    "F7_cells": target["cells"],
                    "comparison_states": len(other_rows),
                    "F7_mean_expression": target["mean_expression"],
                    "median_other_state_mean_expression": median_other_mean,
                    "F7_minus_median_other_state_mean": mean_difference,
                    "F7_minus_best_other_state_mean": (
                        target["mean_expression"] - best_other_mean
                    ),
                    "F7_percent_cells_expressing": target[
                        "percent_cells_expressing"
                    ],
                    "median_other_state_percent_expressing": median_other_prevalence,
                    "F7_minus_median_other_state_percent_expressing": prevalence_difference,
                    "F7_highest_mean_across_states": highest_mean,
                    "source_specific_descriptive_support": descriptive_support,
                }
            )

    contrast_path = os.path.join(
        output_dir,
        f"skin_fibroblast_atlas_source_specific_F7_gene_contrasts_{output_version}.csv",
    )
    _write_rows(
        contrast_path,
        contrast_rows,
        [
            "GSE",
            "gene",
            "F7_cells",
            "comparison_states",
            "F7_mean_expression",
            "median_other_state_mean_expression",
            "F7_minus_median_other_state_mean",
            "F7_minus_best_other_state_mean",
            "F7_percent_cells_expressing",
            "median_other_state_percent_expressing",
            "F7_minus_median_other_state_percent_expressing",
            "F7_highest_mean_across_states",
            "source_specific_descriptive_support",
        ],
    )

    gene_inventory_rows = []
    for gene in panel_genes:
        gene_inventory_rows.append(
            {
                "gene": gene,
                "found_in_atlas": gene in symbol_to_index,
                "atlas_gene_index_zero_based": symbol_to_index.get(gene, ""),
            }
        )
    gene_inventory_path = os.path.join(
        output_dir,
        f"skin_fibroblast_atlas_gene_panel_inventory_{output_version}.csv",
    )
    _write_rows(
        gene_inventory_path,
        gene_inventory_rows,
        ["gene", "found_in_atlas", "atlas_gene_index_zero_based"],
    )

    return {
        "metadata_counts": metadata_path,
        "chunk_manifest": chunk_manifest_path,
        "gene_inventory": gene_inventory_path,
        "expression_summary": expression_path,
        "source_specific_contrasts": contrast_path,
        "found_genes": found_genes,
        "missing_genes": missing_genes,
        "state_label_field": state_label_field,
        "target_state": target_state,
        "total_remote_expression_bytes": int(
            sum(
                row["size_bytes"]
                for row in chunk_rows
                if row["size_bytes"] is not None
            )
        ),
    }
