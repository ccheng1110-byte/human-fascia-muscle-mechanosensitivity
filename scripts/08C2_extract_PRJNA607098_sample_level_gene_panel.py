"""Step 08C2: extract sample-level F7/non-F7 gene-panel summaries.

Only observation metadata, the gene index, and expression chunks containing
the frozen panel are read from the public skin-fibroblast atlas Zarr.  The
full H5AD is not downloaded.  This helper performs aggregation only; all
inferential statistics are performed in the paired-sample R script.
"""

from __future__ import annotations

import csv
import os
import re
import time
from collections.abc import Callable

import fsspec
import numpy as np
import zarr


def _retry_read(reader: Callable, label: str, attempts: int = 5):
    last_error = None
    for attempt in range(attempts):
        try:
            return reader()
        except Exception as exc:  # Network libraries expose several transient types.
            last_error = exc
            if attempt + 1 < attempts:
                time.sleep(min(2**attempt, 8))
    raise RuntimeError(f"Remote read failed after {attempts} attempts: {label}") from last_error


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
        categories = _as_text(
            _retry_read(lambda: np.asarray(node["categories"][:]), path + "/categories")
        )
        codes = np.asarray(
            _retry_read(lambda: np.asarray(node["codes"][:]), path + "/codes"),
            dtype=np.int64,
        )
        values = np.full(codes.shape[0], "__MISSING__", dtype=object)
        valid = (codes >= 0) & (codes < categories.shape[0])
        values[valid] = categories[codes[valid]]
        return values.astype(str)
    return _as_text(_retry_read(lambda: np.asarray(node[:]), path))


def _write_rows(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _recover_source_samples(
    observation_index: np.ndarray,
    source_mask: np.ndarray,
    expected_samples: int,
) -> np.ndarray:
    sample_id = np.full(observation_index.shape[0], "__MISSING__", dtype=object)
    source_indices = np.flatnonzero(source_mask)
    for index in source_indices:
        value = str(observation_index[index])
        if "_" not in value:
            continue
        barcode, candidate_sample = value.rsplit("_", 1)
        valid_barcode = re.fullmatch(r"[ACGTN]{12,}(?:-\d+)?", barcode) is not None
        valid_sample = re.fullmatch(r"SRS\d+", candidate_sample) is not None
        if valid_barcode and valid_sample:
            sample_id[index] = candidate_sample

    source_samples = sample_id[source_mask].astype(str)
    if np.any(source_samples == "__MISSING__"):
        raise RuntimeError(
            "At least one PRJNA607098 observation index failed the frozen "
            "barcode_SRS sample reconstruction rule."
        )
    observed_samples = np.unique(source_samples)
    if observed_samples.size != int(expected_samples):
        raise RuntimeError(
            f"Expected {expected_samples} reconstructed samples but found "
            f"{observed_samples.size}."
        )
    return sample_id.astype(str)


def run_step08c2_extract(
    panel_genes,
    output_dir: str,
    zarr_url: str = (
        "https://storage.googleapis.com/haniffalab/skin-fibroblast/"
        "zarr/adata_webportal.zarr"
    ),
    source_id: str = "PRJNA607098",
    target_state: str = "F7: Fascia-like myofibroblast",
    expected_samples: int = 12,
    allow_missing_genes: bool = False,
):
    panel_genes = list(dict.fromkeys(str(value) for value in panel_genes))
    os.makedirs(output_dir, exist_ok=True)

    filesystem = fsspec.filesystem("http")
    mapper = fsspec.mapping.FSMap(
        root=zarr_url,
        fs=filesystem,
        check=False,
        create=False,
    )
    root = _retry_read(
        lambda: zarr.open_consolidated(mapper, mode="r"),
        "consolidated Zarr metadata",
    )

    gse = _decode_obs(root, "obs/GSE")
    celltype = _decode_obs(root, "obs/celltype")
    observation_index = _decode_obs(root, "obs/_index")
    if not (gse.shape[0] == celltype.shape[0] == observation_index.shape[0]):
        raise RuntimeError("Observation metadata arrays have inconsistent lengths.")

    source_mask = gse == source_id
    if not np.any(source_mask):
        raise RuntimeError("No cells were found for source: " + source_id)
    sample_id = _recover_source_samples(
        observation_index,
        source_mask,
        expected_samples,
    )

    gene_symbols = _as_text(
        _retry_read(lambda: np.asarray(root["var/gene_symbol"][:]), "var/gene_symbol")
    )
    symbol_to_index = {}
    for index, symbol in enumerate(gene_symbols):
        if symbol not in symbol_to_index:
            symbol_to_index[symbol] = index
    found_genes = [gene for gene in panel_genes if gene in symbol_to_index]
    missing_genes = [gene for gene in panel_genes if gene not in symbol_to_index]
    if missing_genes and not allow_missing_genes:
        raise RuntimeError(
            "Frozen panel gene(s) missing from the atlas: " + ", ".join(missing_genes)
        )

    x_array = root["X"]
    row_chunk, column_chunk = (int(value) for value in x_array.chunks)
    if row_chunk != int(x_array.shape[0]):
        raise RuntimeError(
            "Unexpected X row chunking; the selective-read rule must be re-audited."
        )

    gene_indices = [symbol_to_index[gene] for gene in found_genes]
    sorted_pairs = sorted(zip(gene_indices, found_genes))
    sorted_indices = [pair[0] for pair in sorted_pairs]
    sorted_genes = [pair[1] for pair in sorted_pairs]
    x_sorted = np.asarray(
        _retry_read(
            lambda: x_array.oindex[:, sorted_indices],
            "frozen gene-panel expression chunks",
        ),
        dtype=np.float32,
    )
    sorted_position = {gene: index for index, gene in enumerate(sorted_genes)}
    x_panel = np.column_stack(
        [x_sorted[:, sorted_position[gene]] for gene in found_genes]
    )
    del x_sorted

    source_samples = sample_id[source_mask]
    source_states = celltype[source_mask]
    source_expression = x_panel[source_mask, :]
    del x_panel

    sample_rows = []
    for sample in sorted(np.unique(source_samples).tolist()):
        sample_mask = source_samples == sample
        for comparison_group, group_mask in (
            ("F7", sample_mask & (source_states == target_state)),
            ("non_F7", sample_mask & (source_states != target_state)),
        ):
            cells = int(np.sum(group_mask))
            if cells == 0:
                raise RuntimeError(
                    f"Sample {sample} has no cells in comparison group {comparison_group}."
                )
            values = source_expression[group_mask, :]
            for gene_position, gene in enumerate(found_genes):
                gene_values = values[:, gene_position]
                sample_rows.append(
                    {
                        "GSE": source_id,
                        "sample_id": sample,
                        "comparison_group": comparison_group,
                        "gene": gene,
                        "cells": cells,
                        "mean_expression": float(np.mean(gene_values, dtype=np.float64)),
                        "median_expression": float(np.median(gene_values)),
                        "percent_cells_expressing": float(100.0 * np.mean(gene_values > 0)),
                    }
                )

    summary_path = os.path.join(
        output_dir,
        "PRJNA607098_sample_F7_nonF7_gene_expression_summary_v1.csv",
    )
    _write_rows(
        summary_path,
        sample_rows,
        [
            "GSE",
            "sample_id",
            "comparison_group",
            "gene",
            "cells",
            "mean_expression",
            "median_expression",
            "percent_cells_expressing",
        ],
    )

    chunk_rows = []
    for chunk_index in sorted({index // column_chunk for index in gene_indices}):
        relative_key = f"X/0.{chunk_index}"
        chunk_rows.append(
            {
                "zarr_key": relative_key,
                "url": zarr_url.rstrip("/") + "/" + relative_key,
                "access_mode": "user_run_R_remote_stream",
            }
        )
    chunk_manifest_path = os.path.join(
        output_dir,
        "PRJNA607098_sample_level_gene_panel_chunk_manifest_v1.csv",
    )
    _write_rows(
        chunk_manifest_path,
        chunk_rows,
        ["zarr_key", "url", "access_mode"],
    )

    return {
        "sample_expression_summary": summary_path,
        "chunk_manifest": chunk_manifest_path,
        "source_cells": int(np.sum(source_mask)),
        "samples": int(np.unique(source_samples).size),
        "found_genes": found_genes,
        "missing_genes": missing_genes,
    }
