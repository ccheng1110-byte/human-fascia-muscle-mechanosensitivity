"""Step 17C: targeted single-cell co-detection summaries for S1.

Only PRJNA607098 observation metadata, the gene index, and the targeted
registry-gene chunks are streamed from the public Zarr. The helper restricts
the primary analysis to target-state cells and returns sample-level summaries;
cells are never treated as independent biological replicates.
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


INTEGRIN_GENES = ["ITGAV", "ITGB1", "ITGA2", "PARVA", "ITGA5", "FERMT2"]
ACTOMYOSIN_GENES = ["CFL1", "LIMK1", "CNN2", "CDC42"]
MECHANOSENSOR_GENES = ["PIEZO2", "TMEM63B", "PANX1"]
HIPPO_GENES = ["NF2", "TEAD1"]


def _retry_read(reader: Callable, label: str, attempts: int = 5):
    last_error = None
    for attempt in range(attempts):
        try:
            return reader()
        except Exception as exc:
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


def _find_optional_obs(root, candidates: list[str]):
    for path in candidates:
        try:
            values = _decode_obs(root, path)
            numeric = np.asarray(
                [float(value) if value not in ("", "NA", "nan", "None", "__MISSING__") else np.nan
                 for value in values],
                dtype=float,
            )
            if np.isfinite(numeric).sum() >= max(20, int(0.5 * numeric.size)):
                return path, numeric
        except (KeyError, AttributeError, ValueError, TypeError):
            continue
    return None, None


def _recover_sample_ids(observation_index: np.ndarray, source_mask: np.ndarray) -> np.ndarray:
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
    if np.any(sample_id[source_mask] == "__MISSING__"):
        raise RuntimeError("The frozen PRJNA607098 barcode/sample reconstruction failed.")
    return sample_id.astype(str)


def _write_rows(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _make_strata(library_values: np.ndarray | None, n_cells: int) -> tuple[list[np.ndarray], bool]:
    if library_values is None:
        return [np.arange(n_cells, dtype=np.int64)], False
    values = np.asarray(library_values, dtype=float)
    finite = np.isfinite(values)
    if finite.mean() < 0.90 or np.unique(values[finite]).size < 4:
        return [np.arange(n_cells, dtype=np.int64)], False
    quantiles = np.nanquantile(values, [0.0, 0.25, 0.5, 0.75, 1.0])
    breaks = np.unique(quantiles)
    if breaks.size < 3:
        return [np.arange(n_cells, dtype=np.int64)], False
    codes = np.digitize(values, breaks[1:-1], right=True)
    strata = [np.flatnonzero(codes == code) for code in np.unique(codes)]
    strata = [index for index in strata if index.size > 0]
    if len(strata) < 2:
        return [np.arange(n_cells, dtype=np.int64)], False
    return strata, True


def _summarize_pair(
    values: np.ndarray,
    library_values: np.ndarray | None,
    gene_a: str,
    gene_b: str,
    source_id: str,
    sample_id: str,
    target_state: str,
    pair_family: str,
    permutation_reps: int,
    rng: np.random.Generator,
) -> dict:
    a = np.asarray(values[:, 0] > 0, dtype=bool)
    b = np.asarray(values[:, 1] > 0, dtype=bool)
    n_cells = int(a.size)
    strata, library_stratified = _make_strata(library_values, n_cells)

    observed_count = int(np.sum(a & b))
    expected_count = 0.0
    for index in strata:
        expected_count += float(index.size) * float(np.mean(a[index])) * float(np.mean(b[index]))

    null_counts = np.zeros(int(permutation_reps), dtype=float)
    if permutation_reps > 0:
        for rep in range(int(permutation_reps)):
            null_count = 0
            for index in strata:
                null_count += int(np.sum(a[index] & rng.permutation(b[index])))
            null_counts[rep] = null_count
        permutation_p = float((1.0 + np.sum(null_counts >= observed_count)) / (permutation_reps + 1.0))
    else:
        permutation_p = np.nan

    observed_rate = observed_count / n_cells if n_cells else np.nan
    expected_rate = expected_count / n_cells if n_cells else np.nan
    excess_rate = observed_rate - expected_rate if n_cells else np.nan
    enrichment = observed_rate / expected_rate if expected_rate > 0 else np.nan
    return {
        "source": source_id,
        "sample_id": sample_id,
        "target_state": target_state,
        "pair_family": pair_family,
        "gene_a": gene_a,
        "gene_b": gene_b,
        "cells": n_cells,
        "gene_a_detection_rate": float(np.mean(a)) if n_cells else np.nan,
        "gene_b_detection_rate": float(np.mean(b)) if n_cells else np.nan,
        "observed_co_detection_rate": observed_rate,
        "expected_independence_rate": expected_rate,
        "excess_co_detection_rate": excess_rate,
        "co_detection_enrichment_ratio": enrichment,
        "permutation_p_one_sided": permutation_p,
        "permutation_reps": int(permutation_reps),
        "library_size_stratified": bool(library_stratified),
    }


def run_step17c_prjna(
    output_dir: str,
    panel_genes,
    zarr_url: str = (
        "https://storage.googleapis.com/haniffalab/skin-fibroblast/"
        "zarr/adata_webportal.zarr"
    ),
    source_id: str = "PRJNA607098",
    target_state: str = "F7: Fascia-like myofibroblast",
    expected_samples: int = 12,
    minimum_cells_per_sample: int = 20,
    permutation_reps: int = 200,
    random_seed: int = 20260826,
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
    sample_id = _recover_sample_ids(observation_index, source_mask)
    source_target_mask = source_mask & (celltype == target_state)
    source_samples = sample_id[source_target_mask]
    sample_ids = sorted(np.unique(source_samples).tolist())
    if len(sample_ids) != int(expected_samples):
        raise RuntimeError(
            f"Expected {expected_samples} target-state samples but found {len(sample_ids)}."
        )

    library_field, library_values_all = _find_optional_obs(
        root,
        [
            "obs/n_genes",
            "obs/nFeature_RNA",
            "obs/total_counts",
            "obs/total_counts_RNA",
            "obs/n_counts",
            "obs/library_size",
            "obs/LibrarySize",
        ],
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

    pair_specs = []
    for family, left_genes, right_genes in (
        ("integrin_x_actomyosin", INTEGRIN_GENES, ACTOMYOSIN_GENES),
        ("integrin_x_mechanosensor", INTEGRIN_GENES, MECHANOSENSOR_GENES),
        ("integrin_x_hippo", INTEGRIN_GENES, HIPPO_GENES),
    ):
        for gene_a in left_genes:
            for gene_b in right_genes:
                pair_specs.append((family, gene_a, gene_b))
    pair_genes = sorted({gene for _, gene_a, gene_b in pair_specs for gene in (gene_a, gene_b)})
    missing_pair_genes = [gene for gene in pair_genes if gene not in symbol_to_index]
    if missing_pair_genes:
        raise RuntimeError("Required candidate pair gene(s) missing: " + ", ".join(missing_pair_genes))

    sorted_pairs = sorted((symbol_to_index[gene], gene) for gene in pair_genes)
    sorted_indices = [index for index, _ in sorted_pairs]
    sorted_genes = [gene for _, gene in sorted_pairs]
    x_array = root["X"]
    if int(x_array.chunks[0]) != int(x_array.shape[0]):
        raise RuntimeError("Unexpected X row chunking; selective-read rule must be re-audited.")
    x_sorted = np.asarray(
        _retry_read(
            lambda: x_array.oindex[:, sorted_indices],
            "S1 targeted coexpression gene chunks",
        ),
        dtype=np.float32,
    )
    sorted_position = {gene: index for index, gene in enumerate(sorted_genes)}
    x_panel = np.column_stack(
        [x_sorted[:, sorted_position[gene]] for gene in pair_genes]
    )
    pair_position = {gene: index for index, gene in enumerate(pair_genes)}
    del x_sorted

    source_target_indices = np.flatnonzero(source_target_mask)
    target_expression = x_panel[source_target_indices, :]
    target_sample_ids = sample_id[source_target_indices]
    target_library = (
        library_values_all[source_target_indices]
        if library_values_all is not None
        else None
    )

    rng = np.random.default_rng(int(random_seed))
    coexpression_rows = []
    detection_rows = []
    for sample in sample_ids:
        sample_mask = target_sample_ids == sample
        sample_indices = np.flatnonzero(sample_mask)
        cells = int(sample_indices.size)
        if cells < int(minimum_cells_per_sample):
            raise RuntimeError(
                f"Target state has only {cells} cells in sample {sample}; "
                f"minimum is {minimum_cells_per_sample}."
            )
        sample_values = target_expression[sample_indices, :]
        sample_library = target_library[sample_indices] if target_library is not None else None
        for gene in pair_genes:
            detected = sample_values[:, pair_position[gene]] > 0
            detection_rows.append(
                {
                    "source": source_id,
                    "sample_id": sample,
                    "target_state": target_state,
                    "gene": gene,
                    "cells": cells,
                    "detected_cells": int(np.sum(detected)),
                    "detection_rate": float(np.mean(detected)),
                    "library_size_field": library_field or "NOT_AVAILABLE",
                }
            )
        for pair_family, gene_a, gene_b in pair_specs:
            values = sample_values[:, [pair_position[gene_a], pair_position[gene_b]]]
            row = _summarize_pair(
                values,
                sample_library,
                gene_a,
                gene_b,
                source_id,
                sample,
                target_state,
                pair_family,
                int(permutation_reps),
                rng,
            )
            row["library_size_field"] = library_field or "NOT_AVAILABLE"
            coexpression_rows.append(row)

    coexpression_path = os.path.join(
        output_dir, "PRJNA607098_S1_cell_level_coexpression_v1.csv"
    )
    _write_rows(
        coexpression_path,
        coexpression_rows,
        [
            "source", "sample_id", "target_state", "pair_family", "gene_a", "gene_b",
            "cells", "gene_a_detection_rate", "gene_b_detection_rate",
            "observed_co_detection_rate", "expected_independence_rate",
            "excess_co_detection_rate", "co_detection_enrichment_ratio",
            "permutation_p_one_sided", "permutation_reps", "library_size_stratified",
            "library_size_field",
        ],
    )
    detection_path = os.path.join(
        output_dir, "PRJNA607098_S1_gene_detection_inventory_v1.csv"
    )
    _write_rows(
        detection_path,
        detection_rows,
        [
            "source", "sample_id", "target_state", "gene", "cells",
            "detected_cells", "detection_rate", "library_size_field",
        ],
    )

    gene_inventory_path = os.path.join(
        output_dir, "PRJNA607098_S1_targeted_gene_coverage_v1.csv"
    )
    _write_rows(
        gene_inventory_path,
        [
            {
                "source": source_id,
                "gene": gene,
                "found_in_atlas": gene in symbol_to_index,
                "atlas_gene_index_zero_based": int(symbol_to_index[gene]) if gene in symbol_to_index else "",
                "used_in_S1_pair_analysis": gene in pair_genes,
            }
            for gene in panel_genes
        ],
        [
            "source", "gene", "found_in_atlas",
            "atlas_gene_index_zero_based", "used_in_S1_pair_analysis",
        ],
    )

    chunk_column_size = int(x_array.chunks[1])
    chunk_rows = []
    for chunk_index in sorted({index // chunk_column_size for index in sorted_indices}):
        relative_key = f"X/0.{chunk_index}"
        chunk_rows.append(
            {
                "zarr_key": relative_key,
                "url": zarr_url.rstrip("/") + "/" + relative_key,
                "access_mode": "user_run_R_remote_stream",
            }
        )
    chunk_path = os.path.join(
        output_dir, "PRJNA607098_S1_targeted_chunk_manifest_v1.csv"
    )
    _write_rows(chunk_path, chunk_rows, ["zarr_key", "url", "access_mode"])

    return {
        "source": source_id,
        "target_state": target_state,
        "source_cells": int(np.sum(source_mask)),
        "target_state_cells": int(np.sum(source_target_mask)),
        "samples": int(len(sample_ids)),
        "found_registry_genes": int(len(found_genes)),
        "registry_genes_requested": int(len(panel_genes)),
        "pair_rows": int(len(coexpression_rows)),
        "library_size_field": library_field or "NOT_AVAILABLE",
        "coexpression": coexpression_path,
        "detection": detection_path,
        "gene_inventory": gene_inventory_path,
        "chunk_manifest": chunk_path,
        "missing_genes": missing_genes,
    }
