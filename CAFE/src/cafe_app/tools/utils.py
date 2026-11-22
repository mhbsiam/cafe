"""Shared CAFE utilities: sparse-matrix handling, cluster-column resolution, CSV-to-AnnData conversion."""
import os
import threading
import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import pyarrow as pa
import pyarrow.csv as pv
from scipy.sparse import issparse
from scipy.stats import mannwhitneyu, kruskal, wasserstein_distance, median_abs_deviation
from statsmodels.stats.multitest import multipletests

import streamlit as st

# Named constants
RANDOM_STATE = 50
# numba's default thread pool aborts the process if two sessions start parallel work at once; hold this around it.
# Leiden holds it too: igraph has one random generator per process, so two runs at once would change each other's labels.
NUMBA_LOCK = threading.Lock()
MARKER_HIGH_THRESHOLD = 1000   # expression above this -> marker annotated with "+"
MARKER_LOW_THRESHOLD = 0       # expression below this -> marker annotated with "-"


def run_umap_parallel(adata, min_dist):
    """UMAP layout on all CPU cores, with sc.tl.umap's defaults; coordinates differ slightly between runs."""
    # sc.tl.umap always seeds the layout, and umap-learn runs a seeded layout on one thread, so call it directly.
    from umap.umap_ import find_ab_params, simplicial_set_embedding

    params = adata.uns["neighbors"]["params"]
    use_rep = params.get("use_rep", "X")
    a, b = find_ab_params(1.0, min_dist)
    with NUMBA_LOCK:
        adata.obsm["X_umap"], _ = simplicial_set_embedding(
            data=adata.X if use_rep == "X" else adata.obsm[use_rep],
            graph=adata.obsp["connectivities"].tocoo(),
            n_components=2,
            initial_alpha=1.0,
            a=a,
            b=b,
            gamma=1.0,
            negative_sample_rate=5,
            n_epochs=500 if adata.n_obs <= 10000 else 200,
            init="spectral",
            random_state=np.random.RandomState(RANDOM_STATE),
            metric=params.get("metric", "euclidean"),
            metric_kwds={},
            densmap=False,
            densmap_kwds={},
            output_dens=False,
            parallel=True,
        )
    adata.uns["umap"] = {"params": {"a": a, "b": b}}


# One Leiden run is a noisy reference: on real data, which resolution looked most stable changed with its seed.
N_REFERENCE = 5


def _neighbor_graph(adata):
    """adata's saved neighbor graph as an igraph graph, the form sc.tl.leiden clusters on."""
    return sc._utils.get_igraph_from_adjacency(adata.obsp["connectivities"], directed=False)


def _subset_graph(adata, idx):
    """Neighbor graph on the cells at idx, built with adata's own neighbor settings."""
    params = adata.uns["neighbors"]["params"]
    use_rep = params.get("use_rep", "X")
    sub = sc.AnnData((adata.X if use_rep == "X" else adata.obsm[use_rep])[idx])
    with NUMBA_LOCK:
        sc.pp.neighbors(
            sub,
            n_neighbors=params["n_neighbors"],
            method="umap",
            use_rep="X",
            metric=params.get("metric", "euclidean"),
            random_state=RANDOM_STATE,
        )
    return _neighbor_graph(sub)


def _leiden_labels(graph, resolution, seed):
    """Leiden labels on an igraph graph; the same call sc.tl.leiden(flavor="igraph", n_iterations=2) makes."""
    from scanpy._utils.random import set_igraph_random_state

    with NUMBA_LOCK, set_igraph_random_state(seed):
        part = graph.community_leiden(
            objective_function="modularity", weights="weight", resolution=resolution, n_iterations=2
        )
    return np.array(part.membership)


def leiden_stability(adata, resolutions, seeds=0, subsamples=10, max_cells=100_000, progress=None):
    """ARI stability of Leiden at each resolution. Returns (table, number of cells used).

    seeds >= 2 adds seed_ari/seed_sd (pairwise ARI between reruns with different seeds on one graph).
    subsamples >= 1 adds sub_ari/sub_sd (ARI between Leiden on a new graph built from a random 80% of
    the cells and each of N_REFERENCE reference runs). Above max_cells, everything runs on a random sample.
    """
    from itertools import combinations
    from sklearn.metrics import adjusted_rand_score

    rng = np.random.default_rng(RANDOM_STATE)
    cells = adata
    if adata.n_obs > max_cells:  # sample once; the sample's graph and its subsamples come from the same cells
        params = adata.uns["neighbors"]["params"]
        use_rep = params.get("use_rep", "X")
        base = np.sort(rng.choice(adata.n_obs, max_cells, replace=False))
        cells = sc.AnnData((adata.X if use_rep == "X" else adata.obsm[use_rep])[base])
        cells.uns["neighbors"] = {"params": {**params, "use_rep": "X"}}
        graph = _subset_graph(cells, np.arange(max_cells))
    else:
        graph = _neighbor_graph(adata)
    n = cells.n_obs
    # Positions within the n cells in use, so reference labels and subsample labels line up.
    sub_idx = [np.sort(rng.choice(n, int(0.8 * n), replace=False)) for _ in range(subsamples)]
    n_runs = max(seeds, N_REFERENCE)
    done, total = 0, subsamples + len(resolutions) * (n_runs + subsamples)

    def tick():
        nonlocal done
        done += 1
        if progress:
            progress(done, total)

    sub_graphs = []
    for idx in sub_idx:
        sub_graphs.append(_subset_graph(cells, idx))
        tick()

    rows = []
    for res in resolutions:
        runs = []
        for s in range(n_runs):
            runs.append(_leiden_labels(graph, res, RANDOM_STATE + s))
            tick()
        row = {"resolution": res, "clusters": len(set(runs[0]))}
        if seeds >= 2:
            aris = [adjusted_rand_score(a, b) for a, b in combinations(runs[:seeds], 2)]
            row["seed_ari"], row["seed_sd"] = np.mean(aris), np.std(aris)
        if subsamples:
            aris = []
            for idx, g in zip(sub_idx, sub_graphs):
                labels = _leiden_labels(g, res, RANDOM_STATE)
                aris += [adjusted_rand_score(ref[idx], labels) for ref in runs[:N_REFERENCE]]
                tick()
            row["sub_ari"], row["sub_sd"] = np.mean(aris), np.std(aris)
        rows.append(row)
    return pd.DataFrame(rows), n


def suggest_resolution(table, seeds=0):
    """One-standard-error rule on the mean of the measures that ran.

    Returns (suggested, lowest, best_score): the stable range is the run of neighbouring resolutions
    around the best one that score within one standard error of it; the suggestion is its highest.
    """
    # The reference run drives most of the spread in sub_ari, so its standard error counts references, not subsamples.
    parts = [(table[f"{k}_ari"], table[f"{k}_sd"] / np.sqrt(n))
             for k, n in (("seed", seeds), ("sub", N_REFERENCE)) if f"{k}_ari" in table]
    table = table.assign(score=np.mean([p[0] for p in parts], axis=0), se=np.mean([p[1] for p in parts], axis=0))
    table = table.sort_values("resolution").reset_index(drop=True)
    best = int(table["score"].idxmax())
    ok = (table["score"] >= table["score"][best] - table["se"][best]).to_numpy()
    lo = hi = best
    while hi + 1 < len(ok) and ok[hi + 1]:
        hi += 1
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    return float(table["resolution"][hi]), float(table["resolution"][lo]), float(table["score"][best])


def cluster_stability(adata, resolution, n_seeds=5):
    """Per-cluster reproducibility of adata.obs['leiden'].

    Reruns Leiden n_seeds times with other seeds on the same graph and averages each cluster's best
    Jaccard overlap with the rerun's clusters (1 = always the same cells, near 0 = arbitrary split).
    """
    graph = _neighbor_graph(adata)
    ref = adata.obs["leiden"].astype(str).to_numpy()
    scores = []
    for s in range(1, n_seeds + 1):
        ct = pd.crosstab(ref, _leiden_labels(graph, resolution, RANDOM_STATE + s))
        inter = ct.to_numpy()
        union = inter.sum(axis=1, keepdims=True) + inter.sum(axis=0, keepdims=True) - inter
        scores.append((inter / union).max(axis=1))
    return pd.DataFrame(
        {"cluster": ct.index, "cells": inter.sum(axis=1), "stability": np.mean(scores, axis=0)}
    ).sort_values("cells", ascending=False, ignore_index=True)


def get_cluster_col(adata):
    """Return 'cell_type' if present in adata.obs, else 'leiden'."""
    return 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'


def safe_flatten(X):
    """1-D array from a sparse or dense matrix (densifies sparse to drop explicit zeros)."""
    if issparse(X):
        return np.asarray(X.toarray()).flatten()
    return np.asarray(X).flatten()


def safe_expression(X, method='Mean'):
    """Mean or median expression, densifying sparse input first."""
    flat = safe_flatten(X)
    if method == 'Median':
        return float(np.median(flat))
    return float(np.mean(flat))


def safe_sort_clusters(labels):
    """Sort cluster labels by int when all are integers, else by string."""
    labels = list(labels)
    try:
        return sorted(labels, key=int)
    except (ValueError, TypeError):
        return sorted(labels, key=str)


# CSV -> AnnData conversion
def _process_single_table(table, file_name):
    """Cast types, replace SampleID/Group from filename, drop 'sample' cols; (table, None) or (None, msg)."""
    # Cast integer columns to float64
    numeric_columns = [col for col in table.column_names
                       if pa.types.is_integer(table.schema.field(col).type)]
    for col in numeric_columns:
        table = table.set_column(table.column_names.index(col), col,
                                 table.column(col).cast(pa.float64()))

    # Cast existing SampleID/Group to string (if present)
    if 'SampleID' in table.column_names:
        table = table.set_column(table.column_names.index('SampleID'), 'SampleID',
                                 table.column('SampleID').cast(pa.string()))
    if 'Group' in table.column_names:
        table = table.set_column(table.column_names.index('Group'), 'Group',
                                 table.column('Group').cast(pa.string()))

    # B11 fix: drop existing SampleID/Group before appending filename-derived ones
    if 'SampleID' in table.column_names:
        table = table.drop(['SampleID'])
    if 'Group' in table.column_names:
        table = table.drop(['Group'])

    # Parse filename: expected format 'sampleID_group.csv'
    stem = os.path.splitext(file_name)[0] if file_name.lower().endswith('.csv') else file_name
    name_parts = stem.split('_')
    if len(name_parts) != 2:
        return None, f"Skipping file {file_name} as it does not have the expected format 'sampleID_group.csv'"

    table = table.append_column('SampleID', pa.array([name_parts[0]] * len(table)))
    table = table.append_column('Group', pa.array([name_parts[1]] * len(table)))

    # Drop any column named 'sample' (case-insensitive)
    columns_to_drop = [col for col in table.column_names if col.lower() == 'sample']
    if columns_to_drop:
        table = table.drop(columns_to_drop)

    return table, None


def _combine_tables(named_tables):
    """[(file_name, arrow table)] -> (DataFrame or None, skip messages, notices); shared by the app loaders and hpc.py.

    Every file must have the first valid file's columns (reordered to match); markers are cast to float64."""
    tables, skipped, notices = [], [], []
    reference = None
    for name, table in named_tables:
        # All-blank or unnamed columns (e.g. a trailing comma) read as null-typed and would drop every row.
        keep = [i for i, f in enumerate(table.schema) if f.name.strip() and not pa.types.is_null(f.type)]
        if len(keep) < table.num_columns:
            notices.append(f"{name}: dropped {table.num_columns - len(keep)} empty column(s)")
            table = table.select(keep)

        processed, skip_msg = _process_single_table(table, name)
        if processed is None:
            skipped.append(skip_msg)
            continue
        if reference is None:
            # Text columns of the first valid file (e.g. "Time", "Well") are not markers; drop them from every file.
            non_numeric = [f.name for f in processed.schema if f.name not in ('SampleID', 'Group')
                           and not (pa.types.is_integer(f.type) or pa.types.is_floating(f.type))]
        processed = processed.drop([c for c in non_numeric if c in processed.column_names])
        if processed.num_columns <= 2:
            skipped.append(f"Skipping file {name}: it has no marker data.")
            continue
        if reference is not None and set(processed.column_names) != set(reference):
            skipped.append(f"Skipping file {name}: its columns do not match the first file's columns.")
            continue
        try:
            if reference is not None:
                processed = processed.select(reference)
            for i, col in enumerate(processed.column_names):
                if col not in ('SampleID', 'Group'):
                    processed = processed.set_column(i, col, processed.column(col).cast(pa.float64()))
        except pa.ArrowException as exc:
            skipped.append(f"Skipping file {name}: marker values must be numeric ({exc}).")
            continue
        if reference is None:
            reference = processed.column_names
            if non_numeric:
                notices.append(f"Dropped non-numeric column(s): {', '.join(non_numeric)}")
        tables.append(processed)

    if not tables:
        return None, skipped, notices

    df = pa.concat_tables(tables).to_pandas()
    n_rows = len(df)
    df = df.dropna()
    if len(df) < n_rows:
        notices.append(f"Dropped {n_rows - len(df)} of {n_rows} rows with missing values")
    if df.empty:
        notices.append("No rows left after removing rows with missing values")
        return None, skipped, notices
    return df, skipped, notices


def _report_load(skipped, notices, save_notice=False):
    """Show loader messages; with save_notice they also go to st.session_state['_load_notice'] for display after a rerun."""
    msgs = skipped + notices
    for msg in msgs:
        st.warning(msg)
    # Only Data Processing reads the key; other pages' loads must not overwrite or clear its unseen notice.
    if not save_notice:
        return
    # Replace, never append: some pages reload on every rerun, and a clean load must not keep an older notice.
    if msgs:
        st.session_state["_load_notice"] = "; ".join(msgs)
    else:
        st.session_state.pop("_load_notice", None)


def _df_to_anndata(df):
    """Convert a combined DataFrame (markers + SampleID/Group) to AnnData."""
    expr_data = df.drop(columns=['SampleID', 'Group'])
    expr_data = expr_data.select_dtypes(include=[float, int])

    metadata = df[['SampleID', 'Group']]
    expr_data.index = expr_data.index.astype(str)
    metadata.index = metadata.index.astype(str)

    adata = sc.AnnData(expr_data)
    adata.obs = metadata
    adata.var_names = expr_data.columns.astype(str)
    return adata


def process_uploaded_csv_files(uploaded_files):
    """Streamlit-uploaded CSVs into one AnnData, with a progress bar and column-consistency checks."""
    start_time = time.time()
    named_tables = []

    progress_bar = st.progress(0)
    elapsed_time_placeholder = st.empty()
    status_placeholder = st.empty()

    total_steps = len(uploaded_files) + 2
    current_step = 0

    for idx, uploaded_file in enumerate(uploaded_files):
        current_step += 1
        progress_bar.progress(current_step / total_steps)
        status_placeholder.write(f"Processing file {idx + 1}/{len(uploaded_files)}: {uploaded_file.name}")
        named_tables.append((uploaded_file.name, pv.read_csv(uploaded_file)))

    current_step += 1
    progress_bar.progress(current_step / total_steps)
    status_placeholder.write("Combining all files into a single table")

    df, skipped, notices = _combine_tables(named_tables)
    _report_load(skipped, notices, save_notice=True)
    if df is None:
        progress_bar.progress(1.0)
        return None

    adata = _df_to_anndata(df)

    current_step += 1
    progress_bar.progress(1.0)
    elapsed_time = time.time() - start_time
    elapsed_time_placeholder.write(f"Data loading and processing completed in {elapsed_time:.2f} seconds")
    status_placeholder.write("Processing complete.")

    return adata


def inspect_uploaded_csv_files(uploaded_files):
    """Diagnostics dict for the upload-review UI (per-file info, marker mismatches, bad names); builds no AnnData."""
    files = []
    marker_sets = {}
    for uploaded_file in uploaded_files:
        try:
            table = pv.read_csv(uploaded_file)
        finally:
            # Rewind (best-effort) so the buffer can be re-read when AnnData is built.
            try:
                uploaded_file.seek(0)
            except (AttributeError, OSError, ValueError):
                pass

        # Markers: columns that aren't metadata, a stray 'sample' column, or empty (dropped on load).
        markers = [
            f.name for f in table.schema
            if f.name not in ("SampleID", "Group") and f.name.lower() != "sample"
            and f.name.strip() and not pa.types.is_null(f.type)
        ]

        name = uploaded_file.name
        stem = os.path.splitext(name)[0] if name.lower().endswith(".csv") else name
        name_parts = stem.split("_")
        name_ok = len(name_parts) == 2
        sample = name_parts[0] if name_ok else stem
        group = name_parts[1] if name_ok else "N/A"

        files.append({
            "file": name,
            "SampleID": sample,
            "Group": group,
            "n_cells": table.num_rows,
            "markers": markers,
            "name_ok": name_ok,
        })
        marker_sets[name] = set(markers)

    sets = list(marker_sets.values())
    union = set().union(*sets) if sets else set()
    common = set.intersection(*sets) if sets else set()

    # Flag a file that lacks any marker present in another (catches case/spelling diffs).
    mismatches = [
        {"file": name, "missing": sorted(union - s)}
        for name, s in marker_sets.items()
        if union - s
    ]
    bad_names = [f["file"] for f in files if not f["name_ok"]]

    # Duplicate file names or SampleIDs would merge into one sample; block them via bad_names.
    duplicates, seen_files, seen_ids = [], set(), {}
    for f in files:
        if f["file"] in seen_files:
            duplicates.append(f"{f['file']} (same file name uploaded twice)")
        elif f["name_ok"] and f["SampleID"] in seen_ids:
            duplicates.append(f"{f['file']} (SampleID '{f['SampleID']}' also used by {seen_ids[f['SampleID']]})")
        seen_files.add(f["file"])
        if f["name_ok"]:
            seen_ids.setdefault(f["SampleID"], f["file"])

    return {
        "files": files,
        "all_markers": sorted(union),
        "common_markers": sorted(common),
        "mismatches": mismatches,
        "bad_names": bad_names + duplicates,
        "duplicates": duplicates,
    }


def process_uploaded_csv_to_df(uploaded_files):
    """Like process_uploaded_csv_files but returns the raw DataFrame (with SampleID/Group) instead of AnnData."""
    named_tables = []

    progress_bar = st.progress(0)
    status_placeholder = st.empty()

    total_steps = len(uploaded_files) + 2
    current_step = 0

    for idx, uploaded_file in enumerate(uploaded_files):
        current_step += 1
        progress_bar.progress(current_step / total_steps)
        status_placeholder.write(f"Processing file {idx + 1}/{len(uploaded_files)}: {uploaded_file.name}")
        named_tables.append((uploaded_file.name, pv.read_csv(uploaded_file)))

    current_step += 1
    progress_bar.progress(current_step / total_steps)
    status_placeholder.write("Combining all files into a single table")

    df, skipped, notices = _combine_tables(named_tables)
    _report_load(skipped, notices)
    if df is None:
        progress_bar.progress(1.0)
        return None

    progress_bar.progress(1.0)
    status_placeholder.write("Processing complete.")

    return df


def process_directory_csv_files(input_dir):
    """CSV files from a directory into one AnnData; None if no valid CSVs are found."""
    start_time = time.time()
    named_tables = []

    csv_files = [f for f in os.listdir(input_dir) if f.lower().endswith('.csv')]
    if not csv_files:
        st.write("No CSV files found in the input directory.")
        return None

    for idx, csv_file in enumerate(csv_files):
        st.write(f"Processing file {idx + 1}/{len(csv_files)}: {csv_file}")
        file_path = os.path.join(input_dir, csv_file)
        named_tables.append((csv_file, pv.read_csv(file_path)))

    df, skipped, notices = _combine_tables(named_tables)
    _report_load(skipped, notices)
    if df is None:
        st.write("No valid CSV files were processed.")
        return None

    adata = _df_to_anndata(df)
    st.write(f"Data loading and processing completed in {time.time() - start_time:.2f} seconds")
    return adata


# Batch correction (ComBat), guarded against confounding with biology.
def run_combat_correction(adata, batch_key='Batch', covariates=('Group',)):
    """ComBat on a technical Batch; refuses (ok=False) if missing/unassigned/single-level/confounded. Returns (ok, msg)."""
    if batch_key not in adata.obs.columns:
        return False, f"Batch column '{batch_key}' not found; assign batches first."

    if adata.obs[batch_key].isnull().any():
        return False, (
            f"Some cells have no '{batch_key}' assigned. Assign every sample to a batch "
            "before running ComBat."
        )

    n_batches = adata.obs[batch_key].nunique()
    if n_batches < 2:
        return False, "Only one batch level; nothing to correct."

    # Confounding guard: a batch nested within one covariate level is uncorrectable.
    for cov in covariates:
        if cov not in adata.obs.columns:
            continue
        if adata.obs[cov].nunique() < 2:
            continue
        crosstab = pd.crosstab(adata.obs[batch_key], adata.obs[cov])
        batches_spanning = (crosstab > 0).sum(axis=1)  # covariate levels per batch
        if (batches_spanning <= 1).all():
            return False, (
                f"Batch is confounded with {cov}; ComBat would remove biological "
                f"signal. Assign batches that each span multiple {cov} levels."
            )

    protected = [c for c in covariates if c in adata.obs.columns and adata.obs[c].nunique() >= 2]
    # ComBat needs numeric covariates: k-1 dummy columns per covariate (a raw string column is singular).
    dummies = (pd.get_dummies(adata.obs[protected].astype(str), drop_first=True, dtype=float)
               if protected else pd.DataFrame())
    temp_cols = [f"_combat_cov_{i}" for i in range(dummies.shape[1])]
    for name, col in zip(temp_cols, dummies.columns):
        adata.obs[name] = dummies[col].to_numpy()
    try:
        sc.pp.combat(adata, key=batch_key, covariates=temp_cols or None)
    except (np.linalg.LinAlgError, ValueError):
        return False, (
            f"Batch is confounded with {', '.join(protected) or 'the design'} (singular ComBat "
            "design), so batch correction was not applied."
        )
    finally:
        adata.obs.drop(columns=temp_cols, inplace=True)
    return True, (
        f"ComBat applied across {n_batches} batches"
        + (f", with {', '.join(protected)} passed to ComBat as covariate(s) to protect that signal."
           if protected else ".")
    )


# QC: data-driven MAD threshold recommendation.
def recommend_nmads(values, cap_pct=5.0, use_lower=True, use_upper=True,
                    grid=None):
    """Smallest n_MADs on grid keeping each active tail <= cap_pct; grid.max() if none, None if degenerate."""
    if grid is None:
        grid = np.arange(2.0, 8.0001, 0.5)
    grid = np.asarray(grid, dtype=float)

    vals = np.asarray(values, dtype=float)
    n = vals.shape[0]
    if n == 0:
        return None
    med = np.median(vals)
    mad = median_abs_deviation(vals)
    if mad == 0 or np.isnan(mad):
        return None

    def tails_ok(nm):
        if use_lower:
            low_pct = np.count_nonzero(vals < med - nm * mad) / n * 100.0
            if low_pct > cap_pct:
                return False
        if use_upper:
            high_pct = np.count_nonzero(vals > med + nm * mad) / n * 100.0
            if high_pct > cap_pct:
                return False
        return True

    for nm in np.sort(grid):
        if tails_ok(nm):
            return float(nm)
    return float(np.max(grid))


# Batch effect diagnostics: per-marker EMD (1-D Wasserstein) of each batch vs pooled, via scipy.
def _robust_scale(col):
    """Robust per-marker scale (normalised MAD in SD units), falling back to SD then 1.0."""
    scale = median_abs_deviation(col, scale="normal")
    if scale == 0 or np.isnan(scale):
        scale = np.std(col)
    if scale == 0 or np.isnan(scale):
        scale = 1.0
    return scale


def marker_robust_scales(adata):
    """Per-marker robust scales as a 1-D array; compute before correction and reuse for before/after EMD."""
    X = adata.X
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=float)
    return np.array([_robust_scale(X[:, j]) for j in range(X.shape[1])])


def compute_batch_emd(
    adata, batch_key="Batch", max_cells_per_batch=20000, seed=RANDOM_STATE, scales=None
):
    """Per-marker EMD (robust-scaled 1-D Wasserstein) of each batch vs pooled; returns marker/emd_max/emd_mean/worst_batch, subsampled per batch."""
    empty = pd.DataFrame(columns=["marker", "emd_max", "emd_mean", "worst_batch"])
    if batch_key not in adata.obs.columns:
        return empty

    batch_values = pd.Series(adata.obs[batch_key].to_numpy())
    levels = [b for b in pd.unique(batch_values) if pd.notnull(b)]
    if len(levels) < 2:
        return empty

    X = adata.X
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    rng = np.random.default_rng(seed)
    batch_idx = {}
    for b in levels:
        idx = np.where(batch_values.to_numpy() == b)[0]
        if idx.shape[0] > max_cells_per_batch:
            idx = rng.choice(idx, size=max_cells_per_batch, replace=False)
        batch_idx[b] = idx
    ref_idx = np.concatenate([batch_idx[b] for b in levels])

    rows = []
    for j, marker in enumerate(adata.var_names):
        col = X[:, j]
        scale = _robust_scale(col) if scales is None else scales[j]
        scaled = col / scale
        pooled_ref = scaled[ref_idx]
        dists = np.array(
            [wasserstein_distance(scaled[batch_idx[b]], pooled_ref) for b in levels]
        )
        rows.append(
            {
                "marker": str(marker),
                "emd_max": float(dists.max()),
                "emd_mean": float(dists.mean()),
                "worst_batch": str(levels[int(dists.argmax())]),
            }
        )

    return (
        pd.DataFrame(rows)
        .sort_values("emd_max", ascending=False)
        .reset_index(drop=True)
    )


def summarize_batch_effect(emd_df, adata, batch_key="Batch", covariates=("Group",), flag=0.5):
    """EMD table -> (verdict, info_card kind, n_flagged markers over `flag`, confounded covariates)."""
    confounded = []
    for cov in covariates:
        if cov not in adata.obs.columns:
            continue
        if adata.obs[cov].nunique() < 2:
            continue
        crosstab = pd.crosstab(adata.obs[batch_key], adata.obs[cov])
        batches_spanning = (crosstab > 0).sum(axis=1)
        if (batches_spanning <= 1).all():
            confounded.append(cov)

    if emd_df is None or emd_df.empty:
        return (
            "Need at least two assigned batches to check for batch effects.",
            "info",
            0,
            confounded,
        )

    n_markers = len(emd_df)
    n_flagged = int((emd_df["emd_max"] > flag).sum())
    top = emd_df.iloc[0]

    if n_flagged == 0:
        verdict = (
            f"Batches look well mixed. No marker exceeds an EMD of {flag:g} "
            f"(largest: {top['marker']} at {top['emd_max']:.2f}). Batch correction "
            "may be unnecessary and could remove real signal."
        )
        kind = "success"
    else:
        verdict = (
            f"{n_flagged} of {n_markers} markers show substantial batch separation "
            f"(EMD above {flag:g}). The strongest is {top['marker']}, shifted most in "
            f"batch {top['worst_batch']} (EMD {top['emd_max']:.2f}). Batch correction "
            "is likely warranted."
        )
        kind = "warning"

    if confounded:
        verdict += (
            " Note: batch is confounded with "
            + ", ".join(confounded)
            + ", so some of this separation may reflect biology rather than a "
            "technical effect. Correcting it could remove the differences you are "
            "testing for."
        )
        kind = "warning"

    return verdict, kind, n_flagged, confounded


# Differential abundance testing: aggregate to per-sample proportions (the correct unit, not per-cell).
def compute_sample_proportions(adata, cluster_col):
    """Per-sample cluster-proportion matrix (rows sum to 1, zero-cell combos kept as 0) and SampleID->Group."""
    obs = adata.obs[['SampleID', 'Group', cluster_col]].copy()
    obs[cluster_col] = obs[cluster_col].astype(str)

    # Per-sample x cluster counts; unstack fills unseen combos with 0.
    counts = (
        obs.groupby(['SampleID', cluster_col], observed=True)
           .size()
           .unstack(fill_value=0)
    )
    totals = counts.sum(axis=1)
    prop_df = counts.div(totals, axis=0).fillna(0.0)

    sample_group = (
        obs.drop_duplicates('SampleID')
           .set_index('SampleID')['Group']
           .reindex(prop_df.index)
    )
    return prop_df, sample_group


def compute_sample_frequencies_long(adata, cluster_col, as_percent=True):
    """Long-form frequency table (SampleID, Group, <cluster_col>, frequency); percent when as_percent, else fraction."""
    prop_df, sample_group = compute_sample_proportions(adata, cluster_col)
    if as_percent:
        prop_df = prop_df * 100
    long = (
        prop_df.reset_index()
               .melt(id_vars='SampleID', var_name=cluster_col,
                     value_name='frequency')
    )
    long['Group'] = long['SampleID'].map(sample_group)
    return long


def differential_abundance_test(prop_df, sample_group, alpha=0.05, min_samples=2):
    """Per-cluster abundance test across groups (Mann-Whitney/Kruskal + BH-FDR); NaN when a group has < min_samples."""
    groups = sorted(sample_group.dropna().unique(), key=str)
    rows = []
    for cluster in prop_df.columns:
        row = {'cluster': cluster}
        group_vectors = {}
        for g in groups:
            samples = sample_group.index[sample_group == g]
            vals = prop_df.loc[samples, cluster].dropna().to_numpy()
            group_vectors[g] = vals
            row[f'mean_prop_{g}'] = float(np.mean(vals)) if len(vals) else np.nan

        enough = len(groups) >= 2 and all(len(v) >= min_samples for v in group_vectors.values())
        stat, p = np.nan, np.nan
        if enough:
            try:
                if len(groups) == 2:
                    stat, p = mannwhitneyu(group_vectors[groups[0]],
                                           group_vectors[groups[1]],
                                           alternative='two-sided')
                else:
                    stat, p = kruskal(*group_vectors.values())
            except ValueError:
                stat, p = np.nan, np.nan  # e.g. all proportions identical

        means = {g: row[f'mean_prop_{g}'] for g in groups}
        if all(m == m for m in means.values()):  # no NaNs
            row['direction'] = f'enriched in {max(means, key=means.get)}'
        else:
            row['direction'] = 'n/a'
        row['statistic'] = float(stat) if stat == stat else np.nan
        row['p_value'] = float(p) if p == p else np.nan
        rows.append(row)

    results = pd.DataFrame(rows)
    results['p_adj'] = np.nan
    valid = results['p_value'].notna()
    if valid.any():
        _, p_adj, _, _ = multipletests(results.loc[valid, 'p_value'].to_numpy(),
                                       alpha=alpha, method='fdr_bh')
        results.loc[valid, 'p_adj'] = p_adj
    results['significant'] = results['p_adj'] < alpha
    return results.sort_values('p_adj', na_position='last').reset_index(drop=True)


# Path validation (traversal prevention).
def validate_path(path, base):
    """Return abs path if it resolves inside base, else raise ValueError."""
    abs_path = os.path.abspath(path)
    if not abs_path.startswith(base + os.sep) and abs_path != base:
        raise ValueError(f"Path '{path}' is outside the allowed directory.")
    return abs_path


def user_dir(path):
    """Abs path for a directory the user typed, "" if blank; any location is fine (local desktop/HPC app)."""
    # Windows "Copy as path" wraps the path in quotes.
    path = path.strip().strip('"\'').strip()
    # abspath("") is the cwd, so keep blank input blank for the caller's empty check.
    return os.path.abspath(os.path.expanduser(path)) if path else ""


# Standardised plot controls.
@dataclass
class PlotSettings:
    """Values returned by plot_controls; unrequested fields stay None."""
    width: Optional[float] = None
    height: Optional[float] = None
    dot_size: Optional[float] = None
    cmap: Optional[str] = None
    file_format: Optional[str] = None


def plot_controls(
    key_prefix,
    *,
    include=("width", "height", "dot_size", "cmap", "file_format"),
    width=(5, 20, 8),
    height=(5, 20, 6),
    dot_size=(1, 100, 5),
    dot_size_step=1,
    size_unit="inches",
    colormaps=None,
    default_cmap=None,
    file_formats=("PNG", "PDF", "SVG", "JPEG"),
):
    """Render plot controls -> PlotSettings; key_prefix must be unique, include selects which controls to show."""
    s = PlotSettings()
    unit = f" ({size_unit})" if size_unit else ""

    size_row = [c for c in ("width", "height", "dot_size") if c in include]
    if size_row:
        cols = st.columns(len(size_row))
        for col, name in zip(cols, size_row):
            with col:
                if name == "width":
                    s.width = st.slider(f"Plot width{unit}", width[0], width[1], width[2],
                                        key=f"{key_prefix}_width")
                elif name == "height":
                    s.height = st.slider(f"Plot height{unit}", height[0], height[1], height[2],
                                         key=f"{key_prefix}_height")
                elif name == "dot_size":
                    s.dot_size = st.slider("Dot size", dot_size[0], dot_size[1], dot_size[2],
                                           step=dot_size_step, key=f"{key_prefix}_dotsize")

    opt_row = [c for c in ("cmap", "file_format") if c in include]
    if opt_row:
        cols = st.columns(len(opt_row))
        for col, name in zip(cols, opt_row):
            with col:
                if name == "cmap":
                    opts = list(colormaps or [])
                    idx = opts.index(default_cmap) if default_cmap in opts else 0
                    s.cmap = st.selectbox("Colormap", options=opts, index=idx,
                                          key=f"{key_prefix}_cmap")
                elif name == "file_format":
                    s.file_format = st.radio("File format", tuple(file_formats),
                                             key=f"{key_prefix}_format")
    return s
