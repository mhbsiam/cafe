"""Shared CAFE utilities: sparse-matrix handling, cluster-column resolution, CSV-to-AnnData conversion."""
import os
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
MARKER_HIGH_THRESHOLD = 1000   # expression above this -> marker annotated with "+"
MARKER_LOW_THRESHOLD = 0       # expression below this -> marker annotated with "-"


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
    name_parts = file_name.replace('.csv', '').split('_')
    if len(name_parts) != 2:
        return None, f"Skipping file {file_name} as it does not have the expected format 'sampleID_group.csv'"

    table = table.append_column('SampleID', pa.array([name_parts[0]] * len(table)))
    table = table.append_column('Group', pa.array([name_parts[1]] * len(table)))

    # Drop any column named 'sample' (case-insensitive)
    columns_to_drop = [col for col in table.column_names if col.lower() == 'sample']
    if columns_to_drop:
        table = table.drop(columns_to_drop)

    return table, None


def _tables_to_anndata(tables):
    """Combine Arrow tables and convert to AnnData."""
    combined_table = pa.concat_tables(tables)
    df = combined_table.to_pandas()

    if df.isnull().values.any():
        st.error("Missing values detected. Dropping missing rows.")
        df = df.dropna()

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
    tables = []
    reference_columns = None

    progress_bar = st.progress(0)
    elapsed_time_placeholder = st.empty()
    status_placeholder = st.empty()

    total_steps = len(uploaded_files) + 2
    current_step = 0

    for idx, uploaded_file in enumerate(uploaded_files):
        current_step += 1
        progress_bar.progress(current_step / total_steps)
        status_placeholder.write(f"Processing file {idx + 1}/{len(uploaded_files)}: {uploaded_file.name}")
        table = pv.read_csv(uploaded_file)

        if reference_columns is None:
            reference_columns = set(table.column_names)
        else:
            current_columns = set(table.column_names)
            if current_columns != reference_columns:
                st.warning(f"Warning: Column names in {uploaded_file.name} do not match the reference columns!")
                st.write(f"Found columns: {current_columns}")
                continue

        processed, skip_msg = _process_single_table(table, uploaded_file.name)
        if processed is None:
            st.write(skip_msg)
            continue
        tables.append(processed)

    if not tables:
        progress_bar.progress(1.0)
        return None

    current_step += 1
    progress_bar.progress(current_step / total_steps)
    status_placeholder.write("Combining all files into a single table")

    adata = _tables_to_anndata(tables)

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

        # Markers: columns that aren't metadata or a stray 'sample' column.
        markers = [
            c for c in table.column_names
            if c not in ("SampleID", "Group") and c.lower() != "sample"
        ]

        name = uploaded_file.name
        name_parts = name.replace(".csv", "").split("_")
        name_ok = len(name_parts) == 2
        sample = name_parts[0] if name_ok else name.replace(".csv", "")
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

    return {
        "files": files,
        "all_markers": sorted(union),
        "common_markers": sorted(common),
        "mismatches": mismatches,
        "bad_names": bad_names,
    }


def process_uploaded_csv_to_df(uploaded_files):
    """Like process_uploaded_csv_files but returns the raw DataFrame (with SampleID/Group) instead of AnnData."""
    tables = []
    reference_columns = None

    progress_bar = st.progress(0)
    status_placeholder = st.empty()

    total_steps = len(uploaded_files) + 2
    current_step = 0

    for idx, uploaded_file in enumerate(uploaded_files):
        current_step += 1
        progress_bar.progress(current_step / total_steps)
        status_placeholder.write(f"Processing file {idx + 1}/{len(uploaded_files)}: {uploaded_file.name}")
        table = pv.read_csv(uploaded_file)

        if reference_columns is None:
            reference_columns = set(table.column_names)
        else:
            current_columns = set(table.column_names)
            if current_columns != reference_columns:
                st.warning(f"Warning: Column names in {uploaded_file.name} do not match the reference columns!")
                st.write(f"Found columns: {current_columns}")
                continue

        processed, skip_msg = _process_single_table(table, uploaded_file.name)
        if processed is None:
            st.write(skip_msg)
            continue
        tables.append(processed)

    if not tables:
        progress_bar.progress(1.0)
        return None

    current_step += 1
    progress_bar.progress(current_step / total_steps)
    status_placeholder.write("Combining all files into a single table")

    combined_table = pa.concat_tables(tables)
    df = combined_table.to_pandas()

    if df.isnull().values.any():
        st.error("Missing values detected. Dropping missing rows.")
        df = df.dropna()

    progress_bar.progress(1.0)
    status_placeholder.write("Processing complete.")

    return df


def process_directory_csv_files(input_dir):
    """CSV files from a directory into one AnnData; None if no valid CSVs are found."""
    start_time = time.time()
    tables = []

    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    if not csv_files:
        st.write("No CSV files found in the input directory.")
        return None

    for idx, csv_file in enumerate(csv_files):
        st.write(f"Processing file {idx + 1}/{len(csv_files)}: {csv_file}")
        file_path = os.path.join(input_dir, csv_file)
        table = pv.read_csv(file_path)

        processed, skip_msg = _process_single_table(table, csv_file)
        if processed is None:
            st.write(skip_msg)
            continue
        tables.append(processed)

    if not tables:
        st.write("No valid CSV files were processed.")
        return None

    adata = _tables_to_anndata(tables)
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

    protected = [c for c in covariates if c in adata.obs.columns]
    sc.pp.combat(adata, key=batch_key)
    return True, (
        f"ComBat applied across {n_batches} batches"
        + (f" (crossed with {', '.join(protected)}, so that signal is retained)."
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
def validate_path(path, base=None):
    """Return abs path if it resolves inside base (default cwd), else raise ValueError."""
    if base is None:
        base = os.path.abspath(os.getcwd())
    abs_path = os.path.abspath(path)
    if not abs_path.startswith(base + os.sep) and abs_path != base:
        raise ValueError(f"Path '{path}' is outside the allowed directory.")
    return abs_path


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
