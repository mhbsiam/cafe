"""Shared utilities for CAFE.

Centralises repeated patterns so that every page and module uses the same
sparse-matrix handling, cluster-column resolution, and CSV-to-AnnData
conversion logic.
"""
import os
import time
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import pyarrow as pa
import pyarrow.csv as pv
from scipy.sparse import issparse
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests

import streamlit as st

# ---------------------------------------------------------------------------
# Named constants (replaces magic numbers scattered across the codebase)
# ---------------------------------------------------------------------------
RANDOM_STATE = 50
MARKER_HIGH_THRESHOLD = 1000   # expression above this -> marker annotated with "+"
MARKER_LOW_THRESHOLD = 0       # expression below this -> marker annotated with "-"

# ---------------------------------------------------------------------------
# Cluster-column resolution
# ---------------------------------------------------------------------------
def get_cluster_col(adata):
    """Return ``'cell_type'`` if present in ``adata.obs``, else ``'leiden'``."""
    return 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'


# ---------------------------------------------------------------------------
# Sparse / dense safe helpers (B7 / B8)
# ---------------------------------------------------------------------------
def safe_flatten(X):
    """Return a 1-D ``np.ndarray`` from either a sparse or dense matrix.

    ``.X.flatten()`` on a sparse matrix returns a ``np.matrix`` with wrong
    semantics (includes explicit zeros).  ``.X.toarray()`` on a dense matrix
    raises ``AttributeError``.  This helper handles both cases.
    """
    if issparse(X):
        return np.asarray(X.toarray()).flatten()
    return np.asarray(X).flatten()


def safe_expression(X, method='Mean'):
    """Compute mean or median expression, working for sparse **and** dense.

    On a sparse matrix ``X.mean()`` returns a 1x1 ``np.matrix`` and
    ``np.median(X)`` operates on the raw sparse data (including explicit
    zeros), both of which are wrong.  This densifies first.
    """
    flat = safe_flatten(X)
    if method == 'Median':
        return float(np.median(flat))
    return float(np.mean(flat))


# ---------------------------------------------------------------------------
# Cluster-label sorting (B9)
# ---------------------------------------------------------------------------
def safe_sort_clusters(labels):
    """Sort cluster labels, using ``key=int`` only when all labels are integers.

    Falls back to natural (string) sort when labels are non-integer (e.g.
    after cluster merging or cell-type annotation).
    """
    labels = list(labels)
    try:
        return sorted(labels, key=int)
    except (ValueError, TypeError):
        return sorted(labels, key=str)


# ---------------------------------------------------------------------------
# CSV -> AnnData conversion (A4 — replaces 6 copy-pasted implementations)
# ---------------------------------------------------------------------------
def _process_single_table(table, file_name):
    """Process one Arrow table: cast types, fix duplicate SampleID/Group (B11),
    append SampleID/Group from filename, drop 'sample' columns.

    Returns the processed Arrow table or ``None`` if the filename is invalid.
    """
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
    """Process Streamlit-uploaded CSV files into a single AnnData.

    Used by Data_Processing.py and Selective_Clustering.py.
    Shows a progress bar and validates column consistency across files.
    """
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


def process_uploaded_csv_to_df(uploaded_files):
    """Process Streamlit-uploaded CSV files into a single pandas DataFrame.

    Like :func:`process_uploaded_csv_files` but returns the raw DataFrame
    (with SampleID/Group columns) instead of an AnnData.  Used by
    Selective_Clustering.py which needs the DataFrame for marker subsetting.
    """
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
    """Process CSV files from a filesystem directory into a single AnnData.

    Used by Cluster_Evaluation.py and cafe_hpc.py.
    Returns ``None`` if no valid CSV files are found.
    """
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


# ---------------------------------------------------------------------------
# Batch correction (ComBat) — guarded against confounding with biology
# ---------------------------------------------------------------------------
def run_combat_correction(adata, batch_key='Batch', covariates=('Group',)):
    """Run ComBat keyed on a technical batch without destroying biological signal.

    ComBat centers each ``batch_key`` level to a common mean.  Keying it on the
    biological condition (``Group``) would remove the very signal downstream
    comparisons test, so this helper requires a dedicated technical ``Batch``
    column and *refuses* when that batch is confounded with any of ``covariates``
    (i.e. every batch lies entirely within a single covariate level).  When the
    batch is genuinely crossed with the biological groups, centering the batches
    leaves the between-group differences intact.

    Note: scanpy's ``covariates=`` argument builds a rank-deficient design for a
    categorical covariate (batch and covariate one-hots are collinear), so it is
    not passed through here; ``covariates`` is used only for the confounding
    guard.

    Returns ``(ok, message)``.  ``adata`` is only mutated when ``ok`` is True.
    Correction is refused (``ok=False``) when the batch column is missing, has
    unassigned cells, has fewer than two levels, or is confounded.
    """
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

    # Confounding guard: a batch nested entirely within one covariate level means
    # ComBat cannot separate technical from biological variation.
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


# ---------------------------------------------------------------------------
# Differential abundance testing (fixes the pseudoreplicated per-cell χ²)
# ---------------------------------------------------------------------------
# The old tab34 test ran chi2_contingency on a per-*cell* contingency table.
# With hundreds of thousands of cells the χ² statistic is enormous and p≈0 for
# essentially any dataset (pseudoreplication — the unit of analysis must be the
# *sample*, not the cell).  These helpers instead aggregate cells to per-sample
# cluster proportions and compare those proportions across groups, which is the
# statistically correct unit for compositional/abundance testing.
def compute_sample_proportions(adata, cluster_col):
    """Aggregate cells to a per-sample cluster-proportion matrix.

    Returns ``(prop_df, sample_group)`` where ``prop_df`` is indexed by
    ``SampleID`` with one column per cluster and rows summing to 1 (each value
    is the fraction of that sample's cells falling in that cluster), and
    ``sample_group`` maps ``SampleID`` -> ``Group``.

    Samples with zero cells in a cluster contribute a proportion of **0** (they
    are not dropped) — dropping them would bias the downstream abundance test.
    """
    obs = adata.obs[['SampleID', 'Group', cluster_col]].copy()
    obs[cluster_col] = obs[cluster_col].astype(str)

    # Per-sample x cluster cell counts; unstack fills unseen combos with 0.
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
    """Long-form per-sample cluster frequency table.

    Columns: SampleID, Group, <cluster_col>, frequency.
    Built on compute_sample_proportions, so zero-cell (sample, cluster)
    combinations are retained as 0 and the denominator is each sample's
    total cells across all clusters.  frequency is a percentage (0-100)
    when as_percent (default), else a fraction (0-1).
    """
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
    """Per-cluster differential abundance across groups (sample-level).

    For each cluster, the per-sample proportions are compared across groups
    with Mann-Whitney U (exactly 2 groups) or Kruskal-Wallis (>2 groups), then
    p-values are corrected across clusters with Benjamini-Hochberg FDR.

    A cluster is only tested when every group has at least ``min_samples``
    samples with data; otherwise its p-value is ``NaN`` (and excluded from the
    FDR correction).  Returns a DataFrame sorted by adjusted p-value with
    columns: ``cluster``, ``mean_prop_<group>`` (one per group), ``statistic``,
    ``p_value``, ``p_adj``, ``direction``, ``significant``.
    """
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
                # e.g. all proportions identical -> test undefined
                stat, p = np.nan, np.nan

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


# ---------------------------------------------------------------------------
# Path validation (S1 / S2 — path traversal prevention)
# ---------------------------------------------------------------------------
def validate_path(path, base=None):
    """Validate that *path* resolves inside *base* (default: cwd).

    Returns the absolute path if valid, otherwise raises ``ValueError``.
    """
    if base is None:
        base = os.path.abspath(os.getcwd())
    abs_path = os.path.abspath(path)
    if not abs_path.startswith(base + os.sep) and abs_path != base:
        raise ValueError(f"Path '{path}' is outside the allowed directory.")
    return abs_path


# ---------------------------------------------------------------------------
# Standardised plot controls (de-duplicates ~30 repeated slider/selectbox blocks)
# ---------------------------------------------------------------------------
@dataclass
class PlotSettings:
    """Values returned by :func:`plot_controls`. Unrequested fields stay ``None``."""
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
    """Render a consistent block of plot controls and return a :class:`PlotSettings`.

    ``key_prefix`` must be unique per form; every widget key is derived from it.
    ``width``/``height``/``dot_size`` are ``(min, max, default)`` tuples. Pass a
    subset via ``include`` to omit controls a given plot does not need (e.g.
    ``include=("width", "height")`` for a plot with no dots or colormap).
    """
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
