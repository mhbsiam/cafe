import datetime
import io
import os
import random
import re
import sys
import tempfile
import time
import zipfile

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import streamlit as st
from scipy.stats import median_abs_deviation

from theme import (
    TOKENS,
    apply_theme,
    info_card,
    loading_status,
    page_header,
    section_header,
    stepper,
)
from utils import (
    inspect_uploaded_csv_files,
    process_uploaded_csv_files,
    run_combat_correction,
)

sc.settings.n_jobs = -1
apply_theme()

# ---------------------------------------------------------------------------
# Pipeline step labels
# ---------------------------------------------------------------------------
_STEPS = ["Upload", "QC", "PCA", "Batch", "UMAP + Leiden", "Export"]

# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------
def mad_outlier_tails(values, nmads):
    """Robust MAD-based outlier detection on a per-cell metric.

    Returns (lower_mask, upper_mask, med, lower_cut, upper_cut). If the metric
    is degenerate (MAD == 0) nothing is flagged.
    """
    values = np.asarray(values, dtype=float)
    med = np.median(values)
    mad = median_abs_deviation(values)
    if mad == 0 or np.isnan(mad):
        empty = np.zeros(values.shape[0], dtype=bool)
        return empty, empty, med, np.nan, np.nan
    lower_cut = med - nmads * mad
    upper_cut = med + nmads * mad
    return values < lower_cut, values > upper_cut, med, lower_cut, upper_cut


def get_total_signal(adata):
    """Per-cell total marker signal, cached against the exact X matrix."""
    X = adata.X
    cached = st.session_state.get("_qc_total_signal")
    if cached is not None and cached[0] is X:
        return cached[1]
    ts = np.asarray(X.sum(axis=1)).ravel().astype(float)
    st.session_state["_qc_total_signal"] = (X, ts)
    return ts


def qc_preview_tables(
    total_signal, groups, nmads_map, default_nmads, use_lower, use_upper
):
    """Build summary and sensitivity tables for the QC preview."""
    grid = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    summary_rows = []
    sens_rows = []
    for label, idx in groups:
        vals = np.asarray(total_signal[idx], dtype=float)
        n = int(vals.shape[0])
        med = float(np.median(vals)) if n else float("nan")
        mad = float(median_abs_deviation(vals)) if n else float("nan")
        degenerate = (n == 0) or (mad == 0) or np.isnan(mad)

        def pct_removed(nm):
            if degenerate:
                return 0.0
            removed = 0
            if use_lower:
                removed += int(np.count_nonzero(vals < med - nm * mad))
            if use_upper:
                removed += int(np.count_nonzero(vals > med + nm * mad))
            return round(removed / n * 100, 2)

        this_nmads = nmads_map.get(label, default_nmads)
        lo = float("nan") if degenerate else med - this_nmads * mad
        hi = float("nan") if degenerate else med + this_nmads * mad
        summary_rows.append(
            {
                "SampleID": label,
                "n cells": n,
                "median": round(med, 1),
                "MAD": round(mad, 2),
                "n_MADs": this_nmads,
                "lower cut": None if np.isnan(lo) else round(lo, 1),
                "upper cut": None if np.isnan(hi) else round(hi, 1),
                "% removed": pct_removed(this_nmads),
            }
        )
        sens = {"SampleID": label}
        for g in grid:
            sens[f"{g:g}"] = pct_removed(g)
        sens_rows.append(sens)
    summary_df = pd.DataFrame(summary_rows).set_index("SampleID")
    sens_df = pd.DataFrame(sens_rows).set_index("SampleID")
    return summary_df, sens_df


# ---------------------------------------------------------------------------
# Reproducibility report
#
# Parameters chosen at each pipeline step are stashed into st.session_state as
# they are applied (keys prefixed ``_params_``). At export time they are
# assembled into a plain-text methods/parameters file bundled in the ZIP.
# ---------------------------------------------------------------------------
def _safe_filename(name):
    """Filesystem-safe token for a sample id used as a CSV filename."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "sample"


def _methods_paragraph(adata, qc, pca, batch, umap):
    """A draft methods sentence users can paste into a manuscript."""
    parts = [
        f"Spectral flow cytometry data from {adata.obs['SampleID'].nunique()} "
        f"samples ({adata.n_obs:,} cells, {adata.n_vars} markers) were analyzed in CAFE."
    ]
    if qc.get("performed"):
        if qc.get("per_sample"):
            parts.append(
                "Low-quality events were removed per sample using median absolute "
                "deviation (MAD) outlier detection on each cell's total marker signal."
            )
        else:
            parts.append(
                f"Low-quality events were removed using median absolute deviation (MAD) "
                f"outlier detection (±{qc.get('n_mads')} MADs) on each cell's total "
                "marker signal."
            )
    if pca.get("performed"):
        parts.append(
            f"Principal component analysis retained {pca.get('n_components')} components "
            f"({pca.get('variance_threshold')}% variance)."
        )
    if batch.get("method", "None") != "None":
        parts.append(
            f"{batch.get('method')} batch correction was applied across batches while "
            f"preserving the {', '.join(batch.get('covariates', []))} covariate(s)."
        )
    parts.append(
        f"A neighborhood graph (n_neighbors={umap.get('n_neighbors')}, "
        f"metric={umap.get('metric')}) was computed on "
        f"{umap.get('use_rep', 'the expression matrix')}, embedded with UMAP "
        f"(min_dist={umap.get('min_dist')}), and clustered with the Leiden algorithm "
        f"(resolution={umap.get('resolution')}), yielding {umap.get('n_clusters')} "
        "clusters. A fixed random seed (50) was used throughout."
    )
    return " ".join(parts)


def build_parameters_report(adata):
    """Assemble a human-readable summary of every parameter used in this run."""
    s = st.session_state
    qc = s.get("_params_qc", {})
    pca = s.get("_params_pca", {})
    batch = s.get("_params_batch", {})
    umap = s.get("_params_umap", {})

    L = []

    def add(line=""):
        L.append(line)

    add("CAFE — Data Processing: Methods & Parameters")
    add("=" * 52)
    add(f"Generated: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    add("Random seed: 50 (neighbors, UMAP, Leiden)")
    add("")

    add("INPUT DATA")
    add("-" * 52)
    add(f"Samples          : {adata.obs['SampleID'].nunique()}")
    add(f"Groups           : {', '.join(map(str, pd.unique(adata.obs['Group'])))}")
    add(f"Cells (final)    : {adata.n_obs:,}")
    add(f"Markers ({adata.n_vars:>3})     : {', '.join(map(str, adata.var_names))}")
    add("")

    add("QUALITY CONTROL")
    add("-" * 52)
    if not qc.get("performed"):
        add("Not performed — no cells removed.")
    else:
        mode = "per-sample" if qc.get("per_sample") else "global"
        add(f"Method           : MAD outlier removal on per-cell total signal ({mode})")
        add(f"Tails removed    : {qc.get('tails', 'Both')} (debris + doublets)")
        nm = qc.get("n_mads")
        if isinstance(nm, dict):
            add("n_MADs           :")
            for sid, v in nm.items():
                add(f"    {sid}: {v}")
        else:
            add(f"n_MADs           : {nm}")
        add(f"Cells before     : {qc.get('n_before', 0):,}")
        add(f"Cells removed    : {qc.get('n_removed', 0):,} ({qc.get('pct', 0):.2f}%)")
        add(f"Cells retained   : {qc.get('retained', adata.n_obs):,}")
    add("")

    add("PRINCIPAL COMPONENT ANALYSIS")
    add("-" * 52)
    if not pca.get("performed"):
        add("Not performed.")
    else:
        add(f"SVD solver       : {pca.get('solver')}")
        add(f"Variance kept    : {pca.get('variance_threshold')}%")
        add(f"Components used  : {pca.get('n_components')}")
    add("")

    add("BATCH CORRECTION")
    add("-" * 52)
    if batch.get("method", "None") == "None":
        add("Not performed.")
    else:
        add(f"Method           : {batch.get('method')}")
        add(f"Batch key        : {batch.get('batch_key')}")
        add(f"Covariates       : {', '.join(batch.get('covariates', []))}")
        mapping = batch.get("mapping", {})
        if mapping:
            add("Sample -> Batch  :")
            for sid, b in mapping.items():
                add(f"    {sid}: {b}")
    add("")

    add("NEIGHBORS / UMAP / LEIDEN")
    add("-" * 52)
    add(f"Representation   : {umap.get('use_rep', 'raw expression (X)')}")
    add(f"n_neighbors      : {umap.get('n_neighbors')}")
    add(f"Distance metric  : {umap.get('metric')}")
    add(f"UMAP min_dist    : {umap.get('min_dist')}")
    add(f"Leiden flavor    : {umap.get('leiden_flavor')}")
    add(f"Leiden resolution: {umap.get('resolution')}")
    add(f"Clusters found   : {umap.get('n_clusters')}")
    add("")

    add("SOFTWARE VERSIONS")
    add("-" * 52)
    add(f"Python           : {sys.version.split()[0]}")
    add(f"scanpy           : {sc.__version__}")
    add(f"anndata          : {ad.__version__}")
    add(f"numpy            : {np.__version__}")
    add(f"pandas           : {pd.__version__}")
    add("")

    add("METHODS (draft for manuscript)")
    add("-" * 52)
    add(_methods_paragraph(adata, qc, pca, batch, umap))
    add("")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# Session state initialization
# ---------------------------------------------------------------------------
def _init_state():
    st.session_state.setdefault("adata", None)
    st.session_state.setdefault("uploaded_files", None)
    st.session_state.setdefault("batch_correction_done", False)
    st.session_state.setdefault("umap_computed", False)
    st.session_state.setdefault("pca_done", False)
    st.session_state.setdefault("leiden_computed", False)
    st.session_state.setdefault("qc_done", False)
    st.session_state.setdefault("pca_selected", False)


_init_state()

# ---------------------------------------------------------------------------
# Page header
# ---------------------------------------------------------------------------
page_header(
    "Data Processing",
    subtitle="Upload CSVs, clean the data, run dimensionality reduction and clustering, then export a ZIP of results.",
)

# ---------------------------------------------------------------------------
# Stepper state
# ---------------------------------------------------------------------------
def _current_step() -> int:
    if st.session_state.adata is None:
        return 0
    if not st.session_state.qc_done:
        return 1
    if not st.session_state.pca_done:
        return 2
    if not st.session_state.batch_correction_done:
        return 3
    if not st.session_state.umap_computed:
        return 4
    return 5


def _done_steps() -> set[int]:
    done = set()
    if st.session_state.adata is not None:
        done.add(0)
    if st.session_state.qc_done:
        done.add(1)
    if st.session_state.pca_done:
        done.add(2)
    if st.session_state.batch_correction_done:
        done.add(3)
    if st.session_state.umap_computed:
        done.add(4)
    return done


stepper(_STEPS, _current_step(), _done_steps())

# ---------------------------------------------------------------------------
# Step 0: Upload CSVs
# ---------------------------------------------------------------------------
if st.session_state.adata is None:
    section_header(
        "Upload your CSV files",
        subtitle="Files must follow the SampleID_Group.csv naming convention and contain the same markers.",
        step=1,
    )

    info_card(
        title="CSV format checklist",
        body=(
            "- Name each file <code>SampleID_Group.csv</code>, e.g. <code>ABC001_Aged.csv</code><br>"
            "- Include the same markers in every file<br>"
            "- Use identical marker names (case-sensitive)<br>"
            "- Use unique sample names<br>"
            "- Two biological groups are recommended for statistics"
        ),
        kind="info",
    )

    # Phase A — pick files. Once loaded we stash the raw files and a
    # diagnostics summary in session, then rerun into the review phase (B) so
    # the user can confirm samples, groups, cell counts and markers before
    # committing to QC.
    if st.session_state.get("_upload_review") is None:
        with st.form(key="upload_csv_form"):
            uploaded_files = st.file_uploader(
                "Choose CSV files", type="csv", accept_multiple_files=True
            )
            submit_upload = st.form_submit_button("Load CSV files", type="primary")

        if uploaded_files and submit_upload:
            st.session_state["_upload_files"] = uploaded_files
            st.session_state["_upload_review"] = inspect_uploaded_csv_files(uploaded_files)
            st.rerun()

        st.stop()

    # Phase B — review what was uploaded.
    review = st.session_state["_upload_review"]
    file_rows = review["files"]

    sample_df = pd.DataFrame(
        [
            {
                "SampleID": r["SampleID"],
                "Group": r["Group"],
                "Cells": r["n_cells"],
                "File": r["file"],
            }
            for r in file_rows
        ]
    )
    total_cells = int(sample_df["Cells"].sum())

    st.markdown("**Samples detected**")
    st.caption(
        f"{len(sample_df)} file(s), {sample_df['Group'].nunique()} group(s), "
        f"{total_cells:,} cells total."
    )
    st.dataframe(
        sample_df.style.format({"Cells": "{:,}"}),
        hide_index=True,
        width="stretch",
    )

    st.markdown("**Markers available**")
    common = review["common_markers"]
    if common:
        st.caption(f"{len(common)} marker(s) shared across all files.")
        st.write(", ".join(common))
    else:
        st.caption("No marker is shared across every file — see the mismatch report below.")

    mismatches = review["mismatches"]
    bad_names = review["bad_names"]
    can_proceed = not mismatches and not bad_names

    if bad_names:
        info_card(
            title="Filename convention",
            body=(
                "These file(s) don't match <code>SampleID_Group.csv</code> and would be "
                "skipped:<br>" + "<br>".join(f"• <code>{n}</code>" for n in bad_names) +
                "<br><br>Rename them and re-upload."
            ),
            kind="error",
        )

    if mismatches:
        rows = "".join(
            f"• <code>{m['file']}</code> is missing: "
            f"{', '.join('<code>%s</code>' % x for x in m['missing'])}<br>"
            for m in mismatches
        )
        info_card(
            title="Marker names don't match across files",
            body=(
                "Marker names are compared case-sensitively. The following file(s) are "
                "missing markers that appear in other files (a mismatch is often just a "
                "case or spelling difference):<br><br>" + rows +
                "<br>Fix the marker names so every file matches, then re-upload."
            ),
            kind="error",
        )

    col_proceed, col_reset = st.columns([1, 1])
    with col_proceed:
        proceed = st.button(
            "Proceed to QC",
            type="primary",
            disabled=not can_proceed,
            key="dp_upload_proceed",
        )
    with col_reset:
        reset = st.button("Upload different files", key="dp_upload_reset")

    if reset:
        st.session_state.pop("_upload_review", None)
        st.session_state.pop("_upload_files", None)
        st.rerun()

    if proceed and can_proceed:
        # Rewind buffers so process_uploaded_csv_files can re-read them.
        for f in st.session_state["_upload_files"]:
            try:
                f.seek(0)
            except Exception:
                pass
        adata = process_uploaded_csv_files(st.session_state["_upload_files"])
        if adata is not None:
            st.session_state.adata = adata
            st.session_state.pop("_upload_review", None)
            st.session_state.pop("_upload_files", None)
            st.rerun()
        else:
            st.error("Could not build an AnnData object from the uploaded files. Check the CSV format and try again.")

    st.stop()

adata = st.session_state.adata

# ---------------------------------------------------------------------------
# Step 1: Quality Control
# ---------------------------------------------------------------------------
if not st.session_state.qc_done:
    section_header(
        "Quality Control",
        subtitle="Remove doublets, debris and unmixing outliers using MAD outlier detection on each cell's total marker signal.",
        step=2,
    )

    sample_ids = list(pd.unique(adata.obs["SampleID"])) if adata is not None else []

    col_controls, col_preview = st.columns([2, 3])

    with col_controls:
        info_card(
            title="QC settings",
            body="Choose whether to run QC globally, per sample, or skip it. Both the lower and upper tails of the total-signal distribution are removed (debris and doublets).",
            kind="info",
        )

        qc_option = st.radio(
            "Select QC option",
            ("None", "Perform QC", "Perform QC per sample"),
            key="qc_option",
        )
        qc_per_sample = qc_option == "Perform QC per sample"

        if qc_option == "Perform QC":
            nmads = st.slider(
                "Number of MADs (higher = fewer cells removed)",
                min_value=2.0,
                max_value=8.0,
                value=5.0,
                step=0.5,
                key="qc_nmads",
            )
        else:
            nmads = 5.0

        qc_tails = "Both"

        nmads_map = {}
        if qc_per_sample and sample_ids:
            with st.expander("Per-sample MAD thresholds", expanded=True):
                st.caption(
                    "Each sample gets its own MAD threshold, computed within that sample. "
                    "Lower a sample's value to filter it more aggressively."
                )
                for sid in sample_ids:
                    nmads_map[sid] = st.number_input(
                        f"n_MADs for {sid}",
                        min_value=2.0,
                        max_value=8.0,
                        value=5.0,
                        step=0.5,
                        key=f"qc_nmads_{sid}",
                    )

    with col_preview:
        if qc_option != "None":
            total_signal_preview = get_total_signal(adata)
            use_lower_pv = qc_tails in ("Both", "Lower only (debris)")
            use_upper_pv = qc_tails in ("Both", "Upper only (doublets)")
            if qc_per_sample and sample_ids:
                preview_groups = [
                    (sid, np.where(adata.obs["SampleID"].values == sid)[0])
                    for sid in sample_ids
                ]
            else:
                preview_groups = [("All samples", np.arange(adata.n_obs))]

            summary_df, sens_df = qc_preview_tables(
                total_signal_preview,
                preview_groups,
                nmads_map,
                nmads,
                use_lower_pv,
                use_upper_pv,
            )

            st.caption(
                "Preview of what the current thresholds would remove. Nothing is removed until you click Apply QC."
            )
            st.dataframe(
                summary_df.style.background_gradient(
                    subset=["% removed"], cmap="YlGnBu"
                ),
                width='stretch',
            )
            st.caption(
                "Sensitivity — % of cells removed at each n_MADs. A row that climbs steeply as n_MADs drops is threshold-sensitive."
            )
            st.dataframe(
                sens_df.style.background_gradient(cmap="YlGnBu", axis=1),
                width='stretch',
            )

    apply_qc_button = st.button(
        label="Apply QC", type="primary", key="apply_qc_button"
    )

    if apply_qc_button:
        if qc_option == "None":
            st.session_state.qc_done = True
            st.session_state["_params_qc"] = {"performed": False}
            st.info("QC not selected. Proceeding without removing cells.")
        else:
            start_time = time.time()
            X = adata.X
            total_signal = np.asarray(X.sum(axis=1)).ravel().astype(float)
            adata.obs["total_signal"] = total_signal

            use_lower = qc_tails in ("Both", "Lower only (debris)")
            use_upper = qc_tails in ("Both", "Upper only (doublets)")
            remove_mask = np.zeros(adata.n_obs, dtype=bool)

            if qc_per_sample:
                for sample_id in adata.obs["SampleID"].unique():
                    idx = np.where(adata.obs["SampleID"].values == sample_id)[0]
                    sample_nmads = nmads_map.get(sample_id, nmads)
                    lower_mask, upper_mask, _, _, _ = mad_outlier_tails(
                        total_signal[idx], sample_nmads
                    )
                    if use_lower:
                        remove_mask[idx[lower_mask]] = True
                    if use_upper:
                        remove_mask[idx[upper_mask]] = True
                lower_cut = upper_cut = None
            else:
                lower_mask, upper_mask, med, lower_cut, upper_cut = mad_outlier_tails(
                    total_signal, nmads
                )
                if use_lower:
                    remove_mask |= lower_mask
                if use_upper:
                    remove_mask |= upper_mask

            n_before = adata.n_obs
            n_removed = int(remove_mask.sum())
            pct_removed = (n_removed / n_before * 100) if n_before else 0.0

            qc_breakdown = pd.DataFrame(
                {
                    "SampleID": adata.obs["SampleID"].values,
                    "Group": adata.obs["Group"].values,
                    "removed": remove_mask,
                }
            )
            removed_summary = qc_breakdown.groupby(["Group", "SampleID"]).agg(
                cells_before=("removed", "size"),
                cells_removed=("removed", "sum"),
            )
            removed_summary["percent_removed"] = (
                removed_summary["cells_removed"] / removed_summary["cells_before"] * 100
            ).round(2)
            if qc_per_sample:
                removed_summary["n_mads"] = [
                    nmads_map.get(sid, nmads) for (_, sid) in removed_summary.index
                ]
            else:
                removed_summary["n_mads"] = nmads

            fig, ax = plt.subplots()
            ax.hist(total_signal, bins=100, color=TOKENS.primary_light, edgecolor=TOKENS.primary)
            if not qc_per_sample and lower_cut is not None:
                if use_lower:
                    ax.axvline(
                        lower_cut, color=TOKENS.primary, linestyle="--", label="Lower cutoff"
                    )
                if use_upper:
                    ax.axvline(
                        upper_cut, color=TOKENS.primary, linestyle="--", label="Upper cutoff"
                    )
                ax.legend()
            ax.set_xlabel("Total marker signal per cell")
            ax.set_ylabel("Number of cells")
            ax.set_title("Per-cell total signal (MAD outlier cutoffs)")
            st.pyplot(fig, width='stretch')
            plt.close(fig)

            adata = adata[~remove_mask].copy()
            st.session_state.adata = adata
            st.session_state.qc_done = True
            st.session_state.pop("_qc_total_signal", None)

            st.session_state["_params_qc"] = {
                "performed": True,
                "per_sample": qc_per_sample,
                "tails": qc_tails,
                "n_mads": dict(nmads_map) if qc_per_sample else nmads,
                "n_before": n_before,
                "n_removed": n_removed,
                "pct": round(pct_removed, 2),
                "retained": adata.n_obs,
            }

            st.success(
                f"Removed {n_removed:,} of {n_before:,} cells ({pct_removed:.2f}%). "
                f"Retained {adata.n_obs:,} cells."
            )
            st.write("Cells removed per sample:")
            st.dataframe(removed_summary, width='stretch')
            st.caption(f"QC completed in {time.time() - start_time:.2f} seconds")
        st.rerun()

    st.stop()

# ---------------------------------------------------------------------------
# Step 2: PCA (optional)
# ---------------------------------------------------------------------------
if not st.session_state.pca_done:
    section_header(
        "Principal Component Analysis",
        subtitle="Optionally reduce dimensionality before computing neighbors and UMAP.",
        step=3,
    )

    if "pca_selected" not in st.session_state:
        st.session_state.pca_selected = False

    if not st.session_state.pca_selected:
        with st.form(key="pca_selection_form"):
            pca_option = st.radio(
                "Select PCA option", ("Skip PCA", "Perform PCA"), key="pca_option"
            )
            proceed_button = st.form_submit_button(label="Proceed", type="primary")

        if proceed_button:
            if pca_option == "Skip PCA":
                st.session_state.pca_done = True
                st.session_state.pca_selected = False
                st.session_state["_params_pca"] = {"performed": False}
            else:
                st.session_state.pca_selected = True
            st.rerun()
        st.stop()

    with st.form(key="pca_options_form"):
        pca_solver = st.selectbox(
            "Choose the SVD solver for PCA",
            ("Auto", "Full", "Arpack", "Randomized"),
            key="pca_solver",
        )
        variance_threshold = st.slider(
            "Select the explained variance threshold (%) to retain",
            min_value=70,
            max_value=99,
            value=95,
            key="variance_threshold",
        )
        apply_pca_button = st.form_submit_button(label="Apply PCA", type="primary")

    if apply_pca_button:
        start_time = time.time()
        loading_status("Running PCA", steps=["Computing PCA", "Selecting components", "Generating plots"], current_index=0)
        progress_bar = st.progress(0)
        progress = 0

        sc.tl.pca(adata, svd_solver=pca_solver.lower())
        progress += 25
        progress_bar.progress(progress)

        pca_variance_ratio = adata.uns["pca"]["variance_ratio"].cumsum()
        threshold_mask = pca_variance_ratio >= variance_threshold / 100
        if not threshold_mask.any():
            st.warning(
                f"No principal component reaches {variance_threshold}% variance. Using all {len(pca_variance_ratio)} components."
            )
            num_components = len(pca_variance_ratio)
        else:
            num_components = threshold_mask.argmax() + 1

        progress += 25
        progress_bar.progress(progress)

        adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :num_components]
        adata.varm["PCs"] = adata.varm["PCs"][:, :num_components]
        adata.uns["pca"]["variance"] = adata.uns["pca"]["variance"][:num_components]
        adata.uns["pca"]["variance_ratio"] = adata.uns["pca"]["variance_ratio"][:num_components]
        pca_variance = adata.uns["pca"]["variance_ratio"]

        progress += 25
        progress_bar.progress(progress)
        st.session_state.adata = adata
        st.session_state.pca_done = True
        st.session_state["_params_pca"] = {
            "performed": True,
            "solver": pca_solver,
            "variance_threshold": variance_threshold,
            "n_components": int(num_components),
        }

        fig, ax = plt.subplots()
        ax.bar(range(1, num_components + 1), pca_variance, color=TOKENS.primary)
        ax.set_xlabel("Principal Component")
        ax.set_ylabel("Explained Variance Ratio")
        ax.set_title("Explained Variance by PCA Components")
        ax.set_xticks(range(1, num_components + 1))
        st.pyplot(fig, width='stretch')
        plt.close(fig)

        st.write("PCA Results by Group:")
        fig, ax = plt.subplots()
        sc.pl.pca(st.session_state.adata, color="Group", ax=ax, show=False)
        st.pyplot(fig, width='stretch')
        plt.close(fig)

        st.success(f"PCA completed in {time.time() - start_time:.2f} seconds. Retained {num_components} components.")
        st.rerun()

    st.stop()

# ---------------------------------------------------------------------------
# Step 3: Batch correction (optional)
# ---------------------------------------------------------------------------
if not st.session_state.batch_correction_done:
    section_header(
        "Batch Correction",
        subtitle="Correct for technical batch effects. Batches must be distinct from biological groups.",
        step=4,
    )

    sample_ids = (
        sorted(pd.unique(adata.obs["SampleID"])) if adata is not None else []
    )
    sample_group = (
        dict(zip(adata.obs["SampleID"], adata.obs["Group"]))
        if adata is not None
        else {}
    )

    info_card(
        title="Use a technical batch, not a biological group",
        body="Correcting on <code>Group</code> would erase the differences you are testing for. Assign each sample to a technical batch such as acquisition day or instrument.",
        kind="warning",
    )

    with st.form(key="batch_correction_form"):
        batch_correction_method = st.radio(
            "Select batch correction method",
            ("None", "ComBat"),
            key="batch_correction_method",
        )

        st.caption("Batch assignment (used only when ComBat is selected):")
        batch_assignment_df = st.data_editor(
            pd.DataFrame(
                {
                    "SampleID": sample_ids,
                    "Group": [sample_group.get(sid, "") for sid in sample_ids],
                    "Batch": ["1"] * len(sample_ids),
                }
            ),
            disabled=["SampleID", "Group"],
            hide_index=True,
            width='stretch',
            key="batch_assignment_editor",
        )

        submit_button = st.form_submit_button(label="Apply", type="primary")

    if submit_button:
        if "Group" not in adata.obs.columns:
            st.error("The 'Group' column does not exist in adata.obs.")
            st.session_state.batch_correction_done = True
            st.session_state["_params_batch"] = {"method": "None"}
            st.rerun()
        elif batch_correction_method == "None":
            st.session_state.batch_correction_done = True
            st.session_state["_params_batch"] = {"method": "None"}
            st.info("No batch correction applied.")
            st.rerun()
        else:
            batch_map = dict(
                zip(
                    batch_assignment_df["SampleID"].astype(str),
                    batch_assignment_df["Batch"].astype(str),
                )
            )
            adata.obs["Batch"] = adata.obs["SampleID"].astype(str).map(batch_map)

            start_time = time.time()
            ok, msg = run_combat_correction(adata, batch_key="Batch", covariates=("Group",))
            elapsed_time = time.time() - start_time

            if not ok:
                st.error(msg)
                adata.obs.drop(columns=["Batch"], inplace=True, errors="ignore")
            else:
                st.success(f"{msg} (in {elapsed_time:.2f} seconds)")
                st.subheader("PCA After Batch Correction")
                fig, ax = plt.subplots()
                sc.pl.pca(adata, color="Group", ax=ax, show=False)
                st.pyplot(fig, width='stretch')
                plt.close(fig)

                st.session_state.adata = adata
                st.session_state.batch_correction_done = True
                st.session_state["_params_batch"] = {
                    "method": "ComBat",
                    "batch_key": "Batch",
                    "covariates": ["Group"],
                    "mapping": batch_map,
                }
                st.rerun()

    st.stop()

# ---------------------------------------------------------------------------
# Step 4: UMAP + Leiden clustering
# ---------------------------------------------------------------------------
if not st.session_state.umap_computed:
    section_header(
        "UMAP + Leiden Clustering",
        subtitle="Compute the nearest-neighbor graph, UMAP embedding, and Leiden clusters.",
        step=5,
    )

    with st.form(key="umap_leiden_form"):
        col1, col2 = st.columns(2)
        with col1:
            resolution = st.slider(
                "Leiden Resolution",
                min_value=0.01,
                max_value=2.5,
                value=0.5,
                step=0.01,
                help="Lower values yield fewer, larger clusters; higher values yield more, smaller clusters.",
            )
            n_neighbors = st.slider(
                "n_neighbors",
                min_value=5,
                max_value=50,
                value=15,
                step=1,
                help="Controls the local neighborhood size.",
            )
        with col2:
            min_dist = st.slider(
                "UMAP min_dist",
                min_value=0.0,
                max_value=1.0,
                value=0.1,
                step=0.01,
                help="Controls how tightly UMAP packs points together.",
            )
            metric_options = [
                "euclidean",
                "manhattan",
                "cosine",
                "correlation",
                "minkowski",
            ]
            selected_metric = st.selectbox(
                "Distance metric",
                options=metric_options,
                index=0,
                help="Determines how distances between data points are calculated.",
            )

        flavor_options = ["igraph", "leidenalg"]
        selected_flavor = st.selectbox(
            "Leiden flavor",
            options=flavor_options,
            index=0,
        )

        dim_option = ["None", "X_pca", "X_umap"]
        pca_default_index = 1 if st.session_state.get("pca_selected", False) else 0
        selected_dim = st.selectbox(
            "Dimension reduction metric to use",
            options=dim_option,
            index=pca_default_index,
            help="Select None if you have not reduced the data using PCA.",
        )

        submit_button = st.form_submit_button(
            label="Apply and Compute", type="primary"
        )

    if submit_button:
        adata = st.session_state.get("adata")
        if not st.session_state.get("leiden_computed"):
            progress_bar = st.progress(0)
            elapsed_time_placeholder = st.empty()
            status_placeholder = st.empty()

            progress = 0
            start_time = time.time()

            def update_elapsed_time():
                elapsed = time.time() - start_time
                elapsed_time_placeholder.markdown(
                    f"**Elapsed Time:** {elapsed:.2f} seconds"
                )

            with st.spinner("Working…"):
                update_elapsed_time()
                status_placeholder.markdown("**Step 1/4: Computing neighbors…**")
                progress_bar.progress(progress)

                use_rep = None if selected_dim == "None" else selected_dim
                if use_rep is not None and use_rep not in adata.obsm:
                    st.error(
                        f"{selected_dim} not found in adata.obsm. Please run the corresponding dimensionality reduction first."
                    )
                    st.stop()
                sc.pp.neighbors(
                    adata,
                    n_neighbors=n_neighbors,
                    method="umap",
                    use_rep=use_rep,
                    metric=selected_metric,
                    random_state=50,
                )
                progress += 25
                progress_bar.progress(progress)

                update_elapsed_time()
                status_placeholder.markdown(
                    f"**Step 2/4: Running UMAP with n_neighbors={n_neighbors}, min_dist={min_dist}, metric={selected_metric}…**"
                )
                sc.tl.umap(adata, min_dist=min_dist, random_state=50)
                progress += 25
                progress_bar.progress(progress)

                update_elapsed_time()
                status_placeholder.markdown("**Step 3/4: Computing Leiden clustering…**")
                sc.tl.leiden(
                    adata,
                    resolution=resolution,
                    random_state=50,
                    flavor=selected_flavor,
                    n_iterations=2,
                    directed=False,
                )
                progress += 25
                progress_bar.progress(progress)

                update_elapsed_time()
                status_placeholder.markdown("**Step 4/4: Saving outputs…**")

                st.session_state["_params_umap"] = {
                    "resolution": resolution,
                    "n_neighbors": n_neighbors,
                    "min_dist": min_dist,
                    "metric": selected_metric,
                    "leiden_flavor": selected_flavor,
                    "use_rep": use_rep if use_rep is not None else "raw expression (X)",
                    "n_clusters": int(adata.obs["leiden"].nunique()),
                    "random_state": 50,
                }

                random_number = random.randint(1000, 9999)
                zip_file_name = f"analysis_outputs_{random_number}.zip"

                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    with tempfile.NamedTemporaryFile(
                        delete=False, suffix=".h5ad"
                    ) as temp_file:
                        temp_path = temp_file.name
                        adata.write(temp_path)

                    zip_file.write(
                        temp_path,
                        arcname=f"adata_{resolution}_{random_number}.h5ad",
                    )
                    os.remove(temp_path)

                    cluster_sample_counts = (
                        adata.obs.groupby(["leiden", "SampleID"])
                        .size()
                        .unstack(fill_value=0)
                    )
                    cluster_sample_counts_csv = io.BytesIO()
                    cluster_sample_counts.to_csv(cluster_sample_counts_csv)
                    cluster_sample_counts_csv.seek(0)
                    zip_file.writestr(
                        f"cluster_sample_counts_{random_number}.csv",
                        cluster_sample_counts_csv.getvalue(),
                    )

                    sample_totals = cluster_sample_counts.sum()
                    cluster_sample_percentages = (
                        cluster_sample_counts.div(sample_totals) * 100
                    )
                    cluster_sample_frequencies_csv = io.BytesIO()
                    cluster_sample_percentages.to_csv(
                        cluster_sample_frequencies_csv
                    )
                    cluster_sample_frequencies_csv.seek(0)
                    zip_file.writestr(
                        f"cluster_sample_frequencies_{random_number}.csv",
                        cluster_sample_frequencies_csv.getvalue(),
                    )

                    median_df = adata.to_df()
                    median_df["leiden"] = adata.obs["leiden"]
                    median_df["SampleID"] = adata.obs["SampleID"]
                    medians = median_df.groupby(["leiden", "SampleID"]).median()
                    median_csv = io.BytesIO()
                    medians.to_csv(median_csv)
                    median_csv.seek(0)
                    zip_file.writestr(
                        f"median_expression_{random_number}.csv",
                        median_csv.getvalue(),
                    )

                    # Reproducibility: plain-text methods & parameters report.
                    zip_file.writestr(
                        f"analysis_parameters_{random_number}.txt",
                        build_parameters_report(adata),
                    )

                    # Round-trip: one CSV per sample carrying the cluster label
                    # (and UMAP coordinates) as extra columns, so users can
                    # re-import into FlowJo and overlay clusters on their gates.
                    expr_df = adata.to_df()
                    umap_coords = adata.obsm.get("X_umap")
                    leiden_num = pd.to_numeric(adata.obs["leiden"], errors="coerce").to_numpy()
                    sample_ids_arr = adata.obs["SampleID"].to_numpy()
                    group_arr = adata.obs["Group"].to_numpy()
                    for sid in pd.unique(sample_ids_arr):
                        mask = sample_ids_arr == sid
                        sub = expr_df.loc[mask].copy()
                        sub["Leiden_Cluster"] = leiden_num[mask]
                        sub["Group"] = group_arr[mask]
                        if umap_coords is not None:
                            sub["UMAP_1"] = umap_coords[mask, 0]
                            sub["UMAP_2"] = umap_coords[mask, 1]
                        sample_csv = io.BytesIO()
                        sub.to_csv(sample_csv, index=False)
                        sample_csv.seek(0)
                        zip_file.writestr(
                            f"per_sample_clustered/{_safe_filename(sid)}.csv",
                            sample_csv.getvalue(),
                        )

                zip_buffer.seek(0)

                st.session_state["_last_zip_buffer"] = zip_buffer.getvalue()
                st.session_state["_last_zip_name"] = zip_file_name
                st.session_state["leiden_computed"] = True
                st.session_state["umap_computed"] = True

                progress += 25
                progress_bar.progress(progress)
                update_elapsed_time()
                status_placeholder.markdown("**Computation completed!**")
                st.success("Clustering complete. Download your results below.")
                st.rerun()

    st.stop()

# ---------------------------------------------------------------------------
# Step 5: Export
# ---------------------------------------------------------------------------
adata = st.session_state.get("adata")
if adata is None:
    st.error("No processed data found. Please complete the pipeline first.")
    st.stop()
else:
    section_header(
        "Export Results",
        subtitle="Download the processed AnnData object and per-cluster tables.",
        step=6,
    )

    info_card(
        title="Pipeline summary",
        body=(
            f"- <strong>{adata.n_obs:,}</strong> cells after QC<br>"
            f"- <strong>{adata.obs['SampleID'].nunique()}</strong> samples<br>"
            f"- <strong>{adata.obs['leiden'].nunique()}</strong> Leiden clusters<br>"
            f"- PCA: {'yes' if st.session_state.get('pca_selected') else 'no'}<br>"
            f"- Batch correction: {'yes' if 'Batch' in adata.obs.columns else 'no'}"
        ),
        kind="success",
    )

    zip_buffer = st.session_state.get("_last_zip_buffer")
    zip_file_name = st.session_state.get("_last_zip_name", "analysis_outputs.zip")
    if zip_buffer is not None:
        st.download_button(
            label="Download All Outputs as ZIP",
            data=zip_buffer,
            file_name=zip_file_name,
            mime="application/zip",
            type="primary",
        )
    else:
        st.warning("No ZIP output found in session. Return to UMAP + Leiden and recompute.")

    info_card(
        title="What's in the ZIP",
        body=(
            "- <code>adata_*.h5ad</code> — processed AnnData for the Visualization page<br>"
            "- <code>analysis_parameters_*.txt</code> — every parameter used, plus a draft methods paragraph<br>"
            "- <code>per_sample_clustered/</code> — one CSV per sample with a <code>Leiden_Cluster</code> "
            "column (and UMAP coordinates) to re-import into FlowJo and overlay clusters on your gates<br>"
            "- cluster count, frequency, and median-expression tables"
        ),
        kind="info",
    )

    info_card(
        title="Next step",
        body="Upload the downloaded <code>.h5ad</code> file to the Visualization page to explore clusters and produce publication figures.",
        kind="info",
    )

    if st.button("Start over", type="secondary", key="dp_start_over"):
        for key in [
            "adata",
            "uploaded_files",
            "batch_correction_done",
            "umap_computed",
            "pca_done",
            "leiden_computed",
            "qc_done",
            "pca_selected",
            "_qc_total_signal",
            "_last_zip_buffer",
            "_last_zip_name",
            "_upload_review",
            "_upload_files",
        ]:
            st.session_state.pop(key, None)
        st.rerun()
