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
import matplotlib.colors  # noqa: F401  (ensures pandas Styler.background_gradient can resolve matplotlib.colors)
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
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
    compute_batch_emd,
    inspect_uploaded_csv_files,
    marker_robust_scales,
    process_uploaded_csv_files,
    recommend_nmads,
    run_combat_correction,
    summarize_batch_effect,
)

sc.settings.n_jobs = -1
apply_theme()

_STEPS = ["Upload", "QC", "PCA", "Batch", "UMAP + Leiden"]


# Core helpers
def mad_outlier_tails(values, nmads):
    """MAD outlier tails: (lower_mask, upper_mask, med, lower_cut, upper_cut); nothing flagged if MAD == 0."""
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
    total_signal, groups, nmads_map, default_nmads, use_lower, use_upper,
    cap_pct=5.0,
):
    """Summary and sensitivity tables for the QC preview (cap_pct = target % removed for the recommendation)."""
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
        rec = recommend_nmads(
            vals, cap_pct=cap_pct, use_lower=use_lower, use_upper=use_upper
        )
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
                "recommended": rec,
                "% at rec": None if rec is None else pct_removed(rec),
            }
        )
        sens = {"SampleID": label}
        for g in grid:
            sens[f"{g:g}"] = pct_removed(g)
        sens_rows.append(sens)
    summary_df = pd.DataFrame(summary_rows).set_index("SampleID")
    sens_df = pd.DataFrame(sens_rows).set_index("SampleID")
    return summary_df, sens_df


def _seed_per_sample_nmads(recs, cap):
    """Seed each per-sample n_MADs from its recommendation; reseed when the cap changes, else keep overrides."""
    prev = st.session_state.get("_qc_reco_cap_prev")
    cap_changed = prev is not None and prev != cap
    for sid, rec in recs.items():
        key = f"qc_nmads_{sid}"
        if key not in st.session_state or cap_changed:
            st.session_state[key] = 5.0 if pd.isna(rec) else float(rec)
    st.session_state["_qc_reco_cap_prev"] = cap


# Review gate: each step stashes result items under a per-step key and reruns;
# the gate renders them with one Continue button so results don't flash past.
def _snapshot_fig(fig, tight=True):
    """PNG bytes of a figure, then close it; tight=False keeps the exact canvas (for matched pixel sizes)."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight" if tight else None)
    plt.close(fig)
    return buf.getvalue()


def _render_review_items(items):
    """Render a stored list of ``(kind, payload)`` result items."""
    renderers = {
        "success": lambda p: st.success(p),
        "info": lambda p: st.info(p),
        "warning": lambda p: st.warning(p),
        "write": lambda p: st.write(p),
        "caption": lambda p: st.caption(p),
        "subheader": lambda p: st.subheader(p),
        "image": lambda p: st.image(p, width="stretch"),
        "image_row": lambda p: [
            col.image(png, caption=cap, width="stretch")
            for col, (png, cap) in zip(st.columns(len(p)), p)
        ],
        "dataframe": lambda p: st.dataframe(p, width="stretch"),
    }
    for kind, payload in items:
        renderers[kind](payload)


def _review_gate(pending_key, done_key, next_label):
    """Render a step's stashed results with one Continue button that sets done_key and advances; else no-op."""
    items = st.session_state.get(pending_key)
    if items is None:
        return
    _render_review_items(items)
    st.divider()
    if st.button(
        f"Continue to {next_label} →",
        type="primary",
        key=f"continue_{done_key}",
    ):
        st.session_state[done_key] = True
        st.session_state.pop(pending_key, None)
        st.rerun()
    st.stop()


# Reproducibility report: per-step params (keys prefixed _params_) assembled
# into a plain-text methods file bundled in the ZIP.
def _safe_filename(name):
    """Filesystem-safe token for a sample id used as a CSV filename."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_") or "sample"


# Canonical CAFE reference — cite this in any publication using CAFE outputs.
_CAFE_CITATION = (
    "Md Hasanul Banna Siam, Md Akkas Ali, Donald Vardaman, Satwik Acharyya, "
    "Mallikarjun Patil, Daniel J Tyrrell, CAFE: An Integrated Web App for "
    "High-Dimensional Analysis and Visualization in Spectral Flow Cytometry, "
    "Bioinformatics, 2025, btaf176, https://doi.org/10.1093/bioinformatics/btaf176"
)


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
    if batch.get("method", "None") != "None":
        parts.append(
            f"{batch.get('method')} batch correction was applied across batches while "
            f"preserving the {', '.join(batch.get('covariates', []))} covariate(s)."
        )
    if pca.get("performed"):
        parts.append(
            f"Principal component analysis retained {pca.get('n_components')} components "
            f"({pca.get('variance_threshold')}% variance)."
        )
    parts.append(
        f"A neighborhood graph (n_neighbors={umap.get('n_neighbors')}, "
        f"metric={umap.get('metric')}) was computed on "
        f"{umap.get('use_rep', 'the expression matrix')}, embedded with UMAP "
        f"(min_dist={umap.get('min_dist')}), and clustered with the Leiden algorithm "
        f"(resolution={umap.get('resolution')}), yielding {umap.get('n_clusters')} "
        "clusters. A fixed random seed (50) was used throughout."
    )
    parts.append(f"When using CAFE, please cite: {_CAFE_CITATION}")
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

    add("BATCH CORRECTION")
    add("-" * 52)
    if batch.get("method", "None") == "None":
        add("Not performed.")
    else:
        add(f"Method           : {batch.get('method')}")
        add(f"Batch key        : {batch.get('batch_key')}")
        add(f"Covariates       : {', '.join(batch.get('covariates', []))}")
        if batch.get("emd_before") is not None and batch.get("emd_after") is not None:
            add(
                f"Batch effect     : mean EMD {batch['emd_before']:.3f} -> "
                f"{batch['emd_after']:.3f} (per-marker worst-batch EMD, after ComBat)"
            )
        mapping = batch.get("mapping", {})
        if mapping:
            add("Sample -> Batch  :")
            for sid, b in mapping.items():
                add(f"    {sid}: {b}")
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

    add("CITATION")
    add("-" * 52)
    add("If you use CAFE in your research, please cite:")
    add(_CAFE_CITATION)
    add("")
    return "\n".join(L)


def build_outputs_zip(adata, resolution, umap_pngs=None):
    """Package all analysis outputs into a ZIP and return (bytes, filename); umap_pngs is {name: png_bytes}."""
    random_number = random.randint(1000, 9999)
    zip_file_name = f"analysis_outputs_{random_number}.zip"

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as zip_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_file:
            temp_path = temp_file.name
            adata.write(temp_path)
        zip_file.write(
            temp_path, arcname=f"adata_{resolution}_{random_number}.h5ad"
        )
        os.remove(temp_path)

        cluster_sample_counts = (
            adata.obs.groupby(["leiden", "SampleID"], observed=False).size().unstack(fill_value=0)
        )
        cluster_sample_counts_csv = io.BytesIO()
        cluster_sample_counts.to_csv(cluster_sample_counts_csv)
        cluster_sample_counts_csv.seek(0)
        zip_file.writestr(
            f"cluster_sample_counts_{random_number}.csv",
            cluster_sample_counts_csv.getvalue(),
        )

        sample_totals = cluster_sample_counts.sum()
        cluster_sample_percentages = cluster_sample_counts.div(sample_totals) * 100
        cluster_sample_frequencies_csv = io.BytesIO()
        cluster_sample_percentages.to_csv(cluster_sample_frequencies_csv)
        cluster_sample_frequencies_csv.seek(0)
        zip_file.writestr(
            f"cluster_sample_frequencies_{random_number}.csv",
            cluster_sample_frequencies_csv.getvalue(),
        )

        median_df = adata.to_df()
        median_df["leiden"] = adata.obs["leiden"]
        median_df["SampleID"] = adata.obs["SampleID"]
        medians = median_df.groupby(["leiden", "SampleID"], observed=False).median()
        median_csv = io.BytesIO()
        medians.to_csv(median_csv)
        median_csv.seek(0)
        zip_file.writestr(
            f"median_expression_{random_number}.csv", median_csv.getvalue()
        )

        # Reproducibility: plain-text methods & parameters report.
        zip_file.writestr(
            f"analysis_parameters_{random_number}.txt",
            build_parameters_report(adata),
        )

        # One CSV per sample with cluster label + UMAP coords, for FlowJo re-import.
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

        # The UMAP figures shown on the results page, as they were rendered.
        for name, png_bytes in (umap_pngs or {}).items():
            if png_bytes is not None:
                zip_file.writestr(f"umap_plots/{name}_{random_number}.png", png_bytes)

    zip_buffer.seek(0)
    return zip_buffer.getvalue(), zip_file_name


# Batch effect diagnostics (read-only), shown before the correction choice.
def _sample_median_zscore(adata):
    """Per-sample median expression, z-scored per marker across samples."""
    df = adata.to_df()
    df["SampleID"] = adata.obs["SampleID"].astype(str).to_numpy()
    med = df.groupby("SampleID").median(numeric_only=True)
    z = (med - med.mean(axis=0)) / med.std(axis=0)
    return z.fillna(0.0)


def render_batch_diagnostics(adata, batch_map):
    """Render batch-effect diagnostics for batch_map (writes adata.obs['Batch'], not X); returns per-marker EMD."""
    adata.obs["Batch"] = adata.obs["SampleID"].astype(str).map(batch_map)

    if adata.obs["Batch"].isnull().any() or (adata.obs["Batch"].astype(str) == "").any():
        info_card(
            title="Assign every sample to a batch",
            body="Some samples have no batch assigned. Fill in the Batch column above before checking.",
            kind="warning",
        )
        return pd.DataFrame()
    if adata.obs["Batch"].nunique() < 2:
        info_card(
            title="Only one batch",
            body="All samples are in a single batch, so there is nothing to compare. Assign at least two batches to check for batch effects.",
            kind="info",
        )
        return pd.DataFrame()

    emd_df = compute_batch_emd(adata, batch_key="Batch")
    verdict, kind, _, _ = summarize_batch_effect(emd_df, adata, batch_key="Batch")
    info_card(title="Batch effect check", body=verdict, kind=kind)

    st.markdown("**Per-marker batch effect (Earth Mover's Distance)**")
    st.caption(
        "EMD between each batch and the pooled data, per marker (robustly scaled so "
        "markers are comparable). Higher = more batch separation. `emd_max` is the "
        "worst batch for that marker."
    )
    with st.expander("How the batch effect is measured"):
        st.markdown(
            "**Why per marker?** Batch effects in cytometry are channel-specific: "
            "each marker can drift on its own (staining efficiency, detector gain, "
            "spillover, day-to-day acquisition). Correction methods (CytoNorm, ComBat) "
            "also work marker by marker, so it helps to see *which* markers are "
            "affected, not just whether some overall effect exists.\n\n"
            "**What is EMD?** Earth Mover's Distance (a.k.a. Wasserstein distance) "
            "measures how different two distributions are. Picture a marker's values "
            "as a pile of sand: EMD is the least amount of sand you would have to "
            "shovel, and how far, to reshape one pile into the other. For each "
            "marker we compare one batch's values against all batches pooled together. "
            "If a batch is shifted or stretched relative to the rest (a technical "
            "effect), the piles don't line up and EMD is large; if the batch looks "
            "like everyone else, almost no sand moves and EMD stays near 0.\n\n"
            "Values are divided by the pooled spread (median absolute deviation) so "
            "markers with naturally large ranges don't dominate. Read the number as "
            "roughly *how many robust standard deviations a batch is shifted*. As a "
            "rough guide, EMD above ~0.5 flags a marker worth a closer look."
        )
    st.dataframe(
        emd_df.style.format({"emd_max": "{:.3f}", "emd_mean": "{:.3f}"})
        .background_gradient(subset=["emd_max"], cmap="YlOrRd"),
        hide_index=True,
        width="stretch",
    )

    ranked = emd_df.sort_values("emd_max")
    fig, ax = plt.subplots(figsize=(6, max(2.0, 0.3 * len(ranked))))
    ax.barh(ranked["marker"], ranked["emd_max"], color=TOKENS.primary)
    ax.set_xlabel("EMD (worst batch vs pooled)")
    ax.set_title("Per-marker batch effect")
    st.pyplot(fig, width="stretch")
    plt.close(fig)

    # Sample-level view: do samples group by Batch (technical) or Group (biology)?
    n_samples = adata.obs["SampleID"].nunique()
    if n_samples >= 2 and adata.n_vars >= 2:
        st.markdown("**Sample similarity (median expression)**")
        st.caption(
            "Each row is a sample. The colour strips mark Batch and Group. If samples "
            "cluster by Batch, that points to a technical effect; if they cluster by "
            "Group, that is your biology, so do not correct it away."
        )
        z = _sample_median_zscore(adata)
        sample_batch = {sid: batch_map.get(str(sid), "?") for sid in z.index}
        sample_group = dict(
            zip(adata.obs["SampleID"].astype(str), adata.obs["Group"].astype(str))
        )
        batch_levels = sorted(set(sample_batch.values()))
        group_levels = sorted(set(sample_group.get(str(s), "?") for s in z.index))
        batch_pal = dict(zip(batch_levels, sns.color_palette("tab10", len(batch_levels))))
        group_pal = dict(zip(group_levels, sns.color_palette("Set2", len(group_levels))))
        row_colors = pd.DataFrame(
            {
                "Batch": [batch_pal[sample_batch[s]] for s in z.index],
                "Group": [group_pal[sample_group.get(str(s), "?")] for s in z.index],
            },
            index=z.index,
        )
        try:
            g = sns.clustermap(
                z,
                row_colors=row_colors,
                cmap="RdBu_r",
                center=0,
                col_cluster=adata.n_vars >= 3,
                figsize=(min(1.0 + 0.4 * adata.n_vars, 12), min(1.5 + 0.35 * n_samples, 12)),
                cbar_kws={"label": "z-score"},
            )
            g.ax_heatmap.set_xlabel("Marker")
            g.ax_heatmap.set_ylabel("Sample")
            st.pyplot(g.fig, width="stretch")
            plt.close(g.fig)
        except Exception as exc:  # clustering can fail on degenerate input
            st.caption(f"Could not draw the sample heatmap ({exc}).")

    if "X_pca" in adata.obsm:
        st.markdown("**PCA coloured by batch**")
        fig, ax = plt.subplots()
        sc.pl.pca(adata, color="Batch", ax=ax, show=False)
        st.pyplot(fig, width="stretch")
        plt.close(fig)
    else:
        st.caption("Run PCA (previous step) to also see a PCA scatter coloured by batch.")

    return emd_df


# Session state initialization
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

page_header(
    "Data Processing",
    subtitle="Upload CSVs, clean the data, run dimensionality reduction and clustering, then export a ZIP of results.",
)


# Stepper state
def _current_step() -> int:
    if st.session_state.adata is None:
        return 0
    if not st.session_state.qc_done:
        return 1
    if not st.session_state.pca_done:
        return 2
    if not st.session_state.batch_correction_done:
        return 3
    return 4


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
    if st.session_state.leiden_computed:
        done.add(4)
    return done


stepper(_STEPS, _current_step(), _done_steps())

# Step 0: Upload CSVs
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

    # Phase A: pick files, then rerun into the review phase (B).
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

# Step 1: Quality Control
if not st.session_state.qc_done:
    section_header(
        "Quality Control",
        subtitle="Remove doublets, debris and unmixing outliers using MAD outlier detection on each cell's total marker signal.",
        step=2,
    )

    _review_gate("_qc_pending", "qc_done", "Batch correction")

    sample_ids = list(pd.unique(adata.obs["SampleID"])) if adata is not None else []

    info_card(
        title="QC settings",
        body="Choose whether to run QC globally, per sample, or skip it. Both the lower and upper tails of the total-signal distribution are removed (debris and doublets).",
        kind="info",
    )

    with st.expander("How MAD-based filtering works, and how to tweak it"):
        st.markdown(
            "**What it looks at.** For every cell we add up all its marker intensities "
            "into one number, the *total marker signal*. Dead cells and debris tend to "
            "have an abnormally low total; doublets and aggregates tend to have an "
            "abnormally high one. Both are technical junk we want to drop.\n\n"
            "**How the cutoffs are set.** We take the median of that total-signal "
            "distribution and its MAD (median absolute deviation). MAD is a robust "
            "measure of spread: because it is built from medians, a handful of extreme "
            "cells cannot inflate it the way they would inflate a standard deviation. A "
            "cell is flagged when its total signal falls below "
            "`median - n x MAD` (debris) or above `median + n x MAD` (doublets). Both "
            "tails are removed.\n\n"
            "**How to tweak it.** The *Number of MADs* (`n`) sets how far from the "
            "median the cutoffs sit. A higher `n` widens the window and removes fewer "
            "cells (more lenient); a lower `n` tightens it and removes more (stricter). "
            "The default of 5 is a common moderate choice. Choose *Perform QC per "
            "sample* to compute the median and MAD within each sample, which is useful "
            "when samples differ in overall brightness or were acquired on different "
            "days.\n\n"
            "**Before you commit.** Nothing is removed until you click *Apply QC*. The "
            "summary table below shows each cutoff and the resulting % removed, and the "
            "sensitivity table shows how that % changes across a range of `n` values, so "
            "a sample whose removal climbs steeply as `n` drops is threshold-sensitive."
        )

    with st.expander("Recommended n_MADs — quick reference"):
        st.markdown(
            "There is no universally correct MAD value — it trades off removing "
            "technical junk against discarding real cells. Use this as a starting "
            "point:"
        )
        st.table(
            pd.DataFrame(
                [
                    {
                        "n_MADs": "2–3",
                        "Stringency": "Strict",
                        "When to use": "Clean, well-controlled samples; aggressively trims tails (risk: removes real cells).",
                    },
                    {
                        "n_MADs": "4–5",
                        "Stringency": "Moderate (default 5)",
                        "When to use": "Robust general-purpose choice for most spectral flow data.",
                    },
                    {
                        "n_MADs": "6–8",
                        "Stringency": "Lenient",
                        "When to use": "Noisy samples or rare populations you want to preserve.",
                    },
                ]
            ).set_index("n_MADs")
        )
        st.caption(
            "In per-sample mode each sample's n_MADs is auto-set from your *target "
            "% removed per tail*: the tightest n that keeps each tail (debris and "
            "doublets) at or under that target. Open *Fine-tune per sample* to "
            "override individual samples. It is a transparent heuristic, not a "
            "biological ground truth."
        )

    qc_option = st.radio(
        "Select QC option",
        ("None", "Perform QC", "Perform QC per sample"),
        key="qc_option",
    )
    qc_per_sample = qc_option == "Perform QC per sample"

    qc_tails = "Both"
    nmads = 5.0
    nmads_map = {}
    reco_cap = 5.0

    if qc_option != "None":
        total_signal_preview = get_total_signal(adata)
        use_lower_pv = qc_tails in ("Both", "Lower only (debris)")
        use_upper_pv = qc_tails in ("Both", "Upper only (doublets)")

        if qc_per_sample and sample_ids:
            preview_groups = [
                (sid, np.where(adata.obs["SampleID"].values == sid)[0])
                for sid in sample_ids
            ]
            # Target cap auto-fills every sample's n_MADs; fine-tuning overrides it.
            reco_cap = st.number_input(
                "Target % removed per tail",
                min_value=1.0,
                max_value=20.0,
                value=5.0,
                step=0.5,
                key="qc_reco_cap",
                help="Each sample's n_MADs is auto-set to the tightest value that "
                "keeps each tail (debris and doublets) at or under this %. Changing "
                "it re-fills the per-sample values below.",
            )
            recs = {
                sid: recommend_nmads(
                    total_signal_preview[idx],
                    cap_pct=reco_cap,
                    use_lower=use_lower_pv,
                    use_upper=use_upper_pv,
                )
                for sid, idx in preview_groups
            }
            _seed_per_sample_nmads(recs, reco_cap)
            with st.expander("Fine-tune per sample (optional)", expanded=False):
                st.caption(
                    "Each sample starts from the recommended value for the target "
                    "above. Override any you like — changing the target resets these."
                )
                for sid in sample_ids:
                    nmads_map[sid] = st.number_input(
                        f"n_MADs for {sid}",
                        min_value=2.0,
                        max_value=8.0,
                        step=0.5,
                        key=f"qc_nmads_{sid}",
                    )
        else:
            nmads = st.slider(
                "Number of MADs (higher = fewer cells removed)",
                min_value=2.0,
                max_value=8.0,
                value=5.0,
                step=0.5,
                key="qc_nmads",
            )
            preview_groups = [("All samples", np.arange(adata.n_obs))]

        summary_df, sens_df = qc_preview_tables(
            total_signal_preview,
            preview_groups,
            nmads_map,
            nmads,
            use_lower_pv,
            use_upper_pv,
            cap_pct=reco_cap,
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
            hist_png = _snapshot_fig(fig)

            adata = adata[~remove_mask].copy()
            st.session_state.adata = adata
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

            # Stash results for the review gate.
            st.session_state["_qc_pending"] = [
                ("image", hist_png),
                (
                    "success",
                    f"Removed {n_removed:,} of {n_before:,} cells ({pct_removed:.2f}%). "
                    f"Retained {adata.n_obs:,} cells.",
                ),
                ("write", "Cells removed per sample:"),
                ("dataframe", removed_summary),
                ("caption", f"QC completed in {time.time() - start_time:.2f} seconds"),
            ]
        st.rerun()

    st.stop()

# Step 3: PCA (optional). Runs after batch correction so X_pca inherits it;
# guarded to stay dormant until the batch step (below in the file) completes.
if st.session_state.batch_correction_done and not st.session_state.pca_done:
    section_header(
        "Principal Component Analysis",
        subtitle="Optionally reduce dimensionality before computing neighbors and UMAP.",
        step=4,
    )

    _review_gate("_pca_pending", "pca_done", "UMAP + Leiden")

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
                # No PCA means nothing to compare against; drop any pre-ComBat plots.
                st.session_state.pop("_pre_combat_pca_batch", None)
                st.session_state.pop("_pre_combat_pca_group", None)
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
        threshold_warning = None
        if not threshold_mask.any():
            threshold_warning = (
                f"No principal component reaches {variance_threshold}% variance. "
                f"Using all {len(pca_variance_ratio)} components."
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
        variance_png = _snapshot_fig(fig)

        # Cumulative variance over all components (pca_variance_ratio was built before X_pca was sliced).
        cum_pct = np.asarray(pca_variance_ratio) * 100.0
        n_total = len(cum_pct)
        retained_pct = float(cum_pct[num_components - 1])
        fig, ax = plt.subplots()
        ax.plot(
            range(1, n_total + 1), cum_pct, marker="o", markersize=3,
            color=TOKENS.primary,
        )
        ax.axhline(
            variance_threshold, color="#c0392b", linestyle="--", linewidth=1,
            label=f"{variance_threshold}% threshold",
        )
        ax.axvline(num_components, color="#c0392b", linestyle=":", linewidth=1)
        ax.scatter([num_components], [retained_pct], color="#c0392b", zorder=5)
        ax.annotate(
            f"{num_components} components → {retained_pct:.1f}%",
            xy=(num_components, retained_pct),
            xytext=(0.5, min(retained_pct - 12, 88)),
            textcoords=("data", "data"),
            fontsize=9,
        )
        ax.set_xlabel("Number of components")
        ax.set_ylabel("Cumulative explained variance (%)")
        ax.set_title("Components retained at the variance threshold")
        ax.set_ylim(0, 100)
        ax.legend(loc="lower right")
        cumulative_png = _snapshot_fig(fig)

        # If ComBat ran, it stashed pre-correction PCA plots — show a before/after grid.
        before_batch = st.session_state.get("_pre_combat_pca_batch")
        before_group = st.session_state.get("_pre_combat_pca_group")
        if before_batch is not None and before_group is not None:
            fig, ax = plt.subplots()
            sc.pl.pca(st.session_state.adata, color="Batch", ax=ax, show=False)
            after_batch_png = _snapshot_fig(fig)
            fig, ax = plt.subplots()
            sc.pl.pca(st.session_state.adata, color="Group", ax=ax, show=False)
            after_group_png = _snapshot_fig(fig)
            scatter_items = [
                ("subheader", "PCA before vs after ComBat — colored by Batch"),
                ("image_row", [(before_batch, "Before ComBat"), (after_batch_png, "After ComBat")]),
                ("subheader", "PCA before vs after ComBat — colored by Group"),
                ("image_row", [(before_group, "Before ComBat"), (after_group_png, "After ComBat")]),
            ]
        else:
            fig, ax = plt.subplots()
            sc.pl.pca(st.session_state.adata, color="Group", ax=ax, show=False)
            pca_scatter_png = _snapshot_fig(fig)
            scatter_items = [
                ("write", "PCA Results by Group:"),
                ("image", pca_scatter_png),
            ]

        # Stash results for the review gate.
        pca_items = []
        if threshold_warning:
            pca_items.append(("warning", threshold_warning))
        pca_items += [
            ("image", variance_png),
            (
                "write",
                f"**{num_components} of {n_total}** components reach the "
                f"{variance_threshold}% cumulative-variance threshold "
                f"(retaining {retained_pct:.1f}%):",
            ),
            ("image", cumulative_png),
            *scatter_items,
            (
                "success",
                f"PCA completed in {time.time() - start_time:.2f} seconds. "
                f"Retained {num_components} components.",
            ),
        ]
        st.session_state["_pca_pending"] = pca_items
        st.rerun()

    st.stop()

# Step 2: Batch correction (optional). Runs before PCA so ComBat corrects X first.
if not st.session_state.batch_correction_done:
    section_header(
        "Batch Correction",
        subtitle="Correct for technical batch effects. Batches must be distinct from biological groups.",
        step=3,
    )

    _review_gate("_batch_pending", "batch_correction_done", "PCA")

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

    # Outside the correction form so edits persist across reruns and feed the diagnostic.
    st.markdown("**Batch assignment**")
    st.caption(
        "Assign each sample to a technical batch (e.g. acquisition day or instrument). "
        "Used both for the batch-effect check and for ComBat."
    )
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
        width="stretch",
        key="batch_assignment_editor",
    )
    batch_map = dict(
        zip(
            batch_assignment_df["SampleID"].astype(str),
            batch_assignment_df["Batch"].astype(str),
        )
    )

    # Diagnostic: is there actually a batch effect?
    if st.button("Check for batch effects", key="dp_check_batch"):
        st.session_state["_show_batch_diag"] = True
    if st.session_state.get("_show_batch_diag"):
        with st.spinner("Measuring batch effects…"):
            render_batch_diagnostics(adata, batch_map)

    # Correction choice.
    with st.form(key="batch_correction_form"):
        batch_correction_method = st.radio(
            "Select batch correction method",
            ("None", "ComBat"),
            key="batch_correction_method",
        )
        submit_button = st.form_submit_button(label="Apply", type="primary")

    if submit_button:
        if "Group" not in adata.obs.columns:
            st.error("The 'Group' column does not exist in adata.obs.")
            st.session_state.batch_correction_done = True
            st.session_state["_params_batch"] = {"method": "None"}
            st.session_state["_show_batch_diag"] = False
            # Drop the Batch column the read-only diagnostic left behind.
            adata.obs.drop(columns=["Batch"], inplace=True, errors="ignore")
            st.rerun()
        elif batch_correction_method == "None":
            st.session_state.batch_correction_done = True
            st.session_state["_params_batch"] = {"method": "None"}
            st.session_state["_show_batch_diag"] = False
            # No correction means no Batch annotation downstream.
            adata.obs.drop(columns=["Batch"], inplace=True, errors="ignore")
            st.info("No batch correction applied.")
            st.rerun()
        else:
            adata.obs["Batch"] = adata.obs["SampleID"].astype(str).map(batch_map)

            # Fix the per-marker scale pre-correction and reuse it after, so before/after EMD share a denominator.
            emd_scales = marker_robust_scales(adata)
            emd_before = compute_batch_emd(adata, batch_key="Batch", scales=emd_scales)

            # Capture pre-ComBat PCA now (before ComBat mutates adata.X in place), so the
            # PCA step can show a before/after comparison. Store only PNG bytes — bdata is
            # discarded — keeping memory flat on large datasets.
            bdata = adata.copy()
            sc.tl.pca(bdata)
            fig, ax = plt.subplots()
            sc.pl.pca(bdata, color="Batch", ax=ax, show=False)
            st.session_state["_pre_combat_pca_batch"] = _snapshot_fig(fig)
            fig, ax = plt.subplots()
            sc.pl.pca(bdata, color="Group", ax=ax, show=False)
            st.session_state["_pre_combat_pca_group"] = _snapshot_fig(fig)
            del bdata

            start_time = time.time()
            ok, msg = run_combat_correction(adata, batch_key="Batch", covariates=("Group",))
            elapsed_time = time.time() - start_time

            if not ok:
                st.error(msg)
                adata.obs.drop(columns=["Batch"], inplace=True, errors="ignore")
            else:
                emd_after = compute_batch_emd(adata, batch_key="Batch", scales=emd_scales)
                batch_items = [("success", f"{msg} (in {elapsed_time:.2f} seconds)")]

                # Before/after batch effect comparison.
                if not emd_before.empty and not emd_after.empty:
                    comp = emd_before[["marker", "emd_max"]].merge(
                        emd_after[["marker", "emd_max"]],
                        on="marker",
                        suffixes=("_before", "_after"),
                    )
                    comp["change"] = comp["emd_max_after"] - comp["emd_max_before"]
                    comp_styled = (
                        comp.rename(
                            columns={
                                "emd_max_before": "EMD before",
                                "emd_max_after": "EMD after",
                            }
                        )
                        .style.format(
                            {"EMD before": "{:.3f}", "EMD after": "{:.3f}", "change": "{:+.3f}"}
                        )
                        .hide(axis="index")
                    )
                    lim = max(comp["emd_max_before"].max(), comp["emd_max_after"].max(), 0.01) * 1.05
                    fig, ax = plt.subplots()
                    ax.scatter(comp["emd_max_before"], comp["emd_max_after"], color=TOKENS.primary)
                    ax.plot([0, lim], [0, lim], color=TOKENS.primary, linestyle="--", linewidth=1)
                    ax.set_xlim(0, lim)
                    ax.set_ylim(0, lim)
                    ax.set_xlabel("EMD before")
                    ax.set_ylabel("EMD after")
                    ax.set_title("Batch effect before vs after ComBat")
                    scatter_png = _snapshot_fig(fig)
                    n_improved = int((comp["change"] < 0).sum())
                    mean_before = float(comp["emd_max_before"].mean())
                    mean_after = float(comp["emd_max_after"].mean())
                    verdict = (
                        f"Batch effect reduced on {n_improved}/{len(comp)} markers "
                        f"(mean worst-batch EMD {mean_before:.3f} → {mean_after:.3f})."
                    )
                    batch_items += [
                        (
                            "success" if mean_after < mean_before else "warning",
                            verdict,
                        ),
                        ("subheader", "Batch effect before vs after"),
                        (
                            "caption",
                            "EMD per marker before and after ComBat. Negative change = "
                            "the batch effect shrank. Points below the diagonal improved.",
                        ),
                        ("dataframe", comp_styled),
                        ("image", scatter_png),
                    ]

                st.session_state.adata = adata
                st.session_state["_show_batch_diag"] = False
                st.session_state["_params_batch"] = {
                    "method": "ComBat",
                    "batch_key": "Batch",
                    "covariates": ["Group"],
                    "mapping": batch_map,
                    "emd_before": (
                        float(emd_before["emd_max"].mean()) if not emd_before.empty else None
                    ),
                    "emd_after": (
                        float(emd_after["emd_max"].mean()) if not emd_after.empty else None
                    ),
                }
                # Stash results for the review gate.
                st.session_state["_batch_pending"] = batch_items
                st.rerun()

    st.stop()

# Step 4: UMAP + Leiden. Final step: compute form, then the terminal
# results-and-download page. umap_computed stays False (no later step).
if not st.session_state.umap_computed:
    section_header(
        "UMAP + Leiden Clustering",
        subtitle="Compute the nearest-neighbor graph, UMAP embedding, and Leiden clusters.",
        step=5,
    )

    # Once clustering has run (_umap_pending set), this is the terminal page:
    # UMAPs, package the ZIP, and download inline.
    if st.session_state.get("_umap_pending") is not None:
        adata = st.session_state.get("adata")
        if adata is None:
            st.error("No processed data found. Please complete the pipeline first.")
            st.stop()

        st.success("Clustering complete.")

        umap_fig_leiden = st.session_state.get("_umap_fig_leiden")
        umap_fig_group = st.session_state.get("_umap_fig_group")
        if umap_fig_leiden is not None and umap_fig_group is not None:
            umap_col1, umap_col2 = st.columns(2)
            with umap_col1:
                st.image(umap_fig_leiden, width="stretch")
            with umap_col2:
                st.image(umap_fig_group, width="stretch")

        # One-line run summary.
        _pca = "PCA" if st.session_state.get("pca_selected") else "no PCA"
        _batch = st.session_state.get("_params_batch", {}).get("method", "None")
        _batch = "no batch correction" if _batch == "None" else f"{_batch} correction"
        st.caption(
            f"{adata.n_obs:,} cells  ·  {adata.obs['SampleID'].nunique()} samples  ·  "
            f"{adata.obs['leiden'].nunique()} Leiden clusters  ·  {_pca}  ·  {_batch}"
        )

        # UMAP PNGs to bundle, named here so build_outputs_zip stays state-agnostic.
        group_key = "Group" if "Group" in adata.obs else "SampleID"
        umap_pngs = {
            "umap_leiden_clusters": st.session_state.get("_umap_fig_leiden"),
            f"umap_{_safe_filename(group_key)}": st.session_state.get("_umap_fig_group"),
        }

        # Let the user take in the UMAPs first, then package. The 2s beat runs once
        # (guarded by _packaging_started) and we rerun into the packaging phase so the
        # spinner shows on a clean render, after the plots are on screen. Bytes are
        # cached so the download button is instant. Non-fatal on failure: the fallback
        # build button covers it.
        if st.session_state.get("_last_zip_buffer") is None:
            if not st.session_state.get("_packaging_started"):
                st.session_state["_packaging_started"] = True
                time.sleep(2)
                st.rerun()
            resolution = st.session_state.get("_umap_resolution", 0.5)
            try:
                with st.spinner("Packaging your files for download…"):
                    zip_bytes, zip_name = build_outputs_zip(adata, resolution, umap_pngs)
                st.session_state["_last_zip_buffer"] = zip_bytes
                st.session_state["_last_zip_name"] = zip_name
            except Exception as exc:  # noqa: BLE001 — surface, don't crash the page
                st.warning(
                    f"Couldn't package the download automatically ({exc}). "
                    "Use the button below to build it."
                )

        zip_buffer = st.session_state.get("_last_zip_buffer")
        zip_file_name = st.session_state.get("_last_zip_name", "analysis_outputs.zip")

        if zip_buffer is not None:
            st.download_button(
                label="Download all outputs (.zip)",
                data=zip_buffer,
                file_name=zip_file_name,
                mime="application/zip",
                type="primary",
                width="stretch",
            )
        elif st.button(
            "Download all outputs (.zip)",
            type="primary",
            key="dp_build_zip",
            width="stretch",
        ):
            # Fallback: build on demand, then rerun to render the ready button.
            resolution = st.session_state.get("_umap_resolution", 0.5)
            with st.spinner(
                "Preparing your download — this can take a moment for large datasets…"
            ):
                zip_bytes, zip_name = build_outputs_zip(adata, resolution, umap_pngs)
            st.session_state["_last_zip_buffer"] = zip_bytes
            st.session_state["_last_zip_name"] = zip_name
            st.rerun()

        with st.expander("What's in the download"):
            st.markdown(
                "- **`adata_*.h5ad`** — processed AnnData for the Visualization page\n"
                "- **`analysis_parameters_*.txt`** — every parameter used, a draft "
                "methods paragraph, and the CAFE citation\n"
                "- **`per_sample_clustered/`** — one CSV per sample with a "
                "`Leiden_Cluster` column and UMAP coordinates, ready to re-import "
                "into FlowJo\n"
                "- **`umap_plots/`** — the two UMAPs above (Leiden clusters and "
                "group) as PNGs\n"
                "- cluster count, frequency, and median-expression tables"
            )

        st.caption(
            "Next: upload the `.h5ad` to the Visualization page to explore clusters "
            "and build publication figures."
        )

        st.divider()

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
                "_packaging_started",
                "_umap_resolution",
                "_upload_review",
                "_upload_files",
                "_show_batch_diag",
                "_qc_pending",
                "_pca_pending",
                "_batch_pending",
                "_umap_pending",
                "_pre_combat_pca_batch",
                "_pre_combat_pca_group",
            ]:
                st.session_state.pop(key, None)
            st.rerun()

        st.stop()

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

        col3, col4 = st.columns(2)
        with col3:
            flavor_options = ["igraph", "leidenalg"]
            selected_flavor = st.selectbox(
                "Leiden flavor",
                options=flavor_options,
                index=0,
            )
        with col4:
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
                status_placeholder.markdown("**Step 4/4: Finalizing…**")

                # Two UMAPs (Leiden, group) snapshot to PNG. Fixed canvas + axes
                # position gives both identical pixel size regardless of legend width.
                group_key = "Group" if "Group" in adata.obs else "SampleID"

                def _umap_panel(color, title):
                    fig = plt.figure(figsize=(6, 5))
                    ax = fig.add_axes([0.12, 0.12, 0.60, 0.78])
                    sc.pl.umap(
                        adata,
                        color=color,
                        ax=ax,
                        show=False,
                        title=title,
                        frameon=True,
                    )
                    return _snapshot_fig(fig, tight=False)

                st.session_state["_umap_fig_leiden"] = _umap_panel(
                    "leiden", "Leiden clusters"
                )
                st.session_state["_umap_fig_group"] = _umap_panel(
                    group_key, group_key
                )

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

                # Drop the now-stale ZIP; stash resolution for output filenames.
                st.session_state.pop("_last_zip_buffer", None)
                st.session_state.pop("_last_zip_name", None)
                st.session_state["_umap_resolution"] = resolution
                st.session_state["leiden_computed"] = True

                progress += 25
                progress_bar.progress(progress)
                update_elapsed_time()
                status_placeholder.markdown("**Computation completed!**")
                # _umap_pending flips this page into its terminal state on rerun.
                st.session_state["_umap_pending"] = [
                    ("success", "Clustering complete."),
                ]
                st.rerun()

    st.stop()
