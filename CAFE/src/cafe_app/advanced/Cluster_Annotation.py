import io
import os
import random
import tempfile
import zipfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import streamlit as st

from theme import IMG_DIR, apply_theme, page_header, section_header
from utils import safe_flatten, safe_sort_clusters
from celltype_gmm import (
    compute_positive_posteriors,
    score_cell_types,
    assign_cell_types,
    suggest_thresholds,
)

st.set_page_config(layout="centered")
apply_theme()

st.logo(os.path.join(IMG_DIR, 's_logo.png'))

page_header(
    "Cluster Annotation",
    subtitle="Annotate Leiden clusters into cell types and save the AnnData file.",
)

# Confidence-level colours reused across the probabilistic UMAP and legends.
_CONF_COLORS = {"High": "#0f5070", "Low": "#f0a500", "Ambiguous": "#c0392b"}

# High-contrast, colourblind-safe qualitative palette (Okabe–Ito order tuned so
# the first two groups read as blue vs vermillion) for the Group overlay.
_GROUP_COLORS = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#000000",
]


def _styled_umap(adata, color, palette, legend_loc="right margin"):
    """Spine-less, tick-less UMAP coloured by `color`, legend in the right margin; returns fig."""
    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.umap(
        adata, color=color, palette=palette, ax=ax, show=False, legend_loc=legend_loc
    )
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    return fig


# ── Method selection (first choice on the page) ───────────────────────────────

method = st.radio(
    "Annotation method",
    ["Manual", "Probabilistic (GMM)"],
    horizontal=True,
    help=(
        "Manual: name each Leiden cluster yourself. "
        "Probabilistic: define marker signatures and let CAFE gate every "
        "marker with per-sample Gaussian mixtures and assign cells by confidence."
    ),
)

# ── Upload ────────────────────────────────────────────────────────────────────

with st.form(key='upload_form'):
    uploaded_file = st.file_uploader("Upload AnnData (.h5ad) file", type="h5ad")
    submit_upload = st.form_submit_button("Load AnnData", type="primary")

if 'adata' not in st.session_state:
    st.session_state.adata = None

if uploaded_file and submit_upload:
    st.session_state.adata = sc.read_h5ad(uploaded_file)
    st.write(f"AnnData loaded with shape: {st.session_state.adata.shape}")

    if 'X_umap' not in st.session_state.adata.obsm.keys():
        if 'neighbors' not in st.session_state.adata.uns:
            sc.pp.neighbors(st.session_state.adata)
        sc.tl.umap(st.session_state.adata)

    # Clear any stale per-cluster annotation values from a previous session
    for key in [k for k in st.session_state if k.startswith("annotation_")]:
        del st.session_state[key]

# ── Annotation workspace ──────────────────────────────────────────────────────

if st.session_state.adata is not None:
    adata = st.session_state.adata

    # Manual mode (unchanged behaviour)
    if method == "Manual":
        # Detect clusters and sort numerically where possible
        leiden_clusters = sorted(
            adata.obs['leiden'].astype(str).unique(),
            key=lambda x: int(x) if x.isdigit() else x
        )
        n_total = len(leiden_clusters)

        # Ensure session_state keys exist (first render after load)
        for c in leiden_clusters:
            if f"annotation_{c}" not in st.session_state:
                st.session_state[f"annotation_{c}"] = ""

        # ── Two-column workspace: UMAP left, annotation table right ───────────
        col_map, col_ann = st.columns([1, 1], gap="large")

        with col_map:
            st.markdown("**Cluster map**")
            st.caption("Use cluster numbers on the map to fill in names on the right.")
            fig, ax = plt.subplots(figsize=(5, 5))
            sc.pl.umap(
                adata,
                color='leiden',
                palette="tab20c",
                ax=ax,
                show=False,
                legend_loc='on data'
            )
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title('')
            st.pyplot(fig, width='stretch')
            plt.close(fig)

        with col_ann:
            st.markdown("**Cluster names**")

            # Live coverage indicator
            n_filled = sum(
                1 for c in leiden_clusters
                if st.session_state.get(f"annotation_{c}", "").strip()
            )
            coverage_fraction = n_filled / n_total if n_total > 0 else 0
            st.progress(coverage_fraction)
            st.caption(f"{n_filled} of {n_total} clusters annotated")

            st.markdown("---")

            # Per-cluster text inputs — keyed into session_state so values persist
            # across reruns triggered by typing in other rows
            for c in leiden_clusters:
                st.text_input(
                    f"Cluster {c}",
                    key=f"annotation_{c}",
                    placeholder="e.g. CD4 T cell",
                )

        st.markdown("---")

        # ── Apply button ──────────────────────────────────────────────────────
        if st.button("Apply Cell Type Annotations", type="primary"):
            # Build mapping from per-cluster inputs
            cell_type_mapping = {
                c: st.session_state[f"annotation_{c}"].strip()
                for c in leiden_clusters
            }

            # Validate: detect empty (unmapped) clusters
            unmapped_clusters = sorted(
                [c for c, name in cell_type_mapping.items() if not name],
                key=lambda x: int(x) if x.isdigit() else x
            )

            if unmapped_clusters:
                st.error(
                    "Every Leiden cluster must be annotated before continuing — "
                    "unmapped cells would be silently dropped from all downstream "
                    "counts and frequencies. Missing annotation for cluster(s): "
                    + ", ".join(unmapped_clusters)
                )
            else:
                adata.obs['cell_type'] = (
                    adata.obs['leiden'].astype(str).map(cell_type_mapping)
                )
                st.success("Cell type annotations applied successfully!")

                random_number = random.randint(1000, 9999)
                zip_buffer = io.BytesIO()

                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_file:
                        temp_path = temp_file.name
                        adata.write_h5ad(temp_path)
                        zip_file.write(temp_path, arcname=f"annotated_adata_{random_number}.h5ad")
                    os.remove(temp_path)

                zip_buffer.seek(0)

                st.download_button(
                    label="Download Annotated Data as ZIP",
                    data=zip_buffer,
                    file_name=f"annotated_data_{random_number}.zip",
                    mime="application/zip"
                )

                st.write("UMAP visualization with 'cell_type' annotations:")
                fig, ax = plt.subplots()

                sc.pl.umap(
                    adata,
                    color='cell_type',
                    palette="tab20c",
                    ax=ax,
                    show=False,
                    legend_loc='on data'
                )

                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.set_xticks([])
                ax.set_yticks([])

                st.pyplot(fig)

                umap_buffer = io.BytesIO()
                fig.savefig(umap_buffer, dpi=300, format='png', bbox_inches='tight')
                umap_buffer.seek(0)
                plt.close(fig)

                st.download_button(
                    label="Download UMAP as PNG",
                    data=umap_buffer,
                    file_name=f"UMAP_cell_type_{random_number}.png",
                    mime="image/png"
                )

    # Probabilistic mode (per-sample GMM gating)
    else:
        st.markdown(
            "Define each cell type by its marker signature. CAFE fits a "
            "per-sample 2-component Gaussian mixture to every marker, turns the "
            "positive-component posterior into a per-marker confidence, and "
            "assigns each cell to its best-matching type."
        )

        markers = list(adata.var_names)
        if "SampleID" not in adata.obs.columns:
            st.warning(
                "No 'SampleID' column found — marker gates will be fit globally "
                "across all cells instead of per sample."
            )

        # ── Marker signature editor ───────────────────────────────────────────
        section_header(
            "Marker signatures",
            subtitle="Mark each defining marker + (positive) or − (negative). Leave blank to ignore.",
        )

        template = pd.DataFrame({"Cell type": ["", "", ""]})
        for m in markers:
            template[m] = ""

        column_config = {
            "Cell type": st.column_config.TextColumn("Cell type", width="medium"),
        }
        for m in markers:
            column_config[m] = st.column_config.SelectboxColumn(
                m, options=["", "+", "-"], width="small"
            )

        edited = st.data_editor(
            template,
            column_config=column_config,
            num_rows="dynamic",
            hide_index=True,
            width='stretch',
            key="celltype_def_editor",
        )

        # Build {cell_type: {marker: +/-}} from non-empty rows.
        definitions = {}
        for _, row in edited.iterrows():
            name = str(row["Cell type"]).strip()
            if not name:
                continue
            states = {m: row[m] for m in markers if row[m] in ("+", "-")}
            if states:
                definitions[name] = states

        # ── Confidence thresholds ─────────────────────────────────────────────
        section_header("Confidence thresholds")

        auto_tune = st.checkbox(
            "Auto-tune thresholds from the data (recommended)",
            value=True,
            help=(
                "Place each threshold at the natural valley in the score / margin "
                "distributions your data produces, instead of guessing. Falls back "
                "to sensible defaults when a distribution has no clear split."
            ),
        )
        st.caption(
            "The sliders below apply **only when Auto-tune is off**. With Auto-tune "
            "on, CAFE picks these from the data and reports the values it used. "
            "For all three, **higher = stricter** (more Unassigned/Ambiguous); "
            "**lower = more permissive**."
        )
        c1, c2, c3 = st.columns(3)
        with c1:
            min_score = st.slider(
                "Min match score", 0.0, 1.0, 0.5, 0.05,
                help=(
                    "Minimum score a cell needs to match a type at all.\n\n"
                    "**Higher = stricter** → more cells become 'Unassigned'.\n"
                    "**Lower** → more cells get a call, even weak ones."
                ),
            )
        with c2:
            high_margin = st.slider(
                "High-confidence margin", 0.0, 1.0, 0.5, 0.05,
                help=(
                    "How far the best type must beat the runner-up to be labelled 'High'.\n\n"
                    "**Higher = stricter** → fewer 'High' cells, more 'Low'.\n"
                    "**Lower** → more cells qualify as 'High'."
                ),
            )
        with c3:
            low_margin = st.slider(
                "Ambiguous margin", 0.0, 1.0, 0.15, 0.05,
                help=(
                    "If the best type barely beats the runner-up (margin below this), "
                    "the call is 'Ambiguous' — e.g. double-positive cells.\n\n"
                    "**Higher = stricter** → more cells flagged 'Ambiguous'.\n"
                    "**Lower** → fewer 'Ambiguous'."
                ),
            )
        seed = int(st.number_input("Random seed", min_value=0, value=0, step=1))

        if not definitions:
            st.info("Add at least one cell type with one or more +/− markers to enable the run.")

        if st.button("Run probabilistic annotation", type="primary", disabled=not definitions):
            progress_bar = st.progress(0.0, text="Fitting per-sample marker gates…")
            posteriors, diagnostics = compute_positive_posteriors(
                adata,
                random_state=seed,
                progress=lambda f: progress_bar.progress(min(f, 1.0)),
            )
            progress_bar.empty()

            scores = score_cell_types(posteriors, definitions)

            if auto_tune:
                tuned = suggest_thresholds(scores)
                min_score = tuned["min_score"]
                high_margin = tuned["high_margin"]
                low_margin = tuned["low_margin"]
                st.info(
                    f"Auto-tuned thresholds — min match score **{min_score:.2f}**, "
                    f"high-confidence margin **{high_margin:.2f}**, "
                    f"ambiguous margin **{low_margin:.2f}**."
                )
            elif high_margin < low_margin:
                st.error(
                    "The high-confidence margin must be greater than or equal to the "
                    f"ambiguous margin (currently high = {high_margin:.2f} < "
                    f"ambiguous = {low_margin:.2f}). Adjust the sliders so the "
                    "'High' band sits above the 'Ambiguous' band."
                )
                st.stop()

            calls = assign_cell_types(
                scores, min_score=min_score, high_margin=high_margin, low_margin=low_margin
            )

            adata.obs["cell_type"] = pd.Categorical(calls["cell_type"].to_numpy())
            adata.obs["cell_type_confidence"] = calls["confidence_score"].to_numpy()
            adata.obs["confidence_level"] = pd.Categorical(
                calls["confidence_level"].to_numpy(),
                categories=["High", "Low", "Ambiguous"],
            )
            if scores.shape[1]:
                adata.obsm["celltype_scores"] = scores.to_numpy()
                adata.uns["celltype_score_names"] = list(scores.columns)

            st.session_state["prob_diagnostics"] = diagnostics
            st.success("Probabilistic annotation applied.")

        # ── Results (persist across reruns via adata.obs) ─────────────────────
        if {"cell_type", "confidence_level"} <= set(adata.obs.columns):
            section_header("Results")

            counts = (
                adata.obs.groupby(["cell_type", "confidence_level"], observed=True)
                .size()
                .unstack(fill_value=0)
            )
            st.markdown("**Cells per type × confidence level**")
            st.dataframe(counts, width='stretch')

            umap_col1, umap_col2 = st.columns(2)
            with umap_col1:
                st.markdown("**UMAP — cell type**")
                fig_ct = _styled_umap(adata, "cell_type", "tab20c")
                st.pyplot(fig_ct, width='stretch')
                plt.close(fig_ct)
            with umap_col2:
                st.markdown("**UMAP — confidence**")
                fig_cf = _styled_umap(
                    adata, "confidence_level", _CONF_COLORS, legend_loc="right margin"
                )
                st.pyplot(fig_cf, width='stretch')
                plt.close(fig_cf)

            # ── Confident-subset UMAPs: High + Low only (Ambiguous excluded) ──
            st.markdown("---")
            level = adata.obs["confidence_level"].astype(str)
            subset = adata[level.isin(["High", "Low"])].copy()
            n_sub = subset.n_obs
            if n_sub == 0:
                st.info(
                    "No High- or Low-confidence cells to plot — every cell was "
                    "called Ambiguous."
                )
            else:
                sub_col1, sub_col2 = st.columns(2)
                with sub_col1:
                    st.markdown(
                        f"**UMAP — cell type** · High + Low subset (n = {n_sub:,})"
                    )
                    st.caption("Ambiguous-confidence cells excluded.")
                    fig_sub_ct = _styled_umap(subset, "cell_type", "tab20c")
                    st.pyplot(fig_sub_ct, width='stretch')
                    plt.close(fig_sub_ct)
                with sub_col2:
                    if "Group" in subset.obs.columns:
                        st.markdown(
                            f"**UMAP — group** · High + Low subset (n = {n_sub:,})"
                        )
                        st.caption("Ambiguous-confidence cells excluded.")
                        # Drop carried-over Group colours so scanpy uses the palette below.
                        subset.uns.pop("Group_colors", None)
                        n_groups = subset.obs["Group"].nunique()
                        grp_palette = (
                            _GROUP_COLORS if n_groups <= len(_GROUP_COLORS) else "tab20"
                        )
                        fig_sub_grp = _styled_umap(subset, "Group", grp_palette)
                        st.pyplot(fig_sub_grp, width='stretch')
                        plt.close(fig_sub_grp)
                    else:
                        st.markdown("**UMAP — group**")
                        st.info("No 'Group' column found in this AnnData.")

            # Per-cluster purity — validates gates against the Leiden clustering.
            if "leiden" in adata.obs.columns:
                st.markdown("**Per-cluster purity** — dominant cell type per Leiden cluster")
                leiden_str = adata.obs["leiden"].astype(str)
                purity_rows = []
                for cl in safe_sort_clusters(leiden_str.unique()):
                    sub = adata.obs.loc[leiden_str == str(cl), "cell_type"]
                    if len(sub) == 0:
                        continue
                    vc = sub.value_counts()
                    purity_rows.append({
                        "Leiden": cl,
                        "Dominant cell type": vc.index[0],
                        "Purity": round(float(vc.iloc[0]) / len(sub), 3),
                        "n cells": int(len(sub)),
                    })
                st.dataframe(pd.DataFrame(purity_rows), hide_index=True, width='stretch')

            # ── GMM gate diagnostics ──────────────────────────────────────────
            diag = st.session_state.get("prob_diagnostics")
            with st.expander("GMM gate diagnostics"):
                if diag is not None and not diag.empty:
                    show = diag.copy()
                    show["cutoff"] = show["cutoff"].round(2)
                    show["separation"] = show["separation"].round(2)
                    cols = ["SampleID", "marker", "cutoff", "separation", "informative"]
                    if "source" in show.columns:
                        cols.append("source")
                    st.dataframe(show[cols], hide_index=True, width='stretch')
                    st.caption(
                        "`source` shows which gate was used per sample: **per-sample** "
                        "(own bimodal fit), **pooled** (per-sample fit wasn't bimodal, so "
                        "the all-sample gate was applied), or **neutral** (unimodal even "
                        "pooled → posterior 0.5, marker ignored for that sample)."
                    )

                    if "SampleID" in adata.obs.columns:
                        sample_opts = sorted(adata.obs["SampleID"].astype(str).unique())
                    else:
                        sample_opts = ["__all__"]
                    d1, d2 = st.columns(2)
                    sel_sample = d1.selectbox("Sample", sample_opts)
                    sel_marker = d2.selectbox("Marker", markers)

                    if "SampleID" in adata.obs.columns:
                        mask = adata.obs["SampleID"].astype(str).to_numpy() == sel_sample
                    else:
                        mask = np.ones(adata.n_obs, dtype=bool)
                    vals = safe_flatten(adata[mask, sel_marker].X)

                    row = show[(show["SampleID"] == sel_sample) & (show["marker"] == sel_marker)]
                    cutoff = float(row["cutoff"].iloc[0]) if not row.empty else float("nan")

                    fig_d, ax_d = plt.subplots(figsize=(5, 3))
                    ax_d.hist(vals, bins=80, color="#9cc3d5", edgecolor="#0f5070")
                    if np.isfinite(cutoff):
                        ax_d.axvline(cutoff, color="#c0392b", ls="--", label=f"gate ≈ {cutoff:.1f}")
                        ax_d.legend()
                    ax_d.set_xlabel(f"{sel_marker} (scaled expression)")
                    ax_d.set_ylabel("cells")
                    st.pyplot(fig_d, width='stretch')
                    plt.close(fig_d)
                else:
                    st.caption("Run the annotation to see per-sample marker gates.")

            # ── Downloads: full + high-conf + high-and-low-conf ────────────────
            section_header("Download")
            st.caption(
                "One ZIP with the full annotated AnnData, a high-confidence subset "
                "(High only), and a high-and-low-confidence subset (High + Low, "
                "with Ambiguous cells left out)."
            )
            if st.button("Prepare download (full + high-conf + high&low-conf)"):
                random_number = random.randint(1000, 9999)
                level = adata.obs["confidence_level"].astype(str)
                subsets = {
                    "annotated_full": adata,
                    "high_confidence": adata[level == "High"],
                    "high_and_low_confidence": adata[level.isin(["High", "Low"])],
                }
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for name, ad_obj in subsets.items():
                        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_file:
                            temp_path = temp_file.name
                            ad_obj.copy().write_h5ad(temp_path)
                            zip_file.write(temp_path, arcname=f"{name}_{random_number}.h5ad")
                        os.remove(temp_path)
                zip_buffer.seek(0)
                st.download_button(
                    label="Download annotated data as ZIP",
                    data=zip_buffer,
                    file_name=f"probabilistic_annotation_{random_number}.zip",
                    mime="application/zip",
                )
