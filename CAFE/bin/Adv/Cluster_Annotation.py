import io
import os
import random
import re
import tempfile
import zipfile
from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import scanpy as sc
import scipy.cluster.hierarchy as sch
import streamlit as st
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.sparse import issparse

from theme import apply_theme, page_header

st.set_page_config(layout="centered")
apply_theme()

st.logo(os.path.join('bin', 'img', 's_logo.png'))

page_header(
    "Advanced Annotation",
    subtitle="Annotate Leiden clusters into cell types and save the AnnData file.",
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
        sc.tl.umap(st.session_state.adata)

    # Clear any stale per-cluster annotation values from a previous session
    for key in [k for k in st.session_state if k.startswith("annotation_")]:
        del st.session_state[key]

# ── Annotation workspace ──────────────────────────────────────────────────────

if st.session_state.adata is not None:
    adata = st.session_state.adata

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

    # ── Two-column workspace: UMAP left, annotation table right ───────────────
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

    # ── Apply button ──────────────────────────────────────────────────────────
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
