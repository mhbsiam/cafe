import os
import time
import gc

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc
import streamlit as st
from plotly.subplots import make_subplots
from scipy.spatial.distance import cdist
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score

from theme import IMG_DIR, apply_theme, page_header
from utils import process_directory_csv_files, validate_path

st.set_page_config(layout="centered")
apply_theme()

sc.settings.n_jobs = -1

st.logo(os.path.join(IMG_DIR, 's_logo.png'))

page_header(
    "Evaluate Leiden Clustering",
    subtitle="Generate Leiden clusterings across a range of parameters and compare clustering-quality metrics.",
)

# ── Parameter sections ────────────────────────────────────────────────────────

with st.expander("Data I/O", expanded=True):
    st.caption("Paths to the folder containing your per-sample CSV files and to the folder where processed AnnData objects and UMAP PNGs will be saved.")
    input_dir = st.text_input("Input CSV Directory Path")
    output_dir = st.text_input("Output Directory Path")

with st.expander("PCA settings", expanded=True):
    st.caption("Control how many principal components are carried forward. A higher variance threshold retains more PCs and preserves more biological signal at the cost of compute time.")
    pca_solver = st.selectbox("PCA SVD solver", ['auto', 'full', 'arpack', 'randomized'])
    variance_threshold = st.slider("Variance Threshold for PCA (%)", min_value=50, max_value=100, value=95)

with st.expander("UMAP & Leiden settings", expanded=True):
    st.caption("These parameters govern both the embedding and the clustering sweep. Each combination of resolution produces a separate AnnData and UMAP image.")
    leiden_resolutions = st.text_input("Leiden Resolutions (comma-separated)", value="0.1,0.3,0.5,0.6,0.7")
    col1, col2 = st.columns(2)
    with col1:
        n_neighbors = st.selectbox("UMAP n_neighbors", [10, 15, 20, 30, 50], index=1)
        min_dist = st.selectbox("UMAP min_dist", [0.1, 0.2, 0.3, 0.4, 0.5])
    with col2:
        distance_metric = st.selectbox("Distance Metric", ['euclidean', 'manhattan', 'cosine'])
        leiden_flavor = st.selectbox("Leiden Clustering Flavor", ['igraph', 'leidenalg'])

# ── Helper functions ──────────────────────────────────────────────────────────

def process_csv_files(input_dir, output_dir):
    # Delegate to utility loader
    return process_directory_csv_files(input_dir)

def perform_pca(adata, pca_solver, variance_threshold, output_dir):
    st.write("Running PCA...")
    start_time = time.time()
    sc.tl.pca(adata, svd_solver=pca_solver.lower())
    pca_variance_ratio = adata.uns['pca']['variance_ratio'].cumsum()
    threshold_mask = pca_variance_ratio >= variance_threshold / 100
    if not threshold_mask.any():
        st.warning(f"No principal component reaches {variance_threshold}% variance. Using all {len(pca_variance_ratio)} components.")
        num_components = len(pca_variance_ratio)
    else:
        num_components = threshold_mask.argmax() + 1

    adata.obsm['X_pca'] = adata.obsm['X_pca'][:, :num_components]
    adata.varm['PCs'] = adata.varm['PCs'][:, :num_components]
    adata.uns['pca']['variance'] = adata.uns['pca']['variance'][:num_components]
    adata.uns['pca']['variance_ratio'] = adata.uns['pca']['variance_ratio'][:num_components]

    st.write(f"PCA applied using the {pca_solver} solver and retaining {num_components} components.")

    fig, ax = plt.subplots()
    ax.plot(range(1, len(pca_variance_ratio) + 1), pca_variance_ratio * 100, marker='o', color='skyblue')
    ax.axhline(y=variance_threshold, color='r', linestyle='--', label=f'{variance_threshold}% Variance')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Cumulative Explained Variance (%)')
    ax.set_title('Cumulative Explained Variance by PCA Components')
    ax.set_xticks(range(1, len(pca_variance_ratio) + 1))
    ax.legend()
    st.pyplot(fig, width='stretch')
    plt.close(fig)

def run_umap_leiden(adata, resolutions, n_neighbors_list, min_dist, metric, output_dir, flavor):

    total_steps = len(resolutions) * len(n_neighbors_list)
    progress_bar = st.progress(0)
    status_text = st.empty()
    current_step = 0

    silhouette_scores, ch_scores, db_scores, inertia_values, num_clusters = [], [], [], [], []
    resolutions_used, neighbors_used = [], []

    for n_neighbors in n_neighbors_list:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, random_state=50)
        sc.tl.umap(adata, min_dist=min_dist, random_state=50)

        for res_idx, resolution in enumerate(resolutions):
            status_text.caption(
                f"Step {current_step + 1} of {total_steps}: "
                fr"n\_neighbors={n_neighbors}, resolution={resolution}"
            )
            sc.tl.leiden(adata, resolution=resolution, random_state=50, flavor=flavor, n_iterations=2, directed=False)

            X = adata.obsm['X_pca']
            cluster_labels = adata.obs['leiden'].astype(int)
            unique_labels = np.unique(cluster_labels)

            if len(unique_labels) >= 2:
                silhouette_avg = silhouette_score(X, cluster_labels)
                ch_score = calinski_harabasz_score(X, cluster_labels)
                db_score = davies_bouldin_score(X, cluster_labels)
            else:
                silhouette_avg = np.nan
                ch_score = np.nan
                db_score = np.nan

            silhouette_scores.append(silhouette_avg)
            ch_scores.append(ch_score)
            db_scores.append(db_score)

            centroids = np.array([X[cluster_labels == label].mean(axis=0) for label in unique_labels])
            inertia = np.sum(np.min(cdist(X, centroids, 'euclidean'), axis=1) ** 2)
            inertia_values.append(inertia)
            num_clusters.append(len(unique_labels))
            resolutions_used.append(resolution)
            neighbors_used.append(n_neighbors)

            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.umap(adata, color='leiden', size=10, ax=ax, show=False)
            umap_plot_path = os.path.join(output_dir, f'umap_n{n_neighbors}_r{resolution}.png')
            plt.savefig(umap_plot_path, dpi=300, bbox_inches='tight')
            st.pyplot(fig, width='stretch')
            st.caption(f"Saved: umap_n{n_neighbors}_r{resolution}.png")
            plt.close(fig)

            adata_output_path = os.path.join(output_dir, f'adata_n{n_neighbors}_r{resolution}.h5ad')
            adata.write(adata_output_path)

            current_step += 1
            progress_bar.progress(current_step / total_steps)

    status_text.empty()

    # ── Results table ─────────────────────────────────────────────────────────
    result_df = pd.DataFrame({
        'n_neighbors': neighbors_used,
        'Resolution': resolutions_used,
        'Silhouette Score': silhouette_scores,
        'Calinski-Harabasz Score': ch_scores,
        'Davies-Bouldin Score': db_scores,
        'Number of Clusters': num_clusters,
        'Inertia': inertia_values
    })
    st.subheader("Clustering Evaluation Metrics")
    st.dataframe(result_df, width='stretch')

    # ── Metric interpretation guide ───────────────────────────────────────────
    with st.expander("How to read these metrics", expanded=False):
        st.markdown(
            """
**Silhouette Score** (range −1 to 1, **higher is better**)
Measures how similar each cell is to its own cluster compared to neighbouring clusters.
Values above 0.5 indicate well-separated clusters; values near 0 suggest overlapping clusters.

**Calinski-Harabasz Score** (**higher is better**)
Ratio of between-cluster dispersion to within-cluster dispersion.
A higher value means clusters are dense and well-separated relative to each other.

**Davies-Bouldin Score** (≥ 0, **lower is better**)
Average similarity of each cluster to its most similar neighbour.
Values close to 0 indicate compact, well-separated clusters.

**Inertia: Elbow Method** (**look for the "elbow"**)
Sum of squared distances from each cell to its cluster centroid.
Plot against number of clusters and look for the point where the curve bends sharply. Adding more clusters beyond that point yields diminishing returns.
            """
        )

    # ── Unified Plotly metric dashboard ───────────────────────────────────────
    st.subheader("Metric Overview")
    st.caption("Interactive plot. Hover for exact values. All four metrics are plotted together for quick comparison.")

    # Build hover text from resolution labels
    hover_labels = [
        f"n_neighbors={nb}<br>resolution={r}<br>clusters={nc}"
        for nb, r, nc in zip(neighbors_used, resolutions_used, num_clusters)
    ]

    # CAFE palette: red accent + two complementary tones + neutral
    _RED   = "#ff4b4b"
    _TEAL  = "#2a9d8f"
    _AMBER = "#e76f51"
    _SLATE = "#457b9d"

    _LINE_OPTS = dict(width=2.5)
    _MARKER_OPTS = dict(size=8, line=dict(width=1.5, color="white"))

    dashboard = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "Silhouette Score (higher ↑ better)",
            "Calinski-Harabasz Score (higher ↑ better)",
            "Davies-Bouldin Score (lower ↓ better)",
            "Inertia: Elbow Method",
        ),
        vertical_spacing=0.16,
        horizontal_spacing=0.12,
    )

    dashboard.add_trace(
        go.Scatter(
            x=resolutions_used, y=silhouette_scores,
            mode="lines+markers",
            name="Silhouette",
            line=dict(color=_RED, **_LINE_OPTS),
            marker=dict(color=_RED, **_MARKER_OPTS),
            hovertemplate="%{text}<br>Silhouette: %{y:.4f}<extra></extra>",
            text=hover_labels,
        ),
        row=1, col=1,
    )

    dashboard.add_trace(
        go.Scatter(
            x=resolutions_used, y=ch_scores,
            mode="lines+markers",
            name="Calinski-Harabasz",
            line=dict(color=_TEAL, **_LINE_OPTS),
            marker=dict(color=_TEAL, **_MARKER_OPTS),
            hovertemplate="%{text}<br>CH Score: %{y:,.1f}<extra></extra>",
            text=hover_labels,
        ),
        row=1, col=2,
    )

    dashboard.add_trace(
        go.Scatter(
            x=resolutions_used, y=db_scores,
            mode="lines+markers",
            name="Davies-Bouldin",
            line=dict(color=_AMBER, **_LINE_OPTS),
            marker=dict(color=_AMBER, **_MARKER_OPTS),
            hovertemplate="%{text}<br>DB Score: %{y:.4f}<extra></extra>",
            text=hover_labels,
        ),
        row=2, col=1,
    )

    dashboard.add_trace(
        go.Scatter(
            x=num_clusters, y=inertia_values,
            mode="lines+markers",
            name="Inertia",
            line=dict(color=_SLATE, **_LINE_OPTS),
            marker=dict(color=_SLATE, **_MARKER_OPTS),
            hovertemplate="%{text}<br>Clusters: %{x}<br>Inertia: %{y:,.0f}<extra></extra>",
            text=hover_labels,
        ),
        row=2, col=2,
    )

    dashboard.update_layout(
        height=560,
        showlegend=False,
        font=dict(family="Arial, sans-serif", size=12, color="#31333F"),
        paper_bgcolor="#FFFFFF",
        plot_bgcolor="#F8F9FA",
        margin=dict(l=10, r=10, t=48, b=10),
    )

    # Axis labels
    dashboard.update_xaxes(title_text="Leiden Resolution", row=1, col=1, gridcolor="#E8E8E8")
    dashboard.update_xaxes(title_text="Leiden Resolution", row=1, col=2, gridcolor="#E8E8E8")
    dashboard.update_xaxes(title_text="Leiden Resolution", row=2, col=1, gridcolor="#E8E8E8")
    dashboard.update_xaxes(title_text="Number of Clusters", row=2, col=2, gridcolor="#E8E8E8")
    dashboard.update_yaxes(title_text="Silhouette Score", row=1, col=1, gridcolor="#E8E8E8")
    dashboard.update_yaxes(title_text="CH Score", row=1, col=2, gridcolor="#E8E8E8")
    dashboard.update_yaxes(title_text="DB Score", row=2, col=1, gridcolor="#E8E8E8")
    dashboard.update_yaxes(title_text="Inertia", row=2, col=2, gridcolor="#E8E8E8")

    st.plotly_chart(dashboard, width='stretch')


if st.button('Run Analysis', type="primary"):
    with st.spinner("Working..."):
        if not input_dir or not output_dir:
            st.error("Please specify both input and output directories.")
        else:
            # Validate paths
            try:
                input_dir = validate_path(input_dir)
                output_dir = validate_path(output_dir)
            except ValueError as e:
                st.error(str(e))
                st.stop()
            adata = process_csv_files(input_dir, output_dir)

            if adata is not None:
                resolutions = [float(r.strip()) for r in leiden_resolutions.split(',') if r.strip()]
                perform_pca(adata, pca_solver, variance_threshold, output_dir)
                run_umap_leiden(adata, resolutions, [n_neighbors], min_dist, distance_metric, output_dir, leiden_flavor)
                gc.collect()
