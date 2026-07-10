import io
import os
import random
import tempfile
import time
import zipfile

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pv
import scanpy as sc
import streamlit as st

from theme import IMG_DIR, apply_theme, page_header
from utils import process_uploaded_csv_to_df, run_combat_correction

st.set_page_config(layout="centered")
apply_theme()

# A2: deployment-mode flag (set by run_w.py for web deployment)
IS_WEB = st.session_state.get('cafe_deployment') == 'web'

sc.settings.n_jobs = -1

st.logo(os.path.join(IMG_DIR, 's_logo.png'))

def process_csv_files(uploaded_files):
    # A4: delegate to shared utils implementation (returns DataFrame)
    return process_uploaded_csv_to_df(uploaded_files)

def perform_batch_correction(adata):
    with st.expander("Batch correction", expanded=True):
        st.caption(
            "Use a *technical* batch distinct from the biological Group "
            "(e.g. acquisition day, instrument) — correcting on Group erases "
            "the differences you're testing for. Assign each sample below; "
            "batches must span multiple groups."
        )

        sample_ids = sorted(pd.unique(adata.obs['SampleID']))
        sample_group = dict(zip(adata.obs['SampleID'], adata.obs['Group']))

        with st.form(key='batch_correction_form'):
            batch_correction_method = st.radio(
                "Select batch correction method",
                ('None', 'ComBat')
            )

            st.caption("Batch assignment (used only when ComBat is selected):")
            batch_assignment_df = st.data_editor(
                pd.DataFrame({
                    'SampleID': sample_ids,
                    'Group': [sample_group.get(sid, '') for sid in sample_ids],
                    'Batch': ['1'] * len(sample_ids),
                }),
                disabled=['SampleID', 'Group'],
                hide_index=True,
                width='stretch',
                key='selective_batch_assignment_editor',
            )

            submit_button = st.form_submit_button(label='Apply', type="primary")

    if submit_button:
        if 'Group' not in adata.obs.columns:
            st.write("The 'Group' column does not exist in adata.obs.")
            st.session_state.batch_correction_done = True
        elif batch_correction_method == 'None':
            st.write("No batch correction applied.")
            st.session_state.batch_correction_done = True
        else:
            batch_map = dict(zip(
                batch_assignment_df['SampleID'].astype(str),
                batch_assignment_df['Batch'].astype(str),
            ))
            adata.obs['Batch'] = adata.obs['SampleID'].astype(str).map(batch_map)

            start_time = time.time()
            ok, msg = run_combat_correction(adata, batch_key='Batch', covariates=('Group',))
            elapsed_time = time.time() - start_time

            if not ok:
                st.error(msg)
                adata.obs.drop(columns=['Batch'], inplace=True, errors='ignore')
            else:
                st.write(f"{msg} (in {elapsed_time:.2f} seconds)")
                st.session_state.batch_correction_done = True

    return adata

def compute_umap_leiden(adata):
    with st.expander("UMAP & Leiden settings", expanded=True):
        st.caption(
            "These parameters control the embedding and clustering. "
            "Lower resolution yields fewer, larger clusters; higher resolution yields more, smaller ones."
        )

        with st.form(key='umap_leiden_form'):
            resolution = st.slider(
                "Leiden Resolution",
                min_value=0.01,
                max_value=2.5,
                value=0.5,
                step=0.01,
                help="Lower values yield fewer, larger clusters; higher values yield more, smaller clusters."
            )

            col1, col2 = st.columns(2)
            with col1:
                n_neighbors = st.slider(
                    "n_neighbors",
                    min_value=5,
                    max_value=50,
                    value=15,
                    step=1,
                    help="Controls the local neighborhood size."
                )
                min_dist = st.slider(
                    "UMAP min_dist",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.1,
                    step=0.01,
                    help="Controls how tightly UMAP packs points together."
                )
            with col2:
                metric_options = [
                    'euclidean', 'manhattan', 'cosine', 'correlation',
                    'chebyshev', 'canberra', 'minkowski', 'hamming'
                ]
                selected_metric = st.selectbox(
                    "Distance metric",
                    options=metric_options,
                    index=0,
                    help="Determines how distances between data points are calculated."
                )

            submit_button = st.form_submit_button(label='Apply and Compute', type="primary")

    if submit_button:

        with st.spinner("Working....."):

            random_state = 50
            if st.session_state.get('pca_done', False) and 'X_pca' in adata.obsm:
                st.write("Using PCA results for clustering.")
                sc.pp.neighbors(
                    adata,
                    n_neighbors=n_neighbors,
                    method='umap',
                    use_rep='X_pca',
                    metric=selected_metric,
                    random_state=random_state
                )
            else:
                st.write("Using raw data for clustering.")
                sc.pp.neighbors(
                    adata,
                    n_neighbors=n_neighbors,
                    method='umap',
                    use_rep=None,
                    metric=selected_metric,
                    random_state=random_state
                )

            sc.tl.umap(adata, min_dist=min_dist, random_state=random_state)

            sc.tl.leiden(
                adata,
                resolution=resolution,
                flavor="igraph",
                n_iterations=2,
                directed=False,
                random_state=random_state
            )

            cluster_sample_counts = adata.obs.groupby(['leiden', 'SampleID'], observed=False).size().unstack(fill_value=0)

            sample_totals = cluster_sample_counts.sum()
            cluster_sample_percentages = cluster_sample_counts.div(sample_totals) * 100

            median_df = adata.to_df()
            median_df['leiden'] = adata.obs['leiden']
            median_df['SampleID'] = adata.obs['SampleID']
            medians = median_df.groupby(['leiden', 'SampleID'], observed=False).median()

            def create_zip_with_outputs(selected_files):
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for file_path, arcname in selected_files:
                        zip_file.write(file_path, arcname=arcname)
                        os.remove(file_path)
                zip_buffer.seek(0)
                return zip_buffer

            random_number = random.randint(1000, 9999)

            selected_files = []

            # Use temporary files
            with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_h5ad:
                adata.write(temp_h5ad.name)
                selected_files.append((temp_h5ad.name, f"adata_flowanalysis_{random_number}.h5ad"))

            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as temp_csv1:
                cluster_sample_counts.to_csv(temp_csv1.name)
                selected_files.append((temp_csv1.name, f"cluster_sample_counts_{random_number}.csv"))

            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as temp_csv2:
                cluster_sample_percentages.to_csv(temp_csv2.name)
                selected_files.append((temp_csv2.name, f"cluster_sample_frequencies_{random_number}.csv"))

            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as temp_csv3:
                medians.to_csv(temp_csv3.name)
                selected_files.append((temp_csv3.name, f"median_expression_{random_number}.csv"))

            zip_buffer = create_zip_with_outputs(selected_files)

            st.download_button(
                label="Download Selected Results as ZIP",
                data=zip_buffer,
                file_name=f"Selected_clustering_{random_number}.zip",
                mime="application/zip"
            )

            del cluster_sample_counts
            del sample_totals
            del cluster_sample_percentages
            del median_df
            del medians
            st.session_state['leiden_computed'] = True
            st.session_state['umap_computed'] = True

    return adata

# ── Session state initialisation ─────────────────────────────────────────────

if 'adata' not in st.session_state:
    st.session_state['adata'] = None
if 'pca_done' not in st.session_state:
    st.session_state['pca_done'] = False
if 'batch_correction_done' not in st.session_state:
    st.session_state['batch_correction_done'] = False
if 'leiden_computed' not in st.session_state:
    st.session_state['leiden_computed'] = False
if 'umap_computed' not in st.session_state:
    st.session_state['umap_computed'] = False
if 'perform_pca' not in st.session_state:
    st.session_state['perform_pca'] = False


page_header(
    "Semi Supervised Clustering",
    subtitle="Build Leiden clusters from a chosen subset of markers (or a marker-defined cell type) to isolate a specific cell population.",
)

# A2: In web mode, show offline-only warning and halt
if IS_WEB:
    st.warning("This module is only available for offline use. Please download the App from https://github.com/mhbsiam/cafe/releases and run it locally.")
    st.stop()

option = st.radio(
    "Choose your option:",
    ["Subset by removing markers", "Subset by cell type"]
)

if st.button("Apply", type="primary"):
    # B10 fix: only clear selective-clustering-specific keys, not all session state.
    # Wiping everything destroys adata and other pages' state.
    for key in list(st.session_state.keys()):
        if key.startswith('selective_') or key.startswith('subset_'):
            del st.session_state[key]
    st.session_state['adata'] = None
    st.session_state['pca_done'] = False
    st.session_state['batch_correction_done'] = False
    st.session_state['leiden_computed'] = False
    st.session_state['umap_computed'] = False
    st.session_state['perform_pca'] = False
    st.session_state['option_confirmed'] = option
    st.rerun()

# Only show the pipeline after Apply has been pressed at least once
_confirmed_option = st.session_state.get('option_confirmed')
if _confirmed_option in ['Subset by removing markers', 'Subset by cell type']:
    option = _confirmed_option  # use the confirmed option for all downstream logic

if _confirmed_option in ['Subset by removing markers', 'Subset by cell type']:

    # ── Pipeline progress indicator ───────────────────────────────────────────
    # Determine current step (1–5) from session state flags
    _steps = [
        "Upload & subset",
        "PCA",
        "Batch correction",
        "UMAP & Leiden",
        "Complete",
    ]
    if st.session_state.get('umap_computed', False):
        _current_step = 5
    elif st.session_state.get('batch_correction_done', False):
        _current_step = 4
    elif st.session_state.get('pca_done', False):
        _current_step = 3
    elif st.session_state['adata'] is not None:
        _current_step = 2
    else:
        _current_step = 1

    _step_label = _steps[_current_step - 1]
    st.progress(_current_step / len(_steps))
    st.caption(f"**Step {_current_step} of {len(_steps)} — {_step_label}**")

    st.markdown("---")

    # ── Step 1: Upload & subset ───────────────────────────────────────────────
    if st.session_state['adata'] is None:
        with st.expander("Upload & subset configuration", expanded=True):
            if option == 'Subset by removing markers':
                st.caption(
                    "Upload your per-sample CSV files. Then select the markers to **retain** "
                    "for AnnData creation — deselect any markers you want to exclude from clustering."
                )
            else:
                st.caption(
                    "Upload your per-sample CSV files. Then choose a marker whose positive "
                    "expression (> 0) defines the cell type subset you want to re-cluster."
                )

            uploaded_files = st.file_uploader("Choose CSV files", type="csv", accept_multiple_files=True)

            if not uploaded_files:
                st.warning("Please upload CSV files to proceed.")
            else:
                df = process_csv_files(uploaded_files)

                if df is None:
                    st.error(
                        "No valid CSV files could be processed. Check that filenames follow "
                        "'sampleID_group.csv' and that columns match across files."
                    )
                    st.stop()

                if option == 'Subset by removing markers':
                    all_columns = df.columns.tolist()
                    marker_columns = [col for col in all_columns if col not in ['SampleID', 'Group', 'Sample']]

                    selected_markers = st.multiselect(
                        "Select markers to retain for AnnData creation",
                        marker_columns,
                        default=marker_columns
                    )

                    if st.button("Submit Marker Selection", type="primary"):

                        expr_data = df[selected_markers]
                        expr_data = expr_data.select_dtypes(include=[float, int])

                        st.write("Selected Markers:")
                        st.write(selected_markers)
                        st.write(f"Expression data shape: {expr_data.shape}")

                        metadata = df[['SampleID', 'Group']]
                        expr_data.index = expr_data.index.astype(str)
                        metadata.index = metadata.index.astype(str)

                        adata = sc.AnnData(expr_data)
                        adata.obs = metadata
                        adata.var_names = expr_data.columns.astype(str)
                        st.session_state['adata'] = adata
                        st.write("AnnData object created.")

                elif option == 'Subset by cell type':

                    expr_data = df.drop(columns=['SampleID', 'Group'])
                    expr_data = expr_data.select_dtypes(include=[float, int])

                    available_markers = expr_data.columns.tolist()
                    selected_marker = st.selectbox(
                        "Select a marker to subset the data based on its expression:",
                        available_markers,
                        help="Cells with zero or negative expression of this marker will be excluded."
                    )

                    # Live cell count preview — updates reactively as marker changes
                    if selected_marker:
                        n_total_cells = len(expr_data)
                        n_passing = int((expr_data[selected_marker] > 0).sum())
                        pct_passing = n_passing / n_total_cells * 100 if n_total_cells > 0 else 0
                        col_a, col_b = st.columns(2)
                        with col_a:
                            st.metric("Cells passing filter", f"{n_passing:,}")
                        with col_b:
                            st.metric("% of total", f"{pct_passing:.1f}%")
                        st.caption(f"Cells with {selected_marker} > 0 out of {n_total_cells:,} total")

                    if st.button("Apply Subsetting", type="primary"):
                        st.write(f"Original data shape: {expr_data.shape}")

                        valid_cells = expr_data[selected_marker] > 0
                        expr_data = expr_data[valid_cells]
                        df_subset = df.loc[expr_data.index]

                        st.write(f"Data shape after subsetting based on '{selected_marker}': {expr_data.shape}")

                        metadata = df_subset[['SampleID', 'Group']]
                        expr_data.index = expr_data.index.astype(str)
                        metadata.index = metadata.index.astype(str)

                        adata = sc.AnnData(expr_data)
                        adata.obs = metadata
                        adata.var_names = expr_data.columns.astype(str)
                        st.session_state['adata'] = adata
                        st.write("AnnData object created after subsetting.")
    else:
        adata = st.session_state['adata']

    # ── Step 2: PCA ───────────────────────────────────────────────────────────
    if st.session_state['adata'] is not None and not st.session_state.get('pca_done', False):

        with st.expander("PCA", expanded=True):
            st.caption("Dimensionality reduction before clustering. PCA is optional — skip it if your panel is small or you prefer to cluster in marker space directly.")

            with st.form(key='pca_selection_form'):
                pca_option = st.radio(
                    "Select PCA option",
                    ('None', 'Perform PCA'),
                    index=1 if st.session_state.get('perform_pca', False) else 0,
                    key='pca_option_radio'
                )

                submit_pca_selection = st.form_submit_button(label='Proceed', type="primary")

        if submit_pca_selection:
            if pca_option == 'Perform PCA':
                st.session_state['perform_pca'] = True
            else:
                st.session_state['perform_pca'] = False
                st.session_state['pca_done'] = True
                st.write("PCA not selected. Proceeding without PCA.")

    if st.session_state.get('perform_pca', False) and not st.session_state.get('pca_done', False):

        with st.expander("PCA settings", expanded=True):
            st.caption("Choose the SVD solver and the variance threshold that determines how many principal components to retain.")

            with st.form(key='pca_options_form'):

                pca_solver = st.selectbox(
                    "Choose the SVD solver for PCA",
                    ('auto', 'full', 'arpack', 'randomized'),
                    key='pca_solver'
                )

                variance_threshold = st.slider(
                    "Select the explained variance threshold (%) to retain",
                    min_value=70, max_value=99,
                    value=95,
                    key='variance_threshold'
                )

                apply_pca_button = st.form_submit_button(label='Apply PCA', type="primary")

        if apply_pca_button:
            adata = st.session_state['adata']
            st.write("Running PCA...")
            sc.tl.pca(adata, svd_solver=pca_solver.lower())

            pca_variance_ratio = adata.uns['pca']['variance_ratio'].cumsum()

            # B12 fix: handle case where no component reaches the variance threshold
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

            st.write(f"PCA applied using the {pca_solver} solver and retaining {num_components} components to meet the {variance_threshold}% explained variance threshold.")

            fig, ax = plt.subplots()
            ax.plot(range(1, len(pca_variance_ratio) + 1), pca_variance_ratio * 100, marker='o', color='skyblue')
            ax.axhline(y=variance_threshold, color='r', linestyle='--', label=f'{variance_threshold}% Variance')
            ax.set_xlabel('Principal Component')
            ax.set_ylabel('Cumulative Explained Variance (%)')
            ax.set_title('Cumulative Explained Variance by PCA Components')
            ax.set_xticks(range(1, len(pca_variance_ratio) + 1))
            ax.legend()
            st.pyplot(fig)
            plt.close(fig)

            st.write("PCA Results by Group:")
            sc.pl.pca(adata, color='Group', show=True)

            st.session_state.pca_done = True
            st.session_state.adata = adata

    # ── Step 3: Batch correction ───────────────────────────────────────────────
    if st.session_state.get('pca_done', False) and not st.session_state.get('batch_correction_done', False):
        adata = perform_batch_correction(st.session_state['adata'])
        st.session_state['adata'] = adata

    # ── Step 4: UMAP & Leiden ─────────────────────────────────────────────────
    if st.session_state.get('batch_correction_done', False) and not st.session_state.get('umap_computed', False):
        adata = compute_umap_leiden(st.session_state['adata'])
        st.session_state['adata'] = adata
else:
    st.caption("Select an option above and press **Apply** to begin.")
