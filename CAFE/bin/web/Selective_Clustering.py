import streamlit as st
import os
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import pyarrow.csv as pv
import pyarrow as pa
import pandas as pd
import time
import gc
import io
import zipfile
import tempfile
import random

st.set_page_config(layout="centered")

sc.settings.n_jobs = -1

image_path = os.path.join('bin', 'img', 's_logo.png')
st.logo(image_path)

# Display an image
image_path = os.path.join('bin', 'img', 'logo_v2.png')
st.image(image_path, caption='', use_container_width=True)

def clear_memory():
    gc.collect()

def process_csv_files(uploaded_files):
    start_time = time.time()
    tables = []
    reference_columns = None

    progress_bar = st.progress(0)
    elapsed_time_placeholder = st.empty()
    status_placeholder = st.empty()

    total_steps = len(uploaded_files) + 5
    current_step = 0

    for idx, uploaded_file in enumerate(uploaded_files):
        current_step += 1
        progress_percentage = current_step / total_steps
        progress_bar.progress(progress_percentage)
        status_placeholder.write(f"Processing file {idx + 1}/{len(uploaded_files)}: {uploaded_file.name}")
        table = pv.read_csv(uploaded_file)

        numeric_columns = [col for col in table.column_names if pa.types.is_integer(table.schema.field(col).type)]
        for col in numeric_columns:
            table = table.set_column(table.column_names.index(col), col, table.column(col).cast(pa.float64()))

        if 'SampleID' in table.column_names:
            table = table.set_column(table.column_names.index('SampleID'), 'SampleID', table.column('SampleID').cast(pa.string()))
        if 'Group' in table.column_names:
            table = table.set_column(table.column_names.index('Group'), 'Group', table.column('Group').cast(pa.string()))

        if reference_columns is None:
            reference_columns = set(table.column_names)
        else:
            current_columns = set(table.column_names)
            if current_columns != reference_columns:
                st.warning(f"Warning: Column names in {uploaded_file.name} do not match the reference columns!")
                st.write(f"Found columns: {current_columns}")
                continue

        file_name = uploaded_file.name
        name_parts = file_name.replace('.csv', '').split('_')

        if len(name_parts) != 2:
            st.write(f"Skipping file {file_name} as it does not have the expected format 'sampleID_group.csv'")
            continue

        table = table.append_column('SampleID', pa.array([name_parts[0]] * len(table)))
        table = table.append_column('Group', pa.array([name_parts[1]] * len(table)))

        columns_to_drop = [col for col in table.column_names if col.lower() == 'sample']
        if columns_to_drop:
            table = table.drop(columns_to_drop)

        tables.append(table)

    combined_table = pa.concat_tables(tables)
    current_step += 1
    progress_bar.progress(current_step / total_steps)
    status_placeholder.write("Combining all files into a single table")

    df = combined_table.to_pandas()

    if df.isnull().values.any():
        st.markdown(
            f"<span style='color:red; font-weight:bold;'>Missing values detected. Dropping missing rows.</span>",
            unsafe_allow_html=True
        )
        df = df.dropna()

    current_step += 1
    progress_bar.progress(1.0)
    elapsed_time = time.time() - start_time
    elapsed_time_placeholder.write(f"Data loading and processing completed in {elapsed_time:.2f} seconds")
    status_placeholder.write("Processing complete.")

    return df

def perform_batch_correction(adata):
    st.subheader("Batch correction")

    with st.form(key='batch_correction_form'):
        batch_correction_method = st.radio(
            "Select batch correction method",
            ('None', 'ComBat')
        )
        submit_button = st.form_submit_button(label='Apply')

    if submit_button:
        if 'Group' in adata.obs.columns:
            if batch_correction_method == 'ComBat':
                start_time = time.time()
                sc.pp.combat(adata, key='Group')
                elapsed_time = time.time() - start_time
                st.write(f"Batch correction completed using ComBat in {elapsed_time:.2f} seconds")
            else:
                st.write("No batch correction applied.")
        else:
            st.write("The 'Group' column does not exist in adata.obs.")

        st.session_state.batch_correction_done = True

    return adata

def compute_umap_leiden(adata):
    st.write("Proceeding with UMAP and Leiden clustering.")

    with st.form(key='umap_leiden_form'):
        resolution = st.slider(
            "Leiden Resolution",
            min_value=0.1,
            max_value=2.5,
            value=0.5,
            step=0.1,
            help="Lower values yield fewer, larger clusters; higher values yield more, smaller clusters."
        )

        n_neighbors = st.slider(
            "n_neighbors for neighbors computation",
            min_value=5,
            max_value=50,
            value=30,
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

        metric_options = [
            'euclidean', 'manhattan', 'cosine', 'correlation',
            'chebyshev', 'canberra', 'minkowski', 'hamming'
        ]
        selected_metric = st.selectbox(
            "Select distance metric for neighbors computation",
            options=metric_options,
            index=0,
            help="Determines how distances between data points are calculated."
        )

        submit_button = st.form_submit_button(label='Apply and Compute')


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

            output_h5ad_path = os.path.join(result_dir, 'adata_flowanalysis.h5ad')
            adata.write(output_h5ad_path)

            cluster_sample_counts = adata.obs.groupby(['leiden', 'SampleID']).size().unstack(fill_value=0)

            sample_totals = cluster_sample_counts.sum()
            cluster_sample_percentages = cluster_sample_counts.div(sample_totals) * 100

            median_df = adata.to_df()
            median_df['leiden'] = adata.obs['leiden']
            median_df['SampleID'] = adata.obs['SampleID']
            medians = median_df.groupby(['leiden', 'SampleID']).median()

            def create_zip_with_outputs(selected_files):
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                    for file_path, arcname in selected_files:
                        zip_file.write(file_path, arcname=arcname)
                        os.remove(file_path)  # Clean up the temporary file immediately
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
            clear_memory()
            st.session_state['leiden_computed'] = True
            st.session_state['umap_computed'] = True

    return adata

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

st.warning("This module is only available for offline use. Please download the App from https://github.com/mhbsiam/cafe/releases and run it locally.")

st.title("Semi Supervised Clustering")
st.write("*Semi supervised clustering allows a user to select a number of markers or a cell type based on marker selection to construct Leiden cluster. This is useful to find out a specific cell population.*")

option = st.radio(
    "Choose your option:",
    ["None", "Subset by removing markers", "Subset by cell type"]
)
