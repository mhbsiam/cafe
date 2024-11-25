import streamlit as st
import os
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
import cmasher as cmr
import pyarrow.pandas_compat
import pandas as pd
import time
import gc
import matplotlib.cm as cm
import pyarrow.csv as pv
import pyarrow as pa
from matplotlib.colors import ListedColormap
from scipy.stats import ttest_ind, mannwhitneyu
import io
import random
import tempfile
import zipfile

sc.settings.n_jobs = -1

def clear_memory():
    gc.collect()

#result_dir = os.path.join(os.getcwd(), 'result')
#os.makedirs(result_dir, exist_ok=True)

image_path = os.path.join('bin', 'img', 's_logo.png')
st.logo(image_path)

# Display an image
image_path = os.path.join('bin', 'img', 'logo_v2.png')
st.image(image_path, caption='', use_column_width=True)

st.title("Data Processing")

st.markdown(
    "*In this section, the app will read and process the uploaded files, combine them into a single DataFrame, and create an AnnData object for further analysis. It includes functionality for performing PCA, batch correction, UMAP, and unsupervised Leiden clustering. Results such as cluster counts, frequencies, and median fluorescence intensities are saved to CSV files, and the processed AnnData object is saved to the current directory/result folder. The resulting AnnData (h5ad file) can be used to create figures and run statistics.*")

option = st.radio(
    "Choose your option:",
    ["None", "Load CSV files"], help="Reload the page to start from the beginning"
)

if "adata" not in st.session_state:
    st.session_state.adata = None
if "uploaded_files" not in st.session_state:
    st.session_state.uploaded_files = None
if "batch_correction_done" not in st.session_state:
    st.session_state.batch_correction_done = False
if "umap_computed" not in st.session_state:
    st.session_state.umap_computed = False
if "pca_done" not in st.session_state:
    st.session_state.pca_done = False
if "leiden_computed" not in st.session_state:
    st.session_state.leiden_computed = False

def process_csv_files(uploaded_files):
    start_time = time.time()
    tables = []
    reference_columns = None

    #Progress bar and placeholders
    progress_bar = st.progress(0)
    elapsed_time_placeholder = st.empty()
    status_placeholder = st.empty()

    total_steps = len(uploaded_files) + 5
    current_step = 0

    #Read and process each CSV file
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

    #Combine all Arrow tables into a single table
    combined_table = pa.concat_tables(tables)
    current_step += 1
    progress_bar.progress(current_step / total_steps)
    status_placeholder.write("Combining all files into a single table")

    df = combined_table.to_pandas()

    if df.isnull().values.any():
        st.markdown(f"<span style='color:red; font-weight:bold;'>Missing values detected. Dropping missing rows.</span>", unsafe_allow_html=True)
        df = df.dropna()

    #Save the combined DataFrame
    #combined_output_path = os.path.join(result_dir, 'combined_updated_files.csv')
    #df.to_csv(combined_output_path, index=False)
    #st.write(f"Combined DataFrame saved to {combined_output_path}")

    expr_data = df.drop(columns=['SampleID', 'Group'])
    expr_data = expr_data.select_dtypes(include=[float, int])

    metadata = df[['SampleID', 'Group']]
    expr_data.index = expr_data.index.astype(str)
    metadata.index = metadata.index.astype(str)

    adata = sc.AnnData(expr_data)
    adata.obs = metadata
    adata.var_names = expr_data.columns.astype(str)
    st.write("Markers in adata:", adata.var_names[:50])

    st.session_state.adata = adata

    current_step += 1
    progress_bar.progress(1.0)
    elapsed_time = time.time() - start_time
    elapsed_time_placeholder.write(f"Data loading and processing completed in {elapsed_time:.2f} seconds")
    status_placeholder.write("Processing complete.")


if option == "Load CSV files":
    if st.session_state.adata is None:
        st.write("Please upload your CSV files.")
        uploaded_files = st.file_uploader("Choose CSV files", type="csv", accept_multiple_files=True)

        if uploaded_files:
            process_csv_files(uploaded_files)


    if 'adata' not in st.session_state:
        st.warning("Failed to initialise AnnData.")
    else:
        adata = st.session_state.adata

    if st.session_state.adata is not None and not st.session_state.get('pca_done', False):
        st.write("Choose whether to perform PCA on the data.")

        with st.form(key='pca_selection_form'):
            pca_option = st.radio(
                "Select PCA option",
                ('None', 'Perform PCA'), key='pca_option'
            )

            proceed_button = st.form_submit_button(label='Proceed')

        if proceed_button:
            st.session_state.pca_selected = (pca_option == 'Perform PCA')

            if not st.session_state.pca_selected:
                st.session_state.pca_done = True
                st.write("PCA not selected. Proceeding with batch correction.")
            else:
                st.write("PCA selected. Now choose PCA options below.")

    if st.session_state.get('pca_selected', False) and not st.session_state.get('pca_done', False):

        with st.form(key='pca_options_form'):
            pca_solver = st.selectbox(
                "Choose the SVD solver for PCA",
                ('Auto', 'Full', 'Arpack', 'Randomized'), key='pca_solver'
            )

            variance_threshold = st.slider(
                "Select the explained variance threshold (%) to retain",
                min_value=70, max_value=99, value=95, key='variance_threshold'
            )

            apply_pca_button = st.form_submit_button(label='Apply PCA')

        if apply_pca_button:
            st.write("Running PCA........")
            progress_bar = st.progress(0)
            progress = 0
            #max_progress = 100

            start_time = time.time()
            adata = st.session_state.adata
            progress += 25
            progress_bar.progress(progress)

            sc.tl.pca(adata, svd_solver=pca_solver.lower())

            progress += 25
            progress_bar.progress(progress)

            pca_variance_ratio = adata.uns['pca']['variance_ratio'].cumsum()
            num_components = (pca_variance_ratio >= variance_threshold / 100).argmax() + 1

            progress += 25
            progress_bar.progress(progress)

            adata.obsm['X_pca'] = adata.obsm['X_pca'][:, :num_components]
            adata.varm['PCs'] = adata.varm['PCs'][:, :num_components]
            adata.uns['pca']['variance'] = adata.uns['pca']['variance'][:num_components]
            adata.uns['pca']['variance_ratio'] = adata.uns['pca']['variance_ratio'][:num_components]

            st.session_state.pca_done = True
            st.session_state.pca_choice = True

            progress += 25
            progress_bar.progress(progress)

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

            pca_variance = adata.uns['pca']['variance_ratio']
            num_components = len(pca_variance)

            fig, ax = plt.subplots()
            ax.bar(range(1, num_components + 1), pca_variance, color='skyblue')
            ax.set_xlabel('Principal Component')
            ax.set_ylabel('Explained Variance Ratio')
            ax.set_title('Explained Variance by PCA Components')
            ax.set_xticks(range(1, num_components + 1))

            st.pyplot(fig)
            plt.close(fig)
            st.write("Generating Figure........")

            st.write("PCA Results by Group:")
            fig, ax = plt.subplots()
            sc.pl.pca(st.session_state.adata, color='Group', ax=ax, show=False)
            st.pyplot(fig)

            elapsed_time = time.time() - start_time
            st.write(f"Completed in {elapsed_time:.2f} seconds")

            clear_memory()

    st.session_state.adata = adata

    if st.session_state.get('pca_done', False) and not st.session_state.batch_correction_done:
        st.subheader("Batch correction")

        with st.form(key='batch_correction_form'):
            batch_correction_method = st.radio(
                "Select batch correction method",
                ('None', 'ComBat'), key='batch_correction_method'
            )
            submit_button = st.form_submit_button(label='Apply')

        if submit_button:
            adata = st.session_state.adata

            if 'Group' in adata.obs.columns:
                if batch_correction_method == 'ComBat':
                    start_time = time.time()
                    sc.pp.combat(adata, key='Group')
                    elapsed_time = time.time() - start_time
                    st.write(f"Batch correction completed using ComBat in {elapsed_time:.2f} seconds")

                    st.subheader("PCA After Batch Correction")

                    # Plot PCA After Correction
                    st.write("PCA After Batch Correction")
                    fig, ax = plt.subplots()
                    sc.pl.pca(st.session_state.adata, color='Group', ax=ax, show=False)
                    st.pyplot(fig)

                elif batch_correction_method == 'None':
                    st.write("No batch correction applied.")
            else:
                st.write("The 'Group' column does not exist in adata.obs.")

            st.session_state.adata = adata
            st.session_state.batch_correction_done = True

            clear_memory()

    if st.session_state.batch_correction_done and not st.session_state.umap_computed:
        st.write("Batch correction has been completed.")

        with st.form(key='umap_leiden_form'):
            resolution = st.slider("Leiden Resolution", min_value=0.0, max_value=2.5, value=0.5, step=0.01,
                                help="Lower values yield fewer, larger clusters; higher values yield more, smaller clusters.")

            n_neighbors = st.slider("n_neighbors for neighbors computation", min_value=5, max_value=50, value=30, step=1,
                                    help="Controls the local neighborhood size.")
            min_dist = st.slider("UMAP min_dist", min_value=0.0, max_value=1.0, value=0.1, step=0.01,
                                help="Controls how tightly UMAP packs points together.")

            flavor_options = ['igraph', 'leidenalg']
            selected_flavor = st.selectbox("Select a flavor for leiden community detection", options=flavor_options, index=0)

            dim_option = ['None', 'X_pca', 'X_umap']
            selected_dim = st.selectbox("Select a dimension reduction metric to use", options=dim_option, index=1, help="Select None if you have not reduced the data using PCA")

            metric_options = ['euclidean', 'manhattan', 'cosine', 'correlation', 'minkowski']
            selected_metric = st.selectbox("Select distance metric for neighbors computation",
                                        options=metric_options,
                                        index=0,  # Default to 'euclidean'
                                        help="Determines how distances between data points are calculated.")

            submit_button = st.form_submit_button(label='Apply and Compute')

        if submit_button:
            if not st.session_state.get('leiden_computed'):
                adata = st.session_state.get('adata')

                progress_bar = st.progress(0)
                elapsed_time_placeholder = st.empty()
                status_placeholder = st.empty()

                progress = 0
                max_progress = 100
                start_time = time.time()

                def update_elapsed_time():
                    elapsed = time.time() - start_time
                    elapsed_time_placeholder.markdown(f"**Elapsed Time:** {elapsed:.2f} seconds")

                with st.spinner("Working....."):
                    update_elapsed_time()
                    status_placeholder.markdown("**Starting computations...**")
                    progress_bar.progress(progress)

                    update_elapsed_time()
                    status_placeholder.markdown("**Step 1/4: Computing neighbors...**")
                    if st.session_state.get('pca_choice'):
                        st.write("Started clustering.")
                        sc.pp.neighbors(adata, n_neighbors=n_neighbors, method='umap', use_rep='X_pca', metric=selected_metric, random_state=50)
                    else:
                        st.write("Started clustering.")
                        sc.pp.neighbors(adata, n_neighbors=n_neighbors, method='umap', use_rep=None, metric=selected_metric, random_state=50)
                    progress += 25
                    progress_bar.progress(progress)

                    update_elapsed_time()
                    status_placeholder.markdown(f"**Step 2/4: Running UMAP with n_neighbors={n_neighbors}, min_dist={min_dist}, and metric={selected_metric}...**")
                    sc.tl.umap(adata, min_dist=min_dist, random_state=50)
                    progress += 25
                    progress_bar.progress(progress)

                    update_elapsed_time()
                    status_placeholder.markdown("**Step 3/4: Computing Leiden clustering...**")
                    sc.tl.leiden(adata, resolution=resolution, random_state=50, flavor="igraph", n_iterations=2, directed=False)
                    progress += 25
                    progress_bar.progress(progress)

                    update_elapsed_time()
                    status_placeholder.markdown("**Step 4/4: Saving outputs...**")

                    random_number = random.randint(1000, 9999)
                    zip_file_name = f"analysis_outputs_{random_number}.zip"

                    zip_buffer = io.BytesIO()

                    with zipfile.ZipFile(zip_buffer, "w") as zip_file:
                        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as temp_file:
                            temp_path = temp_file.name
                            adata.write(temp_path)

                        zip_file.write(temp_path, arcname=f"adata_flowanalysis_{random_number}.h5ad")
                        os.remove(temp_path)

                        cluster_sample_counts = adata.obs.groupby(['leiden', 'SampleID']).size().unstack(fill_value=0)
                        cluster_sample_counts_csv = io.BytesIO()
                        cluster_sample_counts.to_csv(cluster_sample_counts_csv)
                        cluster_sample_counts_csv.seek(0)
                        zip_file.writestr(f"cluster_sample_counts_{random_number}.csv", cluster_sample_counts_csv.getvalue())

                        sample_totals = cluster_sample_counts.sum()
                        cluster_sample_percentages = cluster_sample_counts.div(sample_totals) * 100
                        cluster_sample_frequencies_csv = io.BytesIO()
                        cluster_sample_percentages.to_csv(cluster_sample_frequencies_csv)
                        cluster_sample_frequencies_csv.seek(0)
                        zip_file.writestr(f"cluster_sample_frequencies_{random_number}.csv", cluster_sample_frequencies_csv.getvalue())

                        median_df = adata.to_df()
                        median_df['leiden'] = adata.obs['leiden']
                        median_df['SampleID'] = adata.obs['SampleID']
                        medians = median_df.groupby(['leiden', 'SampleID']).median()
                        median_csv = io.BytesIO()
                        medians.to_csv(median_csv)
                        median_csv.seek(0)
                        zip_file.writestr(f"median_expression_{random_number}.csv", median_csv.getvalue())

                    zip_buffer.seek(0)

                    st.download_button(
                        label="Download All Outputs as ZIP",
                        data=zip_buffer,
                        file_name=zip_file_name,
                        mime="application/zip",
                    )

                    progress += 25
                    progress_bar.progress(progress)

                    update_elapsed_time()
                    status_placeholder.markdown("**Computation completed!**")
                    st.session_state['leiden_computed'] = True
                    st.session_state['umap_computed'] = True

                    st.balloons()

                    st.write("**⚡ Pro Tip: To visualize the results, load the AnnData (h5ad file) in the Visualization tab**")

else:
    st.write("Please upload CSV files.")
