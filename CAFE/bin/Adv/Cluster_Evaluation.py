import os
import time
import gc
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import pyarrow.csv as pv
import pyarrow as pa
import matplotlib.pyplot as plt
import streamlit as st
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from scipy.spatial.distance import cdist

# Set Scanpy settings
sc.settings.n_jobs = -1

#image_path = os.path.join('bin', 'CAFE_logo.png')
#st.image(image_path, caption='', use_column_width=True)

image_path = os.path.join('bin', 'img', 'small_logo.png')
st.logo(image_path)

# Define the result directory at the beginning
result_dir = os.path.join(os.getcwd(), 'result')
os.makedirs(result_dir, exist_ok=True)

# Streamlit page title
st.title("Evaluate Leiden Clustering")
st.write('*This module allows creating multiple Leiden clustering files with varied paremeters. Each adata is saved and can be used to evaluate clustering quality.*')


# Input fields for file directories
input_dir = st.text_input("Input CSV Directory Path")
output_dir = st.text_input("Output Directory Path")

# Dropdowns and sliders for PCA options, UMAP, and Leiden
pca_solver = st.selectbox("PCA SVD solver", ['auto', 'full', 'arpack', 'randomized'])
variance_threshold = st.slider("Variance Threshold for PCA", min_value=50, max_value=100, value=95)

# Input for Leiden resolutions
leiden_resolutions = st.text_input("Leiden Resolutions (comma-separated)", value="0.1,0.3,0.5,0.6,0.7")

# Dropdown for UMAP parameters
n_neighbors = st.selectbox("UMAP n_neighbors", [10, 15, 20, 30, 50], index=1)
min_dist = st.selectbox("UMAP min_dist", [0.1, 0.2, 0.3, 0.4, 0.5])
distance_metric = st.selectbox("Distance Metric", ['euclidean', 'manhattan', 'cosine'])

# Add an option to choose the flavor for Leiden resolution
leiden_flavor = st.selectbox("Leiden Clustering Flavor", ['igraph', 'leidenalg'])

def process_csv_files(input_dir, output_dir):
    start_time = time.time()
    tables = []

    # Get list of CSV files
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    if not csv_files:
        st.write("No CSV files found in the input directory.")
        return None

    # Process each CSV file
    for idx, csv_file in enumerate(csv_files):
        st.write(f"Processing file {idx + 1}/{len(csv_files)}: {csv_file}")
        file_path = os.path.join(input_dir, csv_file)
        table = pv.read_csv(file_path)

        #Convert all into a consistent data type
        numeric_columns = [col for col in table.column_names if pa.types.is_integer(table.schema.field(col).type)]
        for col in numeric_columns:
            # Convert integer columns to float (double) type
            table = table.set_column(table.column_names.index(col), col, table.column(col).cast(pa.float64()))

        #Convert SampleID to string**
        if 'SampleID' in table.column_names:
            table = table.set_column(table.column_names.index('SampleID'), 'SampleID', table.column('SampleID').cast(pa.string()))
        if 'Group' in table.column_names:
            table = table.set_column(table.column_names.index('Group'), 'Group', table.column('Group').cast(pa.string()))


        # Extract SampleID and Group from filename
        name_parts = csv_file.replace('.csv', '').split('_')
        if len(name_parts) != 2:
            st.write(f"Skipping file {csv_file} as it does not have the expected format 'sampleID_group.csv'")
            continue

        # Add SampleID and Group columns
        table = table.append_column('SampleID', pa.array([name_parts[0]] * len(table)))
        table = table.append_column('Group', pa.array([name_parts[1]] * len(table)))
        tables.append(table)

    if not tables:
        st.write("No valid CSV files were processed.")
        return None

    combined_table = pa.concat_tables(tables)
    df = combined_table.to_pandas()
    df.dropna(inplace=True)

    # Create AnnData object
    expr_data = df.drop(columns=['SampleID', 'Group']).select_dtypes(include=[np.number])
    metadata = df[['SampleID', 'Group']]
    adata = sc.AnnData(expr_data)
    adata.obs = metadata
    adata.var_names = expr_data.columns.astype(str)

    st.write(f"Data loading and processing completed in {time.time() - start_time:.2f} seconds")
    return adata

def perform_pca(adata, pca_solver, variance_threshold, output_dir):
    st.write("Running PCA...")
    start_time = time.time()
    sc.tl.pca(adata, svd_solver=pca_solver.lower())
    pca_variance_ratio = adata.uns['pca']['variance_ratio'].cumsum()
    num_components = (pca_variance_ratio >= variance_threshold / 100).argmax() + 1

    adata.obsm['X_pca'] = adata.obsm['X_pca'][:, :num_components]
    adata.varm['PCs'] = adata.varm['PCs'][:, :num_components]
    adata.uns['pca']['variance'] = adata.uns['pca']['variance'][:num_components]
    adata.uns['pca']['variance_ratio'] = adata.uns['pca']['variance_ratio'][:num_components]

    st.write(f"PCA applied using the {pca_solver} solver and retaining {num_components} components.")

    # Plot cumulative explained variance
    fig, ax = plt.subplots()
    ax.plot(range(1, len(pca_variance_ratio) + 1), pca_variance_ratio * 100, marker='o', color='skyblue')
    ax.axhline(y=variance_threshold, color='r', linestyle='--', label=f'{variance_threshold}% Variance')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Cumulative Explained Variance (%)')
    ax.set_title('Cumulative Explained Variance by PCA Components')
    ax.set_xticks(range(1, len(pca_variance_ratio) + 1))
    ax.legend()
    st.pyplot(fig)

def run_umap_leiden(adata, resolutions, n_neighbors_list, min_dist, metric, output_dir, flavor):
    # Precompute neighbors for all resolutions and neighbors
    total_steps = len(resolutions) * len(n_neighbors_list)
    progress_bar = st.progress(0)
    current_step = 0

    silhouette_scores, ch_scores, db_scores, inertia_values, num_clusters = [], [], [], [], []

    for n_neighbors in n_neighbors_list:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, metric=metric, random_state=50)
        sc.tl.umap(adata, min_dist=min_dist, random_state=50)

        for resolution in resolutions:
            st.write(f"Computing Leiden clustering with n_neighbors={n_neighbors}, resolution={resolution}, flavor={flavor}")
            sc.tl.leiden(adata, resolution=resolution, random_state=50, flavor=flavor, n_iterations=2, directed=False)

            # Retrieve PCA embeddings
            X = adata.obsm['X_pca']
            cluster_labels = adata.obs['leiden'].astype(int)

            # Compute silhouette, Calinski-Harabasz, Davies-Bouldin scores
            silhouette_avg = silhouette_score(X, cluster_labels)
            silhouette_scores.append(silhouette_avg)

            ch_score = calinski_harabasz_score(X, cluster_labels)
            ch_scores.append(ch_score)

            db_score = davies_bouldin_score(X, cluster_labels)
            db_scores.append(db_score)

            # Compute cluster centroids and inertia
            unique_labels = np.unique(cluster_labels)
            centroids = np.array([X[cluster_labels == label].mean(axis=0) for label in unique_labels])
            inertia = np.sum(np.min(cdist(X, centroids, 'euclidean'), axis=1) ** 2)
            inertia_values.append(inertia)
            num_clusters.append(len(unique_labels))

            # Save UMAP plots as before
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.umap(adata, color='leiden', size = 10, ax=ax, show=False)
            umap_plot_path = os.path.join(output_dir, f'umap_n{n_neighbors}_r{resolution}.png')
            plt.savefig(umap_plot_path, dpi=300, bbox_inches='tight')
            st.pyplot(fig)
            plt.close(fig)

            # Save the processed AnnData object
            adata_output_path = os.path.join(output_dir, f'adata_n{n_neighbors}_r{resolution}.h5ad')
            adata.write(adata_output_path)
            st.write(f"Processed AnnData object saved as: {adata_output_path}")

            # Update progress bar
            current_step += 1
            progress_bar.progress(current_step / total_steps)

    # --- Step 4: Plot Evaluation Metrics ---
    st.write("Plotting evaluation metrics...")

    # Plot Silhouette Score
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(resolutions, silhouette_scores, marker='o')
    ax1.set_title('Silhouette Score vs. Leiden Resolution')
    ax1.set_xlabel('Leiden Resolution')
    ax1.set_ylabel('Silhouette Score')
    st.pyplot(fig1)

    # Plot Calinski-Harabasz Score
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    ax2.plot(resolutions, ch_scores, marker='o', color='green')
    ax2.set_title('Calinski-Harabasz Score vs. Leiden Resolution')
    ax2.set_xlabel('Leiden Resolution')
    ax2.set_ylabel('Calinski-Harabasz Score')
    st.pyplot(fig2)

    # Plot Davies-Bouldin Score
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    ax3.plot(resolutions, db_scores, marker='o', color='red')
    ax3.set_title('Davies-Bouldin Score vs. Leiden Resolution')
    ax3.set_xlabel('Leiden Resolution')
    ax3.set_ylabel('Davies-Bouldin Score')
    st.pyplot(fig3)

    # Plot Elbow Method: Inertia vs. Number of Clusters
    fig4, ax4 = plt.subplots(figsize=(10, 6))
    ax4.plot(num_clusters, inertia_values, marker='o')
    ax4.set_title('Elbow Method: Inertia vs. Number of Clusters')
    ax4.set_xlabel('Number of Clusters')
    ax4.set_ylabel('Sum of Squared Distances (Inertia)')
    st.pyplot(fig4)

    # Display clustering evaluation metrics
    result_df = pd.DataFrame({
        'Resolution': resolutions,
        'Silhouette Score': silhouette_scores,
        'Calinski-Harabasz Score': ch_scores,
        'Davies-Bouldin Score': db_scores,
        'Number of Clusters': num_clusters,
        'Inertia': inertia_values
    })
    st.write("**Clustering Evaluation Metrics for Each Resolution:**")
    st.dataframe(result_df)


if st.button('Run Analysis'):
    with st.spinner("Working..."):
        if not input_dir or not output_dir:
            st.error("Please specify both input and output directories.")
        else:
            adata = process_csv_files(input_dir, output_dir)

            if adata is not None:
                if pca_solver != 'none':
                    perform_pca(adata, pca_solver, variance_threshold, output_dir)

                # Prepare Leiden resolutions and neighbors
                leiden_resolutions_list = [float(res) for res in leiden_resolutions.split(',')]
                run_umap_leiden(adata, leiden_resolutions_list, [n_neighbors], min_dist, distance_metric, output_dir, leiden_flavor)