#!/usr/bin/env python3

import argparse
import os
import time
import gc
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import pyarrow.csv as pv
import pyarrow as pa
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



# Set Scanpy settings
sc.settings.n_jobs = -1

# Argument parser
parser = argparse.ArgumentParser(description='Process spectral flow cytometry data.')

parser.add_argument('--input', required=True, help='Directory path for CSV input files')
parser.add_argument('--output', required=True, help='Directory path for output files')
parser.add_argument('--pca', choices=['auto', 'full', 'arpack', 'randomized', 'none'], default='none', help='PCA solver to use')
parser.add_argument('--cutoff', type=int, default=95, help='Explained variance cutoff for PCA')
parser.add_argument('--leiden', type=float, nargs='+', required=True, help='Leiden resolution(s)')
parser.add_argument('--nneighbor', type=int, nargs='+', required=True, help='n_neighbors parameter(s) for UMAP and neighbors computation')
parser.add_argument('--distance', type=str, default='euclidean', help='Distance metric to use')
parser.add_argument('--min_dist', type=float, default=0.1, help='UMAP min_dist parameter')

args = parser.parse_args()

def process_csv_files(input_dir, output_dir):
    start_time = time.time()
    tables = []

    # Get list of CSV files
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    if not csv_files:
        print("No CSV files found in the input directory.")
        return None

    # Process each CSV file
    for idx, csv_file in enumerate(csv_files):
        print(f"Processing file {idx + 1}/{len(csv_files)}: {csv_file}")
        file_path = os.path.join(input_dir, csv_file)
        table = pv.read_csv(file_path)

        # Extract SampleID and Group from filename
        name_parts = csv_file.replace('.csv', '').split('_')
        if len(name_parts) != 2:
            print(f"Skipping file {csv_file} as it does not have the expected format 'sampleID_group.csv'")
            continue

        # Add SampleID and Group columns
        table = table.append_column('SampleID', pa.array([name_parts[0]] * len(table)))
        table = table.append_column('Group', pa.array([name_parts[1]] * len(table)))
        tables.append(table)

    # Combine tables
    if not tables:
        print("No valid CSV files were processed.")
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

    # Save combined data
    combined_output_path = os.path.join(output_dir, 'combined_updated_files.csv')
    df.to_csv(combined_output_path, index=False)
    print(f"Combined DataFrame saved to {combined_output_path}")

    print(f"Data loading and processing completed in {time.time() - start_time:.2f} seconds")
    return adata

def perform_pca(adata, pca_solver, variance_threshold, output_dir):
    print("Running PCA...")
    start_time = time.time()
    sc.tl.pca(adata, svd_solver=pca_solver.lower())
    pca_variance_ratio = adata.uns['pca']['variance_ratio'].cumsum()
    num_components = (pca_variance_ratio >= variance_threshold / 100).argmax() + 1

    # Subset the PCA results
    adata.obsm['X_pca'] = adata.obsm['X_pca'][:, :num_components]
    adata.varm['PCs'] = adata.varm['PCs'][:, :num_components]
    adata.uns['pca']['variance'] = adata.uns['pca']['variance'][:num_components]
    adata.uns['pca']['variance_ratio'] = adata.uns['pca']['variance_ratio'][:num_components]

    print(f"PCA applied using the {pca_solver} solver and retaining {num_components} components.")

    # Save PCA variance ratios
    variance_ratio_path = os.path.join(output_dir, 'pca_variance_ratio.csv')
    pd.Series(adata.uns['pca']['variance_ratio']).to_csv(variance_ratio_path, index=False)
    print(f"PCA variance ratios saved to {variance_ratio_path}")

    print(f"PCA completed in {time.time() - start_time:.2f} seconds")

def batch_correction(adata):
    if 'Group' in adata.obs.columns:
        start_time = time.time()
        sc.pp.combat(adata, key='Group')
        print(f"Batch correction completed in {time.time() - start_time:.2f} seconds")
    else:
        print("The 'Group' column does not exist in adata.obs.")

def run_umap_leiden(adata, resolutions, n_neighbors_list, min_dist, metric, output_dir):
    for n_neighbors in n_neighbors_list:
        print(f"\nComputing neighbors with n_neighbors={n_neighbors}, metric={metric}")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, method='umap', metric=metric)
        sc.tl.umap(adata, min_dist=min_dist)
        print(f"UMAP computed with min_dist={min_dist}")

        for resolution in resolutions:
            print(f"\nComputing Leiden clustering with resolution={resolution}")
            sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)
            adata.obs['leiden'] = adata.obs['leiden'].astype('category')

            # Save AnnData object
            output_h5ad_path = os.path.join(output_dir, f'adata_n{n_neighbors}_r{resolution}.h5ad')
            adata.write(output_h5ad_path)
            print(f"Anndata object saved as: {output_h5ad_path}")

            # Generate and save UMAP plot
            umap_plot_path = os.path.join(output_dir, f'umap_n{n_neighbors}_r{resolution}.png')
            print(f"Generating UMAP plot for n_neighbors={n_neighbors}, resolution={resolution}")

            # Adjust figure size and dot size as needed
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.umap(adata, color='leiden', legend_fontsize=8, size=10,
                        show=False, title=f'UMAP n_neighbors={n_neighbors}, resolution={resolution}',
                        palette='tab20', ax=ax, legend_loc='on data')
            ax.axis('off')

            # Save the figure
            plt.savefig(umap_plot_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"UMAP plot saved as: {umap_plot_path}")

            # Compute and save cluster sample counts
            cluster_sample_counts = adata.obs.groupby(['leiden', 'SampleID']).size().unstack(fill_value=0)
            output_cluster_csv = os.path.join(output_dir, f'cluster_counts_n{n_neighbors}_r{resolution}.csv')
            cluster_sample_counts.to_csv(output_cluster_csv)
            print(f"Cluster sample counts saved as: {output_cluster_csv}")

            # Compute and save cluster sample frequencies (%)
            sample_totals = cluster_sample_counts.sum()
            cluster_sample_percentages = cluster_sample_counts.div(sample_totals) * 100
            output_cluster_freq_csv = os.path.join(output_dir, f'cluster_freq_n{n_neighbors}_r{resolution}.csv')
            cluster_sample_percentages.to_csv(output_cluster_freq_csv)
            print(f"Cluster sample frequencies saved as: {output_cluster_freq_csv}")

            # Compute and save median fluorescence intensities
            median_df = adata.to_df()
            median_df['leiden'] = adata.obs['leiden']
            median_df['SampleID'] = adata.obs['SampleID']
            medians = median_df.groupby(['leiden', 'SampleID']).median()
            output_median_csv = os.path.join(output_dir, f'median_expression_n{n_neighbors}_r{resolution}.csv')
            medians.to_csv(output_median_csv)
            print(f"Median expression saved as: {output_median_csv}")

            # Clean up
            del cluster_sample_counts, sample_totals, cluster_sample_percentages, median_df, medians
            gc.collect()

def main():
    # Ensure output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Process CSV files
    adata = process_csv_files(args.input, args.output)
    if adata is None:
        print("No data to process. Exiting.")
        return

    # Perform PCA if requested
    if args.pca.lower() != 'none':
        perform_pca(adata, args.pca, args.cutoff, args.output)

    # Batch correction
    batch_correction(adata)

    # Run UMAP and Leiden clustering
    run_umap_leiden(
        adata,
        resolutions=args.leiden,
        n_neighbors_list=args.nneighbor,
        min_dist=args.min_dist,
        metric=args.distance,
        output_dir=args.output
    )

if __name__ == '__main__':
    main()
