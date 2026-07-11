import os
import time

import numpy as np
import pandas as pd
import streamlit as st
from joblib import Parallel, delayed
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA

from theme import IMG_DIR, apply_theme, page_header
from utils import validate_path

st.set_page_config(layout="centered")
apply_theme()


# Reset configuration when freshly navigating to the page
if st.session_state.get("_page_changed"):
    st.session_state["downsampling_configured"] = False

st.logo(os.path.join(IMG_DIR, "s_logo.png"))

page_header(
    "Downsample your data",
    subtitle="Reduce each file to a target cell count. Uniform preserves population frequencies; density-weighted enriches dominant populations; inverse-density keeps rare cells for clustering.",
)


def run_pca(data, n_components=5):
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data)
    explained_variance = pca.explained_variance_ratio_
    return pca_result


def kde_fft(pca_data):
    kde = gaussian_kde(pca_data.T, bw_method="silverman")
    density = kde(pca_data.T)
    return density


def downsample_uniform(df, n_samples):
    # Equal probability sampling
    downsampled_idx = np.random.choice(df.index, size=n_samples, replace=False)
    return df.loc[downsampled_idx]


def downsample_density_weighted(df, density, n_samples):
    # Density-proportional sampling
    probabilities = density / density.sum()
    downsampled_idx = np.random.choice(
        df.index, size=n_samples, p=probabilities, replace=False
    )
    return df.loc[downsampled_idx]


def downsample_inverse_density(df, density, n_samples):
    # Inverse-density sampling
    inv = 1.0 / np.clip(density, np.finfo(float).eps, None)
    probabilities = inv / inv.sum()
    downsampled_idx = np.random.choice(
        df.index, size=n_samples, p=probabilities, replace=False
    )
    return df.loc[downsampled_idx]


def process_file(file_path, output_dir, method, n_components, n_markers, n_samples):
    df = pd.read_csv(file_path)

    if len(df) <= n_samples:
        # Keep small files unchanged
        downsampled_df = df
        skipped = True
    elif method == "uniform":
        downsampled_df = downsample_uniform(df, n_samples)
        skipped = False
    else:  # PCA + KDE for density methods
        # Use numeric marker columns only, dropping known metadata, so leading
        # non-marker columns (e.g. SampleID/Group) can't corrupt the density estimate.
        marker_data = df.drop(columns=[c for c in ("SampleID", "Group") if c in df.columns])
        marker_data = marker_data.select_dtypes(include=[np.number]).iloc[:, :n_markers]
        pca_result = run_pca(marker_data, n_components)
        density = kde_fft(pca_result)
        if method == "density_weighted":
            downsampled_df = downsample_density_weighted(df, density, n_samples)
        else:  # inverse_density
            downsampled_df = downsample_inverse_density(df, density, n_samples)
        skipped = False

    file_name = os.path.basename(file_path)
    output_file = os.path.join(output_dir, f"downsampled_{file_name}")
    downsampled_df.to_csv(output_file, index=False)
    return file_name if skipped else None


def process_directory_parallel(
    directory, output_dir, method, n_components, n_markers, n_samples, n_jobs
):
    files = [
        os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".csv")
    ]

    start_time = time.time()

    results = Parallel(n_jobs=n_jobs)(
        delayed(process_file)(
            file, output_dir, method, n_components, n_markers, n_samples
        )
        for file in files
    )

    elapsed_time = time.time() - start_time
    skipped_files = [name for name in results if name is not None]
    return elapsed_time, skipped_files


st.subheader("Downsampling Tools")

# Step 1: Choose downsampling method
st.caption("Choose a method, then press Continue to configure it.")
method_label = st.radio(
    "Downsampling method",
    [
        "Uniform (preserve distribution)",
        "Density-weighted (enrich dominant populations)",
        "Inverse-density (enrich rare populations)",
    ],
    help="Uniform keeps relative frequencies intact and is required before frequency-based group "
    "comparisons. Density-weighted (the original CAFE method) enriches dominant populations and "
    "drops rare cells/outliers. Inverse-density boosts rare populations for clustering/visualization. "
    "Both density methods distort relative frequencies.",
)
method = {
    "Uniform (preserve distribution)": "uniform",
    "Density-weighted (enrich dominant populations)": "density_weighted",
    "Inverse-density (enrich rare populations)": "inverse_density",
}[method_label]


if not st.session_state.get("downsampling_configured"):
    if st.button("Continue", type="primary"):
        st.session_state["downsampling_configured"] = True
        st.rerun()

# Step 2: Configure sampling details
if st.session_state.get("downsampling_configured"):
    with st.expander("Data I/O", expanded=True):
        st.caption("Paths to the folder of CSV files to downsample, and where to write the outputs.")
        input_directory = st.text_input("Enter the input directory:", "")
        output_directory = st.text_input("Enter the output directory:", "")

    with st.expander("Sampling", expanded=True):
        st.caption("Target cell count per file. Files that already have fewer cells are copied unchanged.")
        n_samples = st.number_input(
            "Enter the number of samples to downsample:",
            min_value=100,
            max_value=100000,
            value=50000,
            step=100,
        )

    # Density-only controls
    if method in ("density_weighted", "inverse_density"):
        with st.expander("Density estimation", expanded=True):
            st.caption(
                "PCA reduces the marker space before KDE estimates local density. "
                "More components capture more biological detail but increase compute time. "
                "Set markers to match the number of marker columns in your CSV (excluding metadata)."
            )
            n_components = st.slider("Select number of PCA components:", 2, 10, 5)
            n_markers = st.slider("Select number of markers (columns):", 5, 50, 20)
    else:
        n_components, n_markers = 5, 20

    with st.expander("Parallelism", expanded=True):
        st.caption("Number of CPU cores to use for parallel file processing. –1 uses all available cores.")
        n_jobs = st.slider(
            "Select the number of jobs (parallel processing cores):", -1, 4, -1
        )

    if st.button("Apply", type="primary"):
        with st.spinner("The program is running...."):
            if input_directory and output_directory:
                # Validate paths
                try:
                    input_directory = validate_path(input_directory)
                    output_directory = validate_path(output_directory)
                except ValueError as e:
                    st.error(str(e))
                    st.stop()

                if not os.path.isdir(input_directory):
                    st.error(f"Input directory '{input_directory}' does not exist or is not a directory.")
                    st.stop()

                if not os.path.exists(output_directory):
                    os.makedirs(output_directory)

                elapsed_time, skipped_files = process_directory_parallel(
                    input_directory,
                    output_directory,
                    method,
                    n_components,
                    n_markers,
                    n_samples,
                    n_jobs,
                )

                st.write(
                    f"Processing complete. Elapsed time: {elapsed_time:.2f} seconds"
                )
                if skipped_files:
                    st.info(
                        f"{len(skipped_files)} file(s) had fewer than {n_samples} cells and were copied "
                        f"without downsampling: {', '.join(skipped_files)}"
                    )
            else:
                st.error("Please provide both input and output directories.")
