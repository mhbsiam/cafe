import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import gaussian_kde
from joblib import Parallel, delayed
import streamlit as st
import time  # Import the time module

image_path = os.path.join('bin', 'img', 's_logo.png')
st.logo(image_path)

# Display an image
image_path = os.path.join('bin', 'img', 'logo_v2.png')
st.image(image_path, caption='', use_column_width=True)

result_dir = os.path.join(os.getcwd(), 'result/statistics')
os.makedirs(result_dir, exist_ok=True)

st.title("Downsample your data")
st.write("*CAFE will downsample large datasets using PCA (Principal Component Analysis) and KDE (Kernel Density Estimation). KDE is applied to the PCA-transformed data to estimate the density of data points. Based on these density estimates, the code probabilistically downsamples the data, favoring points in regions of higher density. This method reduces sampling bias and preserves original data distribution. If your CSV files are very large, we recommend using an HPC cluster.*")

def run_pca(data, n_components=5):
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data)
    explained_variance = pca.explained_variance_ratio_
    return pca_result

def kde_fft(pca_data):
    kde = gaussian_kde(pca_data.T, bw_method='silverman')
    density = kde(pca_data.T)
    return density

def downsample_data(df, pca_data, density, n_samples=50000):
    probabilities = density / density.sum()

    downsampled_idx = np.random.choice(df.index, size=n_samples, p=probabilities, replace=False)
    downsampled_df = df.loc[downsampled_idx]
    return downsampled_df

def process_file(file_path, output_dir, n_components, n_markers, n_samples):
    df = pd.read_csv(file_path)

    marker_data = df.iloc[:, :n_markers]
    pca_result = run_pca(marker_data, n_components)
    density = kde_fft(pca_result)
    downsampled_df = downsample_data(df, pca_result, density, n_samples=n_samples)
    file_name = os.path.basename(file_path)
    output_file = os.path.join(output_dir, f'downsampled_{file_name}')
    downsampled_df.to_csv(output_file, index=False)

def process_directory_parallel(directory, output_dir, n_components, n_markers, n_samples, n_jobs):
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv')]

    start_time = time.time()

    Parallel(n_jobs=n_jobs)(delayed(process_file)(file, output_dir, n_components, n_markers, n_samples) for file in files)

    elapsed_time = time.time() - start_time
    return elapsed_time

st.subheader('PCA and KDE-based Downsampling Tool')

input_directory = st.text_input('Enter the input directory:', '')

output_directory = st.text_input('Enter the output directory:', '')

n_components = st.slider('Select number of PCA components:', 2, 10, 5)

n_samples = st.number_input('Enter the number of samples to downsample:', min_value=100, max_value=100000, value=50000, step=100)


n_markers = st.slider('Select number of markers (columns):', 5, 50, 20)

n_jobs = st.slider('Select the number of jobs (parallel processing cores):', -1, 4, -1)

if st.button('Apply'):
    with st.spinner("The program is running...."):
        if input_directory and output_directory:
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)

            elapsed_time = process_directory_parallel(input_directory, output_directory, n_components, n_markers, n_samples, n_jobs)

            st.write(f'Processing complete. Elapsed time: {elapsed_time:.2f} seconds')
        else:
            st.error('Please provide both input and output directories.')
