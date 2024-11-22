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

# Define the result directory at the beginning
result_dir = os.path.join(os.getcwd(), 'result/statistics')
os.makedirs(result_dir, exist_ok=True)


st.title("Downsample your data")
st.write("*CAFE will downsample large datasets using PCA (Principal Component Analysis) and KDE (Kernel Density Estimation). KDE is applied to the PCA-transformed data to estimate the density of data points. Based on these density estimates, the code probabilistically downsamples the data, favoring points in regions of higher density. This method reduces sampling bias and preserves original data distribution. If your CSV files are very large, we recommend using an HPC cluster.*")


# Function to run PCA
def run_pca(data, n_components=5):
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data)
    explained_variance = pca.explained_variance_ratio_
    return pca_result

# Function to run FFT-based KDE
def kde_fft(pca_data):
    kde = gaussian_kde(pca_data.T, bw_method='silverman')  # FFT-based KDE with silverman bandwidth
    density = kde(pca_data.T)
    return density

# Function to downsample based on KDE density
def downsample_data(df, pca_data, density, n_samples=50000):
    # Normalize densities to probabilities
    probabilities = density / density.sum()

    # Select n_samples samples based on probabilities
    downsampled_idx = np.random.choice(df.index, size=n_samples, p=probabilities, replace=False)
    downsampled_df = df.loc[downsampled_idx]
    return downsampled_df

# Function to process each file: Read CSV, run PCA, KDE, downsample, and save
def process_file(file_path, output_dir, n_components, n_markers, n_samples):
    # Read the CSV file
    df = pd.read_csv(file_path)

    # Assuming marker data is in specific columns (based on user selection)
    marker_data = df.iloc[:, :n_markers]  # Adjust this based on your actual marker columns

    # Step 1: Run PCA
    pca_result = run_pca(marker_data, n_components)

    # Step 2: Run FFT-based KDE
    density = kde_fft(pca_result)

    # Step 3: Downsample the data
    downsampled_df = downsample_data(df, pca_result, density, n_samples=n_samples)

    # Step 4: Save the downsampled data with "downsampled_" prefix
    file_name = os.path.basename(file_path)
    output_file = os.path.join(output_dir, f'downsampled_{file_name}')
    downsampled_df.to_csv(output_file, index=False)

# Main processing loop (Parallel Processing): Read all CSV files and process them with elapsed time calculation
def process_directory_parallel(directory, output_dir, n_components, n_markers, n_samples, n_jobs):
    # Get a list of all CSV files in the directory
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv')]

    start_time = time.time()  # Record the start time

    # Run the process in parallel using Joblib with user-selected n_jobs
    Parallel(n_jobs=n_jobs)(delayed(process_file)(file, output_dir, n_components, n_markers, n_samples) for file in files)

    # Calculate the elapsed time
    elapsed_time = time.time() - start_time
    return elapsed_time

# Streamlit UI for input and settings
st.subheader('PCA and KDE-based Downsampling Tool')

# Input directory
input_directory = st.text_input('Enter the input directory:', '')

# Output directory
output_directory = st.text_input('Enter the output directory:', '')

# Slider for PCA n_components
n_components = st.slider('Select number of PCA components:', 2, 10, 5)

n_samples = st.number_input('Enter the number of samples to downsample:', min_value=100, max_value=100000, value=50000, step=100)


# Slider for number of markers
n_markers = st.slider('Select number of markers (columns):', 5, 50, 20)

# Slider for number of jobs (parallel processing)
n_jobs = st.slider('Select the number of jobs (parallel processing cores):', -1, 4, -1)

# Process button
if st.button('Apply'):
    with st.spinner("The program is running...."):
        if input_directory and output_directory:
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)

            # Call the parallel processing function and get the elapsed time
            elapsed_time = process_directory_parallel(input_directory, output_directory, n_components, n_markers, n_samples, n_jobs)

            # Display the elapsed time
            st.write(f'Processing complete. Elapsed time: {elapsed_time:.2f} seconds')
        else:
            st.error('Please provide both input and output directories.')
