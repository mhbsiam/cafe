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

sc.settings.n_jobs = -1

# Clear memory
def clear_memory():
    gc.collect()

# Display an image
image_path = os.path.join('bin', 'img', 'logo.png')
st.image(image_path, caption='', use_column_width=True)
image_path2 = os.path.join('bin', 'img', 's_logo.png')
st.logo(image_path2)

# Define the result directory at the beginning
result_dir = os.path.join(os.getcwd(), 'result')
os.makedirs(result_dir, exist_ok=True)

# Title and subtitle with the CAFE acronym
st.title("Cell Analyzer for Flow Experiment")
st.subheader("Explore and analyze flow cytometry data interactively")

st.markdown(
    "This tool will analyze immune cells from spectral flow cytometry and generate publishable figures. For detailed documentation, please follow the tutorial.")
st.markdown("---")

# Documentation and citation section
st.markdown("## Citation")
st.markdown("""
For detailed documentation on how to use CAFE, including tutorials and examples, please follow the tutorial.

If you use CAFE in your research, please cite our paper:

> Vardaman, Donald, Md Akkas Ali, Md Hasanul Banna Siam, Chase Bolding, Harrison Tidwell, Holly R. Stephens, Mallikarjun Patil, and Daniel J. Tyrrell. "Development of a Spectral Flow Cytometry Analysis Pipeline for High-Dimensional Immune Cell Characterization." The Journal of Immunology (2024) https://doi.org/10.4049/jimmunol.2400370.
""")

st.markdown("---")

image_path = os.path.join('bin', 'img', 'uab.png')
st.image(image_path, caption='', use_column_width=True)

# Display an image
image_path = os.path.join('bin', 'img', 'funding.png')
st.image(image_path, caption='', use_column_width=True)

st.markdown("---")

st.subheader("How to Run the App: Step-by-Step Guide")

image_path = os.path.join('bin', 'img', 'da_cafe.png')
st.image(image_path, caption='', use_column_width=True)

st.subheader("01. Preparing the CSV files")
st.markdown("01. Perform manual gating and export scaled-csv-files from FlowJo or similar software. These CSV files will be the input files for CAFE tool.")
st.markdown("02. For CAFE to work, each CSV file needs to be renamed in the following way:")
st.markdown("SampleID_GroupA.csv, **Example: ABC001_Aged.csv**")
st.markdown("SampleID_GroupB.csv, **Example: ABC002_Young.csv**")
st.markdown("**⚡Pro tip:** For optimal performance, CAFE expects two groups to run properly and perform statistics. Example groups: Aged and Young")
st.markdown("03. Ensure that each CSV file:")
st.markdown("- Contains the same number of columns (markers)")
st.markdown("- Uses identical marker names across files")
st.markdown("- Uses unique sample names across files")
st.markdown("- Example of consistent markers: CD45, TCRY, CD28, etc")
st.markdown("- CSV files with mismatched column names may result in errors.")

st.markdown("**CSV file structure:**")
st.markdown("Here’s an example of how your CSV file should look like, with marker names as column headers and numerical values representing the expression levels:")

import pandas as pd

# Data in dictionary format
data = {
    "CD27": [27257.3, 16917.6, 47378.3, 1769.79, 40472.1, 487.425, 31111.4, -249.461, 18263.4],
    "CD38": [5442.59, 7401.0, 4412.47, 54090.1, 2639.16, 23893.5, 2010.18, 32396.4, 455.373],
    "CD127": [13068.2, 10038.9, 19117.0, 265.372, 20865.8, 2089.89, 3647.74, 1188.07, 2590.04],
    "HLA-DR": [-2432.2, -1056.87, -1147.41, 8136.79, -2138.81, 6163.82, -1939.37, 18938.1, -1000.35],
    "CD1c": [-3011.07, -743.682, -654.677, -2202.24, -1856.0, 555.821, 264.905, 455.37, -2214.69],
    "CD141": [2048.14, -3049.75, -5256.23, 3329.33, -3293.17, 3007.88, 1963.75, -3065.3, 2975.21],
    "CD45ra": [-516.212, -542.594, 68264.9, 49450.1, 81219.1, 69839.0, 12751.9, 17207.5, 421.45],
    "CD16": [10914.9, 3165.02, 17845.1, 3431.71, 10629.2, 7129.61, 8712.62, 1354.81, 5342.13],
    "CCR5": [4355.31, 4260.19, 8269.94, 1735.02, 4249.25, 1701.63, 4819.16, 1981.51, 2328.19],
    "CD4": [87660.1, 90922.7, 97427.5, 908.1, 91438.0, 164.528, 89793.9, -38.6853, 54503.7],
    "CD11c": [2256.65, 2309.36, 5624.15, 15106.8, 3298.7, 19372.5, 5114.97, 9176.99, 2726.56],
    "CD56": [1316.62, -686.368, 2726.43, 6331.52, 960.926, 9389.31, 1778.75, 17000.4, 1764.66],
}

# Create a DataFrame
df = pd.DataFrame(data)
st.dataframe(df)

st.subheader("02. Running CAFE")
st.code('''

# 1. Navigate to the tool’s directory
"Tip: Navigate to the directory where your CAFE tool’s Python file (cafe.py) is located."
cd ./path/to/cafe

# 2. Run the tool
python cafe.py (if you have installed via conda)
pixi run cafe (if you have installed via pixi)

''', language="python")

st.markdown("**⚡Pro tip: If you want to restart the program, a quicker approach is to clear cache from the top right settings menu, and reload the page.**")

st.subheader("03. Load files")
st.markdown("You have the option to select the CSV files, or if you have already generated AnnData containing Leiden Clusters, you can load the AnnData.")
st.markdown("If you need to increase the upload file size limit, edit cafe.py file and save.")
st.markdown("Example:")
st.code('''cli.main_run([cafe_path, "--server.maxUploadSize", "5000"])''')
st.markdown("This will increase the upload size to 5GB")
st.markdown("**⚡Pro tip: Avoid increasing server upload size unless you need to since this will cause higher memory usage.**")

st.subheader("04. Follow on-screen options")
st.markdown("Use Buttons, Sliders and Drop down menus to navigate through the app")
