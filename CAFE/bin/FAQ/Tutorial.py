import streamlit as st
import os

image_path = os.path.join('bin', 'img', 'small_logo.png')
st.logo(image_path)

st.subheader("How to Run the App: Step-by-Step Guide")

st.subheader("01. Preparing the CSV files")
st.markdown("01. Perform manual gating and export scaled-csv-files from FlowJo or similar software. These CSV files will be the input files for CAFE tool.")
st.markdown("02. For CAFE to work, each CSV file needs to be renamed in the following way:")
st.markdown("SampleID_GroupA.csv, **Example: ABC001_Aged.csv**")
st.markdown("SampleID_GroupB.csv, **Example: ABC002_Young.csv**")
st.markdown("**⚡Pro tip:** For optimal performance, CAFE expects two groups to run properly and perform statistics. Example groups: Aged and Young")
st.markdown("03. Ensure that each CSV file:")
st.markdown("- Contains the same number of columns (markers)")
st.markdown("- Uses identical marker names across files")
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
python cafe.py

''', language="python")

st.markdown("**⚡Pro tip: If you want to restart the program, a quicker approach is to clear cache from the top right settings menu, select rerun and select Reset.**")
st.markdown("**⚡Or, you can clear cache and reload the web page.**")

st.subheader("03. Load files")
st.markdown("You have the option to select the CSV files, or if you have already generated AnnData containing Leiden Clusters, you can load the AnnData.")
st.markdown("If you need to increase the upload file size limit, add --server.maxUploadSize flag while starting the app")
st.markdown("Example:")
st.code('''streamlit run cafe.py --server.maxUploadSize 2000 ''')
st.markdown("This will increase the upload size to 2TB")
st.markdown("**⚡Pro tip: Avoid increasing server upload size unless you need to since this will cause higher memory usage.**")

st.subheader("04. Follow on-screen options")
st.markdown("Use Buttons, Sliders and Drop down menus to navigate through the app")
