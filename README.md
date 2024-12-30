![Python Versions](https://img.shields.io/badge/python-3.12.5-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![bioRxiv](https://img.shields.io/badge/bioRxiv-https://doi.org/10.1101/2024.12.03.626714-red)

![Logo](CAFE/bin/img/logo.png)

# CAFE (Cell Analyzer for Flow Experiment)

CAFE is an open-source, free, no-code, web-app platform for high-dimensional spectral flow cytometry data (SFCM) analysis. CAFE has been developed in Python and it can run seamlessly on regular computers operating on either Windows or macOS/Linux. The application will allow the analysis and visualization of SFCM data to produce high-resolution publication-ready figures that support reproducible research of immune cell populations.

A lightweight demo version of CAFE is hosted here: [Website](https://tyrrell-lab.com)

Example data for testing available at FigShare: 
01. [Downsampled-CSV](https://figshare.com/articles/dataset/Downsampled_data_from_FlowRepository_FR-FCM-Z3WR/27940719)
02. [Processed-H5AD](https://figshare.com/articles/dataset/Downsampled_data_from_FlowRepository_FR-FCM-Z3WR/27940752)

##
## CAFE Interface

![Logo](CAFE/bin/img/CAFE_interface.jpg)

##
## How to Download CAFE

### Step 1: Download the tool as a ZIP file
1. Click the [Releases](https://github.com/mhbsiam/cafe/releases) button on the right side of this page.
2. Locate the latest release
3. Select **CAFE_version.zip** file to download.

### Step 2: Extract the ZIP File & Install the tool
1. Extract the files
2. Navigate to the " ./CAFE " folder path where you will find the following files: cafe.py, pixi.toml, cafe.yaml etc.
3. Then, follow the instructions below to run the tool using Pixi or Conda.


##
## Installation

### Method 1: Pixi Way

1. **Make sure you are in the correct directory (./CAFE/)**  
   If unsure, open your terminal app (or Powershell in Windows) and navigate to the folder containing the `cafe.py` file by typing:
   ```bash
   cd ./path/to/CAFE
   ```
   **Note:** Windows uses the backslash ( \ ) for the file system delimiter. For example, if the path to CAFE in Linux/macOS is `./documents/CAFE`, it will be `.\documents\CAFE` in Windows.


2. **Install Pixi**  
   Run the following command to install Pixi for Mac/Linux:
   ```bash
   curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v0.34.0 bash
   ```
   For Windows, (1) first run the following code:
   ```bash
   iwr -useb https://pixi.sh/install.ps1 | iex
   ```
   (2) Then run the code below to restore Pixi to a specific version:
   ```bash
   pixi self-update --version 0.34.0
   ```

3. **Run the CAFE tool**
   - **⚠️ Ensure you are in the correct directory (.CAFE/) where you will find the following files: bin folder, cafe.py, pixi.toml, cafe.yaml etc.**
   - Run the tool by typing:
   ```bash
   pixi run cafe
   ```


### Method 2: Conda Way

1. **Ensure Anaconda is installed**  
   If you don't have Anaconda installed, download it from [here](https://www.anaconda.com/download/success).

2. **Make sure you are in the correct directory (./CAFE/)**  
   If unsure, open your terminal app (or Anaconda Powershell) and navigate to the folder containing the `cafe.py` file by typing:
   ```bash
   cd ./path/to/CAFE
   ```

3. **Create the Conda environment**  
   Run the following command:
   ```bash
   conda env create -f cafe.yaml
   ```
   This will automatically create a new environment called `cafe` and install all required packages.

4. **Activate the Conda environment**  
   Type:
   ```bash
   conda activate cafe
   ```

5. **Run the tool**  
   Execute the tool by typing:
   ```bash
   python cafe.py
   ```

   **Note:** If you already have a `cafe` environment from a previous version, remove it first by deactivating the environment and running:
   ```bash
   conda env remove -n cafe -y
   ```
   Then, start from step 3 to set up a fresh environment.
##

## How to Run the App

![Logo](CAFE/bin/img/workflow.png)

### 01. Preparing the CSV Files

1.1 **Perform manual gating and export scaled CSV files**  
   - Open FlowJo or a similar software and upload your FCS files. Perform manual inspection of flow cytometry data and use manual gating to remove debris, dead cells, and doublets. Then, export the CSV files as scaled CSV files. You can also gate on appropriate cell type (e.g CD45+) and export the data to obtain more focused clustering results. 
   - These CSV files will be the input files for the CAFE tool.

1.2 **Rename each CSV file in the following way:**  
   - `SampleID_GroupA.csv` (e.g., `ABC001_Aged.csv`)
   - `SampleID_GroupB.csv` (e.g., `ABC002_Young.csv`)

   **⚡ Pro tip:** For optimal performance, CAFE expects two groups to run properly and perform statistics. Example groups: Aged and Young.

1.3 **Ensure that each CSV file:**
   - Contains the same number of columns (markers).
   - Uses identical marker names across files.
   - Example of consistent markers: CD45, TCRY, CD28, etc.

   **⚡ Pro tip:** CSV files with mismatched column names may result in errors. SampleID should be unique to avoid incorrect results.

**Scaled CSV file structure:**  
Here’s an example of how your CSV file should look like, with marker names as column headers and numerical values representing the expression levels:

| CD27     | CD38     | CD127    | HLA-DR  | CD1c     | CD141    | CD45ra   | CD16     | CCR5     | CD4      | CD11c    | CD56     |
|----------|----------|----------|---------|----------|----------|----------|----------|----------|----------|----------|----------|
| 27257.3 | 5442.59 | 13068.2 | -2432.2 | -3011.07 | 2048.14 | -516.212 | 10914.9 | 4355.31 | 87660.1 | 2256.65 | 1316.62 |
| 16917.6 | 7401.0 | 10038.9 | -1056.87 | -743.682 | -3049.75 | -542.594 | 3165.02 | 4260.19 | 90922.7 | 2309.36 | -686.368 |
| 47378.3 | 4412.47 | 19117.0 | -1147.41 | -654.677 | -5256.23 | 68264.9 | 17845.1 | 8269.94 | 97427.5 | 5624.15 | 2726.43 |
| 1769.79 | 54090.1 | 265.372 | 8136.79 | -2202.24 | 3329.33 | 49450.1 | 3431.71 | 1735.02 | 908.1 | 15106.8 | 6331.52 |
| 40472.1 | 2639.16 | 20865.8 | -2138.81 | -1856.0 | -3293.17 | 81219.1 | 10629.2 | 4249.25 | 91438.0 | 3298.7 | 960.926 |
| 487.425 | 23893.5 | 2089.89 | 6163.82 | 555.821 | 3007.88 | 69839.0 | 7129.61 | 1701.63 | 164.528 | 19372.5 | 9389.31 |
| 31111.4 | 2010.18 | 3647.74 | -1939.37 | 264.905 | 1963.75 | 12751.9 | 8712.62 | 4819.16 | 89793.9 | 5114.97 | 1778.75 |
| -249.461 | 32396.4 | 1188.07 | 18938.1 | 455.37 | -3065.3 | 17207.5 | 1354.81 | 1981.51 | -38.6853 | 9176.99 | 17000.4 |
| 18263.4 | 455.373 | 2590.04 | -1000.35 | -2214.69 | 2975.21 | 421.45 | 5342.13 | 2328.19 | 54503.7 | 2726.56 | 1764.66 |

Data accessible at FlowRepository: FR-FCM-Z3WR

### 02. Running CAFE

```python

# 2.0 Open terminal app (or Anaconda Powershell or Windows Powershell)

# 2.1 Navigate to the directory where your CAFE tool’s Python file (cafe.py) is located.
cd ./path/to/cafe

# 2.2 Run the tool with pixi
pixi run cafe

# or Run the tool with conda
python cafe.py
```
##
## FAQ

**01. Can I load new data without re-running the CAFE app from the beginning?**  
- Absolutely. As mentioned above, you can select **Clear Cache** from the top right setting menu, then either press **Rerun** or simply reload the web page.  
**⚡ Pro Tip:** Clear cache and reload the web page for a quicker restart!

**02. What should I do if I get an error while loading CSV files?**  
- Make sure all CSV files have the same number of markers (i.e., columns) and that the column names match across files.
- Check that your files are named correctly, e.g., `SampleID_GroupA.csv`.

**03. How long does it take for Leiden clustering?**  
- Depending on the user's specification, UMAP computation and Leiden clustering steps can take up to 30 minutes. For example, with an Apple M3 Pro 18GB system, we ran Leiden resolution of 1.0 with UMAP n_neighbors set to 15 using a dataset of 350,000 cells and 12.25M data points, it took ~12 minutes to complete.

##
## How to run Data Processing step in HPC?

This is optional. If you have access to a High-Performance Computing system, you can process your data in HPC using the provided cafe_hpc.py script to make an AnnData Object. You can also test multiple Leiden resolution values to determine which fits the dataset best. This is particularly useful if you have a large dataset.

### Steps:

- Download the standalone script found in CAFE/cafe_hpc.py
- Assign input and output directory
- Assign parameters. Note that some arguments accept multi-parameters such as Leiden resolution
- Run the app

### 1. Choose the Parameters
Decide on the arguments to pass to the script:
- **`--input`**: Path to the input directory containing the CSV files.
- **`--output`**: Path to the output directory where results will be saved.
- **`--pca`**: PCA solver to use (choose from `auto`, `full`, `arpack`, `randomized`, or `none`).
- **`--cutoff`**: Explained variance cutoff for PCA (e.g., 95).
- **`--leiden`**: One or more resolutions for Leiden clustering (e.g., `0.5 1.0 2.0`).
- **`--nneighbor`**: One or more `n_neighbors` values for UMAP and neighbors computation (e.g., `10 20 30`).
- **`--distance`**: Distance metric to use for neighbors computation (default: `euclidean`).
- **`--min_dist`**: `min_dist` parameter for UMAP (default: `0.1`).

### 2. Run the Script
Run the script in your terminal using the following command. For example:
```bash
python3 cafe_hpc.py \
  --input /path/to/input_dir \
  --output /path/to/output_dir \
  --pca auto \
  --cutoff 95 \
  --leiden 0.5 1.0 \
  --nneighbor 15 30 \
  --distance euclidean \
  --min_dist 0.1
```

## 
## Tool Dependencies

```
python = "==3.12.5"
numpy = "==1.26.4"
pandas = "==2.2.3"
scipy = "==1.14.1"
pyarrow = "==18.0.0"
statsmodels = "==0.14.4"
scanpy = "==1.10.3"
cliffs-delta = "==1.0.0"
matplotlib = "==3.9.2"
seaborn = "==0.13.2"
plotly = "==5.24.1"
cmasher = "==1.9.0"
igraph = "==0.10.8"
streamlit = "==1.39.0" -> "==1.41.0"
watchdog = "==5.0.2" -> "==6.0.0"
kaleido = "==0.2.1"
```

##
## System Requirements 

![Logo](CAFE/bin/img/os_selection.png)

## Known Issues

Windows OS has known issues with Scanpy and thus generates a higher number of clusters with a given Leiden resolution value. A numpy conflict may provide a value error message, ignore the warning. For optimal performance, consider using WSL to download the CAFE folder from Github using git clone and run the app using Pixi.
```bash
git clone https://github.com/mhbsiam/cafe.git
```
[Here](https://learn.microsoft.com/en-us/windows/wsl/install) is a tutorial on how to install and activate WSL on a Windows computer.
Finally, while reporting data generated using CAFE, mentioning the operating system used for data generation is highly recommended for replicability. 

## Citation

If you use CAFE in your research, please cite our paper:

Md Hasanul Banna Siam, Md Akkas Ali, Donald Vardaman, Mallikarjun Patil, Satwik Acharyya, and Daniel Tyrrell. “Cafe: An Integrated Web App for High-Dimensional Analysis and Visualization in Spectral Flow Cytometry,” bioRxiv (2024). [https://doi.org/10.1101/2024.12.03.626714](https://doi.org/10.1101/2024.12.03.626714).

Donald Vardaman, Md Akkas Ali, Md Hasanul Banna Siam, Chase Bolding, Harrison Tidwell, Holly R. Stephens, Mallikarjun Patil, and Daniel J. Tyrrell. "Development of a Spectral Flow Cytometry Analysis Pipeline for High-Dimensional Immune Cell Characterization." The Journal of Immunology (2024). [https://doi.org/10.4049/jimmunol.2400370](https://doi.org/10.4049/jimmunol.2400370).
