![Python Versions](https://img.shields.io/badge/python-3.12.5-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![bioRxiv](https://img.shields.io/badge/bioRxiv-https://doi.org/10.1101/2024.12.03.626714-red)

![Logo](CAFE/bin/img/logo.png)

# CAFE (Cell Analyzer for Flow Experiment)

CAFE is an open-source, free, no-code, web-app platform for high-dimensional spectral flow cytometry data (SFCM) analysis. CAFE has been developed in Python and it can run seamlessly on regular computers operating on either Windows or macOS/Linux. The application will allow the analysis and visualization of SFCM data to produce high-resolution publication-ready figures that support reproducible research of immune cell populations.

Example data for testing available at FigShare: 
01. [Downsampled-CSV](https://figshare.com/articles/dataset/Downsampled_data_from_FlowRepository_FR-FCM-Z3WR/27940719)
02. [Processed-H5AD](https://figshare.com/articles/dataset/Downsampled_data_from_FlowRepository_FR-FCM-Z3WR/27940752)

**Update: CAFE is now accepted for publication in Bioinformatics**

Md Hasanul Banna Siam, Md Akkas Ali, Donald Vardaman, Satwik Acharyya, Mallikarjun Patil, Daniel J Tyrrell, CAFE: An Integrated Web App for High-Dimensional Analysis and Visualization in Spectral Flow Cytometry, Bioinformatics, 2025;, btaf176, https://doi.org/10.1093/bioinformatics/btaf176

##
## Documentation
Click [here](https://mhbsiam.github.io/cafe) to access the guide for installation and running CAFE 


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
## Quick Installation

### Using Pixi Package Manager

**For Mac/Linux**  
Run the following command to install Pixi for Mac/Linux:
```bash
curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v0.34.0 bash
```
**For Windows**,  
(1) first run the following code:
```bash
iwr -useb https://pixi.sh/install.ps1 | iex
```
(2) Then run the code below to restore Pixi to a specific version:
```bash
pixi self-update --version 0.34.0
```

##

## Run CAFE

```python

# Open terminal app (or Windows Powershell)

# Navigate to the directory where your CAFE tool’s file (pixi.toml) is present.
cd ./path/to/cafe

# Run the tool with pixi
pixi run cafe

```
##
## Workflow

![Logo](CAFE/bin/img/workflow.png)

## Known Issues

Windows OS has known issues with Scanpy and thus generates a higher number of clusters with a given Leiden resolution value. A numpy conflict may provide a value error message, ignore the warning. For optimal performance, consider using WSL to download the CAFE folder from Github using git clone and run the app using Pixi.
```bash
git clone https://github.com/mhbsiam/cafe.git
```
[Here](https://learn.microsoft.com/en-us/windows/wsl/install) is a tutorial on how to install and activate WSL on a Windows computer.
Finally, while reporting data generated using CAFE, mentioning the operating system used for data generation is highly recommended for replicability. 

## Citation

If you use CAFE in your research, please cite our paper:

Md Hasanul Banna Siam, Md Akkas Ali, Donald Vardaman, Satwik Acharyya, Mallikarjun Patil, Daniel J Tyrrell, CAFE: An Integrated Web App for High-Dimensional Analysis and Visualization in Spectral Flow Cytometry, Bioinformatics, 2025;, btaf176, https://doi.org/10.1093/bioinformatics/btaf176

Donald Vardaman, Md Akkas Ali, Md Hasanul Banna Siam, Chase Bolding, Harrison Tidwell, Holly R. Stephens, Mallikarjun Patil, and Daniel J. Tyrrell. "Development of a Spectral Flow Cytometry Analysis Pipeline for High-Dimensional Immune Cell Characterization." The Journal of Immunology (2024). [https://doi.org/10.4049/jimmunol.2400370](https://doi.org/10.4049/jimmunol.2400370).
