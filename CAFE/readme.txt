Installation:
-----------------------------------------------
Method 01: Pixi way:

01.	Download the tool, open terminal, and Head over to the folder where you have the tool’s python file (i.e. cafe.py):

cd ./path/to/cafe

02.	Install pixi using the following command

For mac/linux: 

curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v0.34.0 bash

For windows: 

iwr -useb https://pixi.sh/install.ps1 | iex
pixi self-update --version 0.34.0


03.	Type the following command to install cafe

pixi install

04.	Run cafe by typing:

pixi run cafe

------------------------------------------------

Method 02: Conda way: 

01. Assuming you have Anaconda installed.
If not, make a fresh install at: https://www.anaconda.com/download/success

02. Download the tool. Open termonal and head over to the folder where you have the tool's python file (i.e. cafe.py) by typing the following command:

cd ./path/to/cafe

03. Paste the following command:

conda env create -f conda.yaml

[This will create a new environment called 'cafe' and install all required packages automatically]

04. Activate the conda environment by typing:

conda activate cafe


05. Run the tool by typing:

python cafe.py


[Note: If you already have a cafe environment for a previous version, remove the environment first. To remove Conda environment, make sure you have deactivated the environment that you are about to delete. Then, type:

conda env remove -n cafe -y 

Having removed the old environment, Start from step 03.
Activate the environment. Run it using step 05.]




