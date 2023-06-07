# ReactivityQM
A quantum chemistry-based workflow that is completely automated and capable of identifying possible nucleophilic and electrophilic sites and calculating their corresponding MCA or MAA.


## Installation

For the installation, we recommend using `conda` to get all the necessary dependencies:

    conda env create -f environment.yml && conda activate ReactivityQM


Then download the binaries of xtb version 6.5.1:

    mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.5.1/xtb-6.5.1-linux-x86_64.tar.xz; tar -xvf ./xtb-6.5.1-linux-x86_64.tar.xz; cd ..


Furthermore, ORCA version 5.0.1 must be installed following the instructions found here: https://sites.google.com/site/orcainputlibrary/setting-up-orca

OBS! 
  1) The path to ORCA must be modified in "src/ReactivityQM/run_orca.py".
  2) The number of available CPUs and memory must be modified to match your hardware.


## Usage

An example of usage via CLI command (takes ~10 min. on 2x16 CPU cores):

    # Create predictions for an intermolecular Heck reaction:
    python src/ReactivityQM/calculator.py --smi 'C=CC' --name 'test' &
    

The calculations are now saved in a "./calculations" folder and a graphical output of the results (in .html format) are found in a "results" folder.
The graphical output presents the user with the most reactive electrophilic and nucleophilic sites highlighted.

An example of using ReactiviQM on a dataframe (takes ~20 min. on 16 CPU cores):

    # Create predictions for "3d" from the Cabri dataset (data/cabri/data_cabri.pkl):
    python heckQM.py

The calculations are now saved in a "./calculations" folder, and a dataframe containing the results are found in "submitit/*_result.pkl"


## Citation 

Our work is available as a preprint on [ChemRxiv](http://doi.org/10.26434/chemrxiv-2022-9pslv), where more information is available. 
```
@article{ree2023reactivityQM,
  doi = {10.26434/chemrxiv-2022-9pslv},
  url = {https://doi.org/10.26434/chemrxiv-2022-9pslv},
  year = {2023},
  month = jul,
  author = {Nicolai Ree and Andreas H. G\"{o}ller and Jan H. Jensen},
  title = {Quantifying Electrophilicity and Nucleophilicity with Automated Quantum Chemistry and its Applications to Retrosynthesis and Covalent Inhibitors}
}
```