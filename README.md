<p align="center">
  <img src="image/logo.png"/>
</p>

---

ESNUEL (**ES**timating **NU**cleophilicity and **EL**ectrophilicity) is a fully automated quantum chemistry (QM)-based workflow that automatically identifies nucleophilic and electrophilic sites and computes methyl cation affinities (MCAs) and methyl anion affinities (MAAs) to estimate nucleophilicity and electrophilicity, respectively.

[TRY ESNUEL: https://www.esnuel.org](https://www.esnuel.org)

## Installation

For the installation, we recommend using `conda` to get all the necessary dependencies:

    conda env create -f environment.yml && conda activate esnuel


Then download the binaries of xtb version 6.5.1:

    mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.5.1/xtb-6.5.1-linux-x86_64.tar.xz; tar -xvf ./xtb-6.5.1-linux-x86_64.tar.xz; cd ..


Furthermore, ORCA version 5.0.1 must be installed following the instructions found here: https://sites.google.com/site/orcainputlibrary/setting-up-orca

OBS! 
  1) The path to ORCA must be modified in "src/esnuel/run_orca.py".
  2) The number of available CPUs and memory must be modified to match your hardware.


## Usage

An example of usage via CLI command:

    # Create predictions for a test molecule (OBS! Only names without "_" are allowed):
    python src/esnuel/calculator.py --smi 'Cn1c(C(C)(C)N)nc(C(=O)NCc2ccc(F)cc2)c(O)c1=O' --name 'testmol' &
    

The calculations are now saved in a "./calculations" folder along with a graphical output of the results (in .html format).
The graphical output presents the user with the most electrophilic and nucleophilic sites within 3 kcal/mol ≈ 12.6 kJ/mol being highlighted.

An example of using ESNUEL in batch mode:

    # Create predictions for a small dataset (example/testmols.csv):
    python src/esnuel/calculator.py -b example/testmols.csv

The calculations are now saved in a "./calculations" folder, and a dataframe containing the results is found in "submitit_results/testmol/*_result.pkl"


The SLURM commands can be modified via the following command line arguments:
  - '--partition': The SLURM partition that you submit to, default='kemi1'.
  - '--parallel_calcs': The number of parallel molecule calculations (the total number of CPU cores requested for each SLURM job = parallel_calcs*cpus_per_calc), default=2.
  - '--cpus_per_calc': The number of cpus per molecule calculation (the total number of CPU cores requested for each SLURM job = parallel_calcs*cpus_per_calc), default=4.
  - '--mem_gb': The total memory usage in gigabytes (gb) requested for each SLURM job, default=20.
  - '--timeout_min': The total allowed duration in minutes (min) of each SLURM job, default=6000.
  - '--slurm_array_parallelism': The maximum number of parallel SLURM jobs to run simultaneously (taking one molecule at a time in batch mode), default=25.


  For the QM calculations, the molecular charge is defined by the formal charge of the molecule using RDKit, and the spin is hardcoded to S=0 (multiplicity=1), as we focus on closed-shell molecules.
  This can be modified in the "calculateEnergy" function in "src/esnuel/calculator.py".


## Citation 

Our work is available as a preprint on [ChemRxiv](https://doi.org/10.26434/chemrxiv-2023-pd395), where more information is available. 
```
@article{ree2023esnuel,
  doi = {10.26434/chemrxiv-2023-pd395},
  url = {https://doi.org/10.26434/chemrxiv-2023-pd395},
  year = {2023},
  month = oct,
  author = {Nicolai Ree and Andreas H. G\"{o}ller and Jan H. Jensen},
  title = {Automated Quantum Chemistry for Estimating Nucleophilicity and Electrophilicity with Applications to Retrosynthesis and Covalent Inhibitors}
}
```
