# This folder contains:
  1) "from_sources": The raw data from different sources, which have been used to design and evaluate ESNUEL.
  2) "results": The final results of ESNUEL along with the scripts used to generate the figures in our paper.


## Data in "from_sources":

| File | Reference | Accessed |
| ------------- | ------------- | ------------- |
| MAA_paper_fig4.smiles.csv | https://dx.doi.org/10.1021/acs.joc.9b03187 | Jun. 30, 2023 |
| MCA_paper_fig5.smiles.csv | https://dx.doi.org/10.1021/acs.joc.0c02327 | Jun. 30, 2023 |
| covalentinhibtors.csv | https://doi.org/10.1007/s10822-020-00342-w | Jun. 30, 2023 |
| covalentinhibtors_modifiedWarheadB.csv | https://doi.org/10.1007/s10822-020-00342-w | Jun. 30, 2023 |
| electrophilicity_with_names.csv | https://doi.org/10.1021/acs.jcim.1c01400 (https://cdb.ics.uci.edu/cgibin/ReactivitiesDatasetsWeb.html) | Jun. 30, 2023 |
| nucleophilicity_with_names.csv | https://doi.org/10.1021/acs.jcim.1c01400 (https://cdb.ics.uci.edu/cgibin/ReactivitiesDatasetsWeb.html) | Jun. 30, 2023 |
| 100_rxn_mechanisms.csv | https://doi.org/10.1021/acs.jcim.1c01400 (https://cdb.ics.uci.edu/cgibin/ReactivitiesDatasetsWeb.html) | Jun. 30, 2023 |

## Computed data in "results":

### Extract precomputed results from submitit_results.tar.gz

    tar -xvzf submitit_results.tar.gz

Then to generate the figures in our paper, please see the following:
 - Figure 1: analyse_results_for_figure1.ipynb
 - Figure 3: analyse_covalentinhibtors.ipynb
 - Figure S2: analyse_results_for_figure1_xtb.ipynb
