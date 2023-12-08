# This folder contains:
  1) "from_sources": The raw data from different sources used to design and evaluate ESNUEL.
  2) "results": The final results of ESNUEL along with the scripts used to generate the figures in our paper.


## Data in "from_sources":

| File | Description | Reference | Accessed |
| ------------- | ------------- | ------------- | ------------- |
| MCA_paper_fig5.smiles.csv | Used in Fig. 1a | https://dx.doi.org/10.1021/acs.joc.0c02327 (extracted from Tab called "Fig 5" in this [file](https://pubs.acs.org/doi/suppl/10.1021/acs.joc.0c02327/suppl_file/jo0c02327_si_002.xlsx)) | June 2023 |
| nucleophilicity_with_names.csv | Used in Fig. 1b | https://doi.org/10.1021/acs.jcim.1c01400 ([https://cdb.ics.uci.edu/cgibin/ReactivitiesDatasetsWeb.html](https://cdb.ics.uci.edu/cgibin/dataset-reactivities/nucleophilicity.csv)) | June 2023 |
| MAA_paper_fig4.smiles.csv | Used in Fig. 1c | https://dx.doi.org/10.1021/acs.joc.9b03187 (extracted from Table S4 in the [SI](https://pubs.acs.org/doi/suppl/10.1021/acs.joc.9b03187/suppl_file/jo9b03187_si_001.pdf)) | June 2023 |
| electrophilicity_with_names.csv | Used in Fig. 1d | https://doi.org/10.1021/acs.jcim.1c01400 ([https://cdb.ics.uci.edu/cgibin/ReactivitiesDatasetsWeb.html](https://cdb.ics.uci.edu/cgibin/dataset-reactivities/eletrophilicity.csv)) | June 2023 |
| covalentinhibtors.csv | Used in Fig. 3 | https://doi.org/10.1007/s10822-020-00342-w ([link to file](https://static-content.springer.com/esm/art%3A10.1007%2Fs10822-020-00342-w/MediaObjects/10822_2020_342_MOESM2_ESM.csv)) | June 2023 |
| covalentinhibtors_modifiedWarheadB.csv | Used in Fig. S12 | https://doi.org/10.1007/s10822-020-00342-w (a modification of [link to file](https://static-content.springer.com/esm/art%3A10.1007%2Fs10822-020-00342-w/MediaObjects/10822_2020_342_MOESM2_ESM.csv)) | June 2023 |
| 100_rxn_mechanisms.csv |  | https://doi.org/10.1021/acs.jcim.1c01400 ([https://cdb.ics.uci.edu/cgibin/ReactivitiesDatasetsWeb.html](https://cdb.ics.uci.edu/cgibin/dataset-reactivities/100_rxn_mechanisms.csv)) | June 2023 |

## Computed data in "results":

### Extract precomputed results from submitit_results.tar.gz

    tar -xvzf submitit_results.tar.gz

Then to generate the figures in our paper, please see the following:
 - Figure 1: analyse_results_for_figure1.ipynb
 - Figure 3: analyse_covalentinhibtors.ipynb
 - Figure S2: analyse_results_for_figure1_xtb.ipynb
