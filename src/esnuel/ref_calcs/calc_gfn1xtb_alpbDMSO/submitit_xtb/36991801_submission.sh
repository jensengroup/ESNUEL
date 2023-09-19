#!/bin/bash

# Parameters
#SBATCH --cpus-per-task=1
#SBATCH --error=/lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/submitit_xtb/%j_0_log.err
#SBATCH --job-name=methyl_cation_xyz
#SBATCH --mem=2GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --open-mode=append
#SBATCH --output=/lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/submitit_xtb/%j_0_log.out
#SBATCH --partition=kemi1
#SBATCH --signal=USR1@90
#SBATCH --time=180
#SBATCH --wckey=submitit

# command
export SUBMITIT_EXECUTOR=slurm
srun --output /lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/submitit_xtb/%j_%t_log.out --error /lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/submitit_xtb/%j_%t_log.err --unbuffered /groups/kemi/ree/anaconda3/envs/uniRXNpred/bin/python -u -m submitit.core._submit /lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/submitit_xtb