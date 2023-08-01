#!/bin/bash

# Parameters
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=8
#SBATCH --error=/lustre/hpc/kemi/ree/PhD/git/ReactivityQM/submitit_results/testmol/%A_%a_0_log.err
#SBATCH --job-name=reactivityQM
#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --open-mode=append
#SBATCH --output=/lustre/hpc/kemi/ree/PhD/git/ReactivityQM/submitit_results/testmol/%A_%a_0_log.out
#SBATCH --partition=kemi1
#SBATCH --signal=USR1@90
#SBATCH --time=6000
#SBATCH --wckey=submitit

# command
export SUBMITIT_EXECUTOR=slurm
srun --output /lustre/hpc/kemi/ree/PhD/git/ReactivityQM/submitit_results/testmol/%A_%a_%t_log.out --error /lustre/hpc/kemi/ree/PhD/git/ReactivityQM/submitit_results/testmol/%A_%a_%t_log.err --unbuffered /groups/kemi/ree/anaconda3/envs/uniRXNpred/bin/python -u -m submitit.core._submit /lustre/hpc/kemi/ree/PhD/git/ReactivityQM/submitit_results/testmol