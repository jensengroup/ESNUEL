#!/bin/bash

# Parameters
#SBATCH --cpus-per-task=24
#SBATCH --error=/lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/calc_gfn1xtb_alpbDMSO/methyl_anion/submitit_orca/%j_0_log.err
#SBATCH --job-name=xtbopt
#SBATCH --mem=20GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --open-mode=append
#SBATCH --output=/lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/calc_gfn1xtb_alpbDMSO/methyl_anion/submitit_orca/%j_0_log.out
#SBATCH --partition=kemi1
#SBATCH --signal=USR1@90
#SBATCH --time=180
#SBATCH --wckey=submitit

# command
export SUBMITIT_EXECUTOR=slurm
srun --output /lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/calc_gfn1xtb_alpbDMSO/methyl_anion/submitit_orca/%j_%t_log.out --error /lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/calc_gfn1xtb_alpbDMSO/methyl_anion/submitit_orca/%j_%t_log.err --unbuffered /groups/kemi/ree/anaconda3/envs/uniRXNpred/bin/python -u -m submitit.core._submit /lustre/hpc/kemi/ree/PhD/applications/universalRXNpredictor/src/uniRXNpred/ref_calcs/calc_gfn1xtb_alpbDMSO/methyl_anion/submitit_orca