import os
import argparse
import numpy as np
import shutil
import subprocess

from rdkit import Chem
from ase.units import Hartree, mol, kcal, kJ

# xTB path and calc setup
path = os.path.dirname(os.path.realpath(__file__)).replace('/src/uniRXNpred/test', '')
XTBHOME = os.path.join(path, 'dep/xtb-6.5.1')
XTBPATH = os.path.join(path, 'dep/xtb-6.5.1/share/xtb')
MANPATH = os.path.join(path, 'dep/xtb-6.5.1/share/man')
LD_LIBRARY_PATH = os.path.join(path, 'dep/xtb-6.5.1/lib')

OMP_NUM_THREADS = '1'
MKL_NUM_THREADS = '1'

def parse_args():
    """
    Argument parser so this can be run from the command line
    """
    parser = argparse.ArgumentParser(description='Run orca input creator')
    parser.add_argument('-n', '--name', default='test_mol.xyz', help='The name of xyz file')
    parser.add_argument('-c', '--charge', default=0, help='The charge of the molecule')
    parser.add_argument('-s', '--spin', default=0, help='The spin of the molecule')
    parser.add_argument('-m', '--method', default=' 2', help='The method of choice')
    parser.add_argument('-a', '--solvent', default='--alpb DMSO', help='The sovlent of choice')
    return parser.parse_args()


def run_xTB(xyzfile, chrg=0, spin=0, method=' 2', solvent='--alpb DMSO', optimize=True):

    global XTBHOME
    global XTBPATH
    global LD_LIBRARY_PATH

    global OMP_NUM_THREADS
    global MKL_NUM_THREADS

    # Set env parameters for xTB
    os.environ['XTBHOME'] = XTBHOME
    os.environ['XTBPATH'] = XTBPATH
    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
    os.environ["OMP_NUM_THREADS"] = OMP_NUM_THREADS
    os.environ['MKL_NUM_THREADS'] = MKL_NUM_THREADS

    # Create calculation directory
    mol_calc_path = os.path.join(os.getcwd(), 'calc', xyzfile.split('.')[0])
    os.makedirs(mol_calc_path, exist_ok=True)
    
    # Copy file to calculation directory
    start_structure_xyz = os.path.join(os.getcwd(), 'calc', xyzfile.split('.')[0], xyzfile)
    shutil.copy(os.path.join(os.getcwd(), xyzfile), start_structure_xyz)

    # Run xTB calc
    if optimize:
        cmd = f'{XTBHOME}/bin/xtb --gfn{method} {start_structure_xyz} --opt --vfukui --lmo --chrg {chrg} --uhf {spin} {solvent}'
    else:
        cmd = f'{XTBHOME}/bin/xtb --gfn{method} {start_structure_xyz} --vfukui --lmo --chrg {chrg} --uhf {spin} {solvent}'

    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=mol_calc_path)
    output = proc.communicate()[0]

    # Save calc output
    with open(f'{mol_calc_path}/{xyzfile.split(".")[0]}.xtbout', 'w') as f:
        f.write(output)

    # Search for the molecular energy
    for i, line in enumerate(output.split('\n')):
        if 'TOTAL ENERGY' in line:
            energy = line.split()[3]

    try: #check if the molecular energy was found.
        energy = float(energy) * Hartree * mol/kJ #* mol/kJ convert energy from Hartree to kJ/mol  or  #* mol/kcal convert energy from Hartree to kcal/mol
    except Exception as e:
        print(e, xyzfile)
        energy = 60000.0

    print(energy)
    return energy


if __name__ == "__main__":

    import submitit
    args = parse_args()
    
    # Slurm settings
    executor = submitit.AutoExecutor(folder="submitit_xtb")
    executor.update_parameters(
        name=args.name,
        cpus_per_task=1,
        mem_gb=2,
        timeout_min=180,
        slurm_partition="kemi1",
        slurm_array_parallelism=1,
    )
    print(executor)

    jobs = []
    with executor.batch():
        job = executor.submit(run_xTB, str(args.name), int(args.charge), int(args.spin * 2), str(args.method), str(args.solvent), True)
        jobs.append(job)
