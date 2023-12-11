# MIT License
#
# Copyright (c) 2023 Nicolai Ree
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import sys
import copy
import argparse
import numpy as np
from itertools import islice
from operator import itemgetter
from concurrent.futures import ThreadPoolExecutor

from rdkit import Chem
from rdkit.Chem import AllChem
# from rdkit import RDLogger
# lg = RDLogger.logger()
# lg.setLevel(RDLogger.CRITICAL)

execution_path = os.getcwd()
base_dir = os.path.dirname(os.path.realpath(__file__)).replace('/src/esnuel', '')

os.chdir(os.path.join(base_dir, 'src/esnuel')) # change path to make the calculation run
from locate_atom_sites import find_electrophilic_sites_and_generate_MAAproducts, find_nucleophilic_sites_and_generate_MCAproducts
from molecule_drawer import generate_structure, generate_output_tables, html_output
import molecule_formats as molfmt
import run_xTB as run_xTB
import run_orca as run_orca


def parse_args():
    """
    Argument parser so this can be run from the command line
    """
    parser = argparse.ArgumentParser(description='Run ESNUEL from the command line')
    parser.add_argument('-s', '--smiles', default='C[C+:20](C)CC(C)(C)C1=C(C=CC(=C1)Br)[OH:10]',
                        help='SMILES input to ESNUEL (only used, if you do not run in batch mode)')
    parser.add_argument('-n', '--name', default='testmol', help='The name of the molecule or job depending on CLI or batch mode. Only names without "_" are allowed.')
    parser.add_argument('-b', '--batch', default=None, help='Path to .csv file for running batched calculations e.g. --> python src/esnuel/calculator.py -b example/testmols.csv')
    parser.add_argument('--partition', default='kemi1', help='The SLURM partition that you submit to')
    parser.add_argument('--parallel_calcs', default=2, help='The number of parallel molecule calculations (the total number of CPU cores requested for each SLURM job = parallel_calcs*cpus_per_calc)')
    parser.add_argument('--cpus_per_calc', default=4, help='The number of cpus per molecule calculation (the total number of CPU cores requested for each SLURM job = parallel_calcs*cpus_per_calc)')
    parser.add_argument('--mem_gb', default=20, help='The total memory usage in gigabyte (gb) requested for each SLURM job')
    parser.add_argument('--timeout_min', default=6000, help='The total allowed duration in minutes (min) of each SLURM job')
    parser.add_argument('--slurm_array_parallelism', default=25, help='The maximum number of parallel SLURM jobs to run simultaneously (taking one molecule at a time in batch mode)')
    return parser.parse_args()


def confsearch_xTB(conf_complex_mols, conf_names, chrg=0, spin=0, method='ff', solvent='', conf_cutoff=10, precalc_path=None):

    global cpus_per_calc
    
    # Run xTB optimizations  
    confsearch_args = []
    for i in range(len(conf_names)):
        if precalc_path:
            confsearch_args.append((conf_names[i]+'_gfn'+method.replace(' ', '')+'.xyz', conf_complex_mols[i], chrg, spin, method, solvent, True, precalc_path[i]))
        else:
            confsearch_args.append((conf_names[i]+'_gfn'+method.replace(' ', '')+'.xyz', conf_complex_mols[i], chrg, spin, method, solvent, True, None))

    with ThreadPoolExecutor(max_workers=cpus_per_calc) as executor:
        results = executor.map(run_xTB.run_xTB, confsearch_args)

    conf_energies = []
    conf_paths = []
    calc_logs = []
    for result in results:
        conf_energy, path_opt, calc_log = result
        conf_energies.append(conf_energy)
        conf_paths.append(path_opt)
        calc_logs.append(calc_log)

    # Find the conformers below cutoff in kJ/mol (3 kcal/mol = 12.6 kJ/mol)
    rel_conf_energies = np.array(conf_energies) - np.min(conf_energies) #covert to relative energies
    below_cutoff = (rel_conf_energies <= conf_cutoff).sum() #get number of conf below cutoff

    conf_tuble = list(zip(conf_names, conf_complex_mols, conf_paths, conf_energies, rel_conf_energies, calc_logs)) #make a tuble
    conf_tuble = sorted(conf_tuble, key=itemgetter(4))[:below_cutoff] #get only the best conf below cutoff

    conf_names, conf_complex_mols, conf_paths, conf_energies, rel_conf_energies, calc_logs = zip(*conf_tuble) #unzip tuble
    conf_names, conf_complex_mols, conf_paths, conf_energies, calc_logs = list(conf_names), list(conf_complex_mols), list(conf_paths), list(conf_energies), list(calc_logs) #tubles to lists
    mol_files = [os.path.join(conf_path, conf_name+'_gfn'+method.replace(' ', '') + '_opt.sdf') for conf_path, conf_name in zip(conf_paths,conf_names)] #list of paths to optimized structures in .sdf format

    # Find only unique conformers
    conf_names, conf_complex_mols, conf_paths, conf_energies, calc_logs = zip(*molfmt.find_unique_confs(list(zip(conf_names, conf_complex_mols, conf_paths, conf_energies, calc_logs)), mol_files, threshold=0.5)) #find unique conformers
    conf_names, conf_complex_mols, conf_paths, conf_energies, calc_logs = list(conf_names), list(conf_complex_mols), list(conf_paths), list(conf_energies), list(calc_logs) #tubles to lists

    return conf_names, conf_complex_mols, conf_paths, conf_energies, calc_logs


def calculateEnergy(args):
    """ Embed the post-insertion complex and calculate the ground-state free energy 
    return: energy [kcal/mol]
    """
    
    global cpus_per_calc
    global mem_gb

    # Calculation settings
    rdkit_mol, name = args
    method=' 1'  # <--- change the method for accurate calculations ('ff', ' 0', ' 1', ' 2'). Remember to change the references, methyl_anion_ref and methyl_cation_ref, if method this is changed!
    solvent = '--alpb DMSO' # <--- change the solvent ('--gbsa solvent_name', '--alpb solvent_name', or ''). Remember to change the references, methyl_anion_ref and methyl_cation_ref, if method this is changed!
    chrg = Chem.GetFormalCharge(rdkit_mol) # checking and adding formal charge to "chrg"
    spin = 0 # OBS! Spin is hardcoded to zero!

    # RDkit conf generator
    rdkit_mol = Chem.AddHs(rdkit_mol)
    rot_bond = Chem.rdMolDescriptors.CalcNumRotatableBonds(rdkit_mol)
    n_conformers = min(1 + 3 * rot_bond, 20)
        
    AllChem.EmbedMultipleConfs(rdkit_mol, numConfs=n_conformers,
            useExpTorsionAnglePrefs=True,
            useBasicKnowledge=True, ETversion=2, randomSeed=90)

    # Unpack confomers and assign conformer names
    conf_mols = [Chem.Mol(rdkit_mol, False, i) for i in range(rdkit_mol.GetNumConformers())]
    conf_names = [name + f'_conf{str(i+1).zfill(2)}' for i in range(rdkit_mol.GetNumConformers())] #change zfill(2) if more than 99 conformers
    conf_names_copy = copy.deepcopy(conf_names)

    # Run a GFN-FF optimization - Prescreening
    conf_names, conf_mols, conf_paths, conf_energies, calc_logs = confsearch_xTB(conf_mols, conf_names, chrg=chrg, spin=spin, method='ff', solvent=solvent, conf_cutoff=10, precalc_path=None)
    
    # Run a GFN?-xTB optimization
    conf_names, conf_mols, conf_paths, conf_energies, calc_logs = confsearch_xTB(conf_mols, conf_names, chrg=chrg, spin=spin, method=method, solvent=solvent, conf_cutoff=10, precalc_path=conf_paths)
    
    # Run Orca single point calculations
    final_conf_energies = []
    final_conf_mols = []
    final_conf_calc_logs = []
    for conf_name, conf_mol, conf_path, conf_energy, calc_log in zip(conf_names, conf_mols, conf_paths, conf_energies, calc_logs):
        # if conf_energy != 60000.0: # do single point calculations on all unique conformers
        #     conf_energy = run_orca.run_orca('xtbopt.xyz', chrg, spin, conf_path, ncores=cpus_per_calc, mem=(int(mem_gb)/2)*1000, optimize=True)
        final_conf_energies.append(conf_energy)
        final_conf_mols.append(conf_mol)
        final_conf_calc_logs.append(calc_log)
    
    # Get only the lowest energy conformer
    minE_index = np.argmin(final_conf_energies)
    
    best_conf_mol = final_conf_mols[minE_index]
    
    best_conf_energy = final_conf_energies[minE_index]
    if best_conf_energy != 60000.0:
        best_conf_energy = run_orca.run_orca('xtbopt.xyz', chrg, spin, conf_paths[minE_index], ncores=cpus_per_calc, mem=(int(mem_gb)/2)*1000, optimize=False) # comment when doing single point calculations on all unique conformers, otherwise this runs a Orca single point calculation on the lowest xTB energy conformer
    
    best_conf_calc_log = final_conf_calc_logs[minE_index]

    
    ### START - CLEAN UP ###
    for conf_name in conf_names_copy:

        conf_path = os.path.join(base_dir, 'calculations', '/'.join(conf_name.split('_')))
        
        if os.path.isfile(os.path.join(conf_path, conf_name+'_gfnff.xyz')):
            os.remove(os.path.join(conf_path, conf_name+'_gfnff.xyz'))
        
        if os.path.isfile(os.path.join(conf_path, conf_name+'_gfn'+ method.replace(' ', '') + '.xyz')):
            os.remove(os.path.join(conf_path, conf_name+'_gfn'+ method.replace(' ', '') + '.xyz'))
        
        # Remove GFNFF-xTB folder
        folder_path = os.path.join(conf_path, 'gfnff')
        if os.path.exists(folder_path):
            for file_remove in os.listdir(folder_path):
                if os.path.isfile(f'{folder_path}/{file_remove}'):
                    os.remove(f'{folder_path}/{file_remove}')
            # checking whether the folder is empty or not
            if len(os.listdir(folder_path)) == 0:
                os.rmdir(folder_path)
            else:
                print("Folder is not empty")
    
        # Clean GFN?-xTB folder
        folder_path = os.path.join(conf_path, 'gfn' + method.replace(' ', ''))
        if os.path.exists(folder_path):
            for file_remove in os.listdir(folder_path):
                if file_remove.split('.')[-1] in ['sdf', 'xtbout', 'out', 'xyz'] and file_remove != 'xtbscreen.xyz':
                    continue
                elif os.path.isfile(f'{folder_path}/{file_remove}'):
                    os.remove(f'{folder_path}/{file_remove}')
    ### END - CLEAN UP ###

    return best_conf_energy, best_conf_mol, best_conf_calc_log



def calc_MAA_and_MCA(reac_smis: str, name: str):

    global parallel_calcs

    reac_smis = reac_smis.split('.')
    reac_mols = [Chem.MolFromSmiles(smi) for smi in reac_smis]
    reac_names = [f'{name}_reac{str(i+1).zfill(2)}' for i in range(len(reac_smis))]

    ### Generate MAA and MCA products ###
    prod_mols = [] # list with molecule objects of all generated products
    prod_names = [] # list with filenames of all generated products
    prod_amounts = [] # list with the number of generated products for each type

    elec_sites_list = [[] for _ in range(len(reac_mols))] # list with electrophilic site indices seperated for each reactant
    elec_names_list = [[] for _ in range(len(reac_mols))] # list with electrophilic site names seperated for each reactant
    elec_prod_smis_list = [[] for _ in range(len(reac_mols))] # list with electrophilic site indices seperated for each reactant
    nuc_sites_list = [[] for _ in range(len(reac_mols))] # list with nucleophilic site indices seperated for each reactant
    nuc_names_list = [[] for _ in range(len(reac_mols))] # list with nucleophilic site names seperated for each reactant
    nuc_prod_smis_list = [[] for _ in range(len(reac_mols))] # list with nucleophilic site indices seperated for each reactant
    
    for i, (reac_mol, reac_name) in enumerate(zip(reac_mols, reac_names)):
        
        # MAA products
        MAA_prod_mols, MAA_prod_smis, elec_sites, elec_names = find_electrophilic_sites_and_generate_MAAproducts(reac_mol) # MAA_prod_mols: unique product mols, MAA_prod_smis: unique product SMILES, elec_sites: site in reactant, elec_names: name of functional group

        elec_sites_list[i].extend(elec_sites) # electrophilic site indices
        elec_names_list[i].extend(elec_names) # electrophilic site names
        elec_prod_smis_list[i].extend(MAA_prod_smis) # electrophilic product SMILES

        prod_mols.extend(MAA_prod_mols) # molecule objects of the generated MAA products
        prod_names.extend([f'{reac_name}prodE{str(i+1).zfill(2)}' for i in range(len(MAA_prod_mols))]) # filenames of the generated MAA products
        prod_amounts.append(len(MAA_prod_mols))

        # MCA products
        MCA_prod_mols, MCA_prod_smis, nuc_sites, nuc_names = find_nucleophilic_sites_and_generate_MCAproducts(reac_mol) # MCA_prod_mols: unique product mols, MCA_prod_smis: unique product SMILES, nuc_sites: site in reactant, nuc_names: name of functional group

        nuc_sites_list[i].extend(nuc_sites) # nucleophilic site indices
        nuc_names_list[i].extend(nuc_names) # nucleophilic site names
        nuc_prod_smis_list[i].extend(MCA_prod_smis) # nucleophilic product SMILES

        prod_mols.extend(MCA_prod_mols) # molecule objects of the generated MCA products
        prod_names.extend([f'{reac_name}prodN{str(i+1).zfill(2)}' for i in range(len(MCA_prod_mols))]) # filenames of the generated MCA products
        prod_amounts.append(len(MCA_prod_mols))
    ### EMD ###


    # Calculate electronic energies in kJ/mol
    all_mols = reac_mols + prod_mols
    all_names = reac_names + prod_names
    args = [(rdkit_mol, name) for rdkit_mol, name in zip(all_mols, all_names)]
    with ThreadPoolExecutor(max_workers=parallel_calcs) as executor:
        results = executor.map(calculateEnergy, args)

    reac_energies = []
    reac_calc_logs = []
    prod_energies = []
    prod_calc_logs = []
    for i, result in enumerate(results):
        best_conf_energy, best_conf_mol, best_conf_calc_log = result
        if i < len(reac_mols):
            reac_energies.append(best_conf_energy)
            reac_calc_logs.append(best_conf_calc_log)
        else:
            prod_energies.append(best_conf_energy)
            prod_calc_logs.append(best_conf_calc_log)
    
    if len(prod_energies) != sum(prod_amounts):
        raise RuntimeError('WARNING! One or several calculations failed: prod_energies != sum(prod_amounts)')


    # Calculate methyl anion affinities (MAAs) and methyl cation affinities (MCAs)  
    MAA_values = [[] for _ in range(len(reac_mols))]
    MAA_calc_logs = [[] for _ in range(len(reac_mols))]
    MCA_values = [[] for _ in range(len(reac_mols))]
    MCA_calc_logs = [[] for _ in range(len(reac_mols))]

    # methyl_anion_ref = -10299.839324620309  # kJ/mol. GFN1-xTB ALPBsolvent DMSO: -3.923001615610 Eh = -10299.839324620309 kJ/mol.
    # methyl_cation_ref = -8328.781129518116 # kJ/mol. GFN1-xTB ALPBsolvent DMSO: -3.172265197289 Eh = -8328.781129518116 kJ/mol.
    # methyl_anion_ref = -10052.419640018385  # kJ/mol. GFN2-xTB ALPBsolvent DMSO: -3.828764434637 Eh = -10052.419640018385 kJ/mol.
    # methyl_cation_ref = -8233.08622397756 # kJ/mol. GFN2-xTB ALPBsolvent DMSO: -3.135816932689 Eh = -8233.08622397756 kJ/mol.
    methyl_anion_ref = -104803.89889591146  # kJ/mol. Orca OPT NumFreq r2SCAN-3c SMDsolvent DMSO on GFN1-xTB ALPBsolvent DMSO geometry: -39.91769694 Eh = -104803.89889591146 kJ/mol.
    methyl_cation_ref = -103897.65006587801 # kJ/mol. Orca OPT NumFreq r2SCAN-3c SMDsolvent DMSO on GFN1-xTB ALPBsolvent DMSO geometry: -39.57252499 Eh = -103897.65006587801 kJ/mol.

    it_prod_energies = iter(prod_energies)
    sliced_prod_energies_list = [list(islice(it_prod_energies, 0, i)) for i in prod_amounts] # all even slices are the MAA energies and all uneven slices are the MCA energies
    it_prod_calc_logs = iter(prod_calc_logs)
    sliced_prod_calc_logs_list = [list(islice(it_prod_calc_logs, 0, i)) for i in prod_amounts] # all even slices are the MAA calc logs and all uneven slices are the MCA calc logs
    for i, reac_energy in enumerate(reac_energies):
        
        if reac_energy not in [60000.0, 120000.0]:
            MAA_values[i] = [-(prod_energy - (reac_energy + methyl_anion_ref)) if prod_energy not in [60000.0, 120000.0] else -np.inf for prod_energy in sliced_prod_energies_list[i*2]]
            MCA_values[i] = [-(prod_energy - (reac_energy + methyl_cation_ref)) if prod_energy not in [60000.0, 120000.0] else -np.inf for prod_energy in sliced_prod_energies_list[i*2+1]]
        else:
            MAA_values[i] = [-np.inf for _ in sliced_prod_energies_list[i*2]]
            MCA_values[i] = [-np.inf for _ in sliced_prod_energies_list[i*2+1]]
        
        MAA_calc_logs[i] = [[reac_calc_logs[i], prod_calc_log] for prod_calc_log in sliced_prod_calc_logs_list[i*2]]
        MCA_calc_logs[i] = [[reac_calc_logs[i], prod_calc_log] for prod_calc_log in sliced_prod_calc_logs_list[i*2+1]]


    ### Draw the output ###
    result_svg = generate_structure(reac_mols, elec_sites_list, MAA_values, nuc_sites_list, MCA_values, molsPerRow=3)

    df_elec = generate_output_tables(reac_mols, elec_names_list, MAA_values, elec_sites_list, MAA_calc_logs, MAA_or_MCA='MAA')
    df_nuc = generate_output_tables(reac_mols, nuc_names_list, MCA_values, nuc_sites_list, MCA_calc_logs, MAA_or_MCA='MCA')

    result_output = html_output('.'.join(reac_smis), result_svg, df_elec, df_nuc)
    
    fd = open(f'{base_dir}/calculations/{name}.html','w')
    fd.write(result_output)
    fd.close()
    ### END ###

    # Print results
    for i, (elec_sites, elec_names, MAAs, nuc_sites, nuc_names, MCAs) in enumerate(zip(elec_sites_list, elec_names_list, MAA_values, nuc_sites_list, nuc_names_list, MCA_values), start=1):
        print(f'Reactant - #{i}')
        print('Electrophilic sites:', '\n', elec_names, '\n', elec_sites, '\n', MAAs)
        print('Nucleophilic sites:', '\n', nuc_names, '\n', nuc_sites, '\n', MCAs)

    return elec_sites_list, elec_names_list, elec_prod_smis_list, MAA_values, MAA_calc_logs, nuc_sites_list, nuc_names_list, nuc_prod_smis_list, MCA_values, MCA_calc_logs


def control_calcs(df):
    
    df['elec_sites'] = pd.Series(dtype='object')
    df['elec_names'] = pd.Series(dtype='object')
    df['elec_prod_smis'] = pd.Series(dtype='object')
    df['MAA_values'] = pd.Series(dtype='object')
    df['MAA_calc_logs'] = pd.Series(dtype='object')
    df['nuc_sites'] = pd.Series(dtype='object')
    df['nuc_names'] = pd.Series(dtype='object')
    df['nuc_prod_smis'] = pd.Series(dtype='object')
    df['MCA_values'] = pd.Series(dtype='object')
    df['MCA_calc_logs'] = pd.Series(dtype='object')
    for idx, row in df.iterrows():

        name = row['name']
        reac_smis = row['smiles'] #row['canon_smiles'], row['smiles']
        
        ### RUN CALCULATIONS ###
        elec_sites, elec_names, elec_prod_smis, MAA_values, MAA_calc_logs, nuc_sites, nuc_names, nuc_prod_smis, MCA_values, MCA_calc_logs = calc_MAA_and_MCA(reac_smis, name)
        
        df.at[idx, 'elec_sites'] = elec_sites
        df.at[idx, 'elec_names'] = elec_names
        df.at[idx, 'elec_prod_smis'] = elec_prod_smis
        df.at[idx, 'MAA_values'] = MAA_values
        df.at[idx, 'MAA_calc_logs'] = MAA_calc_logs
        df.at[idx, 'nuc_sites'] = nuc_sites
        df.at[idx, 'nuc_names'] = nuc_names
        df.at[idx, 'nuc_prod_smis'] = nuc_prod_smis
        df.at[idx, 'MCA_values'] = MCA_values
        df.at[idx, 'MCA_calc_logs'] = MCA_calc_logs

    return df


if __name__ == "__main__":

    import pandas as pd
    import submitit
    import sys

    args = parse_args()

    parallel_calcs = int(args.parallel_calcs)
    cpus_per_calc = int(args.cpus_per_calc)
    mem_gb = int(args.mem_gb)

    executor = submitit.AutoExecutor(folder=os.path.join(base_dir, f'submitit_results/{args.name}'))
    executor.update_parameters(
        name="esnuel",
        cpus_per_task=int(parallel_calcs*cpus_per_calc),
        mem_gb=mem_gb,
        timeout_min=int(args.timeout_min),
        slurm_partition=str(args.partition),
        slurm_array_parallelism=int(args.slurm_array_parallelism),
    )
    print(executor)


    if args.batch is None:
        
        ### CLI JOB ###
        jobs = []
        with executor.batch():
            job = executor.submit(calc_MAA_and_MCA, args.smiles, args.name)
            jobs.append(job)
        ### END ###
    
    else:

        ### Batch JOB ###
        df = pd.read_csv(os.path.join(execution_path, args.batch), sep=',', encoding='unicode_escape') # read .csv file that contains the columns: "name" and "smiles"
        # df = df.head(1)
        
        jobs = []
        with executor.batch():
            chunk_size = 1
            for start in range(0, df.shape[0], chunk_size):
                df_subset = df.iloc[start:start + chunk_size]
                job = executor.submit(control_calcs, df_subset)
                jobs.append(job)
        ### END ###
