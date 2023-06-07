# MIT License
#
# Copyright (c) 2022 Nicolai Ree
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

import argparse
import numpy as np

import lightgbm as lgb
from rdkit import Chem

from molecule_svg import generate_structure
from locate_atom_sites import find_electrophilic_sites_and_generate_MAAproducts, find_nucleophilic_sites_and_generate_MCAproducts
from smi2gcs.DescriptorCreator.PrepAndCalcDescriptor import Generator


def parse_args():
    """
    Argument parser so this can be run from the command line
    """
    parser = argparse.ArgumentParser(description='Run regioselectivity predictions from the command line')
    parser.add_argument('-s', '--smiles', default='C[C+:20](C)CC(C)(C)C1=C(C=CC(=C1)Br)[OH:10]',
                        help='SMILES input for regioselectivity predictions')
    parser.add_argument('-n', '--name', default='testmol', help='The name of the molecule')
    return parser.parse_args()



def calc_atomic_descriptors(smi: str, name: str):
    """
    From a SMILES, do a single conformer embedding, calculate CM5 atomic charges, 
    locate electrophilic and nucleophilic sites, and generate GCS atomic descriptors.
    Return: Lists with atomic indicies, functional group names, and descriptors for the electrophilic and nucleophilic sites.
    """
    
    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    cm5_list = generator.calc_CM5_charges(smi, name=name, optimize=False, save_output=True)

    MAA_prod_mols, MAA_prod_smis, elec_sites, elec_names = find_electrophilic_sites_and_generate_MAAproducts(rdkit_mol)
    elec_atom_indices, elec_descriptor_vector = generator.create_descriptor_vector(elec_sites, des[0], **des[1])
    sorted_elec_sites = []
    sorted_elec_names = []
    for site in elec_atom_indices:
        i = elec_sites.index(site)
        sorted_elec_sites.append(elec_sites[i])
        sorted_elec_names.append(elec_names[i])
    
    MCA_prod_mols, MCA_prod_smis, nuc_sites, nuc_names = find_nucleophilic_sites_and_generate_MCAproducts(rdkit_mol)
    nuc_atom_indices, nuc_descriptor_vector = generator.create_descriptor_vector(nuc_sites, des[0], **des[1])
    sorted_nuc_sites = []
    sorted_nuc_names = []
    for site in nuc_atom_indices:
        i = nuc_sites.index(site)
        sorted_nuc_sites.append(nuc_sites[i])
        sorted_nuc_names.append(nuc_names[i])
    
    return sorted_elec_sites, sorted_elec_names, elec_descriptor_vector, sorted_nuc_sites, sorted_nuc_names, nuc_descriptor_vector


if __name__ == "__main__":
    
    args = parse_args()

    ### SETUP ###
    # Load GCS descriptor generator
    generator = Generator()
    des =('GraphChargeShell', {'charge_type': 'cm5', 'n_shells': 5, 'use_cip_sort': True})

    # Load LightGBM regression model for electrophilicity predictions
    final_model_elec = lgb.Booster(model_file='../../models/LGBMRegressor_electrophilicity/std_way_elec/final_best_model.txt')
    # Load LightGBM regression model for nucleophilicity predictions
    final_model_nuc = lgb.Booster(model_file='../../models/LGBMRegressor_nucleophilicity/std_way/final_best_model.txt')
    ### END ###


    ### Run electrophilicity and nucleophilicity predictions on reactants ###
    reac_smis = args.smiles.split('.')
    reac_mols = [Chem.MolFromSmiles(smi) for smi in reac_smis]

    elec_sites_list = []
    elec_names_list = []
    elec_preds_list = []
    nuc_sites_list = []
    nuc_names_list = []
    nuc_preds_list = []
    for i, smi in enumerate(reac_smis, start=1):
        print(f'Reactant - #{i}')

        # Locate sites and generate atomic descriptors
        elec_sites, elec_names, elec_descriptor_vector, nuc_sites, nuc_names, nuc_descriptor_vector = calc_atomic_descriptors(smi, args.name)
        elec_sites_list.append(elec_sites)
        elec_names_list.append(elec_names)
        nuc_sites_list.append(nuc_sites)
        nuc_names_list.append(nuc_names)
        
        # Run predictions
        if len(elec_descriptor_vector):
            elec_preds = final_model_elec.predict(elec_descriptor_vector, num_iteration=final_model_elec.best_iteration)
            elec_preds_list.append(elec_preds)
        else:
            elec_preds_list.append([])
        print('Electrophilic sites:', elec_names_list[-1], elec_sites_list[-1], elec_preds_list[-1])
        
        if len(nuc_descriptor_vector):
            nuc_preds = final_model_nuc.predict(nuc_descriptor_vector, num_iteration=final_model_nuc.best_iteration)
            nuc_preds_list.append(nuc_preds)
        else:
            nuc_preds_list.append([])
        print('Nucleophilic sites:', nuc_names_list[-1], nuc_sites_list[-1], nuc_preds_list[-1])
    ### END ###
    

    ### Draw the output ###
    result_svg = generate_structure(reac_mols, elec_sites_list, elec_preds_list, nuc_sites_list, nuc_preds_list, molsPerRow=4)
    
    fd = open(f'./{args.name}.svg','w')
    fd.write(result_svg)
    fd.close()
    ### END ###