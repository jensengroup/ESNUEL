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

import numpy as np

from rdkit import Chem
from rdkit.Chem import rdmolfiles, AllChem, rdDetermineBonds
from rdkit.ML.Cluster import Butina
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def convert_xyz_to_sdf(xyzfile, sdffile, chrg=None):

    rdkit_mol = Chem.MolFromXYZFile(xyzfile)
    rdDetermineBonds.DetermineConnectivity(rdkit_mol, useHueckel=True)
    # rdDetermineBonds.DetermineBondOrders(rdkit_mol, charge=chrg)
    # rdDetermineBonds.DetermineBonds(rdkit_mol, charge=chrg, covFactor=1.3, allowChargedFragments=True, useHueckel=False, embedChiral=False, useAtomMap=False)

    if len(Chem.MolToSmiles(rdkit_mol).split('.')) != 1:
        # print('OBS! Trying to detemine bonds without Hueckel')
        rdkit_mol = Chem.MolFromXYZFile(xyzfile)
        rdDetermineBonds.DetermineConnectivity(rdkit_mol, useHueckel=False)

    if len(Chem.MolToSmiles(rdkit_mol).split('.')) != 1:
        # print('OBS! Trying to detemine bonds without Hueckel and covFactor=1.35')
        rdkit_mol = Chem.MolFromXYZFile(xyzfile)
        rdDetermineBonds.DetermineConnectivity(rdkit_mol, useHueckel=False, covFactor=1.35)

    writer = Chem.rdmolfiles.SDWriter(sdffile)
    writer.write(rdkit_mol)
    writer.close()

    return


def get_bonds(sdf_file):
    """ The count line is; aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv,
    where aaa is atom count and bbb is bond count.
    """

    atoms = 0
    bond_list = []

    searchlines = open(sdf_file, 'r').readlines()

    for i, line in enumerate(searchlines):
        words = line.split() #split line into words
        if len(words) < 1:
            continue
        if i == 3:
            atoms = int(line[0:3])
            bonds = int(line[3:6])
        if 'Pd' in words: #find atom index of Pd
            transistion_metal_idx = i - 3
        else:
            transistion_metal_idx = -1
        if i > atoms+3 and i <= atoms+bonds+3:
            atom_1 = int(line[0:3])
            atom_2 = int(line[3:6])
            if (atom_1 == transistion_metal_idx) or (atom_2 == transistion_metal_idx): #skip bonds to Pd
                continue
            if atom_2 > atom_1:
                bond_list.append(tuple((atom_1,atom_2)))
            else:
                bond_list.append(tuple((atom_2,atom_1)))

    bond_list.sort()

    return bond_list


def get_bonds_molblock(molblock):
    """ The count line is; aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv,
    where aaa is atom count and bbb is bond count.
    """

    atoms = 0
    bond_list = []

    searchlines = molblock.split('\n')

    for i, line in enumerate(searchlines):
        words = line.split() #split line into words
        if len(words) < 1:
            continue
        if i == 3:
            atoms = int(line[0:3])
            bonds = int(line[3:6])
        if 'Pd' in words: #find atom index of Pd
            transistion_metal_idx = i - 3
        else:
            transistion_metal_idx = -1
        if i > atoms+3 and i <= atoms+bonds+3:
            atom_1 = int(line[0:3])
            atom_2 = int(line[3:6])
            if (atom_1 == transistion_metal_idx) or (atom_2 == transistion_metal_idx): #skip bonds to Pd
                continue
            if atom_2 > atom_1:
                bond_list.append(tuple((atom_1,atom_2)))
            else:
                bond_list.append(tuple((atom_2,atom_1)))

    bond_list.sort()

    return bond_list


def compare_sdf_structure(start, end, molblockStart=False, molblockEnd=False):
    """
    Returns True if structures are the same

    Return False if there has been a proton transfer
    """
    
    if molblockStart:
        bond_start = get_bonds_molblock(start)
    else:
        bond_start = get_bonds(start)
    
    if molblockEnd:
        bond_end = get_bonds_molblock(end)
    else:    
        bond_end = get_bonds(end)

    return bond_start == bond_end


def find_unique_confs(best_conformers, mol_files, threshold=0.5):
    """ Clustering conformers with RDKit's Butina algorithm
    to find unique conformer from a list of .sdf files
    using either heavy-atom root mean square deviation (RMSD) 
    or heavy-atom torsion fingerprint deviation (TFD) """

    rdkit_mol = next(rdmolfiles.ForwardSDMolSupplier(mol_files[0], sanitize=False, removeHs=True))
    for mol_file in mol_files[1:]:
        mol = next(rdmolfiles.ForwardSDMolSupplier(mol_file, sanitize=False, removeHs=True))
        rdkit_mol.AddConformer(mol.GetConformer(), assignId=True)

    # calculate difference matrix
    diffmat = AllChem.GetConformerRMSMatrix(rdkit_mol, prealigned=False) #threshold=0.5, sanitize=False, load AllChem
    # diffmat = TorsionFingerprints.GetTFDMatrix(rdkit_mol) #threshold=0.01, sanitize=True, load TorsionFingerprints

    # Cluster conformers
    num_confs = rdkit_mol.GetNumConformers()
    clt = Butina.ClusterData(diffmat, num_confs, threshold,
                             isDistData=True, reordering=True)

    # Get unique conformers
    centroid_idx = [c[0] for c in clt] # centroid indexes
    unique_best_conformers = [best_conformers[i] for i in centroid_idx]
    
    return unique_best_conformers


def getResonanceStructures(smiles, flags = 0):
    mol = Chem.MolFromSmiles(smiles)
    suppl = Chem.ResonanceMolSupplier(mol, flags)
    resMols = [mol for mol in suppl]
    return (mol, resMols)


def embed_organic_mol(rdkit_mol: str):
    """ Generate a single conformer embedding from rdkit_mol """
    
    # Set parameters
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 90
    ps.useSmallRingTorsions=True
    ps.ETversion=2
    ps.useExpTorsionAnglePrefs=True
    ps.useBasicKnowledge=True
    
    # Embed
    if AllChem.EmbedMolecule(rdkit_mol, ps) == -1:
        print(f'1st embed; will try useRandomCoords=True')
        
        # Change parameters and try again
        ps.useRandomCoords=True
        ps.maxIterations=1000
        if AllChem.EmbedMolecule(rdkit_mol, ps) == -1:
            raise Exception(f'2nd embed failed')

    return rdkit_mol


def run_rxn(reactant, smarts, chrg=None):
    """ This functions is adapted to work for reaction SMARTS,
    where the atom with index 1 of the product is the attachment point """
    
    # Run reaction
    rxn = AllChem.ReactionFromSmarts(smarts)
    Chem.Kekulize(reactant,clearAromaticFlags=True)
    ps = rxn.RunReactants([reactant])
  
    # Find possible products
    product_mols = []
    product_smis = []
    product_sites = []
    if ps:
        ps = np.asarray(ps)[:,0] # Keep only first product of reaction
        #ps = np.concatenate(ps) # Concatenate all products
    else:
        return product_mols, product_smis, product_sites

    for mol in ps:
        # Canonicalize SMILES
        try:
            site = mol.GetAtomWithIdx(1).GetPropsAsDict()['react_atom_idx'] # OBS! the connecting atom is always the atom with index 1 in the product # Original
            # site = mol.GetAtomWithIdx(2).GetPropsAsDict()['react_atom_idx'] # OBS! the connecting atom is always the atom with index 1 in the product # FOR "-SMe"
            mol = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
            smi = Chem.MolToSmiles(Chem.RemoveHs(mol))
            smi = Chem.CanonSmiles(smi)
            
            # mol3d = embed_organic_mol(Chem.AddHs(Chem.MolFromSmiles(smi)))
            # AllChem.UFFOptimizeMolecule(mol3d)
            # XYZBlock = Chem.MolToXYZBlock(mol3d)
            # molXYZ = Chem.MolFromXYZBlock(XYZBlock)
            # rdDetermineBonds.DetermineBonds(molXYZ, charge=chrg, covFactor=1.3, allowChargedFragments=True, useHueckel=False, embedChiral=False, useAtomMap=False)
            # smi = Chem.MolToSmiles(Chem.RemoveHs(molXYZ))
            
            # if not reactant.HasSubstructMatch(Chem.MolFromSmarts('O=[N+][O-]')) and not reactant.HasSubstructMatch(Chem.MolFromSmarts('[N-]=[N+]=[N-]')) and not reactant.HasSubstructMatch(Chem.MolFromSmarts('N#[N+][CH2-]')) and not reactant.HasSubstructMatch(Chem.MolFromSmarts('C1=CC=CO1')) and not reactant.HasSubstructMatch(Chem.MolFromSmarts('COC1=CC=CC=C1')) and not reactant.HasSubstructMatch(Chem.MolFromSmarts('[O-]C(OCC)=C[N+]1=CC=CC=C1')): # OBS!!! THIS IS A HACK DUE TO ISSUE #4884, [N-]=[N+]=[N-], N#[N+][CH2-], C1=CC=CO1, COC1=CC=CC=C1 and [O-]C(OCC)=C[N+]1=CC=CC=C1 are new issues
            #     _, resMols = getResonanceStructures(smi)
            #     smi = Chem.MolToSmiles(resMols[0])

            # smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        except:
            continue
        
        # Keep only unique products that has the same number of atoms as the reactant as hydrogens are being removed and added
        if smi not in product_smis and Chem.AddHs(Chem.MolFromSmiles(smi)).GetNumAtoms() - 4 == Chem.AddHs(reactant).GetNumAtoms(): # Original
        # if smi not in product_smis and Chem.AddHs(Chem.MolFromSmiles(smi)).GetNumAtoms() - 5 == Chem.AddHs(reactant).GetNumAtoms(): # FOR "-SMe"
            product_mols.append(Chem.MolFromSmiles(smi))
            product_smis.append(smi)
            product_sites.append(site)
    
    return product_mols, product_smis, product_sites


if __name__ == "__main__":
    
    import sys

    print(compare_sdf_structure(sys.argv[1], sys.argv[2]))
