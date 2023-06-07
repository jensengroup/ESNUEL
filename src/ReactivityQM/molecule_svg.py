import matplotlib
import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import PrepareMolForDrawing
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(False)

from collections import defaultdict

# Drawing Options
elec_color = matplotlib.colors.ColorConverter().to_rgb('lightskyblue') # electrophilic sites are highlighted with this color
nuc_color = matplotlib.colors.ColorConverter().to_rgb('gold') # nucleophilic sites are highlighted with this color
arad = 0.25
subImgSize = (300,300)


def find_identical_atoms(rdkit_mol, atom_list):
    len_list = len(atom_list)
    
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(rdkit_mol.GetAtoms()):
        if atom.GetIdx() in atom_list[:len_list]:
            sym_atoms = [int(atom_idx) for atom_idx, ranking in enumerate(atom_rank) if ranking == atom_rank[idx] and atom_idx not in atom_list] 
            atom_list.extend(sym_atoms)
    return atom_list


def draw2d(rdkit_mol, elec_sites, MAA_values, top_score_elec, nuc_sites, MCA_values, top_score_nuc, subImgSize):
    
    global elec_color
    global nuc_color
    global arad

    d2d = rdMolDraw2D.MolDraw2DSVG(subImgSize[0], subImgSize[1])
    dos = d2d.drawOptions()
    dos.atomHighlightsAreCircles = False
    dos.fillHighlights = True

    #label atoms with scores
    atomHighlighs = defaultdict(list)
    highlightRads = {}

    # Locate electrophilic sites within the top_score_elec value
    for idx, val in enumerate(MAA_values):
        if val >= top_score_elec:
            site = elec_sites[idx] # the atomic index of the located site
            identical_sites = find_identical_atoms(rdkit_mol, [site])
            for site in identical_sites:
                atomHighlighs[site].append(elec_color)
                highlightRads[site] = arad

    # Locate nucleophilic sites within the top_score_nuc value
    for idx, val in enumerate(MCA_values):
        if val >= top_score_nuc:
            site = nuc_sites[idx] # the atomic index of the located site
            identical_sites = find_identical_atoms(rdkit_mol, [site])
            for site in identical_sites:
                atomHighlighs[site].append(nuc_color)
                highlightRads[site] = arad
    
    rdkit_mol = PrepareMolForDrawing(rdkit_mol)
    d2d.DrawMoleculeWithHighlights(rdkit_mol, '', dict(atomHighlighs), {}, highlightRads, {})
    d2d.FinishDrawing()

    return d2d.GetDrawingText()


def generate_structure(rdkit_mols, elec_sites_list, MAA_values_list, nuc_sites_list, MCA_values_list, molsPerRow=4):

    global subImgSize

    cutoff = 3.0 # DFT scores are within a range of 3.0 orders of magnitude from experimental results according to https://doi.org/10.1021/acs.jcim.1c01400
    top_score_elec = max(np.concatenate(MAA_values_list), default=np.inf) - cutoff # all electrophilic sites within this value are highlighted
    top_score_nuc = max(np.concatenate(MCA_values_list), default=np.inf) - cutoff # all nucleophilic sites within this value are highlighted

    nRows = len(rdkit_mols) // molsPerRow
    if len(rdkit_mols) % molsPerRow:
        nRows += 1
    if nRows == 1:
        molsPerRow = len(rdkit_mols)
    fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])

    header = """<svg version='1.1' baseProfile='full'
                xmlns='http://www.w3.org/2000/svg'
                        xmlns:rdkit='http://www.rdkit.org/xml'
                        xmlns:xlink='http://www.w3.org/1999/xlink'
                    xml:space='preserve'
    width='{0}px' height='{1}px' viewBox='0 0 {0} {1}'>
    <!-- END OF HEADER -->""".format(fullSize[0],fullSize[1])

    spacer = '<g transform="translate({0},{1})">\n{2}</g>'
    
    cwidth = 0
    cheight = 0
    drawed_mols = []
    for i in range(len(rdkit_mols)):
        res = draw2d(rdkit_mols[i], elec_sites_list[i], MAA_values_list[i], top_score_elec, nuc_sites_list[i], MCA_values_list[i], top_score_nuc, subImgSize)
        res = res.split("\n")
        end_of_header = res.index("<!-- END OF HEADER -->") + 1 
        res = "\n".join(res[end_of_header:-2])
        
        res = "".join(spacer.format(int(cwidth*subImgSize[0]), int(cheight*subImgSize[1]), res))
        drawed_mols.append(res)
        
        if int(i+1) % molsPerRow == 0 and i != 0:
            cheight += 1
            cwidth = 0
        elif molsPerRow == 1:
            cheight += 1
            cwidth = 0
        else:
            cwidth += 1

    svg = header + "\n" + "\n".join(drawed_mols) + "\n</svg>"

    return svg


if __name__ == "__main__":
    
    rdkit_smiles = ['CC(C)(C)[CH:1]=[OH+]']
    rdkit_mols = [Chem.MolFromSmiles(smi) for smi in rdkit_smiles]
    
    elec_sites_list = [[4]]
    MAA_values_list = [[481.16354902635794]]
    nuc_sites_list = [[4]]
    MCA_values_list = [[-2.0]]
    
    result_svg = generate_structure(rdkit_mols, elec_sites_list, MAA_values_list, nuc_sites_list, MCA_values_list, molsPerRow=4)
    
    fd = open('test.svg','w')
    fd.write(result_svg)
    fd.close()
    


