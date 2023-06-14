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

import matplotlib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import PrepareMolForDrawing
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(False)

from collections import defaultdict

# Drawing Options
elec_color = matplotlib.colors.ColorConverter().to_rgb('lightskyblue') + (0.7,) # electrophilic sites are highlighted with this color (before:lightskyblue)
nuc_color = matplotlib.colors.ColorConverter().to_rgb('lightgreen') + (0.7,) # nucleophilic sites are highlighted with this color (before:gold)
arad = 0.25
subImgSize = (350,350)


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
    dos.useBWAtomPalette()
    dos.atomHighlightsAreCircles = False
    dos.fillHighlights = True
    dos.addAtomIndices = True
    dos.minFontSize = 18
    dos.annotationFontScale = 0.85


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

    cutoff = 12.6 # kJ/mol (3 kcal/mol = 12.6 kJ/mol)
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


def generate_output_tables(rdkit_mols, names_list, values_list, sites_list, calc_logs, MAA_or_MCA='MAA'):
    
    values_list_new = []
    sites_list_new = []
    calc_logs_new = []

    for i in range(len(rdkit_mols)):
            
            for idx, val in enumerate(values_list[i]):
                
                site = sites_list[i][idx] # the atomic index of the located site
                identical_sites = find_identical_atoms(rdkit_mols[i], [site])
                for site in identical_sites:
                    
                    values_list_new.append(val)
                    
                    if len(rdkit_mols) > 1: # if more than one molecule add molecule ID to Atom ID
                        sites_list_new.append(f'{i+1}.{site}')
                    else:
                        sites_list_new.append(f'{site}')

                    calc_logs_new.append(",".join([Chem.MolToSmiles(Chem.RemoveHs(Chem.MolFromSmiles(smi))) if smi not in [None, 'xtbopt.xyz was not created'] else str(smi) for smi in calc_logs[i][idx]]))
    

    if MAA_or_MCA == 'MAA':
    
        dict_table = {'Atom ID': sites_list_new, 'MAA Value [kJ/mol]': values_list_new, 'Error Log (Reactant, Product)': calc_logs_new, 'Type': [n.replace('_', ' ').capitalize() for n in np.concatenate(names_list)]}
        df_table = pd.DataFrame(dict_table).sort_values(by=['MAA Value [kJ/mol]'], ascending=False)
    
    elif MAA_or_MCA == 'MCA':
    
        dict_table = {'Atom ID': sites_list_new, 'MCA Value [kJ/mol]': values_list_new, 'Error Log (Reactant, Product)': calc_logs_new, 'Type': [n.replace('_', ' ').capitalize() for n in np.concatenate(names_list)]}
        df_table = pd.DataFrame(dict_table).sort_values(by=['MCA Value [kJ/mol]'], ascending=False)
    
    return df_table


def html_output(rdkit_smiles, svg_file, df_elec, df_nuc):


    html_table_elec = df_elec.to_html(index=False, decimal='.', float_format='%.2f').replace('<table border="1" class="dataframe">', '<table>\n<caption class="caption">Electrophilicity (MAA)</caption>').replace('<tr style="text-align: right;">', '<tr>')
    html_table_nuc = df_nuc.to_html(index=False, decimal='.', float_format='%.2f').replace('<table border="1" class="dataframe">', '<table>\n<caption class="caption">Nucleophilicity (MCA)</caption>').replace('<tr style="text-align: right;">', '<tr>')

    html_head = """
            <!DOCTYPE html>
            <html>
            <head>
            <title>Website Design Example</title>
            <style>
                /* Add some basic styling to improve the appearance */
                body {
                font-family: Arial, sans-serif;
                background-color: #f7f7f7;
                margin: 0;
                padding: 20px;
                }

                .container {
                max-width: 1000px;
                margin: 0 auto;
                background-color: #fff;
                padding: 20px;
                border-radius: 5px;
                box-shadow: 0 2px 5px rgba(0, 0, 0, 0.1);
                display: flex;
                flex-direction: column;
                align-items: center;
                }

                svg {
                display: block;
                max-width: 100%;
                height: auto;
                margin-bottom: 20px;
                margin: 0 auto;
                }

                .text-holder1 {
                font-size: 20px;
                font-weight: bold;
                margin-bottom: 5px;
                color: black;
                }

                .text-holder2 {
                font-size: 16px;
                font-weight: normal;
                margin-bottom: 5px;
                color: black;
                }

                .tables-container {
                display: flex;
                justify-content: center;
                align-items: flex-start;
                }

                .table {
                flex: 0 0 50%;
                padding: 10px;
                }

                table {
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
                }

                th, td {
                padding: 14px;
                text-align: left;
                border-bottom: 1px solid #ddd;
                color: black;
                }

                th {
                background-color: #f2f2f2;
                }

                .table:nth-child(1) caption {
                font-size: 18px;
                font-weight: bold;
                vertical-align: middle;
                line-height: 35px;
                height: 35px;
                margin-bottom: 0px;
                color: black;
                background: #87cefa80;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
                }

                .table:nth-child(2) caption {
                font-size: 18px;
                font-weight: bold;
                vertical-align: middle;
                line-height: 35px;
                height: 35px;
                margin-bottom: 0px;
                color: black;
                background: #90ee9080;
                border-top-left-radius: 5px;
                border-top-right-radius: 5px;
                }
            </style>
            </head>"""
    
    html_body = """
            <body>
            <div class="container">
                <div class="text-holder1">SMILES:</div>
                <div class="text-holder2">{0}</div>
                
                {1}

                <div class="tables-container">
                    <div class="table">
                        {2}
                    </div>

                    <div class="table">
                        {3}
                    </div>
                </div>
            </div>
            </body>
            </html>""".format(rdkit_smiles, svg_file, html_table_elec, html_table_nuc)
            
    return html_head + html_body


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
    


