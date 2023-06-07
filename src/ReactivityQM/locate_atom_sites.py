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
import molecule_formats as molfmt

# n_smirks_names_add_list = ["Phenol",
#                           "Silyl_ether",  
#                           "Pyridine_like_nitrogen"] 

# n_smirks_add_list = ['[OX2H:1][$(c(c)c),$([#6X3;R](=[#6X3;R])[#6X3;R]):2]>>[CH3][OX3H+1:1][*:2]', 
#                     '[#8X2H0:1][#14X4:2]([!#1:3])([!#1:4])[!#1:5]>>[CH3][#8X3H0+:1][*:2]([*:3])([*:4])[*:5]', 
#                     '[#7X2;$([nX2](:*):*),$([#7X2;R](=[*;R])[*;R]):1]>>[CH3][#7X3+:1]'] 

# n_smirks_names_list = ["Ether", 
#                       "Ketone", 
#                       "Amide", 
#                       "Enolate", 
#                       "Aldehyde", 
#                       "Imine", 
#                       "Nitranion", 
#                       "Carbanion", 
#                       "Nitronate", 
#                       "Ester", 
#                       "Carboxylic acid", 
#                       "Amine", 
#                       "Cyanoalkyl/nitrile anion",
#                       "Nitrile", 
#                       "Isonitrile"] + n_smirks_names_add_list

# n_smirks_list = ['[OX2:1]([#6;!$(C([OX2])[#7,#8,#15,#16,F,Cl,Br,I]);!$([#6]=[#8]):2])[#6;!$(C([OX2])[#7,#8,#15,#16]);!$([#6]=[#8]):3]>>[CH3][OX3+:1]([*:2])[*:3]', 
#                 '[OX1H0:1]=[#6X3:2]([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][OX2H0:1][#6X3+:2]([*:3])[*:4]',
#                 '[OX1:1]=[CX3;$([CX3][#6]),$([CX3H]):2][#7X3;!R:3]>>[CH3][OX2:1][CX3:2]=[#7X3+:3]', 
#                 '[#6;$([#6]=,:[#6]-[#8-]),$([#6-]-[#6]=,:[#8]):1]~[#6:2]~[#8;$([#8-]-[#6]=,:[#6]),$([#8]=,:[#6]-[#6-]):3]>>[CH3][#6+0:1][*:2]=[#8+0:3]',
#                 '[OX1:1]=[$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):2]>>[CH3][OX2:1][#6+:2]', 
#                 '[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):1]=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):2]>>[CH3][NX3+:1]=[*:2]',
#                 '[#7X2-:1]>>[CH3][#7X3+0:1]', 
#                 '[#6-;!$([#6X1-]#[#7,#8,#15,#16]):1]>>[CH3][#6+0:1]', 
#                 '[#6:1]=[#7+:2](-[#8-:3])-[#8-:4]>>[CH3][#6:1][#7+:2](=[#8+0:3])-[*:4]', 
#                 '[OX1:1]=[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]):2][#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][OX2:1][#6X3+:2][*:3][*:4]', 
#                 '[OX1:1]=[CX3;$([R0][#6]),$([H1R0]):2][$([OX2H]),$([OX1-]):3]>>[CH3][OX2:1][CX3+:2][*:3]',
#                 '[#7+0;$([N;R;!$([#7X2]);$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]),$([NX3+0;!$([#7X3][CX3;$([CX3][#6]),$([CX3H])]=[OX1])]),$([NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]):1]>>[CH3][#7+:1]', 
#                 '[C:1]=[C:2]=[#7X1-:3]>>[CH3][C:1][C:2]#[#7X1+0:3]',
#                 '[NX1:1]#[CX2;!$(CC=C=[#7X1-]);!$(CC=C):2]>>[CH3][NX2+:1]#[*:2]', 
#                 '[CX1-:1]#[NX2+:2]>>[CH3][CX2+0:1]#[NX2+:2]'] + n_smirks_add_list


# e_smirks_names_add_list = ["Anhydride", 
#                           "Acyl Halide", 
#                           "Halide leaving group"]

# e_smirks_add_list = ['[CX3:1](=[OX1:2])[OX2:3][CX3:4]=[OX1:5]>>[CH3][CX4:1](-[OX1-:2])[*:3][*:4]=[*:5]', 
#                     '[CX3:1](=[OX1:2])[ClX1,BrX1,IX1:3]>>[CH3][CX4:1](-[OX1-:2])[*:3]', 
#                     '[C;!$([CX3]=[OX1]);!$([CX3;!R]=[CX3][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-])]):1][ClX1,BrX1,$(OS(=O)(=O)[CH3]),$(OS(=O)(=O)c1ccc([CH3])cc1),$(OS(=O)(=O)C(F)(F)F),IX1:2]>>[CH3][*:1].[*:2]']

# e_smirks_names_list = ["Oxonium", 
#                       "Carbocation", 
#                       "Ketone", 
#                       "Amide", 
#                       "Ester",
#                       "Iminium", 
#                       "Michael acceptor", 
#                       "Imine", 
#                       "Aldehyde"] + e_smirks_names_add_list

# e_smirks_list = ['[#6:1]=[O+;!$([O]~[!#6]);!$([S]*~[#7,#8,#15,#16]):2]>>[CH3][#6:1][O+0:2]', 
#                 '[#6+:1]>>[CH3][#6+0:1]', 
#                 '[#6X3:1](=[OX1:2])([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][#6X4:1](-[OX1-:2])([*:3])[*:4]', 
#                 '[CX3;$([CX3][#6]),$([CX3H]):1](=[OX1:2])[#7X3;!R:3]>>[CH3][CX4:1](-[OX1-:2])[*:3]',
#                 '[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]),$([#6X3][OX2H0]):1](=[OX1:2])[#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][#6X4:1](-[OX1-:2])[*:3][*:4]', 
#                 '[CX3:1]=[NX3+;!$(N([#8-])[#8-]):2]>>[CH3][CX4:1]-[NX3+0:2]', 
#                 '[CX3;!R:1]=[CX3:2][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]):3]>>[CH3][CX4:1]-[CX3-:2][*:3]',
#                 '[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):1]=[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):2]>>[CH3][CX4:1]-[NX2-:2]',
#                 '[CX3;$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):1]=[OX1:2]>>[CH3][CX4:1]-[OX1-:2]'] + e_smirks_add_list

# e_smarts_list = pd.read_csv('../../data_Baldi/smartsPatterns/SMARTS_electrophilicity.csv', header=None, dtype=str, delimiter='\n')[0].tolist()
# e_smarts_names_list = pd.read_csv('../../data_Baldi/smartsPatterns/SMARTS_electrophilicity_names.csv', header=None, dtype=str, delimiter='\n')[0].tolist()

# n_smarts_list = pd.read_csv('../../data_Baldi/smartsPatterns/SMARTS_nucleophilicity.csv', header=None, dtype=str, delimiter='\n')[0].tolist()
# n_smarts_names_list = pd.read_csv('../../data_Baldi/smartsPatterns/SMARTS_nucleophilicity_names.csv', header=None, dtype=str, delimiter='\n')[0].tolist()


# n_smirks_names_list = ['anion_with_charge_minus1',
#                       'double_bond',
#                       'double_bond_neighbouratom_with_charge_plus1',
#                       'triple_bond',
#                       'triple_bond_neighbouratom_with_charge_plus1',
#                       'atom_with_lone_pair',
#                       ] 

# n_smirks_list = ['[*-:1]>>[CH3][*+0:1]',
#                 '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+0:2]>>[CH3][*:1]-[*+1:2]',
#                 '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+1:2]>>[CH3][*:1]-[*+2:2]',
#                 '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+0:2]>>[CH3][*:1]=[*+1:2]',
#                 '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+1:2]>>[CH3][*:1]=[*+2:2]',
#                 '[!X4;!#1;!#6:1]>>[CH3][*+1:1]',
#                 ]


# e_smirks_names_list = ['cation_with_charge_plus1',
#                       'double_bond',
#                       'double_bond_neighbouratom_with_charge_plus1',
#                       'triple_bond',
#                       'triple_bond_neighbouratom_with_charge_plus1',
#                       ] 

# e_smirks_list = ['[*+:1]>>[CH3][*+0:1]',
#                 '[*+0:1]=[*+0:2]>>[CH3][*:1]-[*-1:2]',
#                 '[*+0:1]=[*+1:2]>>[CH3][*:1]-[*+0:2]',
#                 '[*+0:1]#[*+0:2]>>[CH3][*:1]=[*-1:2]',
#                 '[*+0:1]#[*+1:2]>>[CH3][*:1]=[*+0:2]',
#                 ]


n_smirks_names_list = ["Ether", 
                      "Ketone", 
                      "Amide", 
                      "Enolate", 
                      "Aldehyde", 
                      "Imine", 
                      "Nitranion", 
                      "Carbanion", 
                      "Nitronate", 
                      "Ester", 
                      "Carboxylic acid", 
                      "Amine", 
                      "Cyanoalkyl/nitrile anion",
                      "Nitrile", 
                      "Isonitrile",
                      "Phenol", # added due to rxn100
                      "Silyl_ether", # added due to rxn100
                      "Pyridine_like_nitrogen", # added due to rxn100
                      'anion_with_charge_minus1', # added to capture additional sites
                      'double_bond', # added to capture additional sites
                      'double_bond_neighbouratom_with_charge_plus1', # added to capture additional sites
                      'triple_bond', # added to capture additional sites
                      'triple_bond_neighbouratom_with_charge_plus1', # added to capture additional sites
                      'atom_with_lone_pair' # added to capture additional sites
                      ]

n_smirks_list = ['[OX2:1]([#6;!$(C([OX2])[#7,#8,#15,#16,F,Cl,Br,I]);!$([#6]=[#8]):2])[#6;!$(C([OX2])[#7,#8,#15,#16]);!$([#6]=[#8]):3]>>[CH3][OX3+:1]([*:2])[*:3]', 
                '[OX1H0:1]=[#6X3:2]([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][OX2H0:1][#6X3+:2]([*:3])[*:4]',
                '[OX1:1]=[CX3;$([CX3][#6]),$([CX3H]):2][#7X3;!R:3]>>[CH3][OX2:1][CX3:2]=[#7X3+:3]', 
                '[#6;$([#6]=,:[#6]-[#8-]),$([#6-]-[#6]=,:[#8]):1]~[#6:2]~[#8;$([#8-]-[#6]=,:[#6]),$([#8]=,:[#6]-[#6-]):3]>>[CH3][#6+0:1][*:2]=[#8+0:3]',
                '[OX1:1]=[$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):2]>>[CH3][OX2:1][#6+:2]', 
                '[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):1]=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):2]>>[CH3][NX3+:1]=[*:2]',
                '[#7X2-:1]>>[CH3][#7X3+0:1]', 
                '[#6-;!$([#6X1-]#[#7,#8,#15,#16]):1]>>[CH3][#6+0:1]', 
                '[#6:1]=[#7+:2](-[#8-:3])-[#8-:4]>>[CH3][#6:1][#7+:2](=[#8+0:3])-[*:4]', 
                '[OX1:1]=[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]):2][#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][OX2:1][#6X3+:2][*:3][*:4]', 
                '[OX1:1]=[CX3;$([R0][#6]),$([H1R0]):2][$([OX2H]),$([OX1-]):3]>>[CH3][OX2:1][CX3+:2][*:3]',
                '[#7+0;$([N;R;!$([#7X2]);$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]),$([NX3+0;!$([#7X3][CX3;$([CX3][#6]),$([CX3H])]=[OX1])]),$([NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]):1]>>[CH3][#7+:1]', 
                '[C:1]=[C:2]=[#7X1-:3]>>[CH3][C:1][C:2]#[#7X1+0:3]',
                '[NX1:1]#[CX2;!$(CC=C=[#7X1-]);!$(CC=C):2]>>[CH3][NX2+:1]#[*:2]', 
                '[CX1-:1]#[NX2+:2]>>[CH3][CX2+0:1]#[NX2+:2]',
                '[OX2H:1][$(c(c)c),$([#6X3;R](=[#6X3;R])[#6X3;R]):2]>>[CH3][OX3+1:1][*:2]', # added due to rxn100
                '[#8X2H0:1][#14X4:2]([!#1:3])([!#1:4])[!#1:5]>>[CH3][#8X3H0+:1][*:2]([*:3])([*:4])[*:5]', # added due to rxn100
                '[#7X2;$([nX2](:*):*),$([#7X2;R](=[*;R])[*;R]):1]>>[CH3][#7X3+:1]', # added due to rxn100
                '[*-:1]>>[CH3][*+0:1]', # added to capture additional sites
                '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+0:2]>>[CH3][*:1]-[*+1:2]', # added to capture additional sites
                '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+1:2]>>[CH3][*:1]-[*+2:2]', # added to capture additional sites
                '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+0:2]>>[CH3][*:1]=[*+1:2]', # added to capture additional sites
                '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+1:2]>>[CH3][*:1]=[*+2:2]', # added to capture additional sites
                '[!X4;!#1;!#6:1]>>[CH3][*+1:1]' # added to capture additional sites
                ]
# ### BEGIN SMe ###
# n_smirks_list = ['[OX2:1]([#6;!$(C([OX2])[#7,#8,#15,#16,F,Cl,Br,I]);!$([#6]=[#8]):2])[#6;!$(C([OX2])[#7,#8,#15,#16]);!$([#6]=[#8]):3]>>[CH3][#16X2][OX3+:1]([*:2])[*:3]', 
#                 '[OX1H0:1]=[#6X3:2]([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][#16X2][OX2H0:1][#6X3+:2]([*:3])[*:4]',
#                 '[OX1:1]=[CX3;$([CX3][#6]),$([CX3H]):2][#7X3;!R:3]>>[CH3][#16X2][OX2:1][CX3:2]=[#7X3+:3]', 
#                 '[#6;$([#6]=,:[#6]-[#8-]),$([#6-]-[#6]=,:[#8]):1]~[#6:2]~[#8;$([#8-]-[#6]=,:[#6]),$([#8]=,:[#6]-[#6-]):3]>>[CH3][#16X2][#6+0:1][*:2]=[#8+0:3]',
#                 '[OX1:1]=[$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):2]>>[CH3][#16X2][OX2:1][#6+:2]', 
#                 '[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):1]=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):2]>>[CH3][#16X2][NX3+:1]=[*:2]',
#                 '[#7X2-:1]>>[CH3][#16X2][#7X3+0:1]', 
#                 '[#6-;!$([#6X1-]#[#7,#8,#15,#16]):1]>>[CH3][#16X2][#6+0:1]', 
#                 '[#6:1]=[#7+:2](-[#8-:3])-[#8-:4]>>[CH3][#16X2][#6:1][#7+:2](=[#8+0:3])-[*:4]', 
#                 '[OX1:1]=[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]):2][#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][#16X2][OX2:1][#6X3+:2][*:3][*:4]', 
#                 '[OX1:1]=[CX3;$([R0][#6]),$([H1R0]):2][$([OX2H]),$([OX1-]):3]>>[CH3][#16X2][OX2:1][CX3+:2][*:3]',
#                 '[#7+0;$([N;R;!$([#7X2]);$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]),$([NX3+0;!$([#7X3][CX3;$([CX3][#6]),$([CX3H])]=[OX1])]),$([NX4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])]):1]>>[CH3][#16X2][#7+:1]', 
#                 '[C:1]=[C:2]=[#7X1-:3]>>[CH3][#16X2][C:1][C:2]#[#7X1+0:3]',
#                 '[NX1:1]#[CX2;!$(CC=C=[#7X1-]);!$(CC=C):2]>>[CH3][#16X2][NX2+:1]#[*:2]', 
#                 '[CX1-:1]#[NX2+:2]>>[CH3][#16X2][CX2+0:1]#[NX2+:2]',
#                 '[OX2H:1][$(c(c)c),$([#6X3;R](=[#6X3;R])[#6X3;R]):2]>>[CH3][#16X2][OX3+1:1][*:2]', # added due to rxn100
#                 '[#8X2H0:1][#14X4:2]([!#1:3])([!#1:4])[!#1:5]>>[CH3][#16X2][#8X3H0+:1][*:2]([*:3])([*:4])[*:5]', # added due to rxn100
#                 '[#7X2;$([nX2](:*):*),$([#7X2;R](=[*;R])[*;R]):1]>>[CH3][#16X2][#7X3+:1]', # added due to rxn100
#                 '[*-:1]>>[CH3][#16X2][*+0:1]', # added to capture additional sites
#                 '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+0:2]>>[CH3][#16X2][*:1]-[*+1:2]', # added to capture additional sites
#                 '[*;!$([!X4;!#1;!#6:1])+0:1]=[*+1:2]>>[CH3][#16X2][*:1]-[*+2:2]', # added to capture additional sites
#                 '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+0:2]>>[CH3][#16X2][*:1]=[*+1:2]', # added to capture additional sites
#                 '[*;!$([!X4;!#1;!#6:1])+0:1]#[*+1:2]>>[CH3][#16X2][*:1]=[*+2:2]', # added to capture additional sites
#                 '[!X4;!#1;!#6:1]>>[CH3][#16X2][*+1:1]' # added to capture additional sites
#                 ]
# ### END SMe ###


e_smirks_names_list = ["Oxonium", 
                      "Carbocation", 
                      "Ketone", 
                      "Amide", 
                      "Ester",
                      "Iminium", 
                      "Michael acceptor", 
                      "Imine", 
                      "Aldehyde",
                      "Anhydride", # added due to rxn100
                      "Acyl Halide", # added due to rxn100
                      'cation_with_charge_plus1', # added to capture additional sites
                      'double_bond', # added to capture additional sites
                      'double_bond_neighbouratom_with_charge_plus1', # added to capture additional sites
                      'triple_bond', # added to capture additional sites
                      'triple_bond_neighbouratom_with_charge_plus1' # added to capture additional sites
                      ]

e_smirks_list = ['[#6:1]=[O+;!$([O]~[!#6]);!$([S]*~[#7,#8,#15,#16]):2]>>[CH3][#6:1][O+0:2]', 
                '[#6+:1]>>[CH3][#6+0:1]', 
                '[#6X3:1](=[OX1:2])([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][#6X4:1](-[OX1-:2])([*:3])[*:4]', 
                '[CX3;$([CX3][#6]),$([CX3H]):1](=[OX1:2])[#7X3;!R:3]>>[CH3][CX4:1](-[OX1-:2])[*:3]',
                '[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]),$([#6X3][OX2H0]):1](=[OX1:2])[#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][#6X4:1](-[OX1-:2])[*:3][*:4]', 
                '[CX3:1]=[NX3+;!$(N([#8-])[#8-]):2]>>[CH3][CX4:1]-[NX3+0:2]', 
                '[CX3;!R:1]=[CX3:2][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]):3]>>[CH3][CX4:1]-[CX3-:2][*:3]',
                '[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):1]=[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):2]>>[CH3][CX4:1]-[NX2-:2]',
                '[CX3;$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):1]=[OX1:2]>>[CH3][CX4:1]-[OX1-:2]',
                '[CX3:1](=[OX1:2])[OX2:3][CX3:4]=[OX1:5]>>[CH3][CX4:1](-[OX1-:2])[*:3][*:4]=[*:5]', # added due to rxn100
                '[CX3:1](=[OX1:2])[ClX1,BrX1,IX1:3]>>[CH3][CX4:1](-[OX1-:2])[*:3]', # added due to rxn100
                '[*+:1]>>[CH3][*+0:1]', # added to capture additional sites
                '[*+0:1]=[*+0:2]>>[CH3][*:1]-[*-1:2]', # added to capture additional sites
                '[*+0:1]=[*+1:2]>>[CH3][*:1]-[*+0:2]', # added to capture additional sites
                '[*+0:1]#[*+0:2]>>[CH3][*:1]=[*-1:2]', # added to capture additional sites
                '[*+0:1]#[*+1:2]>>[CH3][*:1]=[*+0:2]' # added to capture additional sites
                ]
# ### BEGIN SMe ###
# e_smirks_list = ['[#6:1]=[O+;!$([O]~[!#6]);!$([S]*~[#7,#8,#15,#16]):2]>>[CH3][#16X2][#6:1][O+0:2]', 
#                 '[#6+:1]>>[CH3][#16X2][#6+0:1]', 
#                 '[#6X3:1](=[OX1:2])([#6;!$([CX3]=[CX3;!R]):3])[#6;!$([CX3]=[CX3;!R]):4]>>[CH3][#16X2][#6X4:1](-[OX1-:2])([*:3])[*:4]', 
#                 '[CX3;$([CX3][#6]),$([CX3H]):1](=[OX1:2])[#7X3;!R:3]>>[CH3][#16X2][CX4:1](-[OX1-:2])[*:3]',
#                 '[#6X3;!$([#6X3][CX3]=[CX3;!R]);$([#6X3][#6]),$([#6X3H]),$([#6X3][OX2H0]):1](=[OX1:2])[#8X2H0:3][#6;!$(C=[O,N,S]):4]>>[CH3][#16X2][#6X4:1](-[OX1-:2])[*:3][*:4]', 
#                 '[CX3:1]=[NX3+;!$(N([#8-])[#8-]):2]>>[CH3][#16X2][CX4:1]-[NX3+0:2]', 
#                 '[CX3;!R:1]=[CX3:2][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]):3]>>[CH3][#16X2][CX4:1]-[CX3-:2][*:3]',
#                 '[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6]):1]=[NX2;$([N][#6]),$([NH]);!$([N][CX3]=[#7,#8,#15,#16]):2]>>[CH3][#16X2][CX4:1]-[NX2-:2]',
#                 '[CX3;$([CX3H][#6;!$([CX3]=[CX3;!R])]),$([CX3H2]):1]=[OX1:2]>>[CH3][#16X2][CX4:1]-[OX1-:2]',
#                 '[CX3:1](=[OX1:2])[OX2:3][CX3:4]=[OX1:5]>>[CH3][#16X2][CX4:1](-[OX1-:2])[*:3][*:4]=[*:5]', # added due to rxn100
#                 '[CX3:1](=[OX1:2])[ClX1,BrX1,IX1:3]>>[CH3][#16X2][CX4:1](-[OX1-:2])[*:3]', # added due to rxn100
#                 '[*+:1]>>[CH3][#16X2][*+0:1]', # added to capture additional sites
#                 '[*+0:1]=[*+0:2]>>[CH3][#16X2][*:1]-[*-1:2]', # added to capture additional sites
#                 '[*+0:1]=[*+1:2]>>[CH3][#16X2][*:1]-[*+0:2]', # added to capture additional sites
#                 '[*+0:1]#[*+0:2]>>[CH3][#16X2][*:1]=[*-1:2]', # added to capture additional sites
#                 '[*+0:1]#[*+1:2]>>[CH3][#16X2][*:1]=[*+0:2]' # added to capture additional sites
#                 ]
# ### END SMe ###

def remove_identical_atoms(rdkit_mol, atom_list):
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(atom_list):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)
    
    atom_list = np.array(atom_list)[idx_list].tolist()
    return atom_list


def remove_identical_atoms_with_associated_list(rdkit_mol, atom_list, associated_list):
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(atom_list):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)
    
    atom_list = np.array(atom_list)[idx_list].tolist()
    associated_list = np.array(associated_list)[idx_list].tolist()
    return atom_list, associated_list


def find_identical_atoms(rdkit_mol, atom_list):
    len_list = len(atom_list)
    
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(rdkit_mol.GetAtoms()):
        if atom.GetIdx() in atom_list[:len_list]:
            sym_atoms = [int(atom_idx) for atom_idx, ranking in enumerate(atom_rank) if ranking == atom_rank[idx] and atom_idx not in atom_list] 
            atom_list.extend(sym_atoms)
    return atom_list


def find_identical_atoms_with_associated_list(rdkit_mol, atom_list, associated_list):
    len_list = len(atom_list)
    
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(rdkit_mol.GetAtoms()):
        if atom.GetIdx() in atom_list[:len_list]:
            sym_atoms = [int(atom_idx) for atom_idx, ranking in enumerate(atom_rank) if ranking == atom_rank[idx] and atom_idx not in atom_list] 
            atom_list.extend(sym_atoms)
            associated_list.extend([associated_list[atom_list[:len_list].index(atom.GetIdx())]]*len(sym_atoms))
    return atom_list, associated_list


def find_electrophilic_sites(rdkit_mol):
    elec_sites = []
    elec_names = []
    elec_smirks = []
    for i, smirks in enumerate(e_smirks_list):
        smarts = smirks.split('>>')[0]
        if rdkit_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            sites = [x[0] for x in rdkit_mol.GetSubstructMatches(Chem.MolFromSmarts(smarts), uniquify=False)]
            # sites = remove_identical_atoms(rdkit_mol, sites)
            for site in sites:
                if site not in elec_sites:
                    elec_sites.append(site)
                    elec_names.append(e_smirks_names_list[i])
                    elec_smirks.append(smirks)
    
    # Remove sites with same canonical rank
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(elec_sites):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)

    elec_sites = np.array(elec_sites)[idx_list].tolist()
    elec_names = np.array(elec_names)[idx_list].tolist()
    elec_smirks = np.array(elec_smirks)[idx_list].tolist()

    return elec_sites, elec_names, elec_smirks


def find_nucleophilic_sites(rdkit_mol):
    nuc_sites = []
    nuc_names = []
    nuc_smirks = []
    for i, smirks in enumerate(n_smirks_list):
        smarts = smirks.split('>>')[0]
        if rdkit_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            sites = [x[0] for x in rdkit_mol.GetSubstructMatches(Chem.MolFromSmarts(smarts), uniquify=False)]
            # sites = remove_identical_atoms(rdkit_mol, sites)
            for site in sites:
                if site not in nuc_sites:
                    nuc_sites.append(site)
                    nuc_names.append(n_smirks_names_list[i])
                    nuc_smirks.append(smirks)
    
    # Remove sites with same canonical rank
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(nuc_sites):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)

    nuc_sites = np.array(nuc_sites)[idx_list].tolist()
    nuc_names = np.array(nuc_names)[idx_list].tolist()
    nuc_smirks = np.array(nuc_smirks)[idx_list].tolist()

    return nuc_sites, nuc_names, nuc_smirks


def find_nucleophilic_sites_and_generate_MCAproducts(rdkit_mol):
    nuc_prods = []
    nuc_smis = []
    nuc_sites = []
    nuc_names = []
    for i, smirks in enumerate(n_smirks_list):

        product_mols, product_smis, sites = molfmt.run_rxn(Chem.AddHs(rdkit_mol), smirks)
            
        for product_mol, product_smi, site in zip(product_mols, product_smis, sites):
            # if product_smi not in nuc_smis:
            if site not in nuc_sites:
                nuc_prods.append(product_mol)
                nuc_smis.append(product_smi)
                nuc_sites.append(site)
                nuc_names.append(n_smirks_names_list[i])
    
    # Remove sites with same canonical rank
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(nuc_sites):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)

    nuc_prods = np.array(nuc_prods)[idx_list].tolist()
    nuc_smis = np.array(nuc_smis)[idx_list].tolist()
    nuc_sites = np.array(nuc_sites)[idx_list].tolist()
    nuc_names = np.array(nuc_names)[idx_list].tolist()

    return nuc_prods, nuc_smis, nuc_sites, nuc_names


def find_electrophilic_sites_and_generate_MAAproducts(rdkit_mol):
    elec_prods = []
    elec_smis = []
    elec_sites = []
    elec_names = []
    for i, smirks in enumerate(e_smirks_list):
            
        product_mols, product_smis, sites = molfmt.run_rxn(Chem.AddHs(rdkit_mol), smirks)
            
        for product_mol, product_smi, site in zip(product_mols, product_smis, sites):
            # if product_smi not in elec_smis:
            if site not in elec_sites:
                elec_prods.append(product_mol)
                elec_smis.append(product_smi)
                elec_sites.append(site)
                elec_names.append(e_smirks_names_list[i])
    
    # Remove sites with same canonical rank
    idx_list = []
    rank_kept = []
    atom_rank = list(Chem.CanonicalRankAtoms(rdkit_mol, breakTies=False))
    for idx, atom in enumerate(elec_sites):
        if atom_rank[atom] not in rank_kept:
            rank_kept.append(atom_rank[atom])
            idx_list.append(idx)

    elec_prods = np.array(elec_prods)[idx_list].tolist()
    elec_smis = np.array(elec_smis)[idx_list].tolist()
    elec_sites = np.array(elec_sites)[idx_list].tolist()
    elec_names = np.array(elec_names)[idx_list].tolist()

    return elec_prods, elec_smis, elec_sites, elec_names