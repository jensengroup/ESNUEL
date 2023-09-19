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
import numpy as np
import subprocess
from ase.units import Hartree, mol, kJ


def run_orca(xyz_file, chrg, spin, path, ncores=2, mem=10000, optimize=False):

    input_file = np.loadtxt(fname=f'{path}/{xyz_file}', skiprows=2, dtype=str, ndmin=2)
    xyz_coords = []
    for l in input_file:
        xyz_coords.append(f'{l[0]:4} {float(l[1]):11.6f} {float(l[2]):11.6f} {float(l[3]):11.6f}\n')

    input_file = open(f"{path}/orca_calc.inp", "w")
    if optimize:
        # input_file.write(f'# ORCA input file \n! OPT NumFreq r2SCAN-3c CPCM \n\n%maxcore {(mem*0.75)/ncores} \n%pal nprocs {ncores} end \n%cpcm smd true SMDsolvent "DMSO" end \n\n*xyz {chrg} {int(2*spin+1)}\n{"".join(xyz_coords)}*')
        # input_file.write(f'# ORCA input file \n! OPT Freq r2SCAN-3c CPCM \n\n%maxcore {(mem*0.75)/ncores} \n%pal nprocs {ncores} end \n%cpcm smd true SMDsolvent "DMSO" end \n\n*xyz {chrg} {int(2*spin+1)}\n{"".join(xyz_coords)}*')
        input_file.write(f'# ORCA input file \n! OPT r2SCAN-3c CPCM \n\n%maxcore {(mem*0.75)/ncores} \n%pal nprocs {ncores} end \n%cpcm smd true SMDsolvent "DMSO" end \n\n*xyz {chrg} {int(2*spin+1)}\n{"".join(xyz_coords)}*')
    else:
        input_file.write(f'# ORCA input file \n! SP r2SCAN-3c CPCM \n\n%maxcore {(mem*0.75)/ncores} \n%pal nprocs {ncores} end \n%cpcm smd true SMDsolvent "DMSO" end \n\n*xyz {chrg} {int(2*spin+1)}\n{"".join(xyz_coords)}*')
    input_file.close()

    # Run ORCA calc --> (OBS! THE FOLLOWING PATHS SHOULD POINT TO ORCA)
    cmd = f'env - PATH="/software/kemi/Orca/orca_5_0_1_linux_x86-64_openmpi411:/software/kemi/openmpi/openmpi-4.1.1/bin:$PATH" LD_LIBRARY_PATH="/software/kemi/openmpi/openmpi-4.1.1/lib:$LD_LIBRARY_PATH" /bin/bash -c "/software/kemi/Orca/orca_5_0_1_linux_x86-64_openmpi411/orca {path}/orca_calc.inp"'
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=f'{path}')
    output = proc.communicate()[0]

    # Save calc output
    with open(f'{path}/orca_calc.out', 'w') as f:
        f.write(output)

    for i, line in enumerate(output.split('\n')):
        if "Check your MOs and check whether a frozen core calculation is appropriate" in line:
            energy = 60000.0
            break
        if "**** WARNING: LOEWDIN FINDS" in line:
            energy = 60000.0
            break
        if "**** WARNING: MULLIKEN FINDS" in line:
            energy = 60000.0
            break
        if "FINAL SINGLE POINT ENERGY" in line and "(Wavefunction not fully converged!)" not in line:
            energy = float(line.split()[-1]) * Hartree * mol/kJ #convert energy from Hartree to kJ/mol
            break
        if "Final Gibbs free energy" in line and "(Wavefunction not fully converged!)" not in line:
            energy = float(line.split()[-2]) * Hartree * mol/kJ #convert energy from Hartree to kJ/mol
            break
        else:
            energy = 60000.0

    return energy