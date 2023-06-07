import os
import argparse
import numpy as np
import subprocess

def parse_args():
    """
    Argument parser so this can be run from the command line
    """
    parser = argparse.ArgumentParser(description='Run orca input creator')
    parser.add_argument('-n', '--name', default='test_mol', help='The name of the molecule')
    parser.add_argument('-c', '--charge', default=0, help='The charge of the molecule')
    parser.add_argument('-s', '--spin', default=0, help='The spin of the molecule')
    return parser.parse_args()


def run_orca(name, xyz_file, chrg, spin, path, ncores=12, mem=10000, optimize=True):
    
    input_file = np.loadtxt(fname=f'{path}/{xyz_file}', skiprows=2, dtype=str, ndmin=2)
    xyz_coords = []
    for l in input_file:
        xyz_coords.append(f'{l[0]:4} {float(l[1]):11.6f} {float(l[2]):11.6f} {float(l[3]):11.6f}\n')

    input_file = open(f"{path}/{name}.inp", "w")
    if optimize:
        input_file.write(f'# ORCA input file \n! OPT NumFreq r2SCAN-3c CPCM \n\n%maxcore {(mem*0.75)/ncores} \n%pal nprocs {ncores} end \n%cpcm smd true SMDsolvent "DMSO" end \n\n*xyz {chrg} {int(2*spin+1)}\n{"".join(xyz_coords)}*')
        # input_file.write(f'# ORCA input file \n! OPT Freq r2SCAN-3c CPCM \n\n%maxcore {(mem*0.75)/ncores} \n%pal nprocs {ncores} end \n%cpcm smd true SMDsolvent "DMSO" end \n\n*xyz {chrg} {int(2*spin+1)}\n{"".join(xyz_coords)}*')
        # input_file.write(f'# ORCA input file \n! OPT r2SCAN-3c CPCM \n\n%maxcore {(mem*0.75)/ncores} \n%pal nprocs {ncores} end \n%cpcm smd true SMDsolvent "DMSO" end \n\n*xyz {chrg} {int(2*spin+1)}\n{"".join(xyz_coords)}*')
    else:
        input_file.write(f'# ORCA input file \n! SP r2SCAN-3c CPCM \n\n%maxcore {(mem*0.75)/ncores} \n%pal nprocs {ncores} end \n%cpcm smd true SMDsolvent "DMSO" end \n\n*xyz {chrg} {int(2*spin+1)}\n{"".join(xyz_coords)}*')
    input_file.close()

    # Run ORCA calc
    cmd = f'env - PATH="/software/kemi/Orca/orca_5_0_1_linux_x86-64_openmpi411:/software/kemi/openmpi/openmpi-4.1.1/bin:$PATH" LD_LIBRARY_PATH="/software/kemi/openmpi/openmpi-4.1.1/lib:$LD_LIBRARY_PATH" /bin/bash -c "/software/kemi/Orca/orca_5_0_1_linux_x86-64_openmpi411/orca {path}/{name}.inp"'
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=f'{path}')
    output = proc.communicate()[0]

    # Save calc output
    with open(f'{path}/{name}.out', 'w') as f:
        f.write(output)

    # Clean up
    file_remove_list = ['charges', 'coordprot.0', 'lmocent.coord', f'{name}_atom46.densities',
                        f'{name}_atom46.out', f'{name}_atom46_property.txt', f'{name}_atom53.densities',
                        f'{name}_atom53.out', f'{name}_atom53_property.txt', f'{name}.cpcm',
                        f'{name}.densities', f'{name}.gbw', 'wbo']
    for file_remove in file_remove_list:
        if os.path.isfile(f'{path}/{file_remove}'):
            os.remove(f'{path}/{file_remove}')

    return

if __name__ == "__main__":

    import submitit
    args = parse_args()
    
    # Slurm settings
    executor = submitit.AutoExecutor(folder="submitit_orca")
    executor.update_parameters(
        name=args.name,
        cpus_per_task=24,
        mem_gb=20,
        timeout_min=180,
        slurm_partition="kemi1",
        slurm_array_parallelism=2,
    )
    print(executor)

    jobs = []
    with executor.batch():
        job = executor.submit(run_orca, f'{args.name}', f'{args.name}.xyz', int(args.charge), args.spin, os.getcwd(), True)
        jobs.append(job)
