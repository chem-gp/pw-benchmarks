import os
import subprocess

# from ase.io import read, write

from ase.build import bulk
from ase.calculators.espresso import Espresso


import os
from datetime import datetime

import numpy as np
import re

def make_calc_dir():
    dir_name = f'{datetime.now().strftime("%Y-%m-%d_%H %M-%S_%f")}'
    dir_name = 'test/' + dir_name
    os.makedirs(dir_name, exist_ok=True)
    abs_dir_path = os.path.abspath(dir_name)
    return abs_dir_path



# https://sites.google.com/site/dceresoli/pseudopotentials
pseudopotentials = {'Na': 'Na.pbe-tm-gipaw-dc.UPF',
                    'Cl': 'Cl.pbe-tm-gipaw.UPF',
                    'Cu': 'Cu.pbe-tm-new-gipaw.UPF',
                    'Sb': 'sb_pbe_v1.4.uspp.F.UPF', # copied from another folder
                    'N': 'N.pbe-tm-new-gipaw-dc.UPF',
                    'F': 'F.pbe-tm-new-gipaw-dc.UPF',
                    # 'Si': 'Si.pbe-tm-new-gipaw-dc.UPF',
                    # # 'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
                    # 'C': 'C.pbe-tm-new-gipaw-dc.UPF',
                    # 'H': 'H.pbe-tm-new-gipaw-dc.UPF',
                    'Si': 'Si.pbe-tm-gipaw.UPF',
                    # 'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
                    'C': 'C.pbe-tm-gipaw.UPF',
                    'H': 'H.pbe-tm-gipaw.UPF'
                    }


def run_pw_scf(calc_directory, structure, num_proc_pw, pw_params):

    print(pw_params)

    os.environ["ASE_ESPRESSO_COMMAND"] = "mpirun -np " + str(num_proc_pw) + " pw.x -in PREFIX.pwi > PREFIX.pwo"
    os.environ["OMP_NUM_THREADS"] = "1,1"
    os.environ["OMP_STACKSIZE"] = "80G"

    os.chdir(calc_directory)

    pseudo_dir = '../../Pseudopotentials'

    calc = Espresso(pseudopotentials=pseudopotentials, pseudo_dir = pseudo_dir,
    outdir = './outdir', **pw_params
                    # prefix = 'crystal', restart_mode = 'from_scratch',
                    # tstress=True, tprnfor=True, nosym=True, 
                    # ecutwfc=10, 
                    # # kpts=(1, 1, 1),
                    # kpts=None, 
                    # ecutrho = 100,
                    # occupations = 'smearing', smearing = 'gauss', degauss = 1.0e-2
                    )
    structure.calc = calc

    # print("Starting qe calculation...")
    print(structure.get_potential_energy())

    os.chdir('../..')
    # print(rocksalt.get_forces())
    return



def run_gipaw(calc_directory, num_proc_gipaw):
    # &inputgipaw
    #     job = 'nmr'
    #     prefix = 'TMS'
    #     tmp_dir = './outdir/'
    #     q_gipaw = 0.01
    #     spline_ps = .true.
    #     use_nmr_macroscopic_shape = .false.
    # /

    os.chdir(calc_directory)


    # num_proc_gipaw = 6
    gipaw_input = """&inputgipaw
        restart_mode = 'from_scratch'
        job = 'nmr'
        prefix = 'crystal'
        tmp_dir = './outdir/'
        diagonalization = 'cg'
        verbosity = 'high'
        q_gipaw = 0.01
        spline_ps = .true.
        use_nmr_macroscopic_shape = .false.
    /
    """
    # diagonalization = 'cg'
    # verbosity = 'high'

    with open("espresso_gipaw.pwi", "w") as f:
        #    lines = ["Adding lines\n", "writing into it \n", "written successfully\n" ]
        f.writelines(gipaw_input)
        f.close()

    print("Starting qe-gipaw calculation...")

    # f = open("espresso_gipaw.pwo", "a")
    # subprocess.run(["mpirun --oversubscribe -np " + str(num_proc_gipaw) + " gipaw.x -in espresso_gipaw.pwi > espresso_gipaw.pwo"], shell=True)
    # subprocess.run(["conda info"], shell=True)

    subprocess.run(["mpirun -np " + str(num_proc_gipaw) + " gipaw.x -in espresso_gipaw.pwi > espresso_gipaw.pwo"], shell=True)

    # p = subprocess.Popen(['mpirun -np', str(num_proc_gipaw), " gipaw.x -in espresso_gipaw.pwi > espresso_gipaw.pwo"],
    #                     cwd=".",
    #                     stdout=subprocess.PIPE,
    #                     stderr=subprocess.STDOUT)

    output_filename = "espresso_gipaw.pwo"
    output_file = os.path.abspath(output_filename)

    os.chdir('../..')

    return output_file





def parse_gipaw_output(gipaw_file, num_atoms):

    file_gipaw = open(gipaw_file, 'r')
    lines = file_gipaw.readlines()

    chem_shifts_tensors = np.full((num_atoms, 3, 3), np.inf)
    chem_shifts_isotropic = np.zeros(num_atoms)


    for i,line in enumerate(lines):   
        if line.startswith('     Total NMR chemical shifts in ppm:'):  

            for j in range(num_atoms):
                atom_line = lines[i+10*j+3]

                atom_line_numbers = re.findall('-?[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+', atom_line)
                atom_index = int(atom_line_numbers[0])
                atom_coords = [ float(a) for a in atom_line_numbers[1:4] ]
                isotropic_shift = float(atom_line_numbers[4])

                chem_shifts_isotropic[j] = float(isotropic_shift)

                # print(atom_index, atom_coords, isotropic_shift)
                
                tensor_row_1 = lines[i+10*j+4]
                tensor_row_1_numbers = re.findall('-?[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+', tensor_row_1)
                # print(tensor_row_1_numbers)
                chem_shifts_tensors[j, 0, : ] = np.array(tensor_row_1_numbers)

                tensor_row_2 = lines[i+10*j+5]
                tensor_row_2_numbers = re.findall('-?[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+', tensor_row_2)
                chem_shifts_tensors[j, 1, : ] = np.array(tensor_row_2_numbers)

                tensor_row_3 = lines[i+10*j+6]
                tensor_row_3_numbers = re.findall('-?[\d]+[.,\d]+|[\d]*[.][\d]+|[\d]+', tensor_row_3)
                chem_shifts_tensors[j, 2, : ] = np.array(tensor_row_3_numbers)

    # print(chem_shifts_isotropic)    
    # print(chem_shifts_tensors)

    return chem_shifts_isotropic, chem_shifts_tensors