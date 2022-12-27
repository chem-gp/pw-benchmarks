import os
import subprocess

# from ase.io import read, write

from ase.build import bulk
from ase.calculators.espresso import Espresso


import os
from datetime import datetime

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
                    'Si': 'Si.pbe-tm-new-gipaw-dc.UPF',
                    'C': 'C.pbe-tm-new-gipaw-dc.UPF',
                    'H': 'H.pbe-tm-new-gipaw-dc.UPF'
                    }


def run_pw_scf(calc_directory, structure, num_proc_pw):

    os.environ["ASE_ESPRESSO_COMMAND"] = "mpirun --oversubscribe -np " + str(num_proc_pw) + " pw.x -in PREFIX.pwi > PREFIX.pwo"
    os.environ["OMP_NUM_THREADS"] = "1,1"

    os.chdir(calc_directory)

    pseudo_dir = '../../Pseudopotentials'

    calc = Espresso(pseudopotentials=pseudopotentials, pseudo_dir = pseudo_dir,
    outdir = './outdir',  prefix = 'crystal', restart_mode = 'from_scratch',
                    tstress=True, tprnfor=True, nosym=True, 
                    ecutwfc=10, kpts=(1, 1, 1), 
                    ecutrho = 100,
                    occupations = 'smearing', smearing = 'gauss', degauss = 1.0e-2
                    )
    structure.calc = calc

    print("Starting qe calculation...")
    print(structure.get_potential_energy())
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
        job = 'nmr'
        prefix = 'crystal'
        tmp_dir = './outdir/'
        q_gipaw = 0.01
        spline_ps = .true.
        use_nmr_macroscopic_shape = .false.
    /
    """

    with open("espresso_gipaw.pwi", "w") as f:
        #    lines = ["Adding lines\n", "writing into it \n", "written successfully\n" ]
        f.writelines(gipaw_input)
        f.close()

    print("Starting qe-gipaw calculation...")

    # f = open("espresso_gipaw.pwo", "a")
    # subprocess.run(["mpirun --oversubscribe -np " + str(num_proc_gipaw) + " gipaw.x -in espresso_gipaw.pwi > espresso_gipaw.pwo"], shell=True)
    subprocess.run(["mpirun -np " + str(num_proc_gipaw) + " gipaw.x -in espresso_gipaw.pwi > espresso_gipaw.pwo"], shell=True)

    output_filename = "espresso_gipaw.pwo"
    output_file = os.path.abspath(output_filename)

    return output_file





# better adopt parser from cclib: 
# https://github.com/cclib/cclib/blob/master/cclib/parser/orcaparser.py#L1329
def parse_gipaw_output(gipaw_filename):
    raise NotImplementedError