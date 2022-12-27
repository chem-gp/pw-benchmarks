import os
import subprocess


def run_pw_scf(calc_directory):
    raise NotImplementedError


def run_gipaw(calc_directory):
    # &inputgipaw
    #     job = 'nmr'
    #     prefix = 'TMS'
    #     tmp_dir = './outdir/'
    #     q_gipaw = 0.01
    #     spline_ps = .true.
    #     use_nmr_macroscopic_shape = .false.
    # /

    os.chdir(calc_directory)


    num_proc_gipaw = 6
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

    # f = open("espresso_gipaw.pwo", "a")
    # subprocess.run(["mpirun --oversubscribe -np " + str(num_proc_gipaw) + " gipaw.x -in espresso_gipaw.pwi > espresso_gipaw.pwo"], shell=True)
    subprocess.run(["mpirun -np " + str(num_proc_gipaw) + " gipaw.x -in espresso_gipaw.pwi > espresso_gipaw.pwo"], shell=True)



    return


def parse_gipaw_output(gipaw_filename):
    raise NotImplementedError