
from runner import pw_runner
from ase.io import read, write

import time


scf_times = []
gipaw_times = []

for i in range(13,14):

    np_pw = i+1
    np_gipaw = i+1

    pw_params = {
        'prefix':'crystal', 
        'restart_mode' : 'from_scratch',
        'tstress':True, 
        'tprnfor':True, 
        'nosym':True, 
        'ecutwfc':10, 
        'kpts':(1, 1, 1),
        # 'kpts':None, 
        'ecutrho' : 100,
        # 'occupations' : 'smearing', 
        # 'smearing' : 'gauss', 
        # 'degauss' : 1.0e-2,
        # 'spline_ps': True
    }

    tms = read("test/tms.xyz")
    # crystal=read("test/MIN-167-350K-CuCbPyz.cif")
    crystal=read("test/HIK-143 293K-activated.cif")
    # crystal=read("test/HIK-143 MeOH.cif")
    # crystal=read("test/KTU-183_2_auto.cif")
    # a = 19.0
    # crystal = tms
    # crystal.set_cell([a, a, a])
    # crystal.set_cell([(a, 0, 0), (0, a, 0), (0, 0, a)])
    # crystal.set_pbc(True)

    calc_dir = pw_runner.make_calc_dir()

    print("Starting QE/SCF calculation...")
    start_time = time.time()
    pw_runner.run_pw_scf(calc_dir, crystal, num_proc_pw=np_pw, pw_params=pw_params)
    scf_time = time.time() - start_time
    print("QE/SCF calculation finished in ", scf_time, " seconds")

    print("Starting gipaw calculation...")
    gipaw_out = pw_runner.run_gipaw(calc_dir, num_proc_gipaw=np_gipaw)
    gipaw_time = time.time() - scf_time - start_time
    print("QE/GIPAW calculation finished in ", gipaw_time, " seconds")

    chemical_shifts_iso, chemical_shifts_tensors = pw_runner.parse_gipaw_output(calc_dir + '/espresso_gipaw.pwo', num_atoms=len(crystal))

    crystal.info['chemical_shifts_iso'] = chemical_shifts_iso
    crystal.info['chemical_shifts_tensors'] = chemical_shifts_tensors

    write(calc_dir+"/crystal.xyz", crystal, format='extxyz')

    scf_times.append(scf_time)
    gipaw_times.append(gipaw_time)


# print("SCF time: ", scf_time)
# print("GIPAW time: ", gipaw_time)

# import matplotlib.pyplot as plt
# import numpy as np

# x = np.arange(1, len(scf_times)+1, 1, dtype=int)
# plt.plot(x, scf_times, label='SCF time')
# plt.plot(x, gipaw_times, label='GIPAW time')
# plt.xticks(x)
# plt.xlabel('Num. of MPI processes')
# plt.ylabel('Time, s.')
# plt.legend()
# plt.show()