# nohup python chem_shifts.py > "test/$(date +"%Y_%m_%d_%I_%M_%p").log" 2>&1 &
# python chem_shifts.py > "test/$(date +"%Y_%m_%d_%I_%M_%p").log" 2>&1
# python chem_shifts.py > "test/$(date +"%Y_%m_%d_%I_%M_%p").log"


# python ./chem_shifts.py

# pkill -9 -f pw.x; pkill -9 -f gipaw.x

from runner import pw_runner
from ase.io import read, write

import time

import os

# Defining parameters:

np_pw = 10
np_gipaw = 10

pw_params = {
    'prefix':'crystal', 
    'restart_mode' : 'from_scratch',
    'tstress':True, 
    'tprnfor':True, 
    'nosym':True, 
    'ecutwfc':60, 
    'kpts':(4, 2, 2),
    # 'kpts':None, 
    'ecutrho' : 240,
    # 'occupations' : 'smearing', 
    # 'smearing' : 'gauss', 
    # 'degauss' : 1.0e-2
}

# Loading structures:

# crystal=read("../structures/MIN-167-350K-CuCbPyz.cif")
# crystal=read("../structures/HIK-143 293K-activated.cif")
# crystal=read("../structures/HIK-143 MeOH.cif")
# crystal=read("../structures/KTU-183_2_auto.cif")
crystal_264=read("../structures/new_systems/ktu_002.cif")
crystal_258=read("../structures/new_systems/KTU-065-25do-b.cif")
crystal_612=read("../structures/new_systems/KTU-TBS-20do_auto.cif")

crystal = crystal_258
tms = read("../structures/tms.xyz")
# crystal = tms
a = 15.0
tms.set_cell([a, a, a])
tms.set_cell([(a, 0, 0), (0, a, 0), (0, 0, a)])
tms.set_pbc(True)
# Creating directory to store all the data:
calc_dir = pw_runner.make_calc_dir()
print("Starting QE/SCF calculation for TMS...")
pw_runner.run_pw_scf(calc_dir, tms, num_proc_pw=np_pw, pw_params=pw_params)
print("Starting gipaw calculation for TMS...")
gipaw_out = pw_runner.run_gipaw(calc_dir, "espresso_gipaw_tms.pwo", num_proc_gipaw=np_gipaw)
chemical_shifts_iso, chemical_shifts_tensors = pw_runner.parse_gipaw_output(calc_dir + '/espresso_gipaw_tms.pwo', num_atoms=len(tms))
tms.info['chemical_shifts_iso'] = chemical_shifts_iso
tms.info['chemical_shifts_tensors'] = chemical_shifts_tensors
write(calc_dir+"/tms.xyz", tms, format='extxyz')
os.rename(calc_dir+"/espresso.pwo", calc_dir+"/espresso_tms.pwo")
print("Starting QE/SCF calculation for crystal...")
pw_runner.run_pw_scf(calc_dir, crystal, num_proc_pw=np_pw, pw_params=pw_params)
print("Starting gipaw calculation for crystal...")
gipaw_out = pw_runner.run_gipaw(calc_dir, "espresso_gipaw.pwo", num_proc_gipaw=np_gipaw)
chemical_shifts_iso, chemical_shifts_tensors = pw_runner.parse_gipaw_output(calc_dir + '/espresso_gipaw.pwo', num_atoms=len(crystal))
crystal.info['chemical_shifts_iso'] = chemical_shifts_iso
crystal.info['chemical_shifts_tensors'] = chemical_shifts_tensors
write(calc_dir+"/crystal.xyz", crystal, format='extxyz')

crystal = crystal_264
tms = read("../structures/tms.xyz")
# crystal = tms
a = 15.0
tms.set_cell([a, a, a])
tms.set_cell([(a, 0, 0), (0, a, 0), (0, 0, a)])
tms.set_pbc(True)
# Creating directory to store all the data:
calc_dir = pw_runner.make_calc_dir()
print("Starting QE/SCF calculation for TMS...")
pw_runner.run_pw_scf(calc_dir, tms, num_proc_pw=np_pw, pw_params=pw_params)
print("Starting GIPAW calculation for TMS...")
gipaw_out = pw_runner.run_gipaw(calc_dir, "espresso_gipaw_tms.pwo", num_proc_gipaw=np_gipaw)
chemical_shifts_iso, chemical_shifts_tensors = pw_runner.parse_gipaw_output(calc_dir + '/espresso_gipaw_tms.pwo', num_atoms=len(tms))
tms.info['chemical_shifts_iso'] = chemical_shifts_iso
tms.info['chemical_shifts_tensors'] = chemical_shifts_tensors
write(calc_dir+"/tms.xyz", tms, format='extxyz')
os.rename(calc_dir+"/espresso.pwo", calc_dir+"/espresso_tms.pwo")
print("Starting QE/SCF calculation for crystal...")
pw_runner.run_pw_scf(calc_dir, crystal, num_proc_pw=np_pw, pw_params=pw_params)
print("Starting GIPAW calculation for crystal...")
gipaw_out = pw_runner.run_gipaw(calc_dir, "espresso_gipaw.pwo", num_proc_gipaw=np_gipaw)
chemical_shifts_iso, chemical_shifts_tensors = pw_runner.parse_gipaw_output(calc_dir + '/espresso_gipaw.pwo', num_atoms=len(crystal))
crystal.info['chemical_shifts_iso'] = chemical_shifts_iso
crystal.info['chemical_shifts_tensors'] = chemical_shifts_tensors
write(calc_dir+"/crystal.xyz", crystal, format='extxyz')

crystal = crystal_612
tms = read("../structures/tms.xyz")
# crystal = tms
a = 15.0
tms.set_cell([a, a, a])
tms.set_cell([(a, 0, 0), (0, a, 0), (0, 0, a)])
tms.set_pbc(True)
# Creating directory to store all the data:
calc_dir = pw_runner.make_calc_dir()
print("Starting QE/SCF calculation for TMS...")
pw_runner.run_pw_scf(calc_dir, tms, num_proc_pw=np_pw, pw_params=pw_params)
print("Starting gipaw calculation for TMS...")
gipaw_out = pw_runner.run_gipaw(calc_dir, "espresso_gipaw_tms.pwo", num_proc_gipaw=np_gipaw)
chemical_shifts_iso, chemical_shifts_tensors = pw_runner.parse_gipaw_output(calc_dir + '/espresso_gipaw_tms.pwo', num_atoms=len(tms))
tms.info['chemical_shifts_iso'] = chemical_shifts_iso
tms.info['chemical_shifts_tensors'] = chemical_shifts_tensors
write(calc_dir+"/tms.xyz", tms, format='extxyz')
os.rename(calc_dir+"/espresso.pwo", calc_dir+"/espresso_tms.pwo")
print("Starting QE/SCF calculation for crystal...")
pw_runner.run_pw_scf(calc_dir, crystal, num_proc_pw=np_pw, pw_params=pw_params)
print("Starting GIPAW calculation for crystal...")
gipaw_out = pw_runner.run_gipaw(calc_dir, "espresso_gipaw.pwo", num_proc_gipaw=np_gipaw)
chemical_shifts_iso, chemical_shifts_tensors = pw_runner.parse_gipaw_output(calc_dir + '/espresso_gipaw.pwo', num_atoms=len(crystal))
crystal.info['chemical_shifts_iso'] = chemical_shifts_iso
crystal.info['chemical_shifts_tensors'] = chemical_shifts_tensors
write(calc_dir+"/crystal.xyz", crystal, format='extxyz')