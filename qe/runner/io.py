# import seaborn as sns
# import matplotlib.pyplot as plt
import numpy as np

from ase.io import read, write

import os


def get_chemical_shifts_from_xyz(struct):

    # os.chdir(save_plots_path)

    if isinstance(struct, str):
        struct = read(struct, format='extxyz')


    chemical_shifts = struct.info['chemical_shifts_iso']
    symbols = struct.get_chemical_symbols()

    chemical_shifts_Si = []
    chemical_shifts_N = []
    chemical_shifts_C = []
    chemical_shifts_H = []
    chemical_shifts_O = []

    indices_Si = []
    indices_N = []
    indices_C = []
    indices_H = []
    indices_O = []

    for i,shift in enumerate(chemical_shifts):   
        if symbols[i] == 'Si':
            chemical_shifts_Si.append(shift)
            indices_Si.append(i)

        if symbols[i] == 'C':
            chemical_shifts_C.append(shift)
            indices_C.append(i)

        if symbols[i] == 'N':
            chemical_shifts_N.append(shift)
            indices_N.append(i)

        if symbols[i] == 'O':
            chemical_shifts_O.append(shift)
            indices_O.append(i)
        
        if symbols[i] == 'H':
            chemical_shifts_H.append(shift)
            indices_H.append(i)


    chemical_shifts_H_np = np.array(chemical_shifts_H)
    chemical_shifts_C_np = np.array(chemical_shifts_C)
    chemical_shifts_N_np = np.array(chemical_shifts_N)
    chemical_shifts_O_np = np.array(chemical_shifts_O)
    chemical_shifts_Si_np = np.array(chemical_shifts_Si)


    chemical_shifts_np = {
        'H': chemical_shifts_H_np,
        'C': chemical_shifts_C_np,
        'N': chemical_shifts_N_np,
        'O': chemical_shifts_O_np,
        'Si': chemical_shifts_Si_np,
        'H_indices': indices_H,
        'C_indices': indices_C,
        'N_indices': indices_N,
        'O_indices': indices_O,
        'Si_indices': indices_Si
    }

    return chemical_shifts_np