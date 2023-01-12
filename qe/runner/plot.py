import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np

import os


def plot_chemical_shifts(crystal, save_plots_path="."):

    # os.chdir(save_plots_path)

    chemical_shifts = crystal.info['chemical_shifts_iso']
    symbols = crystal.get_chemical_symbols()

    chemical_shifts_Si = []
    chemical_shifts_N = []
    chemical_shifts_C = []
    chemical_shifts_H = []

    for i,shift in enumerate(chemical_shifts):   
        if symbols[i] == 'Si':
            chemical_shifts_Si.append(shift)

        if symbols[i] == 'C':
            chemical_shifts_C.append(shift)

        if symbols[i] == 'N':
            chemical_shifts_N.append(shift)
        
        if symbols[i] == 'H':
            chemical_shifts_H.append(shift)

    sns.set()
    sns.set_style("whitegrid", {"axes.grid": False})
    plt.figure(figsize=(10, 1))
    # plt.xticks(fontsize=18)
    # plt.yticks(fontsize=18)

    ax = sns.rugplot(chemical_shifts_H, color="darkgreen", linewidth=0.5, height=0.5)
    ax.set(yticklabels=[])
    ax.invert_xaxis()
    ax.set_title('1H chemical shifts')
    plt.savefig(save_plots_path + "/1H_chemical_shifts.png", bbox_inches='tight', dpi=400)
    plt.show()

    plt.figure(figsize=(10, 1))
    ax = sns.rugplot(chemical_shifts_C, color="darkgreen", linewidth=0.5, height=0.5)
    ax.set(yticklabels=[])
    ax.invert_xaxis()
    ax.set_title('13C chemical shifts')
    plt.savefig(save_plots_path + "/13C_chemical_shifts.png", bbox_inches='tight', dpi=400)
    plt.show()

    plt.figure(figsize=(10, 1))
    ax = sns.rugplot(chemical_shifts_N, color="darkgreen", linewidth=0.5, height=0.5)
    ax.set(yticklabels=[])
    ax.invert_xaxis()
    ax.set_title('15N chemical shifts')
    plt.savefig(save_plots_path + "/15N_chemical_shifts.png", bbox_inches='tight', dpi=400)
    plt.show()

    plt.figure(figsize=(10, 1))
    ax = sns.rugplot(chemical_shifts_Si, color="darkgreen", linewidth=0.5, height=0.5)
    ax.set(yticklabels=[])
    ax.invert_xaxis()
    ax.set_title('29Si chemical shifts')
    plt.savefig(save_plots_path + "/29Si_chemical_shifts.png", bbox_inches='tight', dpi=400)
    plt.show()

    chemical_shifts_H_np = np.array(chemical_shifts_H)
    chemical_shifts_C_np = np.array(chemical_shifts_C)

    chemical_shifts_np = {
        'H': chemical_shifts_H_np,
        'C': chemical_shifts_C_np
    }

    return chemical_shifts_np