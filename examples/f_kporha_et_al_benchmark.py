from libRustBCA import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'/../scripts')
sys.path.append('scripts')
import time
from materials import *
from tomlkit import parse, dumps

def input_file(ion, target, incident_energy, angle, number_ions=1000, pot="WW", Es=None):

    if ion == tungsten:
        ion_interaction_index = 1
    else:
        ion_interaction_index = 0

    mfp = (target["n"]/10**30)**(-1./3.)

    if not Es:
        Es = target["Es"]

    if pot == "WW":
        pot_name = "WW"
    else:
        pot_name = "MORSE"

    cpr = {'CPR': {'n0': 2, 'nmax': 32, 'epsilon': 5e-4, 'complex_threshold': 1E-9, 'far_from_zero': 1e3, 'interval_limit': 1E-4, 'derivative_free': True}}
    options = {
        'name': f'input_file_{ion["symbol"]}_{target["symbol"]}_{np.round(angle, 1)}_{np.round(incident_energy/1000, 4)}_{pot_name}_{Es}',
        'track_trajectories': False, # whether to track trajectories for plotting; memory intensive
        'track_recoils': True, # whether to track recoils; must enable for sputtering
        'track_recoil_trajectories': False, # whether to track recoil trajectories for plotting
        'track_displacements': False, # whether to track collisions with T > Ed for each species
        'track_energy_losses': False, # whether to track detailed collision energies; memory intensive
        'write_buffer_size': 8192, # how big the buffer is for file writing
        'weak_collision_order': 0, # weak collisions at radii (k + 1)*r; enable only when required
        'suppress_deep_recoils': False, # suppress recoils too deep to ever sputter
        'high_energy_free_flight_paths': False, # SRIM-style high energy free flight distances; use with caution
        'num_threads': 6, # number of threads to run in parallel
        'num_chunks': 10, # code will write to file every nth chunk; for very large simulations, increase num_chunks
        'electronic_stopping_mode': 'LOW_ENERGY_NONLOCAL',
        'mean_free_path_model': 'LIQUID', # liquid is amorphous (constant mean free path); gas is exponentially-distributed mean free paths
        'interaction_potential': [['ZBL', 'ZBL'],
                                  ['ZBL', pot]],
        'scattering_integral': [
            [{'GAUSS_MEHLER': {'n_points': 5}}, {'GAUSS_MEHLER': {'n_points': 5}}],
            [{'GAUSS_MEHLER': {'n_points': 5}}, {'GAUSS_MEHLER': {'n_points': 5}}],
        ],

        'root_finder': [
            ["DEFAULTNEWTON", "DEFAULTNEWTON"],
            ["DEFAULTNEWTON", cpr]
        ],
        'seed': 0 # if <0, will generate a seed from thread-local PRNG; if >0, will be used as seed to PRNG
    }

    # material parameters are per-species
    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        # bulk binding energy; typically zero as a model choice
        'Eb': [
            target["Eb"],
        ],
        # surface binding energy
        'Es': [
            Es
        ],
        # cutoff energy - particles with E < Ec stop
        'Ec': [
            target["Ec"],
        ],
        # atomic number
        'Z': [
            target["Z"],
        ],
        # atomic mass
        'm': [
            target["m"],
        ],
        # used to pick interaction potential from matrix in [options]
        'interaction_index': [1],
        'surface_binding_model': {
            "PLANAR": {'calculation': "INDIVIDUAL"}
        },
        'bulk_binding_model': 'INDIVIDUAL'
    }

    particle_parameters = {
        'interaction_index': [ion_interaction_index],
        'length_unit': 'ANGSTROM',
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        # number of computational ions of this species to run at this energy
        'N': [number_ions],
        # atomic mass
        'm': [ion["m"]],
        # atomic number
        'Z': [ion["Z"]],
        # incident energy 
        'E': [incident_energy],
        # cutoff energy - if E < Ec, particle stops
        'Ec': [ion["Ec"]],
        # surface binding energy
        'Es': [ion["Es"]],
        # initial position - if Es significant and E low, start (n)^(-1/3) above surface
        # otherwise 0, 0, 0 is fine; most geometry modes have surface at x=0 with target x>0
        'pos': [[-7.5, 0.0, 0.0]],
        # initial direction unit vector; most geometry modes have x-axis into the surface
        'dir': [
            [
                np.cos(angle*np.pi/180.0),
                np.sin(angle*np.pi/180.0),
                0.0
            ]
        ],
    }

    geometry_0D = {
        'length_unit': 'ANGSTROM',
        # used to correct nonlocal stopping for known compound discrpancies
        'electronic_stopping_correction_factor': 0.0,
        # number densities of each species; toml densities with length_unit ANGSTROM are in 1/A^3
        'densities': [ target["n"]/10**30 ]
    }

    input_data = {
        'options': options,
        'material_parameters': material_parameters,
        'particle_parameters': particle_parameters,
        'geometry_input': geometry_0D
    }

    return input_data


def tungsten_sputtering(ion, energy, angle, num_ions=1000, run_sim=False, pot="WW", Es=None):
    
    input = input_file(ion, tungsten, energy, angle, number_ions=num_ions, pot=pot, Es=Es)
    if run_sim: rustbca_py(input, geometry_mode="0D")
    sputtered = np.atleast_2d(np.genfromtxt(f'{input["options"]["name"]}sputtered.output', delimiter=','))
    if np.size(sputtered) > 0:
        num_sputtered = np.shape(sputtered)[0]
    else:
        num_sputtered = 0
    return num_sputtered/num_ions

# From F. Kphorha et al.
# picked a representative surface at random
# digitized with WebPlotDigitizer
w_tot= np.array([
    [-0.11725293132328218, 0.12903225806451624],
    [14.891122278056951, 0.2129032258064516],
    [30.01675041876047, 0.5032258064516131],
    [44.79061976549414, 0.9483870967741936],
    [59.916247906197654, 1.2516129032258065],
    [74.69011725293133, 1.1677419354838712],
])
w_sp = np.array([
    [-0.23450586264656437, 0.1322580645161291],
    [14.891122278056951, 0.19354838709677402],
    [29.782244556113906, 0.42258064516129035],
    [44.79061976549414, 0.7161290322580646],
    [59.916247906197654, 0.7516129032258064],
    [74.69011725293133, 0.1870967741935483],
])
ne = np.array([
    [-0.11725293132328218, 0.3225806451612905],
    [15.008375209380233, 0.3193548387096774],
    [30.134003350083756, 0.3322580645161288],
    [45.0251256281407, 0.38064516129032255],
    [60.03350083752095, 0.2193548387096773],
    [75.15912897822446, -0.003225806451613078],
])
ar = np.array([
    [-0.4690117252931323, 0.42903225806451606],
    [14.773869346733669, 0.44838709677419364],
    [29.782244556113906, 0.5064516129032259],
    [45.0251256281407, 0.5806451612903225],
    [60.26800670016752, 0.15483870967741908],
    [75.27638190954774, -0.003225806451613078],
])
energy = 200

ions = [tungsten, neon, argon]
datasets = [w_sp, ne, ar]
colors = []

num_ions = 10000
run_sim = False

tungsten["Es"] = 8.79
tungsten["Eb"] = 0.0
tungsten["Ec"] = 8.79
pot_morse ={'MORSE': {'D': 2.87454*1.602e-19, 'alpha': 13.3682*1e9, 'r0': 0.238631*1e-9}}

angles = np.linspace(0.0, 89.9, 12)
for ion, dataset in zip(ions, datasets):

    # keep computational cost to a minimum 
    # probably suppresses R_N
    ion["Ec"] = 8.0

    print(f'Running {ion["name"]}...')

    Y_RustBCA_WW = [tungsten_sputtering(ion, energy, angle, num_ions, run_sim) for angle in angles]
    Y_RustBCA_Morse = [tungsten_sputtering(ion, energy, angle, num_ions, run_sim,  pot=pot_morse) for angle in angles]
    Y_RustBCA_Morse_1175_eV = [tungsten_sputtering(ion, energy, angle, num_ions, run_sim,  pot=pot_morse, Es=11.75) for angle in angles]
    Y_RustBCA_KRC = [sputtering_yield(ion, tungsten, energy, angle, num_ions) for angle in angles]

    color = plt.plot(angles, Y_RustBCA_WW, label=f'RustBCA {ion["symbol"]} on W - ZBL/WW', linestyle='--')[0].get_color()
    colors.append(color)
    plt.plot(angles, Y_RustBCA_KRC, label=f'RustBCA {ion["symbol"]} on W - Kr-C', linestyle=':', color=color)

    plt.plot(angles, Y_RustBCA_Morse, label=f'RustBCA {ion["symbol"]} on W - ZBL/Morse', linestyle='-.', color=color)
    plt.plot(angles, Y_RustBCA_Morse_1175_eV, label=f'RustBCA {ion["symbol"]} on W - ZBL/Morse (Es=11.75 eV)', linestyle='-.', marker='o', color=color)
    plt.plot(dataset[:, 0], dataset[:, 1], label=f'MD {ion["symbol"]} on W', color=color)

plt.title('Comparison to F Kporha et al., W Sputtering Yields')
plt.scatter([0.0, 0.0], [0.13, 0.2], label='Exp. Ne on W (Laegreid et al., Stuart et al.)', marker='*', s=100, color=colors[1])
plt.scatter([0.0, 0.0, 0.0], [0.29, 0.6, 0.32], label='Exp. Ar on W (Laegreid et al., Stuart et al., Somogyvári et al.)', marker='*', s=100, color=colors[2])
plt.scatter([0.0,], [0.123], label='Exp. W on W', marker='*', s=100, color=colors[0])
plt.gca().set_ylim([0.0, 1.4])
plt.legend()
plt.xlabel('E [eV]')
plt.ylabel('Y [at/ion]')

# Produced by Somogyvari et al. (2012)
# digitized with WebPlotDigitizer
data_Somogyvari_Ar_W = np.array([
    [39.95309128868367, 0.00009332543007969924],
    [39.668850177558554, 0.0001621810097358933],
    [54.3104899786253, 0.004365158322401661],
    [79.29116634226808, 0.046773514128719856],
    [104.00501083088497, 0.10964781961431856],
    [156.2422374575878, 0.20892961308540398],
    [202.03488119514256, 0.32359365692962827],
    [255.71248811339458, 0.35481338923357547],
    [303.5082793629512, 0.4365158322401659],
])

plt.figure()
plt.title('Comparison to Somogyvari et al.: Ar on W sputtering')
energies = np.logspace(np.log10(30), np.log10(300), 12)
rustbca_ar_w_krc = np.array([sputtering_yield(argon, tungsten, energy, 0.0, num_ions) for energy in energies])
rustbca_ar_w_ww = np.array([tungsten_sputtering(argon, energy, 0.0, num_ions, run_sim) for energy in energies])
rustbca_ar_w_morse = np.array([tungsten_sputtering(argon, energy, 0.0, num_ions, run_sim,  pot=pot_morse) for energy in energies])
rustbca_ar_w_morse_1175_eV = np.array([tungsten_sputtering(argon, energy, 0.0, num_ions, run_sim,  pot=pot_morse, Es=11.75) for energy in energies])

plt.loglog(data_Somogyvari_Ar_W[:, 0], data_Somogyvari_Ar_W[:, 1], linestyle='', marker='o', label='Exp. Somogyvari et al.')
color = plt.loglog(energies, rustbca_ar_w_ww, linestyle='--', label='RustBCA, ZBL/WW potential')[0].get_color()
plt.loglog(energies, rustbca_ar_w_krc, linestyle=':', label='RustBCA, Kr-C', color=color)
plt.loglog(energies, rustbca_ar_w_morse, linestyle='-.', label='RustBCA, ZBL/Morse', color=color)
plt.loglog(energies, rustbca_ar_w_morse_1175_eV, linestyle='-.', label='RustBCA, ZBL/Morse (Es=11.75eV)', color=color, marker='o')
plt.scatter([100.0, 200.0], [0.1, 0.425], marker='*', label='F Kporha et al., Li Potential')
plt.legend()
plt.xlabel('E [eV]')
plt.ylabel('Y [at/ion]')

# Compiled/produced by Eckstein and Biersack (1986)
# digitized with WebPlotDigitizer from Fig. 1c
data_exp_Saidoh_Sone_1983 = np.array([
    [246.28270037954246, 0.20727494578828917],
    [495.1576962754319, 0.5491441819199058],
    [993.8892543387418, 0.7166462062713239],
    [1943.3723882799943, 1.1996619797357773],
    [6947.722728427526, 3.106222121244088],
])
data_tridyn = np.array([
    [68.14894695911738, 0.0007716158998594478],
    [97.04978626836672, 0.008084212159224518],
    [143.63389679442, 0.04173038662408122],
    [190.63348477243179, 0.09674248249006288],
    [301.12742785416896, 0.23046349454010032],
    [488.4687632973695, 0.5009442817195837],
    [692.5511880621929, 0.7839253399618061],
    [981.6598083496964, 1.1045904715883672],
    [1919.284303781459, 1.7777507596853939],
    [4905.436808264407, 3.1000711459090384],
    [9847.180039865756, 4.207988934522391],
])

energies = np.logspace(np.log10(50), np.log10(8000), 10)
rustbca_w_w_krc = np.array([sputtering_yield(tungsten, tungsten, energy, 0.0, num_ions) for energy in energies])
rustbca_w_w_ww = np.array([tungsten_sputtering(tungsten, energy, 0.0, num_ions, run_sim) for energy in energies])
rustbca_w_w_morse = np.array([tungsten_sputtering(tungsten, energy, 0.0, num_ions, run_sim, pot=pot_morse) for energy in energies])
rustbca_w_w_morse_1175_eV = np.array([tungsten_sputtering(tungsten, energy, 0.0, num_ions, run_sim,  pot=pot_morse, Es=11.75) for energy in energies])
plt.figure()
plt.title('Comparison to Eckstein and Biersack: W on W')
plt.loglog(energies, rustbca_w_w_krc, label='RustBCA Kr-C')
plt.loglog(energies, rustbca_w_w_ww, label='RustBCA WW/ZBL')
plt.loglog(energies, rustbca_w_w_morse, label='RustBCA ZBL/Morse')
plt.loglog(energies, rustbca_w_w_morse_1175_eV, label='RustBCA ZBL/Morse (Es=11.75 eV)')
plt.scatter(data_exp_Saidoh_Sone_1983[:, 0], data_exp_Saidoh_Sone_1983[:, 1], label='Exp. Saidoh Sone 1983')
plt.loglog(data_tridyn[:, 0], data_tridyn[:, 1], label='TRIDYN, Eckstein and Biersack 1986')
plt.scatter([200.0], [0.2], marker='*', label='F Kporha et al., Li Potential')
plt.legend()
plt.show()