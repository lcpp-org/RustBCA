from libRustBCA import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'/../scripts')
sys.path.append('scripts')
import time
from tomlkit import parse, dumps

'''
This script is a first draft of a comprehensive RustBCA input
file creation script using tomlkit.

It includes two geometry modes, 1D and 0D.

It simulates the following situations:

if mode == '1D': 
    H+ (1 keV)
    |
    V
__________
|        |
|    B   | dx = 100 A
_________|
|        |
|  TiB2  | dx = 100 A
_________|
|        |
|   Ti   | dx = 1000 A


if mode == '0D':
    H+ (1 keV)
    |
    V
__________
|        |
|  TiB2  |
|        |
|        |

And calculates implantation profiles, reflection coefficients,
and sputtering yields. It also uses the ergonomic python functions
to compare the result of using the default values for H on B with
the custom values of this input file.

It creates an input file as a nested dictionary which is written to
a TOML file using tomlkit.

It runs the input file with cargo run --release and reads the output files.
'''

run_sim = True
mode = '1D'
incident_energy = 1000.0 # eV
number_ions = 100000 # at least 10k are typically needed for decent results
angle = 45.0 # degrees; measured from surface normal

'''
For organizational purposes, species are commonly defined in dictionaries.
Additional examples can be found in scripts/materials.py, but values 
should be checked for correctness before use. Values are explained
in the relevant sections below.
'''
hydrogen = {
    'symbol': 'H',
    'name': 'hydrogen',
    'Z': 1.0,
    'm': 1.008, # AMU
    'Ec': 0.95, # eV
    'Es': 1.5, # eV
}

titanium = {
    'symbol': 'Ti',
    'name': 'titanium',
    'Z': 22.0,
    'm': 47.867, # AMU
    'Es': 4.84, # eV
    'Ec': 3.5, # eV
    'Eb': 0., # eV
    'Ed': 19.0, # eV
    'n': 5.67e28, # 1/m^3
}

boron = {
    'symbol': 'B',
    'name': 'boron',
    'Z': 5.0,
    'm': 10.811, # AMU
    'n': 1.309E29, # 1/m^3
    'Es': 5.77, # eV
    'Eb': 0., # eV
    'Ec': 5., # eV
    'Ed': 25.0 # eV
}

# species definitions
ion = hydrogen
target1 = boron
target2 = titanium

# geometry definitions
n_i = 0.0328 # 1 / A^3 from n_i = rho_TiB2 / (mB * 2 + mTi)
layer_thicknesses = [100.0, 100.0, 1000.0] # A
layer_1_densities = [boron["n"]/10**30, 0.0] # 1/A^3
layer_2_densities = [n_i * 2, n_i]  # 1/A^3
layer_3_densities = [0.0, titanium["n"]/10**30]  # 1/A^3

options = {
    'name': 'input_file',
    'track_trajectories': False, # whether to track trajectories for plotting; memory intensive
    'track_recoils': True, # whether to track recoils; must enable for sputtering
    'track_recoil_trajectories': False, # whether to track recoil trajectories for plotting
    'track_displacements': False, # whether to track collisions with T > Ed for each species
    'track_energy_losses': False, # whether to track detailed collision energies; memory intensive
    'write_buffer_size': 8192, # how big the buffer is for file writing
    'weak_collision_order': 0, # weak collisions at radii (k + 1)*r; enable only when required
    'suppress_deep_recoils': False, # suppress recoils too deep to ever sputter
    'high_energy_free_flight_paths': False, # SRIM-style high energy free flight distances; use with caution
    'num_threads': 4, # number of threads to run in parallel
    'num_chunks': 10, # code will write to file every nth chunk; for very large simulations, increase num_chunks
    'electronic_stopping_mode': 'LOW_ENERGY_NONLOCAL',
    'mean_free_path_model': 'LIQUID', # liquid is amorphous (constant mean free path); gas is exponentially-distributed mean free paths
    'interaction_potential': [['KR_C']],
    'scattering_integral': [
        [
            {
                'GAUSS_MEHLER': {'n_points': 6}
            }
        ]
    ],

    'root_finder': [
        [
            {
                'NEWTON': {
                    'max_iterations': 100,
                    'tolerance': 1e-6
                }
            }
        ]
    ],
    'seed': 0 # if <0, will generate a seed from thread-local PRNG; if >0, will be used as seed to PRNG
}

# material parameters are per-species
material_parameters = {
    'energy_unit': 'EV',
    'mass_unit': 'AMU',
    # bulk binding energy; typically zero as a model choice
    'Eb': [
        target1["Eb"],
        target2["Eb"]
    ],
    # surface binding energy
    'Es': [
        target1["Es"],
        target2["Es"]
    ],
    # cutoff energy - particles with E < Ec stop
    'Ec': [
        target1["Ec"],
        target2["Ec"]
    ],
    # displacement energy - only used to track displacements
    'Ed': [
        target1["Ed"],
        target2["Ed"]
    ],
    # atomic number
    'Z': [
        target1["Z"],
        target2["Z"]
    ],
    # atomic mass
    'm': [
        target1["m"],
        target2["m"]
    ],
    # used to pick interaction potential from matrix in [options]
    'interaction_index': [0, 0],
    'surface_binding_model': {
        "PLANAR": {'calculation': "INDIVIDUAL"}
    },
    'bulk_binding_model': 'INDIVIDUAL'
}

particle_parameters = {
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
    'pos': [[0.0, 0.0, 0.0]],
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
    'electronic_stopping_correction_factor': 1.0,
    # number densities of each species
    'densities': [2 * n_i, n_i]
}

geometry_1D = {
    'length_unit': 'ANGSTROM',
    # used to correct nonlocal stopping for known compound discrpancies
    'electronic_stopping_correction_factors': [1.0, 1.0, 1.0],
    # thickness of each layer in order from top (x=0) to bottom 
    'layer_thicknesses': layer_thicknesses,
    # number densitiy of each layer in order from top to bottom
    'densities': [
        layer_1_densities,
        layer_2_densities,
        layer_3_densities,
    ]
}

if mode == '1D':
    input_data = {
        'options': options,
        'material_parameters': material_parameters,
        'particle_parameters': particle_parameters,
        'geometry_input': geometry_1D
    }
elif mode == '0D':
    input_data = {
    'options': options,
    'material_parameters': material_parameters,
    'particle_parameters': particle_parameters,
    'geometry_input': geometry_0D
}

# Attempt to cleanup line endings
input_string = dumps(input_data).replace('\r', '')
with  open('examples/input_file.toml', 'w') as input_file:
    input_file.write(input_string)

if run_sim:
    os.system(f'cargo run --release {mode} examples/input_file.toml')

# Read output files - ensure arrays are at least 2D for indexing
sputtered = np.atleast_2d(np.genfromtxt('input_filesputtered.output', delimiter=','))
reflected = np.atleast_2d(np.genfromtxt('input_filereflected.output', delimiter=','))
implanted = np.atleast_2d(np.genfromtxt('input_filedeposited.output', delimiter=','))

print('H on B-TiB2-Ti')
print(f'R_N, R_E: {np.size(reflected[:, 0])/number_ions}, {np.sum(reflected[:, 2]/incident_energy/number_ions)}')
print(f'Y: {np.size(sputtered[:, 0])/number_ions} [at./ion]')
print()
print('H on B (default settings)')
print(f'Y: {sputtering_yield(hydrogen, boron, incident_energy, angle, 10000)} [at./ion]')
print(f'R_N, R_E: {reflection_coefficient(hydrogen, boron, incident_energy, angle, 10000)}')

x = implanted[:, 2]
num_bins=100
bins = np.linspace(0.0, 500.0, num_bins)
plt.hist(x, bins=bins, histtype='step')
if mode == '1D':
    plt.plot([layer_thicknesses[0], layer_thicknesses[0]], [0.0, number_ions], color='gray')
    plt.plot([layer_thicknesses[0] + layer_thicknesses[1], layer_thicknesses[0] + layer_thicknesses[1]], [0.0, number_ions], linestyle='--', color='gray')

plt.ylim([0.0, 3000])
plt.xlim([0.0, 500.0])
plt.xlabel('x [A]')
plt.ylabel('f(x) [counts]')
plt.title('Implantation Distribution H on B-TiB2-Ti')
plt.show()