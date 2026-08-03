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

def run_test(
        interaction_potential='KR_C',
        high_energy_free_flight_path=False,
        electronic_stopping_mode='LOW_ENERGY_NONLOCAL',
        mfp='LIQUID',
        weak_collision_order=0,
        num_threads=6,
        index=0,
        scattering_integral={'GAUSS_MEHLER': {'n_points': 6}},
        run_sim=True,
        num_samples=100000
        ):
    mode = '1D'
    incident_energy = 1000.0 # eV
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
        'name': f'input_file_{index}',
        'track_trajectories': False, # whether to track trajectories for plotting; memory intensive
        'track_recoils': True, # whether to track recoils; must enable for sputtering
        'track_recoil_trajectories': False, # whether to track recoil trajectories for plotting
        'track_displacements': False, # whether to track collisions with T > Ed for each species
        'track_energy_losses': False, # whether to track detailed collision energies; memory intensive
        'write_buffer_size': 8192, # how big the buffer is for file writing
        'weak_collision_order': weak_collision_order, # weak collisions at radii (k + 1)*r; enable only when required
        'suppress_deep_recoils': False, # suppress recoils too deep to ever sputter
        'high_energy_free_flight_paths': high_energy_free_flight_path, # SRIM-style high energy free flight distances; use with caution
        'num_threads': num_threads, # number of threads to run in parallel
        'num_chunks': 10, # code will write to file every nth chunk; for very large simulations, increase num_chunks
        'electronic_stopping_mode': electronic_stopping_mode,
        'mean_free_path_model': mfp, # liquid is amorphous (constant mean free path); gas is exponentially-distributed mean free paths
        'interaction_potential': [[interaction_potential]],
        'scattering_integral': [
            [
                scattering_integral
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
        'N': [num_samples],
        # atomic mass
        'm': [ion["m"]],
        # atomic number
        'Z': [ion["Z"]],
        # incidenet energy 
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
    with  open(f'examples/input_file_{index}.toml', 'w') as input_file:
        input_file.write(input_string)

    if run_sim:
        os.system(f'cargo run --release {mode} examples/input_file_{index}.toml')

    # Read output files - ensure arrays are at least 2D for indexing
    sputtered = np.atleast_2d(np.genfromtxt(f'input_file_{index}sputtered.output', delimiter=','))
    reflected = np.atleast_2d(np.genfromtxt(f'input_file_{index}reflected.output', delimiter=','))
    implanted = np.atleast_2d(np.genfromtxt(f'input_file_{index}deposited.output', delimiter=','))

    return sputtered, reflected, implanted

num_bins = 75
num_samples = 100000
run_sim = True
show_plots = True

# interaction potentials
interaction_potentials = ['KR_C', 'ZBL', 'MOLIERE', 'LENZ_JENSEN']

sim_index = 0

if not run_sim:
    sim_times = np.genfromtxt('sim_times.txt')
else:
    sim_times = []

Y_test = np.array([
    0.02593, 0.02104, 0.02467, 0.02731, 0.02593, 0.0218, 0.02211, 0.02288, 
    0.03577, 0.02593, 0.02587, 0.0296, 0.0255, 0.02568, 0.02587, 0.02593, 
    0.02905, 0.02593, 0.02593, 0.02593, 0.02593, 0.02593, 0.02593, 0.02593, 
    0.02593
])

R_N_test = np.array([
    0.17442, 0.17032, 0.16802, 0.18216, 0.17442, 0.17437, 0.17449, 0.17777, 
    0.28801, 0.17442, 0.17452, 0.22074, 0.17506, 0.17446, 0.17445, 0.17442, 
    0.17722, 0.17442, 0.17442, 0.17442, 0.17442, 0.17442, 0.17442, 0.17442, 
    0.17442
])

Y = np.zeros_like(Y_test)
R_N = np.zeros_like(R_N_test)

for interaction_potential in interaction_potentials:
    start = time.time()
    s, r, i = run_test(interaction_potential=interaction_potential, index=sim_index, num_samples=num_samples, run_sim=run_sim)
    stop = time.time()
    Y[sim_index] = np.shape(s)[0]/num_samples
    R_N[sim_index] = np.shape(r)[0]/num_samples
    sim_time = (stop - start)/1e-3
    if run_sim: sim_times.append(sim_time)
    x = i[:, 2]
    plt.figure(1)
    plt.title('Interaction Potentials')
    plt.hist(x, bins=num_bins, histtype='step', label=f'{interaction_potential} [{np.round(sim_time)} ms]')
    plt.legend()
    plt.xlabel('x [A]')
    plt.ylabel(f'f(x) [counts]')
    print(f'{sim_index} Y: {np.shape(s)[0]/num_samples} R_N: {np.shape(r)[0]/num_samples}')
    sim_index += 1

for weak_collision_order in [0, 1, 2, 3]:

    start = time.time()
    s, r, i = run_test(weak_collision_order=weak_collision_order, index=sim_index, num_samples=num_samples, run_sim=run_sim)
    stop = time.time()
    Y[sim_index] = np.shape(s)[0]/num_samples
    R_N[sim_index] = np.shape(r)[0]/num_samples
    sim_time = (stop - start)/1e-3
    if run_sim: sim_times.append(sim_time)
    x = i[:, 2]
    plt.figure(2)
    plt.title('Weak Collision Orders')
    plt.hist(x, bins=num_bins, histtype='step', label=f'k={weak_collision_order} [{np.round(sim_time)} ms]')
    plt.legend()
    plt.xlabel('x [A]')
    plt.ylabel(f'f(x) [counts]')
    print(f'{sim_index} Y: {np.shape(s)[0]/num_samples} R_N: {np.shape(r)[0]/num_samples}')
    sim_index += 1

for electronic_stopping_mode in ['LOW_ENERGY_LOCAL', 'LOW_ENERGY_NONLOCAL', 'INTERPOLATED', 'LOW_ENERGY_EQUIPARTITION']:
    start = time.time()
    s, r, i = run_test(electronic_stopping_mode=electronic_stopping_mode, index=sim_index, num_samples=num_samples, run_sim=run_sim)
    stop = time.time()
    Y[sim_index] = np.shape(s)[0]/num_samples
    R_N[sim_index] = np.shape(r)[0]/num_samples
    sim_time = (stop - start)/1e-3
    if run_sim: sim_times.append(sim_time)
    x = i[:, 2]
    plt.figure(3)
    plt.title('Electronic Stopping Modes')
    plt.hist(x, bins=num_bins, histtype='step', label=f'{electronic_stopping_mode} [{np.round(sim_time)} ms]')
    plt.legend()
    plt.xlabel('x [A]')
    plt.ylabel(f'f(x) [counts]')
    print(f'{sim_index} Y: {np.shape(s)[0]/num_samples} R_N: {np.shape(r)[0]/num_samples}')
    sim_index += 1

for scattering_integral in ['MENDENHALL_WELLER', {'GAUSS_MEHLER': {'n_points': 6}}, 'GAUSS_LEGENDRE']:
    start = time.time()
    s, r, i = run_test(scattering_integral=scattering_integral, index=sim_index, num_samples=num_samples, run_sim=run_sim)
    stop = time.time()
    Y[sim_index] = np.shape(s)[0]/num_samples
    R_N[sim_index] = np.shape(r)[0]/num_samples
    sim_time = (stop - start)/1e-3
    if run_sim: sim_times.append(sim_time)
    x = i[:, 2]
    plt.figure(4)
    plt.title('Scattering Integrals')
    plt.hist(x, bins=num_bins, histtype='step', label=f'{scattering_integral} [{np.round(sim_time)} ms]')
    plt.legend()
    plt.xlabel('x [A]')
    plt.ylabel(f'f(x) [counts]')
    print(f'{sim_index} Y: {np.shape(s)[0]/num_samples} R_N: {np.shape(r)[0]/num_samples}')
    sim_index += 1

for mfp in ['LIQUID', 'GASEOUS']:
    start = time.time()
    s, r, i = run_test(mfp=mfp, index=sim_index, num_samples=num_samples, run_sim=run_sim)
    stop = time.time()
    Y[sim_index] = np.shape(s)[0]/num_samples
    R_N[sim_index] = np.shape(r)[0]/num_samples
    sim_time = (stop - start)/1e-3
    if run_sim: sim_times.append(sim_time)
    x = i[:, 2]
    plt.figure(5)
    plt.title('MFP Distribution')
    plt.hist(x, bins=num_bins, histtype='step', label=f'{mfp} [{np.round(sim_time)} ms]')
    plt.legend()
    plt.xlabel('x [A]')
    plt.ylabel(f'f(x) [counts]')
    print(f'{sim_index} Y: {np.shape(s)[0]/num_samples} R_N: {np.shape(r)[0]/num_samples}')
    sim_index += 1

num_threads = [1, 2, 3, 4, 5, 6, 7, 8]
sim_index_threads_start = sim_index
for n in num_threads:
    start = time.time()
    s, r, i = run_test(num_threads=n, index=sim_index, num_samples=num_samples, run_sim=run_sim)
    stop = time.time()
    Y[sim_index] = np.shape(s)[0]/num_samples
    R_N[sim_index] = np.shape(r)[0]/num_samples
    sim_time = (stop - start)/1e-3
    if run_sim: sim_times.append(sim_time)
    print(f'{sim_index} Y: {np.shape(s)[0]/num_samples} R_N: {np.shape(r)[0]/num_samples}')
    sim_index += 1

sim_index_threads_stop = sim_index
plt.figure(6)
plt.plot(num_threads, sim_times[sim_index_threads_start]/np.array(sim_times[sim_index_threads_start:sim_index_threads_stop]), label="Amdahl's law; s=0.06")
s = 0.06
p = 1 - s
plt.plot(num_threads, 1/(s + p/np.array(num_threads)), label='RustBCA (i5-8600k, 6 cores)')
plt.xlabel('n threads')
plt.ylabel('t [ms]')
plt.legend()
plt.plot([6, 6], [0, 10], linestyle='--', color='gray')
plt.gca().set_ylim([1, 5])

if run_sim: np.savetxt('sim_times.txt', sim_times)
if show_plots: plt.show()

np.testing.assert_allclose(Y_test, Y)
np.testing.assert_allclose(R_N_test, R_N)


