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

def input_file(ion, target, incident_energy, angle, number_ions=1000):

    mfp = (target["n"]/10**30)**(-1./3.)

    cpr = {'CPR': {'n0': 2, 'nmax': 32, 'epsilon': 1e-3, 'complex_threshold': 1E-9, 'truncation_threshold': 1E-12, 'far_from_zero': 1e3, 'interval_limit': 1E-3, 'derivative_free': True}}
    options = {
        'name': f'input_file_{ion["symbol"]}_{target["symbol"]}_{np.round(angle, 1)}_{np.round(incident_energy/1000, 4)}',
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
        'electronic_stopping_mode': 'INTERPOLATED',
        'mean_free_path_model': 'LIQUID', # liquid is amorphous (constant mean free path); gas is exponentially-distributed mean free paths
        'interaction_potential': [['KR_C', 'KR_C'],
                                  ['KR_C', 'KR_C']],
        'scattering_integral': [
            [{'GAUSS_MEHLER': {'n_points': 5}}, {'GAUSS_MEHLER': {'n_points': 5}}],
            [{'GAUSS_MEHLER': {'n_points': 5}}, {'GAUSS_MEHLER': {'n_points': 5}}],
        ],

        'root_finder': [
            ["DEFAULTNEWTON", "DEFAULTNEWTON"],
            ["DEFAULTNEWTON", "DEFAULTNEWTON"]
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
            target["Es"]
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
        'interaction_index': [0],
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
        'pos': [[-2.*mfp, 0.0, 0.0]],
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

kolasinski = np.array([
[78.52760736196319, 0.07331378299120273],
[98.15950920245393, 0.11290322580645173],
[149.69325153374234, 0.22580645161290347],
[201.2269938650307, 0.35483870967741926],
[250.3067484662576, 0.46627565982404695],
[299.3865030674846, 0.5835777126099708],
[397.5460122699387, 0.7771260997067451],
[500.61349693251526, 0.9149560117302054],
[598.7730061349691, 1.1348973607038124],
[699.3865030674845, 1.2580645161290325],
[998.7730061349691, 1.7243401759530792],
])
tartz  = np.array([
[74.84662576687117, 0.08797653958944296],
[98.15950920245393, 0.12023460410557174],
[228.22085889570548, 0.30498533724340193],
[426.9938650306747, 0.6598240469208212],
[638.0368098159508, 0.9618768328445748],
[840.4907975460121, 1.2023460410557185],
[1040.490797546012, 1.4428152492668622],
[1240.4907975460121, 1.656891495601173],
[1440.4907975460121, 1.8005865102639296],
])
doerner = np.array([
[123.92638036809814, 0.09384164222873892],
[149.69325153374234, 0.10850439882697938],
[174.2331288343558, 0.12023460410557196],
[200, 0.1304985337243405],
])
yalin = np.array([
[200, 0.20674486803519065],
[250.3067484662576, 0.2903225806451615],
[348.46625766871165, 0.4545454545454546],
[500.61349693251526, 0.7360703812316718],
[748.4662576687115, 0.9237536656891496],
])
blandino = np.array([
[500.61349693251526, 0.43695014662756604],
[750.920245398773, 0.6774193548387097],
])
zalm = np.array([
[198.7730061349693, 0.8211143695014664],
[500.61349693251526, 1.598240469208211],
])
bhattacharjee = np.array([
[99.38650306748468, 0.14369501466275691],
[198.7730061349693, 0.32844574780058666],
[299.3865030674846, 0.6070381231671553],
[402.4539877300612, 0.6832844574780059],
[500.61349693251526, 0.7741935483870968],
[603.6809815950919, 0.7976539589442815],
])
weijsenfeld = np.array([
[200, 0.23313782991202325],
[299.38650306748474, 0.4164222873900294],
[399.9999999999997, 0.6011730205278589],
[500.61349693251526, 0.7302052785923754],
[600, 0.9002932551319645],
[699.3865030674845, 1.0835777126099706],
[799.9999999999998, 1.2492668621700878],
[900.6134969325149, 1.4222873900293256],
[1001.2269938650306, 1.5865102639296187],
])
rosenberg = np.array([
[198.77300613496942, 0.3225806451612898],
[299.3865030674846, 0.5425219941348978],
[398.7730061349693, 0.7258064516129035],
[599.9999999999995, 1.0747800586510257],
])

num_ions = 10000
num_energies = 25
run_sim = True
energies = np.logspace(np.log10(25), np.log10(1600), num_energies)
angle = 0.0

datasets = [rosenberg, weijsenfeld, bhattacharjee, zalm, blandino, yalin, doerner, kolasinski, tartz]
dataset_names = [
    'Rosenberg 1962*',
    'Weijsenfeld 1967*',
    'Bhattacharjee 1997*',
    'Zalm 1983*',
    'Blandino 1996*',
    'Yalin 2007*',
    'Doerner 2003*',
    'Kolasinski 2005*',
    'Tartz 2011',
]

Y_Xe_Mo = np.zeros(num_energies)
molybdenum['n'] = 6.452e28
molybdenum['Eb'] = 0.0
molybdenum['Ec'] = 3.0

for index, energy in enumerate(energies):

    input_data = input_file(xenon, molybdenum, energy, angle, num_ions)
    if run_sim: rustbca_py(input_data, geometry_mode="0D")
    sputtered = np.genfromtxt(f'{input_data["options"]["name"]}sputtered.output', delimiter=',')

    Y_Xe_Mo[index] = np.shape(sputtered)[0]/num_ions

for dataset_name, dataset in zip(dataset_names, datasets):
    plt.scatter(dataset[:, 0], dataset[:, 1], label=dataset_name)

plt.plot(energies, Y_Xe_Mo, label='RustBCA Default')
plt.gca().set_xscale('log')
plt.legend()
plt.xlabel('E [eV]')
plt.ylabel('Y [at/ion]')
plt.title('Xe on Mo Sputtering Yields')
plt.show()