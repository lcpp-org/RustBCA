from libRustBCA import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'/../scripts')
sys.path.append('scripts')
from materials import *
from formulas import *

#Corrections to hydrogen to enable low-E simualtions and to match literature Es for H on Ni
hydrogen['Ec'] = 0.1
hydrogen['Es'] = 1.5

#This function simply contains an entire input file as a multi-line f-string to modify some inputs.
def run_morse_potential(energy, index, num_samples=10000, run_sim=True):

    input_file = f'''
    [options]
    name = "morse_{index}"
    track_recoils = false
    weak_collision_order = 0
    electronic_stopping_mode = "LOW_ENERGY_NONLOCAL"
    mean_free_path_model = "LIQUID"
    interaction_potential = [[{{"MORSE"={{D=5.4971E-20, r0=2.782E-10, alpha=1.4198E10}}}}]]
    scattering_integral = [["GAUSS_LEGENDRE"]]
    root_finder = [[{{"CPR"={{n0=2, nmax=100, epsilon=1E-9, complex_threshold=1E-3, truncation_threshold=1E-9, far_from_zero=1E9, interval_limit=1E-12, derivative_free=true}}}}]]
    num_threads = 4
    num_chunks = 10

    [particle_parameters]
    length_unit = "ANGSTROM"
    energy_unit = "EV"
    mass_unit = "AMU"
    N = [ {num_samples} ]
    m = [ 1.008 ]
    Z = [ 1 ]
    E = [ {energy} ]
    Ec = [ 0.1 ]
    Es = [ 1.5 ]
    interaction_index = [ 0 ]
    pos = [ [ -4.4, 0.0, 0.0,] ]
    dir = [ [ 0.9999999999984769, 1.7453292519934434e-6, 0.0,] ]

    [geometry_input]
    length_unit = "ANGSTROM"
    electronic_stopping_correction_factor = 1.09
    densities = [ 0.0914 ]

    [material_parameters]
    energy_unit = "EV"
    mass_unit = "AMU"
    Eb = [ 0.0 ]
    Es = [ 5.61 ]
    Ec = [ 3.0 ]
    Z = [ 28 ]
    m = [ 58.69 ]
    interaction_index = [ 0 ]
    surface_binding_model = {{"PLANAR"={{calculation="INDIVIDUAL"}}}}
    bulk_binding_model = "AVERAGE"
    '''

    with open(f'morse_{index}.toml', 'w') as file:
        file.write(input_file)

    if run_sim: os.system(f'cargo run --release --features cpr_rootfinder 0D morse_{index}.toml')

    reflected_list = np.atleast_2d(np.genfromtxt(f'morse_{index}reflected.output', delimiter=','))
    if np.size(reflected_list) > 0:
        num_reflected = np.shape(reflected_list)[0]
        energy_reflected = np.sum(reflected_list[:, 2])
    else:
        num_reflected = 0
        energy_reflected = 0.

    return num_reflected/num_samples, energy_reflected/energy/num_samples

def run_krc_morse_potential(energy, index, num_samples=10000, run_sim=True):

    input_file = f'''
    [options]
    name = "krc_morse_{index}"
    track_recoils = false
    weak_collision_order = 0
    electronic_stopping_mode = "LOW_ENERGY_NONLOCAL"
    mean_free_path_model = "LIQUID"
    interaction_potential = [[{{"KRC_MORSE"={{D=5.4971E-20, r0=2.782E-10, alpha=1.4198E10, k=7E10, x0=0.75E-10}}}}]]
    scattering_integral = [["GAUSS_LEGENDRE"]]
    root_finder = [[{{"CPR"={{n0=2, nmax=100, epsilon=1E-9, complex_threshold=1E-3, truncation_threshold=1E-9, far_from_zero=1E9, interval_limit=1E-12, derivative_free=true}}}}]]
    num_threads = 4
    num_chunks = 10

    [particle_parameters]
    length_unit = "ANGSTROM"
    energy_unit = "EV"
    mass_unit = "AMU"
    N = [ {num_samples} ]
    m = [ 1.008 ]
    Z = [ 1 ]
    E = [ {energy} ]
    Ec = [ 0.1 ]
    Es = [ 1.5 ]
    interaction_index = [ 0 ]
    pos = [ [ -4.4, 0.0, 0.0,] ]
    dir = [ [ 0.9999999999984769, 1.7453292519934434e-6, 0.0,] ]

    [geometry_input]
    length_unit = "ANGSTROM"
    electronic_stopping_correction_factor = 1.09
    densities = [ 0.0914 ]

    [material_parameters]
    energy_unit = "EV"
    mass_unit = "AMU"
    Eb = [ 0.0 ]
    Es = [ 5.61 ]
    Ec = [ 3.0 ]
    Z = [ 28 ]
    m = [ 58.69 ]
    interaction_index = [ 0 ]
    surface_binding_model = {{"PLANAR"={{calculation="INDIVIDUAL"}}}}
    bulk_binding_model = "AVERAGE"
    '''

    with open(f'krc_morse_{index}.toml', 'w') as file:
        file.write(input_file)

    if run_sim: os.system(f'cargo run --release --features cpr_rootfinder 0D krc_morse_{index}.toml')

    reflected_list = np.atleast_2d(np.genfromtxt(f'krc_morse_{index}reflected.output', delimiter=','))
    if np.size(reflected_list) > 0:
        num_reflected = np.shape(reflected_list)[0]
        energy_reflected = np.sum(reflected_list[:, 2])
    else:
        num_reflected = 0
        energy_reflected = 0.

    return num_reflected/num_samples, energy_reflected/energy/num_samples


#Data digitized from Fig. 1 of https://doi.org/10.1016/0022-3115(84)90433-1
#Modeled reflection coefficients using EAM at 100 eV and below and experimentally at 100 eV and above
#Used to show that TRIM breaks down at lower incident energies
data = np.array(
    [[0.2962771, 0.18217821],
    [0.96344626, 0.4871287],
    [3.0372448, 0.8930693],
    [9.876638, 0.86534655],
    [30.184526, 0.7841584],
    [96.644066, 0.6217822], #note: this is the last EAM data point
    [102.83226, 0.4950495], #note: this is the first experimental data point
    [203.52855, 0.37623763],
    [516.34265, 0.3980198],
    [809.75903, 0.32277226],
    [1006.2257, 0.2990099],
    [2054.3218, 0.17821783],
    [4129.5522, 0.13069306],
    [6890.884, 0.0990099]]
)

#Plotting the EAM data points.
energies = data[:6, 0]
r_benchmark = data[:6, 1]
plt.semilogx(energies, r_benchmark, marker='o', linestyle='', label='EAM')

#Plotting the experimental data points.
energies = data[6:, 0]
r_benchmark = data[6:, 1]
plt.semilogx(energies, r_benchmark, marker='^', linestyle='', label='Exp.')

#Running and plotting the H-Ni simulations with the Morse potential and updated Es
num_energies = 15
energies = np.logspace(-1, 4, num_energies)
run_sim = True
num_samples = 10000
R_N = np.zeros(num_energies)
R_E = np.zeros(num_energies)
R_N_2 = np.zeros(num_energies)
R_E_2 = np.zeros(num_energies)

for index, energy in enumerate(energies):
    R_N[index], R_E[index] = run_krc_morse_potential(energy, index, num_samples=num_samples, run_sim=True)
    R_N_2[index], R_E_2[index] = run_morse_potential(energy, index, num_samples=num_samples, run_sim=True)

plt.semilogx(energies, R_N, label='R_N Morse-Kr-C H-Ni, Es=1.5eV', color='purple')
plt.semilogx(energies, R_N_2, label='R_N Morse H-Ni, Es=1.5eV', color='green')

#Plotting RustBCA data points, using the ergonomic helper function reflection_coefficient().
energies = np.logspace(-1, 4, 50)
r_rustbca = np.array([reflection_coefficient(hydrogen, nickel, energy, 0.0, 10000) for energy in energies])
r_n = r_rustbca[:, 0]
r_e = r_rustbca[:, 1]
plt.semilogx(energies, r_n, label='R_N, Default Settings', color='black')

r_thomas = [thomas_reflection(hydrogen, nickel, energy) for energy in energies]
plt.semilogx(energies, r_thomas, label='Thomas', color='red', linestyle='--')

plt.legend()
plt.title('Reflection Coefficients H+ on Ni')
plt.xlabel('E [eV]')
plt.ylabel('R')
plt.show()
plt.savefig('Morse.png')
