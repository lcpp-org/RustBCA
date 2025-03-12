from libRustBCA import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'../../scripts')
sys.path.append('../../scripts')
from materials import *
from formulas import *

def four_eight(r, alpha, beta):
    return -alpha/r**4 + beta/r**8

def sigmoid(x, k, x0):
    return 1./(1. + np.exp(-k*(x - x0)))

def lennard_jones_12_6(r, sigma, epsilon):
    return 4.*epsilon*((sigma/r)**12 - (sigma/r)**6)

def morse(r, D, alpha, r0):
    return D*(np.exp(-2.*alpha*(r - r0)) - 2*np.exp(-alpha*(r - r0)))

def morse_krc(r, D, alpha, Za, Zb, r0, k, x0):
    return sigmoid(r, k, x0)*morse(r, D, alpha, r0) + sigmoid(r, -k, x0)*krc(r/ANGSTROM, Za, Zb)*EV

def krc(r, Za, Zb):
    a = 0.885*A0*(np.sqrt(Za) + np.sqrt(Zb))**(-2./3.)
    x = r/a*ANGSTROM
    return Za*Zb*Q**2/4./PI/EPS0/(r*ANGSTROM)*(0.19*np.exp(-0.28*x) + 0.47*np.exp(-0.64*x) + 0.34*np.exp(-1.9*x))/EV

#Corrections to hydrogen to enable low-E simualtions and to match literature Es for H on Ni
hydrogen['Ec'] = 0.1
hydrogen['Es'] = 1.5

D=5.4971E-20
r0=2.782E-10
alpha=1.4198E10
k=7E10
x0=0.75E-10
offset_eV = 1
lw = 3
Za = 1
Zb = nickel["Z"]
EV = 1.602E-19
r = np.logspace(-2, np.log(5.0)/np.log(10), 1000)
plt.figure()
color_krc = plt.plot(r, krc(r, Za, Zb) + offset_eV, label='Kr-C', linewidth=lw)[0].get_color()
color_morse = plt.plot(r, morse(r*ANGSTROM, D, alpha, r0)/EV + offset_eV, label='Morse', linewidth=lw)[0].get_color()
color_krc_morse = plt.plot(r, morse_krc(r*ANGSTROM, D, alpha, Za, Zb, r0, k, x0)/EV + offset_eV, label='Kr-C-Morse', linestyle='--', linewidth=lw/2)[0].get_color()
plt.gca().set_yscale('log')
plt.gca().set_ylim([0.25, 2000])
plt.gca().set_xlim([r[0], r[-1]])
plt.gca().set_xscale('log')
plt.gca().set_xlim([0.1, 5])
plt.plot(r, offset_eV + np.zeros(len(r)), color='black', alpha=0.5)
plt.yticks([offset_eV, 10 + offset_eV, 100 + offset_eV, 1000 + offset_eV], [0, 10, 100, 1000])
plt.title('Interaction potentials for H-Ni')
plt.xlabel('r [A]')
plt.ylabel('E [eV]')
plt.legend()
plt.savefig('interaction_potentials_h_ni.png')

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

data_srim = np.array([
    [0.1, 0.0],
    [1, 0.0],
    [1.5, 1929/11844],
    [2, 4414/12894],
    [5, 13240./32030],
    [10, 19977./50006],
    [50, 6173/19054],
    [100, 10196./35964],
    [500, 3350/17192],
    [1000, 4525./29840],
    [5000, 820/16460],
    [10000, 566./24196],
    [50000, 53/17150.],
    [100000, 29./23084],
])

plt.figure()

#Plotting the EAM data points.
energies = data[:6, 0]
r_benchmark = data[:6, 1]
color = plt.semilogx(energies, r_benchmark, marker='o', linestyle='', label='Molecular Dynamics', color='black')[0].get_color()
plt.errorbar(energies, r_benchmark, yerr=0.1, color=color, linestyle='', capsize=2)

#Plotting the experimental data points.
energies = data[6:, 0]
r_benchmark = data[6:, 1]
plt.semilogx(energies, r_benchmark, marker='^', linestyle='', label='Experiment', color=color, markersize=7)

#energies = data_srim[:, 0]
##r_srim = data_srim[:, 1]
#plt.semilogx(energies, r_srim, linestyle='--', label='SRIM')

#Running and plotting the H-Ni simulations with the Morse potential and updated Es
num_energies = 50
energies = np.logspace(-1, 4, num_energies)
run_sim = False
num_samples = 10000
R_N = np.zeros(num_energies)
R_E = np.zeros(num_energies)
R_N_2 = np.zeros(num_energies)
R_E_2 = np.zeros(num_energies)

for index, energy in enumerate(energies):
    R_N[index], R_E[index] = run_krc_morse_potential(energy, index, num_samples=num_samples, run_sim=run_sim)
    R_N_2[index], R_E_2[index] = run_morse_potential(energy, index, num_samples=num_samples, run_sim=run_sim)

plt.semilogx(energies, R_N, label='Morse-Kr-C Potential', color=color_krc_morse, linestyle='--')
plt.semilogx(energies, R_N_2, label='Morse Potential', color=color_morse)

#Plotting RustBCA data points, using the ergonomic helper function reflection_coefficient().
energies = np.logspace(-1, 4, 50)
r_rustbca = np.array([reflection_coefficient(hydrogen, nickel, energy, 0.0, 10000) for energy in energies])
r_n = r_rustbca[:, 0]
r_e = r_rustbca[:, 1]
plt.semilogx(energies, r_n, label='Kr-C Potential', color=color_krc)

r_thomas = [thomas_reflection(hydrogen, nickel, energy) for energy in energies]
plt.semilogx(energies, r_thomas, label='Thomas Empirical Formula', color='red', linestyle='--')

plt.legend()
plt.title('Reflection Coefficients H+ on Ni')
plt.xlabel('E [eV]')
plt.ylabel('R')
plt.gca().set_xlim([0.1, 1e4])
plt.gca().set_ylim([0.0, 1.2])
plt.show()
plt.savefig('Morse.png')
