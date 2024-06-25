import numpy as np
import matplotlib.pyplot as plt
import os, sys
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'/../scripts')
sys.path.append('scripts')
from materials import *
from formulas import *
from rustbca import *

def run_rustbca(energy, ck, index=0, num_samples=10000, run_sim=True):

    input_file = f'''
    [options]
    name = "estop_{index}"
    track_recoils = false
    weak_collision_order = 0
    electronic_stopping_mode = "INTERPOLATED"
    mean_free_path_model = "LIQUID"
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
    Ec = [ 1.0 ]
    Es = [ 1.5 ]
    interaction_index = [ 0 ]
    pos = [ [ -4.4, 0.0, 0.0,] ]
    dir = [ [ 0.999, 0.001, 0.0,] ]

    [geometry_input]
    length_unit = "ANGSTROM"
    electronic_stopping_correction_factor = {ck}
    densities = [ {aluminum["n"]/10**30} ]

    [material_parameters]
    energy_unit = "EV"
    mass_unit = "AMU"
    Eb = [ {aluminum["Eb"]} ]
    Es = [ {aluminum["Es"]} ]
    Ec = [ {aluminum["Ec"]} ]
    Z = [ {aluminum["Z"]} ]
    m = [ {aluminum["m"]} ]
    interaction_index = [ 0 ]
    surface_binding_model = {{"PLANAR"={{calculation="INDIVIDUAL"}}}}
    bulk_binding_model = "AVERAGE"
    '''

    with open(f'estop_{index}.toml', 'w') as file:
        file.write(input_file)

    if run_sim: os.system(f'cargo run --release 0D estop_{index}.toml')

    deposited_list = np.atleast_2d(np.genfromtxt(f'estop_{index}deposited.output', delimiter=','))

    x = deposited_list[:, 2]
    y = deposited_list[:, 3]

    return x, y


def main():
    energy = 100000
    num_samples = 100000
    ck = np.linspace(0.95, 1.05, 3)
    num_bins = 100
    run_sim = False

    for index, ck_ in enumerate(ck):
        x, y = run_rustbca(energy, ck_, index, num_samples=num_samples, run_sim=run_sim)
        plt.figure(1)
        plt.hist(x, bins=num_bins, histtype='step')
        plt.legend(ck)
        plt.xlabel('x [A}]')
        plt.ylabel('f(x) [counts]')
        plt.figure(2)
        plt.hist(y, bins=num_bins, histtype='step')
        plt.legend(ck)
        plt.xlabel('y [A]')
        plt.ylabel('f(y) [counts]')

    plt.show()

if __name__ == '__main__':
    main()