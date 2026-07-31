from libRustBCA import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'/../scripts')
sys.path.append('scripts')
from materials import *
from rustbca import *
from itertools import *

#This function simply contains an entire input file as a multi-line f-string to modify some inputs.
def run(energy, index, num_samples=10000, run_sim=True, a=1000, x0=0, y0=500, z0=500, ux=0.999, uy=0.01, uz=0.00, track_trajectories=False):

    input_file = f'''
    #Use with TRIMESH geometry option only
    [options]
    name = "cube_{index}"
    track_trajectories = {str(track_trajectories).lower()}
    track_recoils = true
    track_recoil_trajectories = {str(track_trajectories).lower()}
    weak_collision_order = 0
    num_threads = 4
    num_chunks = 5

    [material_parameters]
    energy_unit = "EV"
    mass_unit = "AMU"
    Eb = [ 0.0,]
    Es = [ {copper["Es"]},]
    Ec = [ {copper["Es"]},]
    Z = [ {copper["Z"]} ]
    m = [ {copper["m"]} ]
    interaction_index = [0]
    surface_binding_model = {{"PLANAR"={{calculation="TARGET"}}}}
    bulk_binding_model = "AVERAGE"

    [particle_parameters]
    length_unit = "ANGSTROM"
    energy_unit = "EV"
    mass_unit = "AMU"
    N = [ {num_samples} ]
    m = [ 4.008 ]
    Z = [ 2 ]
    E = [ {energy} ]
    Ec = [ 1.0 ]
    Es = [ 0.0 ]
    pos = [ [ {x0}, {y0}, {z0},] ]
    dir = [ [ {ux}, {uy}, {uz},] ]
    interaction_index = [ 0 ]

    [geometry_input]
    length_unit = "ANGSTROM"
    electronic_stopping_correction_factor = 1.0
    densities = [0.05]
    vertices = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, {a}],
        [0.0, {a}, 0.0],
        [{a}, 0.0, 0.0],
        [0.0, {a}, {a}],
        [{a}, 0.0, {a}],
        [{a}, {a}, 0.0],
        [{a}, {a}, {a}]
    ]
    indices = [
        #front
        [3, 1, 0],
        [3, 5, 1],
        #bottom
        [3, 2, 0],
        [6, 2, 3],
        #right
        [3, 6, 7],
        [3, 7, 5],
        #left
        [0, 1, 2],
        [2, 1, 4],
        #top
        [5, 4, 7],
        [1, 4, 5],
        #back
        [4, 7, 6],
        [6, 2, 4]
    ]
    '''

    with open(f'cube_{index}.toml', 'w') as file:
        file.write(input_file)

    if run_sim: os.system(f'cargo run --release --features parry3d TRIMESH cube_{index}.toml')

    reflected_list = np.atleast_2d(np.genfromtxt(f'cube_{index}reflected.output', delimiter=','))
    if np.size(reflected_list) > 0:
        num_reflected = np.shape(reflected_list)[0]
        energy_reflected = np.sum(reflected_list[:, 2])
    else:
        num_reflected = 0
        energy_reflected = 0.

    sputtered_list = np.atleast_2d(np.genfromtxt(f'cube_{index}sputtered.output', delimiter=','))
    if np.size(sputtered_list) > 0:
        num_sputtered = np.shape(sputtered_list)[0]
        energy_sputtered = np.sum(sputtered_list[:, 2])
    else:
        num_sputtered = 0
        energy_sputtered = 0.

    implanted_list = np.atleast_2d(np.genfromtxt(f'cube_{index}deposited.output', delimiter=','))

    return num_reflected/num_samples, num_sputtered/num_samples, reflected_list, sputtered_list, implanted_list

def main():
    do_plots = True
    run_sim = True
    track_trajectories = False
    a = 1000
    num_samples = 100000
    num_angles = 10
    energy = 200
    angles = np.linspace(0.0, 89.9, num_angles)
    
    sputtering_yields_x = np.zeros(num_angles)
    sputtering_yields_z = np.zeros(num_angles)
    sputtering_yields_y = np.zeros(num_angles)

    R_x = np.zeros(num_angles)
    R_y = np.zeros(num_angles)
    R_z = np.zeros(num_angles)

    sputtering_yields_x_minus = np.zeros(num_angles)
    sputtering_yields_z_minus = np.zeros(num_angles)
    sputtering_yields_y_minus = np.zeros(num_angles)

    R_x_minus = np.zeros(num_angles)
    R_y_minus = np.zeros(num_angles)
    R_z_minus = np.zeros(num_angles)

    range_x = np.zeros(num_angles)
    range_y = np.zeros(num_angles)
    range_z = np.zeros(num_angles)
    range_x_minus = np.zeros(num_angles)
    range_y_minus = np.zeros(num_angles)
    range_z_minus = np.zeros(num_angles)

    sim_index = 0
    impl_index = num_angles*3//4
    bins = np.linspace(0, 500, 100)
    
    for angle_index, angle in enumerate(angles):

        R, Y, reflected, sputtered, implanted = run(energy, sim_index, num_samples, track_trajectories=track_trajectories, run_sim=run_sim, ux=np.cos(angle*np.pi/180.), uy=np.sin(angle*np.pi/180.)/np.sqrt(2), uz=np.sin(angle*np.pi/180.)/np.sqrt(2))
        sim_index += 1
        plt.figure(1)
        if angle_index == impl_index: plt.hist(implanted[:, 2], histtype='step', bins=bins, density=False, label='x')
        sputtering_yields_x[angle_index] = Y
        R_x[angle_index] = R
        range_x[angle_index] = np.mean(implanted[:, 2])

        R, Y, reflected, sputtered, implanted = run(energy, sim_index, num_samples, track_trajectories=track_trajectories, run_sim=run_sim, x0=a, ux=-np.cos(angle*np.pi/180.), uy=np.sin(angle*np.pi/180.)/np.sqrt(2), uz=np.sin(angle*np.pi/180.)/np.sqrt(2))
        sim_index += 1
        if angle_index == impl_index: plt.hist(a - implanted[:, 2], histtype='step', bins=bins, density=False, label='-x')
        sputtering_yields_x_minus[angle_index] = Y
        R_x_minus[angle_index] = R
        range_x_minus[angle_index] = a - np.mean(implanted[:, 2])

        R, Y, reflected, sputtered, implanted = run(energy, sim_index, num_samples, track_trajectories=track_trajectories, run_sim=run_sim, x0=a/2., y0=a/2., z0=a, ux=np.sin(angle*np.pi/180.)/np.sqrt(2), uy=np.sin(angle*np.pi/180.)/np.sqrt(2), uz=-np.cos(angle*np.pi/180.))
        sim_index += 1
        if angle_index == impl_index: plt.hist(a - implanted[:, 4], histtype='step', bins=bins, density=False, label='-z')
        sputtering_yields_z_minus[angle_index] = Y
        R_z_minus[angle_index] = R
        range_z_minus[angle_index] = a - np.mean(implanted[:, 4])

        R, Y, reflected, sputtered, implanted = run(energy, sim_index, num_samples, track_trajectories=track_trajectories, run_sim=run_sim, x0=a/2., y0=a/2., z0=0.0, ux=np.sin(angle*np.pi/180.)/np.sqrt(2), uy=np.sin(angle*np.pi/180.)/np.sqrt(2), uz=np.cos(angle*np.pi/180.))
        sim_index += 1
        if angle_index == impl_index: plt.hist(implanted[:, 4], histtype='step', bins=bins, density=False, label='z')
        sputtering_yields_z[angle_index] = Y
        R_z[angle_index] = R
        range_z[angle_index] = np.mean(implanted[:, 4])

        R, Y, reflected, sputtered, implanted = run(energy, sim_index, num_samples, track_trajectories=track_trajectories, run_sim=run_sim, x0=a/2., y0=0.0, z0=a/2., ux=np.sin(angle*np.pi/180.)/np.sqrt(2), uz=np.sin(angle*np.pi/180.)/np.sqrt(2), uy=np.cos(angle*np.pi/180.))
        sim_index += 1
        if angle_index == impl_index: plt.hist(implanted[:, 3], histtype='step', bins=bins, density=False, label='y')
        sputtering_yields_y[angle_index] = Y
        R_y[angle_index] = R
        range_y[angle_index] = np.mean(implanted[:, 3])

        R, Y, reflected, sputtered, implanted = run(energy, sim_index, num_samples, track_trajectories=track_trajectories, run_sim=run_sim, x0=a/2., y0=a, z0=a/2., ux=np.sin(angle*np.pi/180.)/np.sqrt(2), uz=np.sin(angle*np.pi/180.)/np.sqrt(2), uy=-np.cos(angle*np.pi/180.))
        sim_index += 1
        if angle_index == impl_index: plt.hist(a - implanted[:, 3], histtype='step', bins=bins, density=False, label='-y')
        sputtering_yields_y_minus[angle_index] = Y
        R_y_minus[angle_index] = R
        range_y_minus[angle_index] = a - np.mean(implanted[:, 3])

    plt.figure(1)
    plt.title('Implantation')
    plt.legend()

    yields = [
        sputtering_yields_x,
        sputtering_yields_y,
        sputtering_yields_z,
        sputtering_yields_x_minus,
        sputtering_yields_y_minus,
        sputtering_yields_z_minus
    ]

    # All sputtering yields for all directions should be close (3*sigma). This compares every pair
    for Y1 in yields:
        for Y2 in yields:
            np.testing.assert_allclose(Y1, Y2, atol=3./np.sqrt(num_samples))

    coeffs = [
        R_x,
        R_y,
        R_z,
        R_x_minus,
        R_y_minus,
        R_z_minus
    ]

    # All reflection coefficients for all directions should be close (3*sigma). This compares every pair
    for R1 in coeffs:
        assert(R1[-1] > 0.99) # R ~ 1 at 89.9 degrees
        for R2 in coeffs:
            np.testing.assert_allclose(R1, R2, atol=3./np.sqrt(num_samples))

    ranges = [
        range_x,
        range_y,
        range_z,
        range_x_minus,
        range_y_minus,
        range_z_minus,
    ]

    for range1 in ranges:
        for range2 in ranges: # since R ~ 1 at the highest two angles, there's too much variance in ranges
            np.testing.assert_allclose(range1[:-2], range2[:-2], rtol=0.05)

    breakpoint()

    # All ranges 

    if do_plots:
        plt.figure(2)
        plt.title('Sputtering')
        plt.plot(angles, sputtering_yields_x, label='+x directed')
        plt.plot(angles, sputtering_yields_z, label='+z directed')
        plt.plot(angles, sputtering_yields_y, label='+y directed')
        plt.plot(angles, sputtering_yields_x_minus, label='-x directed')
        plt.plot(angles, sputtering_yields_z_minus, label='-z directed')
        plt.plot(angles, sputtering_yields_y_minus, label='-y directed')
        plt.legend()

        plt.figure(3)
        plt.title('Reflection')
        plt.plot(angles, R_x, label='+x directed')
        plt.plot(angles, R_z, label='+z directed')
        plt.plot(angles, R_y, label='+y directed')
        plt.plot(angles, R_x_minus, label='-x directed')
        plt.plot(angles, R_z_minus, label='-z directed')
        plt.plot(angles, R_y_minus, label='-y directed')
        plt.legend()

        plt.show()
        plt.savefig('test_cube.png')
        plt.clf()

    if track_trajectories and do_plots:
        do_trajectory_plot('cube_19', boundary=[[0.0, 0.0], [0.0, a], [a, a], [a, 0.0], [0.0, 0.0]])
        plt.savefig('test_cube_traj_1.png')
        plt.clf()
        
        do_trajectory_plot('cube_20', boundary=[[0.0, 0.0], [0.0, a], [a, a], [a, 0.0], [0.0, 0.0]])
        plt.savefig('test_cube_traj_2.png')
        plt.clf()

if __name__ == '__main__':
    main()