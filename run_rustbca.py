import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, box
from itertools import chain
import os
import toml

Q = 1.602E-19
PI = 3.14159
AMU = 1.66E-27
ANGSTROM = 1E-10
MICRON = 1E-6
CM = 1E-2
EPS0 = 8.85E-12
A0 = 0.52918E-10
K = 1.11265E-10
ME = 9.11E-31
SQRTPI = 1.77245385
SQRT2PI = 2.506628274631
C = 300000000.

def refraction(E, cosx, cosy, cosz, Es, dx, dy):
    costheta = (cosx*dx + cosy*dy)
    return costheta

def main(Zb, Mb, n, Ec, Es, Eb, Ma, Za, E0, N, theta, thickness, depth, name='test_'):

    options = {
        'name': name,
        'track_trajectories': False,
        'track_recoils': True,
        'track_recoil_trajectories': False,
        'write_files': True
    }

    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'Eb': Eb,
        'Es': Es,
        'Ec': Ec,
        'n': n,
        'Z': Zb,
        'm': Mb
    }
    dx = 2.*n**(-1./3.)/np.sqrt(2.*np.pi)/MICRON

    minx, miny, maxx, maxy = 0.0, -thickness/2., depth, thickness/2.
    surface = box(minx, miny, maxx, maxy)
    energy_surface = surface.buffer(dx, cap_style=2, join_style=2)
    simulation_surface = surface.buffer(2.*dx, cap_style=2, join_style=2)

    geometry = {
        'length_unit': 'MICRON',
        'surface': list(surface.exterior.coords),
        'energy_surface': list(energy_surface.exterior.coords),
        'simulation_surface': list(simulation_surface.exterior.coords),
    }

    cosx = np.cos(theta*180./np.pi)
    sinx = np.sin(theta*180./np.pi)

    particle_parameters = {
        'length_unit': 'MICRON',
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'm': [Ma for _ in range(N)],
        'Z': [Za for _ in range(N)],
        'E': [E0 for _ in range(N)],
        'pos': [(-1.1*dx, np.random.uniform(miny, maxy), 0.) for _ in range(N)],
        'dir': [(cosx, sinx, 0.) for _ in range(N)]
    }

    input_file = {
        'options': options,
        'material_parameters': material_parameters,
        'geometry': geometry,
        'particle_parameters': particle_parameters
    }

    with open('input.toml', 'w') as file:
        toml.dump(input_file, file, encoder=toml.TomlNumpyEncoder())

    os.system('rm *.output')
    os.system('cargo run --release')

    do_plots(surface, energy_surface, simulation_surface, name)


def do_plots(surface, energy_surface, simulation_surface, name):
    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.output', delimiter=','))
    trajectory_data = np.genfromtxt(name+'trajectory_data.output', delimiter=',').transpose().astype(int)

    colors = {
        1: 'red',
        29: 'black',
        2: 'red',
        74: 'blue',
    }

    linewidths = {
        1: 1,
        29: 1,
        2: 1,
        74: 1,
    }

    species_names = {
        1: 'Hydrogen',
        2: 'Helium'
    }

    minx, miny, maxx, maxy = simulation_surface.bounds

    fig1, axis1 = plt.subplots()
    plt.plot(*surface.exterior.xy, color='dimgray')
    plt.plot(*energy_surface.exterior.xy, '--', color='dimgray')
    plt.plot(*simulation_surface.exterior.xy, '--', color='dimgray')

    index = 0
    if np.size(trajectories) > 0:
        for trajectory_length in trajectory_data:

            M = trajectories[index, 0]
            Z = trajectories[index, 1]
            E = trajectories[index:(trajectory_length + index), 2]
            x = trajectories[index:(trajectory_length + index), 3]
            y = trajectories[index:(trajectory_length + index), 4]
            z = trajectories[index:(trajectory_length + index), 5]

            index += trajectory_length

            plt.plot(x, y, color = colors[Z], linewidth = linewidths[Z])
            plt.scatter(x[0], y[0], color = colors[Z], s=10, marker='.')

    if np.size(sputtered) > 0:
        plt.scatter(sputtered[:,3], sputtered[:,4], s=50, color='blue', marker='*')
    plt.xlabel('x [um]')
    plt.ylabel('y [um]')
    plt.title('Helium trajectories in tungsten')

    if np.size(deposited) > 0:
        plt.figure(2)
        num_bins = 50
        bins = np.linspace(miny, maxy, num_bins)
        plt.hist(deposited[:, 3], bins=bins)
        plt.title('Helium y-Deposition')
        plt.xlabel('y [um]')
        plt.ylabel('Helium Deposition (A.U.)')

        plt.figure(3)
        num_bins_x = 50
        num_bins_y = 50
        binx = np.linspace(minx, maxx, num_bins_x)
        biny = np.linspace(miny, maxy, num_bins_y)
        plt.hist2d(deposited[:, 2], deposited[:, 3], bins=(binx, biny))
        plt.plot(*surface.exterior.xy, color='dimgray')
        plt.plot(*energy_surface.exterior.xy, '--', color='dimgray')
        plt.plot(*simulation_surface.exterior.xy, '--', color='dimgray')
        plt.title('Helium Deposition')
        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        #plt.axis('square')

    plt.show()

if __name__ == '__main__':
    Zb = 29
    Mb = 63.54
    n = 8.4E28
    Ec = 3.52
    Es = 3.52
    Eb = 1.0
    Ma = 4
    Za = 2
    E0 = 1E6
    N = 10000
    theta = 0.00001
    thickness = 10
    depth = 10
    main(Zb, Mb, n, Ec, Es, Eb, Ma, Za, E0, N, theta, thickness, depth)
