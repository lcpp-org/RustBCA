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

def main(Zb, Mb, n, Ec, Es, Eb, Ma, Za, E0, N, N_, theta, thickness, depth, track_trajectories=False, track_recoils=False, track_recoil_trajectories=False, write_files=True, name='test_'):

    options = {
        'name': name,
        'track_trajectories': track_trajectories,
        'track_recoils': track_recoils,
        'track_recoil_trajectories': track_recoil_trajectories,
        'write_files': write_files
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
        'N': [N_ for _ in range(N)],
        'm': [Ma for _ in range(N)],
        'Z': [Za for _ in range(N)],
        'E': [E0 for _ in range(N)],
        'pos': [(-1.1*dx, 0., 0.) for _ in range(N)],
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
    os.system('rustBCA.exe')

    Y, R, ion_range = do_plots(surface, energy_surface, simulation_surface, name)
    return Y, R, ion_range


def do_plots(surface, energy_surface, simulation_surface, name):
    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.output', delimiter=','))
    trajectory_data = np.genfromtxt(name+'trajectory_data.output', delimiter=',').transpose().astype(int)

    if False:

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
        plt.title('5 MeV Helium Deposition on Copper')

        if np.size(deposited) > 0:
            plt.figure(2)
            num_bins = 200
            bins = np.linspace(minx, maxx, num_bins)
            plt.hist(deposited[:, 2], bins=bins)
            plt.title('5 MeV Helium x-Deposition on Copper')
            plt.xlabel('y [um]')
            plt.ylabel('Helium Deposition (A.U.)')

            plt.figure(3)
            num_bins = 200
            bins = np.linspace(miny, maxy, num_bins)
            plt.hist(deposited[:, 3], bins=bins)
            plt.title('5 MeV Helium y-Deposition on Copper')
            plt.xlabel('y [um]')
            plt.ylabel('Helium Deposition (A.U.)')

            plt.figure(4)
            num_bins_x = 200
            num_bins_y = 200
            binx = np.linspace(minx, maxx, num_bins_x)
            biny = np.linspace(miny, maxy, num_bins_y)
            plt.hist2d(deposited[:, 2], deposited[:, 3], bins=(binx, biny))
            plt.plot(*surface.exterior.xy, color='dimgray')
            plt.plot(*energy_surface.exterior.xy, '--', color='dimgray')
            plt.plot(*simulation_surface.exterior.xy, '--', color='dimgray')
            plt.title('5 MeV Helium Deposition on Copper')
            plt.xlabel('x [um]')
            plt.ylabel('y [um]')
            #plt.axis('square')

    #plt.show()
    if np.size(sputtered) > 0:
        Y = len(sputtered[:,0])/N/N_
    else:
        Y = 0

    if np.size(reflected) > 0:
        R = len(reflected[:,0])/N/N_
    else:
        R = 0

    if np.size(deposited) > 0:
        ion_range = np.mean(deposited[:, 2])
    else:
        ion_range = 0

    return Y, R, ion_range

if __name__ == '__main__':
    hydrogen = {
        'symbol': 'H',
        'name': 'hydrogen',
        'Z': 1,
        'm': 1.008,
    }
    helium = {
        'symbol': 'He',
        'name': 'helium',
        'Z': 2,
        'm': 4.002602
    }

    beryllium = {
        'symbol': 'Be',
        'name': 'beryllium',
        'Z': 4,
        'm': 9.012182,
        'n': 1.235E29,
        'Es': 3.31,
        'Eb': 3.,
        'Ec': 3.,
    }

    boron = {
        'symbol': 'B',
        'name': 'boron',
        'Z': 5,
        'm': 10.811,
        'n': 1.37E29,
        'Es': 5.76,
        'Eb': 3.,
        'Ec': 3.,
    }

    neon = {
        'symbol': 'Ne',
        'name': 'neon',
        'Z': 10,
        'm': 20.1797,
    }

    silicon = {
        'symbol': 'Si',
        'name': 'silicon',
        'Z': 14,
        'm': 28.08553,
        'n': 4.996E28,
        'Es': 4.72,
        'Eb': 2.,
        'Ec': 2.,
    }

    argon = {
        'symbol': 'Ar',
        'name': 'argon',
        'Z': 18,
        'm': 39.948,
    }

    copper = {
        'symbol': 'Cu',
        'name': 'copper',
        'Z': 29,
        'm': 63.546,
        'n': 8.491E28,
        'Es': 3.52,
        'Eb': 3.,
        'Ec': 3.,
    }

    tungsten = {
        'symbol': 'W',
        'name': 'tungsten',
        'Z': 74,
        'm': 183.84,
        'n': 6.306E28,
        'Es': 11.75,
        'Eb': 3.,
        'Ec': 3.,
    }

    beam_species = [hydrogen, helium]#, beryllium, boron, neon, silicon, argon, copper, tungsten]
    target_species = [beryllium, boron, silicon, copper, tungsten]

    N = 1
    N_ = 10000
    theta = 0.00001
    thickness = 1000
    depth = 1000
    energies = np.round(np.logspace(1, 4, 20))

    name = []


    for beam in beam_species:
        Za = beam['Z']
        Ma = beam['m']
        for target in target_species:
            name.append(beam['symbol']+' on '+target['symbol']+' Y')
            name.append(beam['symbol']+' on '+target['symbol']+' R')
            Y = []
            R = []
            ion_range = []
            Zb = target['Z']
            Ec = target['Ec']
            Es = target['Es']
            Eb = target['Eb']
            Mb = target['m']
            n = target['n']
            for energy in energies:
                Y_, R_, ion_range_ = main(Zb, Mb, n, Ec, Es, Eb, Ma, Za, energy, N, N_, theta, thickness, depth, track_trajectories=False, track_recoils=True, track_recoil_trajectories=False, write_files=True, name=str(energy)+beam['symbol']+'_'+target['symbol'])
                Y.append(Y_)
                R.append(R_)
                ion_range.append(ion_range_)
            plt.figure(1)
            handle = plt.loglog(energies, Y)
            plt.loglog(energies, R, ':', color=handle[0].get_color())
            plt.figure(2)
            plt.loglog(energies, ion_range, color=handle[0].get_color())

    plt.figure(1)
    plt.legend(name)
    plt.axis([0, 1E5, 0, 3])
    plt.show()
    breakpoint()
