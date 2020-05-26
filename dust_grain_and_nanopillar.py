import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
from shapely.geometry import Point, Polygon, box
from itertools import chain
import os
import toml
import time
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
from scipy.ndimage import zoom, gaussian_filter
from matplotlib import rcParams
import matplotlib as mpl
rcParams.update({'figure.autolayout': True})

Q = 1.602E-19
PI = 3.14159
AMU = 1.66E-27
ANGSTROM = 1E-10
MICRON = 1E-6
NM = 1E-9
CM = 1E-2
EPS0 = 8.85E-12
A0 = 0.52918E-10
K = 1.11265E-10
ME = 9.11E-31
SQRTPI = 1.77245385
SQRT2PI = 2.506628274631
C = 299792000.

INTERPOLATED = 0
LOW_ENERGY_NONLOCAL = 1
LOW_ENERGY_LOCAL= 2
LOW_ENERGY_EQUIPARTITION = 3

LIQUID = 0
GASEOUS = 1

MOLIERE = 0
KR_C = 1
ZBL = 2
LENZ_JENSEN = 3
TRIDYN = -1

QUADRATURE = 0
MAGIC = 1

titanium = {
    'symbol': 'Ti',
    'name': 'titanium',
    'Z': 22,
    'm': 47.867,
    'Es': 4.84,
    'Ec': 3.5,
    'Eb': 0.,
    'Q': 0.54,
    'n': 5.67e28,
    'W': 2.57,
    's': 2.5
}

hydrogen = {
    'symbol': 'H',
    'name': 'hydrogen',
    'Z': 1,
    'm': 1.008,
    'Ec': 0.95,
    'Es': 1.5,
    'Eb': 0.,
    'n': 5.4E28
}

nitrogen = {
    'symbol': 'N',
    'name': 'nitrogen',
    'Z': 7,
    'm': 14,
    'n': 5.4E28,
    'Es': 0.,
    'Eb': 0.,
    'Ec': 1.0,
    'Q': 1.
}

deuterium = {
    'symbol': 'D',
    'name': 'deuterium',
    'Z': 1,
    'm': 2,
    'Ec': 0.95,
    'Es': 1.5,
}

helium = {
    'symbol': 'He',
    'name': 'helium',
    'Z': 2,
    'm': 4.002602,
    'Ec': 1.0,
    'Es': 0.
}

beryllium = {
    'symbol': 'Be',
    'name': 'beryllium',
    'Z': 4,
    'm': 9.012182,
    'n': 1.2347E29,
    'Es': 3.31,
    'Eb': 0.,
    'Ec': 3.,
    'Q': 1.66,
    'W': 2.32,
    's': 2.5
}

boron = {
    'symbol': 'B',
    'name': 'boron',
    'Z': 5,
    'm': 10.811,
    'n': 1.309E29,
    'Es': 5.76,
    'Eb': 0.,
    'Ec': 5.,
    'Q': 2.62,
    'W': 4.39,
    's': 2.5
}

arsenic = {
    'symbol': 'As',
    'name': 'arsenic',
    'Z': 33,
    'm': 74.921595,
    'n': 4.603E28,
    'Es': 3.12,
    'Eb': 0.,
    'Ec': 1.,
    'Q': None,
    'W': None,
    's': None
}

neon = {
    'symbol': 'Ne',
    'name': 'neon',
    'Z': 10,
    'm': 20.1797,
    'Ec': 1.0,
    'Es': 0.
}

krypton = {
    'symbol': 'Kr',
    'name': 'krypton',
    'Z': 36,
    'm': 83.80,
    'Ec': 1.0,
    'Es': 0.
}

silicon = {
    'symbol': 'Si',
    'name': 'silicon',
    'Z': 14,
    'm': 28.08553,
    'n': 4.996E28,
    'Es': 4.72,
    'Eb': 0.,
    'Ec': 1.5,
    'Q': 0.66,
    'W': 2.32,
    's': 2.5
}

argon = {
    'symbol': 'Ar',
    'name': 'argon',
    'Z': 18,
    'm': 39.948,
    'Ec': 1.0,
    'Es': 0.
}

oxygen = {
    'symbol': 'O',
    'name': 'oxygen',
    'Z': 8,
    'm': 15.9994,
    'Eb': 0.,
    'Ec': 2.0,
    'Es': 2.58,
    'n': 4.291E28
}

aluminum = {
    'symbol': 'Al',
    'name': 'aluminum',
    'Z': 13,
    'm': 26.98,
    'n': 6.026E28,
    'Es': 3.39,
    'Ec': 3.,
    'Eb': 3.,
}

copper = {
    'symbol': 'Cu',
    'name': 'copper',
    'Z': 29,
    'm': 63.546,
    'n': 8.491E28,
    'Es': 3.52,
    'Eb': 0.,
    'Ec': 3.0,
    'Q': 1.0,
    'W': 2.14,
    's': 2.5
}

tungsten = {
    'symbol': 'W',
    'name': 'tungsten',
    'Z': 74,
    'm': 183.84,
    'n': 6.306E28,
    'Es': 11.75,
    'Eb': 0.,
    'Ec': 6.,
    'Q': 0.72,
    'W': 2.14,
    's': 2.8
}

gold = {
    'symbol': 'Au',
    'name': 'gold',
    'Z': 79,
    'm': 196.97,
    'n': 5.901E28,
    'Es': 3.79,
    'Eb': 0.,
    'Ec': 3.,
    'Q': 1.0,
    'W': 0.,
    's': 0.,
}

nickel = {
    'symbol': 'Ni',
    'name': 'nickel',
    'Z': 28,
    'm': 58.69,
    'n': 9.14E28,
    'Es': 4.44,
    'Eb': 0.,
    'Ec': 3.,
    'Q': 1.0,
}

cesium = {
    'symbol': 'Cs',
    'name': 'cesium',
    'Z': 55,
    'm': 132.905,
    'Ec': 0.8,
    'Es': 0.8,
}

xenon = {
    'symbol': 'Xe',
    'name': 'xenon',
    'Z': 54,
    'm': 131.293,
    'Ec': 1.0,
    'Es': 0.
}

def circle(r, n, x0=0., y0=0.):
    theta = np.linspace(0., 2.*np.pi, n + 1)
    x = x0 + r*np.sin(theta)
    y = y0 + r*np.cos(theta)
    boundary = [(x_, y_) for x_, y_ in zip(x, y)]

    triangles = []
    for (x1, y1, x2, y2) in zip(x, y, x[1:], y[1:]):
        triangles.append([x0, x1, x2, y0, y1, y2])
    return triangles, boundary

def main(radius, beam, targets, E0, theta=0.0001, N=10000, N_=1, N_mesh=10):
    triangles, boundary = circle(radius, N_mesh)

    name = 'boron_dust_grain_'

    options = {
        'name': name,
        'track_trajectories': True,
        'track_recoils': True,
        'track_recoil_trajectories': True,
        'write_files': True,
        'stream_size': 8000,
        'print': True,
        'print_num': np.min((10, int(N_*N))),
        'weak_collision_order': 0,
        'suppress_deep_recoils': False,
        'high_energy_free_flight_paths': False,
        'electronic_stopping_mode': LOW_ENERGY_NONLOCAL,
        'mean_free_path_model': LIQUID,
        'interaction_potential': KR_C,
        'scattering_integral': QUADRATURE,
        'tolerance': 1E-3,
        'max_iterations': 100
    }

    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'Eb': [target['Eb'] for target in targets],
        'Es': [target['Es'] for target in targets],
        'Ec': [target['Ec'] for target in targets],
        'n': [target['n'] for target in targets],
        'Z': [target['Z'] for target in targets],
        'm': [target['m'] for target in targets],
        'electronic_stopping_correction_factor': 1.,
        'energy_barrier_thickness': 2.*(sum([target['n'] for target in targets])/len(targets))**(-1./3.)/np.sqrt(2.*np.pi)
    }

    #dx = 2.*n[0]**(-1./3.)/np.sqrt(2.*np.pi)/MICRON
    #dx = n**(-1./3.)*np.cos(theta)/MICRON
    dx = 10.*ANGSTROM/MICRON

    minx, miny, maxx, maxy = -radius, -radius, radius, radius
    surface = box(minx, miny, maxx, maxy)

    energy_surface = surface.buffer(dx, cap_style=2, join_style=2)
    simulation_surface = surface.buffer(10.*dx, cap_style=2, join_style=2)

    mesh_2d_input = {
        'length_unit': 'MICRON',
        'coordinate_sets': triangles,
        'data': [[target['n']/len(targets) for target in targets]]*N_mesh,
        'boundary_points': boundary,
        'simulation_boundary_points': list(simulation_surface.exterior.coords),
    }

    cosx = np.cos(theta*np.pi/180.)
    sinx = np.sin(theta*np.pi/180.)
    positions = [(-radius -dx, np.random.uniform(miny, maxy), 0.) for _ in range(N)]

    particle_parameters = {
        'length_unit': 'MICRON',
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'N': [N_ for _ in range(N)],
        'm': [beam['m'] for _ in range(N)],
        'Z': [beam['Z'] for _ in range(N)],
        'E': [E0 for _ in range(N)],
        'Ec': [beam['Ec'] for _ in range(N)],
        'Es': [beam['Es'] for _ in range(N)],
        'pos': positions,
        'dir': [(cosx, sinx, 0.) for _ in range(N)]
    }

    input_file = {
        'options': options,
        'material_parameters': material_parameters,
        'particle_parameters': particle_parameters,
        'mesh_2d_input': mesh_2d_input
    }

    with open('input.toml', 'w') as file:
        toml.dump(input_file, file, encoder=toml.TomlNumpyEncoder())

    os.system('rm *.output')
    os.system('rustBCA.exe')

    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.output', delimiter=','))
    trajectory_data = np.atleast_1d(np.genfromtxt(name+'trajectory_data.output', delimiter=',').transpose().astype(int))

    Z_max = np.max((beam['Z'], np.max([target['Z'] for target in targets])))
    Z_min = np.min((beam['Z'], np.min([target['Z'] for target in targets])))
    cmap = matplotlib.cm.get_cmap('magma')

    x_max = 0.
    index = 0
    if np.size(trajectories) > 0:
        plt.figure(num='trajectories')
        x_boundary = [x for x, _ in boundary]
        y_boundary = [y for _, y in boundary]
        plt.fill(x_boundary, y_boundary, color='dimgrey', linewidth=0.)

        for trajectory_length in trajectory_data:

            M = trajectories[index, 0]
            Z = trajectories[index, 1]
            E = trajectories[index:(trajectory_length + index), 2]
            x = trajectories[index:(trajectory_length + index), 3]
            y = trajectories[index:(trajectory_length + index), 4]
            z = trajectories[index:(trajectory_length + index), 5]

            if np.max(x) > x_max:
                x_max = np.max(x)

            #plt.scatter(x[0], y[0], color = cmap(Z/Z_max), marker='o', s=5)
            plt.plot(x, y, color = cmap((Z - Z_min)/(Z_max - Z_min)), linewidth = 1)
            index += trajectory_length

        if np.size(sputtered) > 0:
            sputtered_colors = [cmap((Z - Z_min)/(Z_max - Z_min)) for Z in sputtered[:,1]]
            plt.scatter(sputtered[:,3], sputtered[:,4], s=50, color=sputtered_colors, marker='*')

        if np.size(reflected) > 0:
            reflected_colors = [cmap((Z - Z_min)/(Z_max - Z_min)) for Z in reflected[:,1]]
            plt.scatter(reflected[:,3], reflected[:,4], s=50, color=reflected_colors, marker='x')

        if np.size(deposited) > 0:
            deposited_colors = [cmap((Z - Z_min)/(Z_max - Z_min)) for Z in deposited[:,1]]
            plt.scatter(deposited[:,2], deposited[:,3], s=50, color=deposited_colors, marker='^')

        x_box = [minx - 10.*dx, minx - 10.*dx, maxx + 10.*dx, maxx + 10.*dx, minx - 10.*dx]
        y_box = [miny - 10.*dx, maxy + 10.*dx, maxy + 10.*dx, miny - 10.*dx, miny - 10.*dx]
        plt.plot(x_box, y_box, color='dimgray', linewidth=3)

        plt.title(f'{E0} eV {beam["symbol"]} on {"".join([target["symbol"]+" " for target in targets])}Dust Grain')
        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        plt.axis('square')
        plt.xticks([-radius, 0., radius])
        plt.yticks([-radius, 0., radius])

    if np.size(deposited) > 0:
        plt.figure(num='dep')
        num_bins = 200
        bins_x = np.linspace(-1.1*radius, 1.1*radius, num_bins)
        bins_y = np.linspace(-1.1*radius, 1.1*radius, num_bins)
        plt.hist2d(deposited[:, 2], deposited[:, 3], bins=(bins_x, bins_y))
        x_boundary = [x for x, _ in boundary]
        y_boundary = [y for _, y in boundary]
        plt.fill(x_boundary, y_boundary, color='white', linewidth=1., alpha=0.1)
        plt.title(f'{E0} eV {beam["symbol"]} on {"".join([target["symbol"]+" " for target in targets])}Dust Grain')
        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        plt.axis('square')
        plt.xticks([-radius, 0., radius])
        plt.yticks([-radius, 0., radius])

    if np.size(sputtered) > 0:
        plt.figure(num='sput_orig')
        num_bins = 200
        bins_x = np.linspace(-1.1*radius, 1.1*radius, num_bins)
        bins_y = np.linspace(-1.1*radius, 1.1*radius, num_bins)
        plt.hist2d(sputtered[:, -3], sputtered[:, -2], bins=(bins_x, bins_y))
        x_boundary = [x for x, _ in boundary]
        y_boundary = [y for _, y in boundary]
        plt.fill(x_boundary, y_boundary, color='white', linewidth=1., alpha=0.1)
        plt.title(f'{E0} eV {beam["symbol"]} on {"".join([target["symbol"]+" " for target in targets])}Dust Grain')
        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        plt.axis('square')
        plt.xticks([-radius, 0., radius])
        plt.yticks([-radius, 0., radius])

    plt.show()

    if np.size(sputtered) > 0:
        Y = len(sputtered[:,0])/N/N_
    else:
        Y = 0

    if np.size(reflected) > 0:
        R = len(reflected[:,0])/N/N_
    else:
        R = 0

    if np.size(deposited) > 0:
        deposited_below_surface = deposited[:, 2] > 0.
        ion_range = np.mean(deposited[deposited_below_surface, 2])
        ion_straggle = np.std(deposited[deposited_below_surface,2])
    else:
        ion_range = 0
        ion_straggle = 0

    return Y, R, ion_range, ion_straggle

if __name__ == '__main__':
    print(main(0.01, helium, [hydrogen, beryllium, boron, nitrogen], 10000., N=10, N_mesh=100))
