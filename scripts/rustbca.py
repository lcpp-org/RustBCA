import numpy as np
from shapely.geometry import Point, Polygon, box
import os
import toml
import time
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib import rcParams, cm
import matplotlib as mpl
import matplotlib.colors as colors

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

titanium = {
    'symbol': 'Ti',
    'name': 'titanium',
    'Z': 22,
    'm': 47.867,
    'Es': 4.84,
    'Ec': 3.5,
    'Eb': 3.,
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
}

nitrogen = {
    'symbol': 'N',
    'name': 'nitrogen',
    'Z': 7,
    'm': 14,
    'n': 5.4E25,
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
    'Ec': 0.1,
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
    'Eb': 2.58,
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
    'Es': 8.79,
    'Eb': 3.,
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

INTERPOLATED = "INTERPOLATED"
LOW_ENERGY_NONLOCAL = "LOW_ENERGY_NONLOCAL"
LOW_ENERGY_LOCAL= "LOW_ENERGY_LOCAL"
LOW_ENERGY_EQUIPARTITION = "LOW_ENERGY_EQUIPARTITION"

LIQUID = "LIQUID"
GASEOUS = "GASEOUS"

MOLIERE = "MOLIERE"
KR_C = "KR_C"
ZBL = "ZBL"
LENZ_JENSEN = "LENZ_JENSEN"

QUADRATURE = "MENDENHALL_WELLER"
MAGIC = "MAGIC"

def thomas_reflection(ion, target, energy_eV):
    '''
    Thomas et al. (1991) semi-empirical reflection coefficient.

    Args:
        ion (dict): a dictionary with the fields Z (atomic number), m (mass)
        target (dict): a dictionary with the fields Z (atomic number), m (mass)
        energy_eV (float): energy in electron-volts

    Returns:
        R (float): reflection coefficient of ion on target with energy_eV
    '''
    #Thomas et al. empirical reflection coefficient (1991)
    Z1 = ion['Z']
    Z2 = target['Z']
    M1 = ion['m']
    M2 = target['m']
    energy_keV = energy_eV/1E3

    #Thomas-Fermi reduced energy
    reduced_energy = 32.55*energy_keV*M2/((M1 + M2)*Z1*Z2*(Z1**0.23 + Z2**0.23))

    mu = M2/M1
    if mu < 1:
        print('Warning: Thomas et al. reflection coefficient not defined for M2/M1 < 1.')
        return None

    mu_ranges = [1, 3, 6, 7, 12, 15, 20]
    A1 = [0.02129, 0.3680, 0.5173, 0.5173, 0.6192, 0.6192, 0.8250]
    A2 = [16.39, 2.985, 2.549, 2.549, 20.01, 20.01, 21.41]
    A3 = [26.39, 7.122, 5.325, 5.325, 8.922, 8.922, 8.606]
    A4 = [0.9131, 0.5802, 0.5719, 0.5719, 0.6669, 0.6669, 0.6425]
    A5 = [6.249, 4.211, 1.094, 1.094, 1.864, 1.864, 1.907]
    A6 = [2.550, 1.597, 1.933, 1.933, 1.899, 1.899, 1.927]

    a1 = interp1d(mu_ranges, A1, bounds_error = False, fill_value=(A1[0], A1[-1]))
    a2 = interp1d(mu_ranges, A2, bounds_error = False, fill_value=(A2[0], A2[-1]))
    a3 = interp1d(mu_ranges, A3, bounds_error = False, fill_value=(A3[0], A3[-1]))
    a4 = interp1d(mu_ranges, A4, bounds_error = False, fill_value=(A4[0], A4[-1]))
    a5 = interp1d(mu_ranges, A5, bounds_error = False, fill_value=(A5[0], A5[-1]))
    a6 = interp1d(mu_ranges, A6, bounds_error = False, fill_value=(A6[0], A6[-1]))

    return a1(mu)*np.log(a2(mu)*reduced_energy + 2.718)/(1. + a3(mu)*reduced_energy**a4(mu) + a5(mu)*reduced_energy**a6(mu))

def wierzbicki_biersack(ion, target, energy_eV):
    '''
    Wierzbicki-Biersack empirical reflection coefficient (1994); not as widely
        applicable as Thomas et al.

    Args:
        ion (dict): a dictionary with the fields Z (atomic number), m (mass)
        target (dict): a dictionary with the fields Z (atomic number), m (mass)
        energy_eV (float): energy in electron-volts

    Returns:
        R (float): reflection coefficient of ion on target with energy_eV
    '''
    #Wierzbicki and Biersack empirical reflection coefficient (1994)
    Z1 = ion['Z']
    Z2 = target['Z']
    M1 = ion['m']
    M2 = target['m']
    energy_keV = energy_eV/1E3

    #I've never seen this form of the reduced energy before - it's Thomas-Fermi
    reduced_energy = 32.55*energy_keV*M2/((M1 + M2)*Z1*Z2*(Z1**0.23 + Z2**0.23))
    mu = M2/M1

    #Here are some empirical coefficients
    a1 = 5.9638
    b1 = 0.0646
    c1 = 52.211

    a2 = 1.26E-3
    b2 = -0.9305
    c2 = 1.235

    #Wierzbicki and Biersack found that you can separate the dependence on mu, e
    RN_mu = np.exp(a1*np.sqrt(1. - b1*(np.log(mu/c1))**2.))
    RN_e = a2*np.exp(b2*(np.log(reduced_energy + 1.))**c2)

    if not (1.03 < mu <= 240):
        print("Warning: Wierzbicki-Biersack may not be accurate for this ion-target pair")
        print(f'False: 1.03 < {mu} <= 240')

    if not (1 < reduced_energy < 10):
        print("Warning: Wierzbicki-Biersack may not be accurate at this energy")
        print(f'False: 1 < {reduced_energy} <= 10')

    return RN_mu*RN_e

def bohdansky_light_ion(ion, target, energy_eV):
    '''
    Bohdansky sputtering yield formula in the light ion (M1/M2 < 0.5) limit.
    Returns 0 if the target does not have a surface binding energy.

    Args:
        ion (dict): a dictionary with the fields Z (atomic number), m (mass)
        target (dict): a dictionary with the fields Z (atomic number), m (mass), Es (surface binding energy)
        energy_eV (float): energy in electron-volts

    Returns:
        Y (float): sputtering yield in atoms/ion
    '''
    z1 = ion['Z']
    z2 = target['Z']
    m1 = ion['m']
    m2 = target['m']
    Us = target['Es']

    if Us == 0.: return 0
    alpha = 0.2

    reduced_mass_2 = m2/(m1 + m2)
    reduced_mass_1 = m1/(m1 + m2)

    #Following assumptions are for very light ions (m1/m2<0.5)
    K = 0.4
    R_Rp = K*m2/m1 + 1.

    Eth = (1.9 + 3.8*(m1/m2) + 0.134*(m2/m1)**1.24)*Us

    a0 = 0.529*ANGSTROM
    a = 0.885*a0*(z1**(2./3.) + z2**(2./3.))**(-1./2.)
    reduced_energy = 0.03255/(z1*z2*(z1**(2./3.) + z2**(2./3.))**(1./2.))*reduced_mass_2*energy_eV
    sn = 3.441*np.sqrt(reduced_energy)*np.log(reduced_energy + 2.718)/(1. + 6.355*np.sqrt(reduced_energy) + reduced_energy*(-1.708 + 6.882*np.sqrt(reduced_energy)))
    Sn = 8.478*z1*z2/(z1**(2./3.) + z2**(2./3.))**(1./2.)*reduced_mass_1*sn

    sputtering_yield = 0.042/Us*(R_Rp)*alpha*Sn*(1-(Eth/energy_eV)**(2./3.))*(1-(Eth/energy_eV))**2
    if sputtering_yield > 0:
        return sputtering_yield
    else:
        return 0.

def yamamura(ion, target, energy_eV):
    '''
    Yamamura sputtering yield formula for normal incidence.

    Args:
        ion (dict): a dictionary with the fields Z (atomic number), m (mass)
        target (dict): a dictionary with the fields Z (atomic number), m (mass), Es (surface binding energy), Q (Yamamura coefficient)
        energy_eV (float): energy in electron-volts

    Returns:
        Y (float): sputtering yield in atoms/ion
    '''
    #Yamamura sputtering yield implementation
    z1 = ion['Z']
    z2 = target['Z']
    m1 = ion['m']
    m2 = target['m']
    Us = target['Es']
    Q = target['Q']
    #W = target['W']
    #s = target['s']

    reduced_mass_2 = m2/(m1 + m2)
    reduced_mass_1 = m1/(m1 + m2)

    #Lindhard's reduced energy
    reduced_energy = 0.03255/(z1*z2*(z1**(2./3.) + z2**(2./3.))**(1./2.))*reduced_mass_2*energy_eV

    #Yamamura empirical constants
    K = 8.478*z1*z2/(z1**(2./3.) + z2**(2./3.))**(1./2.)*reduced_mass_1
    a_star = 0.08 + 0.164*(m2/m1)**0.4 + 0.0145*(m2/m1)**1.29
    #Sputtering threshold energy
    Eth = (1.9 + 3.8*(m1/m2) + 0.134*(m2/m1)**1.24)*Us
    #Lindhard-Scharff-Schiott nuclear cross section
    sn = 3.441*np.sqrt(reduced_energy)*np.log(reduced_energy + 2.718)/(1. + 6.355*np.sqrt(reduced_energy) + reduced_energy*(-1.708 + 6.882*np.sqrt(reduced_energy)))
    #Lindhard-Scharff electronic cross section
    k = 0.079*(m1 + m2)**(3./2.)/(m1**(3./2.)*m2**(1./2.))*z1**(2./3.)*z2**(1./2.)/(z1**(2./3.) + z2**(2./3.))**(3./4.)
    se = k*np.sqrt(reduced_energy)

    return 0.42*a_star*Q*K*sn/Us/(1. + 0.35*Us*se)*(1. - np.sqrt(Eth/energy_eV))**2.8

def do_trajectory_plot(name, thickness=None, depth=None, boundary=None, plot_final_positions=True, plot_origins=True, show=True):
    '''
    Plots trajectories of ions and recoils from [name]trajectories.output.
    Optionally marks final positions/origins and draws material geometry.

    Geometry input is in same length_unit as the rustbca particles.

    Args:
        name (string): name of rustbca simulation
        thickness list(float): thickness of target, or None to not draw target
        depth list(float): depth of target, or None to not draw target
        boundary list((float, float)): points that make up boundary, or None to not draw boundary
        plot_final_positions (bool): mark final positions (reflected: X sputtered: * deposited ^)
        plot_origins (bool): mark originating locations of particles (o)
        show (bool): whether or not to show plots

    '''

    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.output', delimiter=','))
    trajectory_data = np.atleast_1d(np.genfromtxt(name+'trajectory_data.output', delimiter=',').transpose().astype(int))

    if np.size(trajectories) > 0:
        min_Z = np.min(trajectories[:, 1])
        max_Z = np.max(trajectories[:, 1])
        colormap = cm.ScalarMappable(norm=colors.Normalize(vmin=min_Z, vmax=max_Z), cmap='tab20b')

    fig1, axis1 = plt.subplots()

    index = 0
    x_max = 0


    if np.size(trajectories) > 0:
        for trajectory_length in trajectory_data:

            M = trajectories[index, 0]
            Z = trajectories[index, 1]
            E = trajectories[index:(trajectory_length + index), 2]
            x = trajectories[index:(trajectory_length + index), 3]
            y = trajectories[index:(trajectory_length + index), 4]
            z = trajectories[index:(trajectory_length + index), 5]

            if np.max(x) > x_max:
                x_max = np.max(x)

            if plot_origins: plt.scatter(x[0], y[0], color=colormap.to_rgba(Z), marker='.', s=5)
            plt.plot(x, y, color = colormap.to_rgba(Z), linewidth = 1)

            index += trajectory_length

        if plot_final_positions:
            if np.size(sputtered) > 0:
                sputtered_colors = [colormap.to_rgba(Z) for Z in sputtered[:,1]]
                plt.scatter(sputtered[:,3], sputtered[:,4], s=50, color=sputtered_colors, marker='*')

            if np.size(reflected) > 0:
                reflected_colors = [colormap.to_rgba(Z) for Z in reflected[:,1]]
                plt.scatter(reflected[:,3], reflected[:,4], s=50, color=reflected_colors, marker='x')

            if np.size(deposited) > 0:
                deposited_colors = [colormap.to_rgba(Z) for Z in deposited[:,1]]
                plt.scatter(deposited[:,2], deposited[:,3], s=50, color=deposited_colors, marker='^')

        if thickness and depth:
            x_box = [0., 0., depth, depth, 0.]
            y_box = [-thickness/2., thickness/2., thickness/2., -thickness/2., -thickness/2.]
            plt.plot(x_box, y_box, color='dimgray', linewidth=3)

        elif boundary:
            for point_1, point_2 in zip(boundary[1:], boundary[:-1]):
                plt.plot([point_1[0], point_2[0]], [point_1[1], point_2[1]], line_width=3, color='dimgray')
            plt.plot([boundary[0][0], boundary[-1][0]], [boundary[0][1], boundary[-1][1]], line_width=3, color='dimgray')

        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        plt.title(name+' Trajectories')
        plt.axis('square')
        if show: plt.show()
        plt.savefig(name+'trajectories_.png')
        plt.close()

def generate_rustbca_input(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, E0, N, N_, theta,
    thickness, depth, track_trajectories=True, track_recoils=True,
    track_recoil_trajectories=True, name='test_',
    track_displacements=True, track_energy_losses=True,
    electronic_stopping_mode=LOW_ENERGY_NONLOCAL,
    weak_collision_order=0, ck=1., mean_free_path_model=LIQUID,
    interaction_potential='"KR_C"', high_energy=False, energy_barrier_thickness=(6.306E10)**(-1/3.),
    initial_particle_position = -1*ANGSTROM/MICRON, integral='"MENDENHALL_WELLER"',
    root_finder = '{"NEWTON"={max_iterations=100, tolerance=1e-3}}',
    delta_x_angstrom=5.):


    '''
    Generates a rustbca input file. Assumes eV, amu, and microns for units.
    '''

    options = {
        'name': name,
        'track_trajectories': track_trajectories,
        'track_recoils': track_recoils,
        'track_recoil_trajectories': track_recoil_trajectories,
        'stream_size': 8000,
        'weak_collision_order': weak_collision_order,
        'suppress_deep_recoils': False,
        'high_energy_free_flight_paths': high_energy,
        'num_threads': 4,
        'num_chunks': 10,
        'use_hdf5': False,
        'electronic_stopping_mode': electronic_stopping_mode,
        'mean_free_path_model': LIQUID,
        'track_displacements': track_displacements,
        'track_energy_losses': track_energy_losses,
    }

    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'Eb': Eb,
        'Es': Esb,
        'Ec': Ecb,
        'Z': Zb,
        'm': Mb,
        'interaction_index': np.zeros(len(n), dtype=int),
        'electronic_stopping_correction_factor': ck,
        'surface_binding_model': "AVERAGE"
    }

    dx = delta_x_angstrom*ANGSTROM/MICRON

    minx, miny, maxx, maxy = 0.0, -thickness/2., depth, thickness/2.
    surface = box(minx, miny, maxx, maxy)

    simulation_surface = surface.buffer(dx, cap_style=2, join_style=2)
    mesh_2d_input = {
        'length_unit': 'MICRON',
        'coordinate_sets': [[0., depth, 0., thickness/2., -thickness/2., -thickness/2.], [0., depth, depth, thickness/2., thickness/2., -thickness/2.]],
        'densities': [np.array(n)*(MICRON)**3, np.array(n)*(MICRON)**3],
        'boundary_points': [[0., thickness/2.], [depth, thickness/2.], [depth, -thickness/2.], [0., -thickness/2.]],
        'simulation_boundary_points':  list(simulation_surface.exterior.coords),
        'energy_barrier_thickness': energy_barrier_thickness
    }

    cosx = np.cos(theta*np.pi/180.)
    sinx = np.sin(theta*np.pi/180.)
    positions = [(initial_particle_position, 0., 0.) for _ in range(N)]

    particle_parameters = {
        'length_unit': 'MICRON',
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'N': [N_ for _ in range(N)],
        'm': [Ma for _ in range(N)],
        'Z': [Za for _ in range(N)],
        'E': [E0 for _ in range(N)],
        'Ec': [Eca for _ in range(N)],
        'Es': [Esa for _ in range(N)],
        'interaction_index': np.zeros(N, dtype=int),
        'pos': positions,
        'dir': [(cosx, sinx, 0.) for _ in range(N)],
        'particle_input_filename': ''
    }

    input_file = {
        'material_parameters': material_parameters,
        'particle_parameters': particle_parameters,
        'mesh_2d_input': mesh_2d_input,
        'options': options,
    }

    with open(f'{name}.toml', 'w') as file:
        toml.dump(input_file, file, encoder=toml.TomlNumpyEncoder())

    with open(f'{name}.toml', 'a') as file:
        file.write(f'root_finder = [[{root_finder}]]\n')
        file.write(f'interaction_potential = [[{interaction_potential}]]\n')
        file.write(f'scattering_integral =  [[{integral}]]\n')

def plot_distributions_rustbca(name, beam, target,
    incident_energy=1, incident_angle=0,
    max_collision_contours=4, plot_2d_reflected_contours=False,
    collision_contour_significance_threshold=0.1, plot_garrison_contours=False,
    plot_reflected_energies_by_number_collisions=False,
    plot_scattering_energy_curve=False):

    '''
    Plots rustbca distributions.

    Args:
        name (string): name of rustbca simulation.
        beam (dict): ions striking target; dictionary with fields symbol and Z
        target (dict): target upon which ions are incdient; dictionary with fields symbol and Z
        incident_energy (float): energy in energy_units; used to scale energy spectra
        incident_angle (float): angle in degrees; currently used to generate titles only
        max_collision_contours (int): number of collision event contours to plot on reflected EAD; default 4
        plot_2d_reflected_contours (bool): whether to plot collision event contours on reflected EAD; default False
        collision_contour_significance_threshold (int): threshold of significance for contour plotting; default 0.1
        plot_reflected_energies_by_number_collisions (bool): separte reflected energy spectrum into collision number distributions; default False
        plot_scattering_energy_curve (bool): plot bold, translucent white curve showing theoretical single-collision reflected energies; default False

    '''
    num_bins = 120

    r = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    s = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    d = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))

    if np.size(r) > 0:
        ux = r[:,6]
        uy = r[:,7]
        uz = r[:,8]
        theta = -np.arctan2(ux, np.sqrt(uy**2 + uz**2))*180./np.pi
        theta = 180. - np.arccos(ux)*180./np.pi
        #theta = (np.arccos(r[:,7]) - np.pi/2.)*180./np.pi
        number_collision_events = r[:, -1]

        fig = plt.figure(num='polar_scatter')
        ax = fig.add_subplot(111, projection='polar')
        normalized_energies_rust = r[:,2]/incident_energy

        c2 = ax.scatter(np.arccos(r[:,7]), normalized_energies_rust, s=1, c=number_collision_events)

        plt.legend(['ftridyn'], loc='upper left')
        ax.set_thetamax(180.)
        ax.set_thetamin(0.)
        ax.set_xlabel('E/E0')
        ax.set_yticks([0., 0.5, 1.])
        ax.set_xticks([0., np.pi/6., 2*np.pi/6., 3*np.pi/6., 4*np.pi/6., 5*np.pi/6., np.pi])
        ax.set_xticklabels(['0', '30', '60', '90', '', '', ''])
        plt.title(f'{beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.savefig(name+'polar_scatter.png')
        plt.close()

        plt.figure(num='rustbca_r2d')
        bin_energy = np.linspace(0., 1.2*np.max(r[:, 2]), num_bins)
        bin_angle = np.linspace(0., 90., num_bins)

        bins = (bin_angle, bin_energy)
        heights, xedges, yedges, image = plt.hist2d(theta, r[:, 2], bins=bins)
        np.savetxt(name+'refl_ead.dat', heights)
        np.savetxt(name+'refl_energies.dat', yedges)
        np.savetxt(name+'refl_angles.dat', xedges)

        if plot_2d_reflected_contours:
            n_max = min(int(np.max(number_collision_events)), max_collision_contours)
            cmap = mpl.cm.get_cmap('Wistia')
            colors = [cmap((n - 1)/n_max) for n in range(1, n_max + 1)]

            for k in range(1, n_max + 1):
                if k < n_max:
                    mask = number_collision_events == k
                else:
                    mask = number_collision_events > k

                heights, xedges, yedges = np.histogram2d(r[mask, 2], theta[mask], density=True, bins=num_bins//2)
                x_centers = (xedges[1:] + xedges[:-1])/2.
                y_centers = (yedges[1:] + yedges[:-1])/2.
                plt.contour(x_centers, y_centers, heights.transpose()/np.max(heights), levels=np.linspace(collision_contour_significance_threshold, 1., 1), colors=[colors[k - 1]], linewidths=1, linestyles='--', alpha=0.75)

            norm = mpl.colors.Normalize(vmin=1, vmax=n_max)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ticks=np.array(range(1, n_max + 1))-0.5, boundaries=np.array(range(0, n_max + 1)))
            cbar.set_label('# of Collision Events')
            cbar_labels = list(range(1, n_max + 1))
            cbar_labels[-1] = '≥'+str(cbar_labels[-1])
            cbar.ax.set_yticklabels(cbar_labels)

        if plot_scattering_energy_curve:
            for mass in [titanium['m'], oxygen['m']]:
                energies = bins[1]
                angles = bins[0]*np.pi/180.
                scattering_angles = -angles + np.pi
                final_energies = incident_energy*((np.cos(scattering_angles) + np.sqrt((mass/beam['m'])**2 - np.sin(scattering_angles)**2))/(1. + mass/beam['m']))**2.
                handle = plt.plot(angles*180./np.pi, final_energies,  linestyle='-', color='white', alpha=0.25, linewidth=7)
                plt.legend(handle, ['Single-Collision Scattering'], loc='lower right', fontsize='x-small')

        plt.ylabel('E [eV]')
        plt.xlabel('angle [deg]')
        plt.title(f'Reflected EAD {beam["symbol"]} on {target["symbol"]}' , fontsize='small')
        plt.savefig(name+'rustbca_r_ead.png')
        plt.close()

    if np.size(s) > 0:
        plt.figure(num='rustbca_s2d')

        bin_energy = np.linspace(0., 1.2*np.max(s[:, 2]), num_bins//2)
        bin_angle = np.linspace(0., 90., num_bins//2)
        bins = (bin_angle, bin_energy)
        ux = s[:,6]
        uy = s[:,7]
        uz = s[:,8]
        #theta = -np.arctan2(ux, np.sqrt(uy**2 + uz**2))
        theta = 180. - np.arccos(ux)*180./np.pi

        heights, xedges, yedges, image = plt.hist2d(theta, s[:, 2], bins=bins)
        np.savetxt(name+'sput_ead.dat', heights)
        np.savetxt(name+'sput_energies.dat', yedges)
        np.savetxt(name+'sput_angles.dat', xedges)

        plt.ylabel('E [eV]')
        plt.xlabel('angle [deg]')
        plt.title(f'Sputtered EAD {beam["symbol"]} on {target["symbol"]}' , fontsize='small')
        plt.savefig(name+'rustbca_s_ead.png')
        plt.close()

    #Deposited ion depth distributions
    if np.size(d) > 0:
        plt.figure(num='d')
        plots = []
        labels = []

        heights, bins, rust = plt.hist(d[d[:,2]>0., 2], histtype='step', bins=num_bins, density=True, color='black')
        np.savetxt(name+'depo_dist.dat', heights)
        np.savetxt(name+'depo_depths.dat', bins)
        plots.append(rust[0])
        labels.append('rustbca')

        if len(plots) > 0: plt.legend(plots, labels, fontsize='small', fancybox=True, shadow=True, loc='upper right')
        plt.title(f'Depth Distributions {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.xlabel('x [um]')
        plt.ylabel('f(x)')
        #axis = plt.gca()
        #axis.set_yscale('log')
        plt.savefig(name+'dep.png')
        plt.close()

    #Reflected ion energy distributions
    if np.size(r) > 0:

        plt.figure(num='re')
        plots = []
        labels = []

        heights, bins, rust = plt.hist(r[:, 2]/incident_energy, histtype='step',  bins=num_bins, density=True, color='black')
        number_collision_events = r[:, -1]

        n_max = min(int(np.max(number_collision_events)), max_collision_contours)
        cmap = mpl.cm.get_cmap('brg')
        colors = [cmap((n - 1)/n_max) for n in range(1, n_max + 1)]

        if plot_reflected_energies_by_number_collisions:
            for k in range(1, n_max + 1):
                if k < n_max:
                    mask = number_collision_events == k
                else:
                    mask = number_collision_events > k

                plt.hist(r[mask, 2]/incident_energy, histtype='step', bins = np.linspace(0.0, 1.0, num_bins), density=True, color=colors[k - 1])
            norm = mpl.colors.Normalize(vmin=1, vmax=n_max)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ticks=np.array(range(1, n_max + 1))-0.5, boundaries=np.array(range(0, n_max + 1)))
            cbar.set_label('# of Collision Events')
            cbar_labels = list(range(1, n_max + 1))
            cbar_labels[-1] = '≥'+str(cbar_labels[-1])
            cbar.ax.set_yticklabels(cbar_labels)

        plots.append(rust[0])
        labels.append('rustbca')

        if len(plots) > 0: plt.legend(plots, labels, fontsize='small', fancybox=True, shadow=True, loc='upper left')
        plt.title(f'Refl. Energies {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.xlabel('E/E0')
        plt.ylabel('f(E)')
        plt.gca().set_xlim(0., 1.)
        plt.gca().set_ylim(bottom=0)

        plt.savefig(name+'ref_e.png')
        plt.close()


    if np.size(s) > 0:
        #Sputtered atom energy distributions
        plt.figure(num='se')
        plots = []
        labels = []

        _, _, rust = plt.hist(s[:,2], histtype='step', bins=num_bins, density=True, color='black')
        plots.append(rust[0])
        labels.append('rustbca')

        energies = np.linspace(0., np.max(s[:,2]), num_bins*10)
        de = energies[1] - energies[0]
        thompson = energies/(energies + target['Es'])**(3)
        thompson /= (np.sum(thompson)*de)
        thompson = plt.plot(energies, thompson, linestyle='-.')
        plots.append(thompson[0])
        labels.append('Thompson (n=2)')

        if len(plots) > 0: plt.legend(plots, labels, fontsize='small', fancybox=True, shadow=True, loc='upper right')
        plt.title(f'Sputtered Energies {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.xlabel('E [eV]')
        plt.ylabel('f(E)')
        plt.savefig(name+'spt_e.png')
        plt.close()

        plt.figure(num='so')
        depth_origin = s[:, 10]
        plt.hist(depth_origin, bins=num_bins, histtype='step', color='black')
        plt.title('Depth of Origin of Sputtered Particles')
        plt.xlabel('x [um]')
        plt.ylabel('f(x)')
        plt.savefig(name+'spt_o.png')
        plt.close()

    #Sputtered atom angular distributions
    plt.figure(num='sa')
    plots = []
    labels = []
    ax = None
    bins = np.linspace(0, np.pi, num_bins)

    if np.size(s) > 0:
        hist, bins = np.histogram(np.arccos(s[:,6]), bins=bins, density=True)
        rust1 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust1[0])
        labels.append('rustbca cosx')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(s[:,7]), bins=bins, density=True)
        rust2 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust2[0])
        labels.append('rustbca cosy')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(s[:,8]), bins=bins, density=True)
        rust3 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust3[0])
        labels.append('rustbca cosz')
        ax = plt.gca()

    if len(plots) > 0: plt.legend(plots, labels, fontsize='small', bbox_to_anchor=(0.0, 1.05), fancybox=True, shadow=True, loc='upper left')
    if ax:
        ax.set_thetamin(0.)
        ax.set_thetamax(180.)
        ax.set_yticks([0., 0.5, 1.])
        ax.set_xticks([0., np.pi/6., 2*np.pi/6., 3*np.pi/6., 4*np.pi/6., 5*np.pi/6., np.pi])
        ax.set_xticklabels(['0', '30', '60', '90', '', '', ''])
    plt.title(f'Sputtered Angles {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
    plt.savefig(name+'spt_ang.png')
    plt.close()

    #Reflected atom angular distributions
    plt.figure(num='ra')
    plots = []
    labels = []
    ax = None
    bins = np.linspace(0, np.pi, num_bins)

    if np.size(r) > 0:
        hist, bins = np.histogram(np.arccos(-r[:,6]), bins=bins, density=True)
        rust1 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust1[0])
        labels.append('rustbca cosx')

        hist, bins = np.histogram(np.arccos(r[:,7]), bins=bins, density=True)
        rust2 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust2[0])
        labels.append('rustbca cosy')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(r[:,8]), bins=bins, density=True)
        rust3 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust3[0])
        labels.append('rustbca cosz')
        ax = plt.gca()

    if len(plots) > 0: plt.legend(plots, labels, fontsize='small', fancybox=True, shadow=True, bbox_to_anchor=(-0.25, 0.5), loc='upper left')
    if ax:
        ax.set_thetamin(0.)
        ax.set_thetamax(180.)
        ax.set_yticks([0., 0.5, 1.])
        ax.set_xticks([0., np.pi/6., 2*np.pi/6., 3*np.pi/6., 4*np.pi/6., 5*np.pi/6., np.pi])
        ax.set_xticklabels(['0', '30', '60', '90', '', '', ''])
    plt.title(f'Refl. Angles {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
    plt.savefig(name+'ref_ang.png')
    plt.close()

def plot_displacements(name, displacement_energy, num_bins=100):
    '''
    Plots displacements given a threshold displacement energy.

    Args:
        name (string): name of rustbca simulation.
        displacement_energy (float): threshold energy needed to produce a permanent frenkel pair in energy_units
        num_bins (int): number of histogram bins to use; default 100
    '''
    displacements = np.genfromtxt(f'{name}displacements.output', delimiter=',')
    M = displacements[:,0]
    Z = displacements[:,1]
    Er = displacements[:,2]
    x = displacements[Er>displacement_energy,3]
    y = displacements[Er>displacement_energy,4]

    plt.figure(1)
    plt.hist(x, bins=num_bins, histtype='step', color='black', density=True)
    plt.xlabel('x [um]')
    plt.ylabel('f(x)')
    plt.title(f'Frenkel Pair Production Distribution {name}')
    plt.savefig(name+'displacements.png')
    plt.close()

    plt.figure(2)
    plt.hist2d(x, y, bins=num_bins)
    plt.xlabel('x [um]')
    plt.ylabel('y [um]')
    plt.title(f'Frenkel Pair Production {name}')

    plt.savefig(name+'displacements2D.png')
    plt.close()

def plot_energy_loss(name, N, num_bins=100, thickness=None, depth=None):
    '''
    Plots energy loss plots separated by electronic and nuclear losses from rustbca.

    Args:
        name (string): name of rustbca simulation
        N (int): number of incident ions
        num_bins (int): number of histogram bins; default 100
    '''

    energy_loss = np.genfromtxt(f'{name}energy_loss.output', delimiter=',')
    M = energy_loss[:,0]
    Z = energy_loss[:,1]
    En = energy_loss[:,2]
    Ee = energy_loss[:,3]
    x = energy_loss[:,4]
    y = energy_loss[:,5]

    plt.figure(1)
    plt.hist(x, weights=En/N, bins=num_bins, histtype='step', color='black')
    plt.hist(x, weights=Ee/N, bins=num_bins, histtype='step', color='red')

    plt.xlabel('x [um]')
    plt.ylabel('E [eV/ion]')
    plt.legend(['Nuclear', 'Electronic'])
    plt.title(f'Energy Losses {name}')
    plt.savefig(name+'energy_loss.png')
    plt.close()

    if depth and thickness:
        num_bins = (np.linspace(0., depth, num_bins), np.linspace(-thickness/2., thickness/2., num_bins))

    plt.figure(2)
    plt.xlabel('x [um]')
    plt.ylabel('y [um]')
    plt.title(f'Nuclear Energy Loss {name}')
    plt.hist2d(x, y, weights=En/N, bins=num_bins, norm=colors.SymLogNorm(0.1))

    plt.savefig(name+'nuclear_energy_loss.png')
    plt.close()

    plt.figure(3)
    plt.xlabel('x [um]')
    plt.ylabel('y [um]')
    plt.title(f'Electronic Energy Loss {name}')
    plt.hist2d(x, y, weights=Ee/N, bins=num_bins, norm=colors.SymLogNorm(0.1))

    plt.savefig(name+'electronic_energy_loss.png')
    plt.close()

def plot_all_depth_distributions(name, displacement_energy, N, num_bins=100):
    '''
    Plots all available depth distributions (energy losses, displacements, implantation profiles) on one plot.

    Args:
        name (string): name of rustbca simulation
        displacement_energy (float): threshold energy needed to produce a permanent frenkel pair in energy_units
        N (int): number of incident ions
        num_bins (int): number of histogram bins to use; default 100
    '''

    energy_loss = np.genfromtxt(f'{name}energy_loss.output', delimiter=',')
    M = energy_loss[:,0]
    Z = energy_loss[:,1]
    En = energy_loss[:,2]
    Ee = energy_loss[:,3]
    x = energy_loss[:,4]
    y = energy_loss[:,5]

    plt.figure(1)
    plt.hist(x, weights=En/N, bins=num_bins, histtype='step', color='black', linewidth=2, linestyle='--', density=True)
    plt.hist(x, weights=Ee/N, bins=num_bins, histtype='step', color='red', linewidth=2, linestyle='--', density=True)
    plt.hist(x, weights=(Ee + En)/N, bins=num_bins, histtype='step', color='green', linewidth=2, density=True)

    displacements = np.genfromtxt(f'{name}displacements.output', delimiter=',')
    M = displacements[:,0]
    Z = displacements[:,1]
    Er = displacements[:,2]
    x = displacements[Er>displacement_energy,3]
    y = displacements[Er>displacement_energy,4]
    plt.hist(x, bins=num_bins, histtype='step', color='blue', linewidth=2, density=True)

    deposited = np.genfromtxt(f'{name}deposited.output', delimiter=',')
    x = deposited[:, 2]
    plt.hist(x, bins=num_bins, histtype='step', color='purple', linewidth=2, density=True)

    plt.legend(['ΔE Nuclear', 'ΔE Electronic', 'ΔE Total', f'Damage Ed = {displacement_energy}', 'Deposited'])
    plt.xlabel('x [um]')
    plt.ylabel('A.U.')
    plt.title(f'Depth Distributions {name}')
    plt.savefig(name+'plot_all_depth_distributions.png')
    plt.close()

def run_iead(ions, target, energies, angles, iead, name="default_", N=1):
    '''
    Given an IEAD in the from of a 2D array of counts, run rustbca for the ions that make up the IEAD.

    Args:
        ions (dict): ion dictionary (see top of file)
        target (dict): target dictionary (see top of file)
        energies list(float): list of energies for each column of IEAD
        angles list(float): list of angles for each row of IEAD
        name (string): name of rustbca simulation; default: 'default_'
        N (int): number of BCA ions to run per PIC particle; default 1
    '''
    import itertools

    energy_angle_pairs = list(itertools.product(energies, angles))

    E0 = np.array([pair[0] for pair in energy_angle_pairs])
    theta = np.array([pair[1] for pair in energy_angle_pairs])

    #skip last row of hPIC input because it's usually garbage
    N_ = np.array([iead[i, j]*N for i in range(len(energies)) for j in range(len(angles))], dtype=int)

    E0 = E0[N_>0]
    theta = theta[N_>0]
    N_ = N_[N_>0]

    Zb = [target['Z']]
    Mb = [target['m']]
    n = [target['n']*(MICRON)**3]
    Esb = [target['Es']] #Surface binding energy
    Ecb = [target['Ec']] #Cutoff energy
    Eb = [target['Eb']] #Bulk binding energy

    Za = ions['Z']
    Ma = ions['m']
    Esa = ions['Es']
    Eca = ions['Ec']

    thickness = 100
    depth = 100

    print(f'Generating input file {name}.toml')

    options = {
        'name': name,
        'track_trajectories': False,
        'track_recoils': True,
        'track_recoil_trajectories': False,
        'stream_size': 8000,
        'weak_collision_order': 3,
        'suppress_deep_recoils': False,
        'high_energy_free_flight_paths': False,
        'num_threads': 8,
        'num_chunks': 100,
        'use_hdf5': False,
        'electronic_stopping_mode': LOW_ENERGY_LOCAL,
        'mean_free_path_model': LIQUID,
        'interaction_potential': [["KR_C"]],
        'scattering_integral': [["MENDENHALL_WELLER"]],
        'track_displacements': True,
        'track_energy_losses': False,
    }

    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'Eb': Eb,
        'Es': Esb,
        'Ec': Ecb,
        'Z': Zb,
        'm': Mb,
        'interaction_index': np.zeros(len(n), dtype=int),
        'electronic_stopping_correction_factor': 1.0,
        'surface_binding_model': "AVERAGE"
    }

    dx = 5.*ANGSTROM/MICRON

    minx, miny, maxx, maxy = 0.0, -thickness/2., depth, thickness/2.
    surface = box(minx, miny, maxx, maxy)

    simulation_surface = surface.buffer(10.*dx, cap_style=2, join_style=2)
    mesh_2d_input = {
        'length_unit': 'MICRON',
        'energy_barrier_thickness': sum(n)**(-1./3.)/np.sqrt(2.*np.pi),
        'coordinate_sets': [[0., depth, 0., thickness/2., -thickness/2., -thickness/2.], [0., depth, depth, thickness/2., thickness/2., -thickness/2.]],
        'densities': [n, n],
        'boundary_points': [[0., thickness/2.], [depth, thickness/2.], [depth, -thickness/2.], [0., -thickness/2.]],
        'simulation_boundary_points':  list(simulation_surface.exterior.coords)
    }

    cosx = np.cos(theta*np.pi/180.)
    sinx = np.sin(theta*np.pi/180.)
    positions = [(-dx, 0., 0.) for _ in range(len(N_))]

    particle_parameters = {
        'length_unit': 'MICRON',
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'N': N_,
        'm': [Ma for _ in range(len(N_))],
        'Z': [Za for _ in range(len(N_))],
        'E': E0,
        'Ec': [Eca for _ in range(len(N_))],
        'Es': [Esa for _ in range(len(N_))],
        'interaction_index': np.zeros(len(N_), dtype=int),
        'pos': positions,
        'dir': [(cx, sx, 0.) for cx, sx in zip(cosx, sinx)],
        'particle_input_filename': ''
    }

    input_file = {
        'material_parameters': material_parameters,
        'particle_parameters': particle_parameters,
        'mesh_2d_input': mesh_2d_input,
        'options': options,
    }

    with open(f'{name}.toml', 'w') as file:
        toml.dump(input_file, file, encoder=toml.TomlNumpyEncoder())
    with open(f'{name}.toml', 'a') as file:
        file.write(r'root_finder = [[{"NEWTON"={max_iterations = 100, tolerance=1E-3}}]]')
    os.system(f'rustBCA.exe {name}.toml')

    plot_distributions_rustbca(name, ions, target, incident_energy=np.max(energies))

    s = np.atleast_2d(np.genfromtxt(f'{name}sputtered.output', delimiter=','))
    r = np.atleast_2d(np.genfromtxt(f'{name}reflected.output', delimiter=','))
    d = np.atleast_2d(np.genfromtxt(f'{name}deposited.output', delimiter=','))

    if np.size(s) > 0:
        Y = np.shape(s)[0]/np.sum(N_)
    if np.size(r) > 0:
        R = np.shape(r)[0]/np.sum(N_)

    return Y, R

def beam_target(ions, target, energy, angle, N_=10000, run_sim=True,
    high_energy=False, integral='"MENDENHALL_WELLER"',
    interaction_potential='"KR_C"',
    root_finder='{"NEWTON"={max_iterations=100, tolerance=1e-3}}', tag="",
    track_trajectories = False,
    plot_trajectories = False,
    plot_distributions = True,
    track_energy_losses = False,
    track_displacements = False,
    plot_depth_distributions = False,
    track_recoils = True,
    do_plots = True,
    thickness=100,
    depth=100, delta_x_angstrom=5.):

    '''
    Simplified generation of a monoenergetic, mono-angular beam on target simulation using rustbca.
    '''

    Zb = [target['Z']]
    Mb = [target['m']]
    n = [target['n']]
    Esb = [target['Es']] #Surface binding energy
    Ecb = [target['Ec']] #Cutoff energy
    Eb = [target['Eb']] #Bulk binding energy

    Za = ions['Z']
    Ma = ions['m']
    Esa = ions['Es']
    Eca = ions['Ec']

    N = 1

    name = ions['name']+'_'+target['name']+'_'
    if tag != '':
        name += tag+'_'

    generate_rustbca_input(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy, N, N_,
        angle, thickness, depth, name=name, track_trajectories=track_trajectories,
        track_recoil_trajectories=track_trajectories, track_energy_losses=track_energy_losses,
        track_displacements=track_displacements, track_recoils=track_recoils,
        electronic_stopping_mode="INTERPOLATED", high_energy=high_energy,
        interaction_potential=interaction_potential, integral=integral,
        root_finder=root_finder, delta_x_angstrom=delta_x_angstrom)

    if run_sim: os.system(f'rustBCA.exe {name}.toml')

    if do_plots:
        if track_energy_losses: plot_energy_loss(name, N*N_, num_bins=200)

        if plot_depth_distributions and track_displacements and track_recoils: plot_all_depth_distributions(name, 38., N_*N, num_bins=200)

        if plot_distributions: plot_distributions_rustbca(name, ions, target, incident_angle=angle, incident_energy=energy)

        if track_displacements and track_recoils: plot_displacements(name, 20., N*N_, num_bins=200)

        if track_trajectories and plot_trajectories:
            do_trajectory_plot(name, thickness=thickness, depth=depth, plot_final_positions=True)

    if track_recoils:
        s = np.atleast_2d(np.genfromtxt(f'{name}sputtered.output', delimiter=','))
    else:
        s = np.array([[]])

    r = np.atleast_2d(np.genfromtxt(f'{name}reflected.output', delimiter=','))
    d = np.atleast_2d(np.genfromtxt(f'{name}deposited.output', delimiter=','))

    return s, r, d

def sputtering(ions, target, energies, angle, N_=10000, run_sim=True):
    '''
    Produces sputtering yield curves from rustbca.
    '''

    Zb = [target['Z']]
    Mb = [target['m']]
    n = [target['n']]
    Esb = [target['Es']] #Surface binding energy
    Ecb = [target['Ec']] #Cutoff energy
    Eb = [target['Eb']] #Bulk binding energy

    Za = ions['Z']
    Ma = ions['m']
    Esa = ions['Es']
    Eca = ions['Ec']

    N = 1

    thickness = 1000
    depth = 1000

    name = ions['name']+'_'+target['name']+'_'

    track_trajectories = False
    plot_trajectories = False
    plot_distributions = False
    track_energy_losses = False
    track_displacements = False
    plot_depth_distributions = False

    integral = "MENDENHALL_WELLER"
    options = ["LOW_ENERGY_LOCAL", "LOW_ENERGY_NONLOCAL", "LOW_ENERGY_EQUIPARTITION", "INTERPOLATED"]


    for option_index, option in enumerate(options):
        s = []
        for index, energy in enumerate(energies):

            generate_rustbca_input(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy, N, N_,
                angle, thickness, depth, name=name+str(index)+'_'+str(option_index), track_trajectories=track_trajectories,
                track_recoil_trajectories=track_trajectories, track_energy_losses=track_energy_losses,
                track_displacements=track_displacements, track_recoils=True,
                electronic_stopping_mode=option, integral=integral)

            if run_sim: os.system(f'rustBCA.exe {name+str(index)+"_"+str(option_index)}.toml')
            sputtered = np.atleast_2d(np.genfromtxt(name+str(index)+'_'+str(option_index)+'sputtered.output', delimiter=','))

            if np.size(sputtered) > 0:
                s.append(np.shape(sputtered)[0]/(N_*N))
            else:
                s.append(0.0)

            if track_energy_losses: plot_energy_loss(name+str(index)+'_'+str(option_index), N*N_, num_bins=200)

            if plot_depth_distributions: plot_all_depth_distributions(name+'_'+str(index)+str(option_index), 38., N_*N, num_bins=200)

            if track_trajectories and plot_trajectories:
                do_trajectory_plot(name, thickness=thickness, depth=depth, num_trajectories=100, plot_final_positions=True)

            if plot_distributions: plot_distributions_rustbca(name+str(index)+'_'+str(option_index), hydrogen, boron, incident_angle=angle, incident_energy=energy)

            if track_displacements: plot_displacements(name+str(index)+'_'+str(option_index), 20., N*N_, num_bins=200)

        plt.loglog(energies, s)

    sy = []
    sb = []
    for energy in energies:
        sy.append(yamamura(ions, target, energy))
        sb.append(bohdansky_light_ion(ions, target, energy))

    plt.loglog(energies, sy)
    plt.loglog(energies, sb)
    plt.legend(options+['Yamamura', 'Bohdansky'])
    plt.show()

def different_depths():
    '''
    Tests reflection coefficients based on depth of Morse potential well.
    '''
    EV = 1.602E-19
    alpha = 2.366/ANGSTROM
    r0 = 2.359*ANGSTROM
    D = 0.057*EV
    root_finder = r'{"CPR"={n0=5, nmax=100, epsilon=0.1, complex_threshold=1E-12, truncation_threshold=1E-12, far_from_zero=1E20, interval_limit=1E-20, derivative_free=false}}'
    integral = r'"GAUSS_LEGENDRE"'

    energy = 400.0
    N = 10000
    depths = [D, D/2, 2*D]
    angles = np.linspace(0.01, 89.9, 20)

    for i, D in enumerate(depths):
        reflection = []
        for j, angle in enumerate(angles):

            interaction_potential = f'{{"MORSE"={{alpha={alpha}, r0={r0}, D={D}}}}}'

            s, r, d = beam_target(deuterium, tungsten, energy, angle, N_=N, run_sim=False, high_energy=False,
                integral=integral, root_finder=root_finder, interaction_potential=interaction_potential,
                tag=f'{i}_{j}')

            if np.size(r) > 0:
                reflection.append(np.shape(r)[0]/N)
            else:
                reflection.append(0.0)

        plt.plot(90 - angles, reflection)

    krc_reflection = []
    for j, angle in enumerate(angles):

        s, r, d = beam_target(deuterium, tungsten, energy, angle, N_=N*10, run_sim=True, high_energy=False,
            tag='krc_'+str(j))

        if np.size(r) > 0:
            krc_reflection.append(np.shape(r)[0]/(N*10))
        else:
            krc_reflection.append(0.)

    plt.plot(90 - angles, krc_reflection)

    plt.title('D on W 400 eV (Morse Potential)')
    plt.ylabel('R')
    plt.xlabel('Angles [deg]')
    plt.legend([f'D = {np.round(depth/EV, 4)} eV' for depth in depths]+['Kr-C'])
    plt.axis([0, 90, 0, 1])
    plt.savefig('refl_morse.png')

def plot_3d_distributions(ions, target, name):
    '''
    Uses Mayavi to plot interactive 3D distributions.

    Args:
        name (string): name of rustbca simulation
        ions (dict): ion dictionary; see top of file
        target (dict): target dictionary; see top of file
    '''

    from mayavi import mlab

    r = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    s = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    d = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))

    #sputtered distribution
    s_v = np.sqrt(s[:, 2])

    s_ux = -s[:, 6]*s_v/np.max(s_v)
    s_uy = s[:, 7]*s_v/np.max(s_v)
    s_uz = s[:, 8]*s_v/np.max(s_v)

    num_bins = 25
    s_dist, _ = np.histogramdd(np.array([s_ux, s_uy, s_uz]).transpose(), bins=num_bins, normed=False)
    s_dist /= np.max(s_dist)

    s_field = mlab.pipeline.scalar_field(s_dist)
    mlab.pipeline.iso_surface(s_field, contours=[0.3, 0.7], opacity=0.2)

    mlab.axes()
    mlab.show()

    #reflected_distribution
    r_v = np.sqrt(r[:,2])

    r_ux = -r[:, 6]*r_v/np.max(r_v)
    r_uy = r[:, 7]*r_v/np.max(r_v)
    r_uz = r[:, 8]*r_v/np.max(r_v)

    r_dist, _ = np.histogramdd(np.array([r_ux, r_uy, r_uz]).transpose(), bins=num_bins, normed=False)
    r_dist /= np.max(r_dist)

    r_field = mlab.pipeline.scalar_field(r_dist)
    mlab.pipeline.iso_surface(r_field, contours=[0.3, 0.7], opacity=0.2)

    mlab.axes()
    mlab.show()

def main():
    '''
    Here an example usage of beam_target is shown. This code runs rustbca and produces plots.

    For helium on copper at 1 keV and 0 degrees angle of incidence,
    all the distributions and trajectories are plotted.
    '''

    energy = 1000.
    angle = 0.001

    beam_target(helium, copper, energy, angle, N_ = 100, do_plots=True, run_sim=True,
        plot_trajectories = True, track_trajectories=True, thickness=0.01, depth=0.01,
        interaction_potential=r"'KR_C'")

if __name__ == '__main__':
    main()
