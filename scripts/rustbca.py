import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from shapely.geometry import Point, Polygon, box
import os
import toml
import time
from scipy.interpolate import interp1d
from matplotlib import rcParams
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import cm
from enum import Enum
from collections import namedtuple
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

def do_trajectory_plot(name, file_num='', symmetric=False, thickness=None, depth=None):
    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.output', delimiter=','))
    trajectory_data = np.atleast_1d(np.genfromtxt(name+'trajectory_data.output', delimiter=',').transpose().astype(int))

    if symmetric:
        if np.size(reflected)>0: reflected[:, 4] = abs(reflected[:, 4])
        if np.size(sputtered)>0: sputtered[:, 4] = abs(sputtered[:, 4])
        if np.size(deposited)>0: deposited[:, 3] = abs(deposited[:, 3])
        if np.size(trajectories)>0: trajectories[:, 4] = abs(trajectories[:, 4])

    colormap = cm.ScalarMappable(norm=colors.PowerNorm(gamma=0.1,vmin=1,vmax=74), cmap='viridis')

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

            #plt.scatter(x[0], y[0], color = colors[Z], marker='o', s=5)
            plt.plot(x, y, color = colormap.to_rgba(Z), linewidth = 1)

            index += trajectory_length

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

        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        plt.title(name+' Trajectories')
        plt.axis('square')
        plt.savefig(name+'trajectories_'+file_num+'.png')
        plt.close()

def generate_rustbca_input(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, E0, N, N_, theta,
    thickness, depth, track_trajectories=True, track_recoils=True,
    track_recoil_trajectories=True, name='test_',
    electronic_stopping_mode=LOW_ENERGY_NONLOCAL,
    weak_collision_order=3, ck=1., mean_free_path_model=LIQUID,
    interaction_potential=KR_C):

    options = {
        'name': name,
        'track_trajectories': track_trajectories,
        'track_recoils': track_recoils,
        'track_recoil_trajectories': track_recoil_trajectories,
        'stream_size': 8000,
        'weak_collision_order': 0,
        'suppress_deep_recoils': False,
        'high_energy_free_flight_paths': True,
        'num_threads': 8,
        'num_chunks': 100,
        'use_hdf5': False,
        'electronic_stopping_mode': "INTERPOLATED",
        'mean_free_path_model': LIQUID,
        'interaction_potential': [[interaction_potential]],
        'scattering_integral': [["MENDENHALL_WELLER"]],
    }

    dx = 20.*ANGSTROM/MICRON

    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'Eb': Eb,
        'Es': Esb,
        'Ec': Ecb,
        'n': n,
        'Z': Zb,
        'm': Mb,
        'interaction_index': np.zeros(len(n), dtype=int),
        'electronic_stopping_correction_factor': ck,
        'energy_barrier_thickness': sum(n)**(-1./3.)/np.sqrt(2.*np.pi)/MICRON
    }

    dx = 5.*ANGSTROM/MICRON

    minx, miny, maxx, maxy = 0.0, -thickness/2., depth, thickness/2.
    surface = box(minx, miny, maxx, maxy)

    simulation_surface = surface.buffer(10.*dx, cap_style=2, join_style=2)
    mesh_2d_input = {
        'length_unit': 'MICRON',
        'coordinate_sets': [[0., depth, 0., thickness/2., -thickness/2., -thickness/2.], [0., depth, depth, thickness/2., thickness/2., -thickness/2.]],
        'densities': [n, n],
        'boundary_points': [[0., thickness/2.], [depth, thickness/2.], [depth, -thickness/2.], [0., -thickness/2.]],
        'simulation_boundary_points':  list(simulation_surface.exterior.coords)
    }

    cosx = np.cos(theta*np.pi/180.)
    sinx = np.sin(theta*np.pi/180.)
    positions = [(-dx*np.sin(theta*np.pi/180.), 0., 0.) for _ in range(N)]

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

    with open('input.toml', 'w') as file:
        toml.dump(input_file, file, encoder=toml.TomlNumpyEncoder())

    with open('input.toml', 'a') as file:
        file.write(r'root_finder = [[{"NEWTON"={max_iterations = 100, tolerance=1E-3}}]]')

def plot_distributions_rustbca(name, beam, target, file_num=1, max_collision_contours=4,
    plot_2d_reflected_contours=False, collision_contour_significance_threshold=0.1, slice_angle=45,
    incident_energy=1., incident_angle=0., plot_garrison_contours=False,
    plot_reflected_energies_by_number_collisions=False,
    plot_scattering_energy_curve=False):

    file_num = str(file_num)
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
            cmap = matplotlib.cm.get_cmap('Wistia')
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
                energies = bins[0]
                angles = bins[1]*np.pi/180.
                scattering_angles = -angles + np.pi
                final_energies = incident_energy*((np.cos(scattering_angles) + np.sqrt((mass/beam['m'])**2 - np.sin(scattering_angles)**2))/(1. + mass/beam['m']))**2.
                handle = plt.plot(final_energies, angles*180./np.pi, linestyle='-', color='white', alpha=0.25, linewidth=7)
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
    plt.figure(num='d')
    plots = []
    labels = []

    if np.size(d) > 0:
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
    plt.figure(num='re')
    plots = []
    labels = []
    #bins = np.linspace(0.0, 1000.0, num_bins)
    if np.size(r) > 0:
        heights, bins, rust = plt.hist(r[:, 2]/incident_energy, histtype='step',  bins=num_bins, density=True, color='black')
        number_collision_events = r[:, -1]

        n_max = min(int(np.max(number_collision_events)), max_collision_contours)
        cmap = matplotlib.cm.get_cmap('brg')
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

    if np.size(r) > 0:
        if len(plots) > 0: plt.legend(plots, labels, fontsize='small', fancybox=True, shadow=True, loc='upper left')
        plt.title(f'Refl. Energies {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.xlabel('E/E0')
        plt.ylabel('f(E)')
        plt.gca().set_xlim(0., 1.)
        plt.gca().set_ylim(bottom=0)

        plt.savefig(name+'ref_e.png')
        plt.close()

    #Sputtered atom energy distributions
    plt.figure(num='se')
    plots = []
    labels = []
    #bins=np.linspace(0.0, 100.0, num_bins)
    if np.size(s) > 0:
        _, _, rust = plt.hist(s[:,2], histtype='step', bins=num_bins, density=True, color='black')
        plots.append(rust[0])
        labels.append('rustbca')

    if np.size(s) > 0:
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

    if np.size(s) > 0:
        plt.figure(num='so')
        depth_origin = s[:, 10]
        plt.hist(depth_origin, bins=num_bins)
        plt.title('Depth of origin of sputtered particles')
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

def main():
    Zb = [5]
    Mb = [10.811]
    n = [1.309E29]
    Esb = [5.76] #Surface binding energy
    Ecb = [1.0] #Cutoff energy
    Eb = [0.0] #Bulk binding energy
    Za = 1
    Ma = 1.008
    Esa = 1.0
    Eca = 0.1
    E0 = 1000.
    N = 1
    N_ = 1000000
    theta = 0.001
    thickness = 0.1
    depth = 0.1

    name = 'test'

    track_trajectories = False

    generate_rustbca_input(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, E0, N, N_, theta,
        thickness, depth, name=name, track_trajectories=track_trajectories,
        track_recoil_trajectories=track_trajectories)

    run_sim = True
    if run_sim: os.system('cargo run --release')

    if track_trajectories: do_trajectory_plots(name, thickness=thickness, depth=depth)

    plot_distributions_rustbca(name, hydrogen, boron, incident_angle=theta, incident_energy=E0)

if __name__ == '__main__':
    main()
