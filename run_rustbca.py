import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from shapely.geometry import Point, Polygon, box
from itertools import chain
import os
import toml
import generate_ftridyn_input as g
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
    'Ec': 0.01,
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
    'Ec': 0.1,
    'Es': 0.
}

krypton = {
    'symbol': 'Kr',
    'name': 'krypton',
    'Z': 36,
    'm': 83.80,
    'Ec': 0.1,
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
    'Ec': 0.1,
    'Es': 0.
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
    'Ec': 5.,
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
    'Ec': 0.1,
    'Es': 0.
}

def main(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, E0, N, N_, theta,
    thickness, depth, track_trajectories=False, track_recoils=False,
    track_recoil_trajectories=False, write_files=True, name='test_',
    random=False, free_flight_path=False,
    electronic_stopping_mode=LOW_ENERGY_NONLOCAL,
    weak_collision_order=3, ck=1., mean_free_path_model=LIQUID,
    interaction_potential=KR_C, scattering_integral=QUADRATURE):

    options = {
        'name': name,
        'track_trajectories': track_trajectories,
        'track_recoils': track_recoils,
        'track_recoil_trajectories': track_recoil_trajectories,
        'write_files': write_files,
        'stream_size': 8000,
        'print': True,
        'print_num': np.min((10, int(N_*N))),
        'weak_collision_order': weak_collision_order,
        'suppress_deep_recoils': False,
        'high_energy_free_flight_paths': free_flight_path,
        'electronic_stopping_mode': electronic_stopping_mode,
        'mean_free_path_model': mean_free_path_model,
        'interaction_potential': interaction_potential,
        'scattering_integral': scattering_integral,
        'tolerance': 1E-6,
        'max_iterations': 100
    }

    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'Eb': Eb,
        'Es': Esb,
        'Ec': Ecb,
        'n': n,
        'Z': Zb,
        'm': Mb,
        'electronic_stopping_correction_factor': ck
    }

    #dx = 2.*n**(-1./3.)/np.sqrt(2.*np.pi)/MICRON
    #dx = n**(-1./3.)*np.cos(theta)/MICRON
    dx = 0.

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

    cosx = np.cos(theta*np.pi/180.)
    sinx = np.sin(theta*np.pi/180.)

    if random:
        positions = [(-dx, np.random.uniform(-thickness/2., thickness/2.), 0.) for _ in range(N)]
    else:
        positions = [(-dx, 0., 0.) for _ in range(N)]

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
        'pos': positions,
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

    os.system('rustBCA.exe')

    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.output', delimiter=','))
    trajectory_data = np.atleast_1d(np.genfromtxt(name+'trajectory_data.output', delimiter=',').transpose().astype(int))

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

def do_plots(name, file_num='', symmetric=False, thickness=None, depth=None):
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

    colors = {
        1: 'red',
        29: 'black',
        2: 'red',
        74: 'blue',
        4: 'black',
        5: 'blue',
        18: 'red',
        14: 'blue',
        54: 'red',
    }

    linewidths = {
        1: 1,
        29: 1,
        2: 1,
        74: 1,
        14: 1,
        54: 1,
    }

    #minx, miny, maxx, maxy = simulation_surface.bounds
    #minx = np.min(deposited[:, 2])
    #maxx = np.max(deposited[:, 2])
    #miny = np.min(deposited[:, 3])
    #maxy = np.max(deposited[:, 3])

    fig1, axis1 = plt.subplots()
    #plt.plot(*surface.exterior.xy, color='dimgray')
    #plt.plot(*energy_surface.exterior.xy, '--', color='dimgray')
    #plt.plot(*simulation_surface.exterior.xy, '--', color='dimgray')

    index = 0
    x_max = 0

    first_atom = []

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

            plt.scatter(x[0], y[0], color = colors[Z], marker='o', s=5)
            plt.plot(x, y, color = colors[Z], linewidth = 1)

            if Z > 10:
                first_atom.append(x[0])

            index += trajectory_length

        if np.size(sputtered) > 0:
            sputtered_colors = [colors[Z] for Z in sputtered[:,1]]
            plt.scatter(sputtered[:,3], sputtered[:,4], s=50, color=sputtered_colors, marker='*')

        if np.size(reflected) > 0:
            reflected_colors = [colors[Z] for Z in reflected[:,1]]
            plt.scatter(reflected[:,3], reflected[:,4], s=50, color=reflected_colors, marker='x')

        if np.size(deposited) > 0:
            deposited_colors = [colors[Z] for Z in deposited[:,1]]
            plt.scatter(deposited[:,2], deposited[:,3], s=50, color=deposited_colors, marker='^')

        if thickness and depth:
            x_box = [0., 0., depth, depth, 0.]
            y_box = [-thickness/2., thickness/2., thickness/2., -thickness/2., -thickness/2.]
            plt.plot(x_box, y_box, color='dimgray', linewidth=3)

        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        plt.title(name+' Trajectories')
        plt.axis('square')
        #plt.axis([0., 1.1*x_max, -1.1*thickness/2., 1.1*thickness/2.])
        breakpoint()
        plt.savefig(name+'trajectories_'+file_num+'.png')
        plt.close()

    if np.size(deposited) > 0:
        plt.figure(2)
        num_bins = 100
        #bins = np.linspace(minx, maxx, num_bins)
        plt.hist(deposited[:, 2], bins=num_bins)
        plt.yscale('log')
        plt.title(name+' X Deposition')
        plt.xlabel('x [um]')
        plt.ylabel('Helium Deposition (A.U.)')
        plt.savefig(name+'x'+file_num+'.png')
        plt.close()

        plt.figure(3)
        num_bins = 100
        #bins = np.linspace(miny, maxy, num_bins)
        plt.hist(deposited[:, 3], bins=num_bins)
        plt.yscale('log')
        plt.title(name+' Y Deposition')
        plt.xlabel('y [um]')
        plt.ylabel('Helium Deposition (A.U.)')
        plt.savefig(name+'y'+file_num+'.png')
        plt.close()

        plt.figure(4)
        num_bins_x = 500
        num_bins_y = 500
        #binx = np.linspace(0.0, 1.1*np.max(deposited[:, 2]), num_bins_x)
        binx = np.linspace(-0.1*depth, 1.1*depth, num_bins_x)
        biny = np.linspace(-1.1*thickness/2., 1.1*thickness/2., num_bins_y)
        plt.hist2d(deposited[:, 2], deposited[:, 3], bins=(binx, biny))
        #plt.plot(*surface.exterior.xy, color='dimgray')
        #plt.plot(*energy_surface.exterior.xy, '--', color='dimgray')
        #plt.plot(*simulation_surface.exterior.xy, '--', color='dimgray')
        if thickness and depth:
            x_box = [0., 0., depth, depth, 0.]
            y_box = [-thickness/2., thickness/2., thickness/2., -thickness/2., -thickness/2.]
            plt.plot(x_box, y_box, color='dimgray', linewidth=3)
        plt.title(name+' 2D Deposition')
        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        #plt.axis('square')
        plt.savefig(name+'2D'+file_num+'.png')
        plt.close()

def plot_distributions(name, ftridyn_name, beam, target, file_num=1, max_collision_contours=4,
    plot_2d_reflected_contours=True, collision_contour_significance_threshold=0.1, slice_angle=45,
    incident_energy=1., incident_angle=0., plot_garrison_contours=True,
    plot_reflected_energies_by_number_collisions=True,
    plot_scattering_energy_curve=False):

    file_num = str(file_num)
    num_bins = 120

    r = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    s = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    d = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))

    rf = np.atleast_2d(np.genfromtxt(ftridyn_name+'RFLST.DAT'))
    sf = np.atleast_2d(np.genfromtxt(ftridyn_name+'SPLST.DAT'))
    df = np.atleast_2d(np.genfromtxt(ftridyn_name+'DUMPPRJ.DAT'))

    if np.size(r) > 0 and np.size(rf) > 0:
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
        normalized_energies_ftridyn = rf[:,2]/incident_energy

        c1 = ax.scatter(np.arccos(rf[:,7]), normalized_energies_ftridyn, s=1, color='red')
        c2 = ax.scatter(np.arccos(r[:,7]), normalized_energies_rust, s=1, c=number_collision_events)

        #Ar on Si 50 eV 80 deg
        if beam['symbol'] == 'Ar' and target['symbol'] == 'Si' and incident_energy==80. and incident_angle==80.:
            data_md = np.genfromtxt('comparison_data/ar_si_50ev_80deg.csv', delimiter=',')
            normalized_energies_md = data_md[:, 0]
            angles_md = np.pi/2. -  data_md[:,1]*np.pi/180.
            c3 = ax.scatter(angles_md, normalized_energies_md, s=1, color='blue')

        plt.legend(['ftridyn', 'rustbca', 'MD*'], loc='upper left')
        ax.set_thetamax(180.)
        ax.set_thetamin(0.)
        ax.set_xlabel('E/E0')
        ax.set_yticks([0., 0.5, 1.])
        ax.set_xticks([0., np.pi/6., 2*np.pi/6., 3*np.pi/6., 4*np.pi/6., 5*np.pi/6., np.pi])
        ax.set_xticklabels(['0', '30', '60', '90', '', '', ''])
        plt.title(f'{beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.savefig(name+'polar_scatter.png')
        plt.close()

        plt.figure(num='slice')
        mask = (theta > slice_angle - 5) & (theta < slice_angle + 5)
        plt.hist(r[mask, 2], bins=num_bins, histtype='step', color='black')
        plt.yscale('log')
        plt.title(f'Refl. Energies at theta={slice_angle} {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.xlabel('E [eV]')
        plt.ylabel(f'f(E, a={slice_angle})')

        plt.figure(num='rustbca_r2d')
        bin_energy = np.linspace(0., 1.2*np.max(r[:, 2]), num_bins)
        bin_angle = np.linspace(0., 90., num_bins)

        bins = (bin_energy, bin_angle)
        heights, xedges, yedges, image = plt.hist2d(r[:, 2], theta, bins=bins)

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
            energies = bins[0]
            angles = bins[1]*np.pi/180.
            scattering_angles = -angles + np.pi
            final_energies = incident_energy*((np.cos(scattering_angles) + np.sqrt((target['m']/beam['m'])**2 - np.sin(scattering_angles)**2))/(1. + target['m']/beam['m']))**2.
            handle = plt.plot(final_energies, angles*180./np.pi, linestyle='-', color='white', alpha=0.25, linewidth=7)
            plt.legend(handle, ['Single-Collision Scattering'], loc='lower right', fontsize='x-small')

        plt.xlabel('E [eV]')
        plt.ylabel('alpha [deg]')
        plt.title(f'Refl. EAD rustbca {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.savefig(name+'rustbca_r_ead.png')
        plt.close()

        plt.figure(num='ftridyn_r2d')
        ux = rf[:,6]
        uy = rf[:,7]
        uz = rf[:,8]
        theta = np.arctan2(ux, np.sqrt(uy**2 + uz**2))*180./np.pi
        theta = np.abs(np.arccos(ux)*180./np.pi)
        #theta = (np.arccos(rf[:,7]) - np.pi/2.)*180./np.pi
        plt.hist2d(rf[:, 2], theta, bins=bins)
        plt.xlabel('E [eV]')
        plt.ylabel('alpha [deg]')
        plt.title(f'Refl. EAD F-TRIDYN {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.savefig(name+'ftridyn_r_ead.png')
        plt.close()

        plt.figure(num='slice')
        mask = (theta > slice_angle - 5) & (theta < slice_angle + 5)
        plt.hist(rf[mask, 2], bins=num_bins, histtype='step', color='black', linestyle='--')
        plt.yscale('log')
        plt.savefig(name+'slice.png')
        plt.legend(['rustbca', 'F-TRIDYN'])
        plt.close()

    if np.size(s) > 0 and np.size(sf) > 0:
        plt.figure(num='rustbca_s2d')

        bin_angle = np.linspace(-90., 90., num_bins//2)
        bin_energy = np.linspace(0., 1.2*np.max(s[:, 2]), num_bins//2)
        bins = (bin_energy, bin_angle)
        ux = s[:,6]
        uy = s[:,7]
        uz = s[:,8]
        #theta = -np.arctan2(ux, np.sqrt(uy**2 + uz**2))
        theta = np.arccos(s[:,7]) - np.pi/2.

        if abs(incident_angle) < 1 or plot_garrison_contours:
            garrison = np.zeros((num_bins//2, num_bins//2))
            for i, energy in enumerate(bin_energy):
                for j, angle in enumerate(bin_angle*np.pi/180.):
                    garrison[j, i] = energy*np.cos(angle)/(energy + target['Es'])**4*(energy*np.cos(angle)**2 + target['Es'])

            plt.figure(num='garrison_2d')
            plt.pcolormesh(bin_energy, bin_angle, garrison)
            plt.xlabel('E [eV]')
            plt.ylabel('beta [deg]')
            plt.title(f'Sputtered EAD B. Garrison (1989) {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
            plt.savefig(name+'garrison_s_ead.png')
            plt.close()

        plt.hist2d(s[:, 2], theta*180./np.pi, bins=bins)
        if abs(incident_angle) < 1 or plot_garrison_contours: plt.contour(bin_energy, bin_angle, garrison, levels = 5, cmap = matplotlib.cm.get_cmap('Wistia'), alpha=0.5)
        plt.xlabel('E [eV]')
        plt.ylabel('beta [deg]')
        plt.title(f'Sputtered EAD rustbca {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.savefig(name+'rustbca_s_ead.png')
        plt.close()

        plt.figure(num='ftridyn_s2d')
        ux = sf[:,6]
        uy = sf[:,7]
        uz = sf[:,8]
        #theta = np.arctan2(ux, np.sqrt(uy**2 + uz**2))
        theta = np.arccos(sf[:,7]) - np.pi/2.
        plt.hist2d(sf[:, 2], theta*180./np.pi, bins=bins)
        if abs(incident_angle) < 1 or plot_garrison_contours: plt.contour(bin_energy, bin_angle, garrison, levels = 5, cmap = matplotlib.cm.get_cmap('Wistia'), alpha=0.5)
        plt.xlabel('E [eV]')
        plt.ylabel('beta [deg]')
        plt.title(f'Sputtered EAD F-TRIDYN {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
        plt.savefig(name+'ftridyn_s_ead.png')
        plt.close()

    #Deposited ion depth distributions
    plt.figure(num='d')
    plots = []
    labels = []

    if np.size(d) > 0:
        _, _, rust = plt.hist(d[d[:,2]>0., 2], histtype='step', bins=num_bins, density=True, color='black')
        plots.append(rust[0])
        labels.append('rustbca')

    if np.size(df) > 0:
        _, _, ftridyn = plt.hist(df[df[:,2]>0.,2]*1E-10*1E6,  linestyle='--', histtype='step', bins=num_bins, density=True, color='black')
        plots.append(ftridyn[0])
        labels.append('F-TRIDYN')

    if len(plots) > 0: plt.legend(plots, labels, fontsize='small', fancybox=True, shadow=True, loc='upper right')
    plt.title(f'Depth Distributions {beam["symbol"]} on {target["symbol"]} E0 = {np.round(incident_energy, 1)} eV {np.round(incident_angle, 1)} deg', fontsize='small')
    plt.xlabel('x [um]')
    plt.ylabel('f(x)')
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

        if beam['symbol'] == 'He':
            if target['symbol'] == 'W':

                #Comparison to MD
                #https://doi.org/10.1016/j.jnucmat.2014.01.021

                if abs(incident_angle) < 1 and incident_energy == 80.:
                #normal incidence 80 ev
                    bin_width = 10/incident_energy #eV
                    x_exp = np.array([5.007064, 15.28719, 25.256065, 35.234875, 45.21699, 55.209045, 64.9064, 75.03532])
                    y_exp = np.array([4.2616224, 0.8636774, 1.2412268, 3.4068558, 6.168511, 10.718246, 22.023607, 51.209103])
                    md = plt.scatter(x_exp/incident_energy, y_exp/np.sum(y_exp)/bin_width)
                    plots.append(md)
                    labels.append('MD*')

                if incident_angle == 60. and incident_energy == 80.:
                #60 degrees 80 ev
                    bin_width = 10/incident_energy #eV
                    x_exp = np.array([4.7257385,14.852321,35.105484,45.232067,55.35865,65.14768,74.93671])
                    y_exp = np.array([1.3513514,0.6756757,1.0135136,2.0270271,4.0540543,7.4324327,85.47298])
                    md = plt.scatter(x_exp/incident_energy, y_exp/np.sum(y_exp)/bin_width)
                    plots.append(md)
                    labels.append('MD*')

                if abs(incident_angle) < 1 and incident_energy == 20.:
                #normal incidence 20 eV
                    bin_width = 2/incident_energy #eV
                    x_exp = np.array([0.18390535,2.2103872,4.2315845,6.2556005,8.278207,10.3015175,12.232877,14.260767,16.323889,18.284842]) + 1
                    y_exp = np.array([0.,2.2690182,1.1046767,1.7723962,1.5240853,1.7337896,1.9443713,5.131175,31.218742,50.665966])
                    md = plt.scatter(x_exp/incident_energy, y_exp/np.sum(y_exp)/bin_width)
                    plots.append(md)
                    labels.append('MD*')

            if target['symbol'] == 'Ni':
            #Further comparison to MD
            #V Rosato, G Maino and A Ventura (1989)

            #He on Ni
            #0.1 eV
            #NORMALIZED ENERGIES
                if incident_energy == 0.1 and abs(incident_angle) < 1:
                    x_exp = np.array([0.847, 0.890, 0.953])
                    bin_width = 0.05
                    y_exp = np.array([11.39, 62.28, 25.74])
                    md = plt.scatter(x_exp + bin_width/2., y_exp/np.sum(y_exp)/bin_width)
                    plots.append(md)
                    labels.append('T=0K MD*')

                    x_high_temp = np.array([0.3979953,0.45450273,0.4969773,0.55361,0.60341424,0.65146434,0.70176977,0.75113547,0.7996868,0.8500548,0.91395456,0.949538,1.0056069,1.0548472,1.097823,1.1473767,1.1973062,1.2474236,1.2971653,1.3539859,1.403477,1.4457635,1.4955677])
                    bin_width = 0.05
                    y_high_temp = np.array([1.3336438,2.060239,2.0788698,2.1037107,0.7219376,9.164726,4.9759355,6.050303,11.686073,7.1464057,6.1217203,5.0846143,8.26735,10.043471,7.2550845,7.27682,5.193293,2.057134,1.0262382,0.0,0.3710604,1.4423226,0.060549606])
                    md = plt.scatter(x_high_temp + bin_width/2., y_high_temp/np.sum(y_high_temp)/bin_width)
                    plots.append(md)
                    labels.append('T=500K MD*')

                #1.0 eV
                if incident_energy == 1.0 and abs(incident_angle) < 1:
                    x_exp = np.array([0.44827586,0.50344825,0.6034483,0.6586207,0.7,0.7586207,0.80344826,0.8586207,0.9137931,0.962069,1.0034482,1.0655172])
                    bin_width = 0.05
                    y_exp = np.array([4.930362,3.091922,5.9331474,6.935933,9.94429,9.025069,17.799442,18.802229,12.785516,4.0947075,1.2534819,4.0947075])
                    md = plt.scatter(x_exp + bin_width/2., y_exp/np.sum(y_exp)/bin_width)
                    plots.append(md)
                    labels.append('MD*')

                #10.0 eV
                if incident_energy == 10. and abs(incident_angle) < 1:
                    x_exp = np.array([0.10666667,0.16333333,0.20333333,0.26333332,0.31,0.35,0.40333334,0.45,0.50333333,0.55,0.6,0.6,0.6533333,0.6933333,0.75333333,0.8066667,0.85333335,0.9])
                    bin_width = 0.05
                    y_exp = np.array([2.20339,5.338983,6.5254235,5.5084743,2.118644,2.118644,3.4745762,6.5254235,2.118644,2.118644,2.118644,7.7118645,5.338983,10.084745,11.101695,21.271187,2.118644,1.101695])
                    md = plt.scatter(x_exp + bin_width/2., y_exp/np.sum(y_exp)/bin_width)
                    plots.append(md)
                    labels.append('MD*')

                #50.0 eV
                if incident_energy == 50. and abs(incident_angle) < 1:
                    x_exp = np.array([0.006821852,0.056072842,0.10084647,0.10599209,0.15076572,0.20248929,0.25374505,0.29785043,0.34646657,0.40079635,0.44556996,0.4948878,0.54858273,0.6010414,0.64644986,0.6989085,0.74338144,0.7926324,0.848566,0.89848524,0.95221364])
                    bin_width = 0.05
                    y_exp = np.array([2.308403,2.3176434,2.3260438,0.075602725,0.084003024,6.8479266,0.10332372,2.3630052,4.5109625,2.4948897,2.5032902,2.2874024,2.4100468,6.6975613,4.5672445,8.854759,24.960659,24.969898,2.4663289,0.22428808,0.23436844])
                    md = plt.scatter(x_exp + bin_width/2., y_exp/np.sum(y_exp)/bin_width)
                    plots.append(md)
                    labels.append('MD*')

    if np.size(rf) > 0:
        _, _, ftridyn = plt.hist(rf[:,2]/incident_energy,  linestyle='--', histtype='step', bins=num_bins, density=True, color='black')
        plots.append(ftridyn[0])
        labels.append('F-TRIDYN')

    if np.size(rf) > 0 or np.size(r) > 0:
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

    if np.size(sf) > 0:
        _, _, ftridyn = plt.hist(sf[:,2], histtype='step',  linestyle='--', bins=num_bins, density=True, color='black')
        plots.append(ftridyn[0])
        labels.append('F-TRIDYN')

    if np.size(s) > 0 and np.size(sf) > 0:
        energies = np.linspace(0., np.max(s[:,2]), num_bins)
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

    #Sputtered atom angular distributions
    plt.figure(num='sa')
    plots = []
    labels = []
    ax = None
    bins = np.linspace(0, np.pi, num_bins)

    if np.size(s) > 0 and np.size(sf) > 0:
        #hist, bins = np.histogram(np.arccos(-s[:,6]), bins=bins, density=True)
        #rust1 = plt.plot(bins[:-1], hist/np.max(hist), linestyle='--')
        #plots.append(rust1[0])
        #labels.append('rustbca cosx')

        #hist, bins = np.histogram(np.arccos(sf[:,6]), bins=bins, density=True)
        #ftridyn1 = plt.plot(bins[:-1], hist/np.max(hist), color=rust1[0].get_color())
        #plots.append(ftridyn1[0])
        #labels.append('ftridyn cosx')

        hist, bins = np.histogram(np.arccos(s[:,7]), bins=bins, density=True)
        rust2 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust2[0])
        labels.append('rustbca cosy')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(sf[:,7]), bins=bins, density=True)
        ftridyn2 = plt.polar(bins[:-1], hist/np.max(hist), color=rust2[0].get_color(), linestyle='--')
        plots.append(ftridyn2[0])
        labels.append('ftridyn cosy')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(s[:,8]), bins=bins, density=True)
        rust3 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust3[0])
        labels.append('rustbca cosz')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(sf[:,8]), bins=bins, density=True)
        ftridyn3 = plt.polar(bins[:-1], hist/np.max(hist), color=rust3[0].get_color(), linestyle='--')
        plots.append(ftridyn3[0])
        labels.append('ftridyn cosz')
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

    if np.size(r) > 0 and np.size(rf) > 0:
        hist, bins = np.histogram(np.arccos(-r[:,6]), bins=bins, density=True)
        rust1 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust1[0])
        labels.append('rustbca cosx')

        #80 eV normal incidence He on W
        if beam['symbol'] == 'He' and target['symbol'] == 'W' and incident_energy == 80. and abs(incident_angle) < 1:
            x_exp = np.array([1.3445379,11.092437,20.840336,30.92437,41.008404,51.092438,60.840336,70.92437]) + 5
            y_exp = np.array([6.835206,14.419476,18.63296,19.662922,17.602997,17.041199,5.243446,0.37453184])
            y_exp/=np.max(y_exp)
            exp1 = plt.polar(x_exp/180.*np.pi, y_exp, linestyle=':', color=rust1[0].get_color())
            plots.append(exp1[0])
            labels.append('MD')

        hist, bins = np.histogram(np.arccos(rf[:,6]), bins=bins, density=True,)
        ftridyn1 = plt.polar(bins[:-1], hist/np.max(hist), color=rust1[0].get_color(), linestyle='--')
        plots.append(ftridyn1[0])
        labels.append('ftridyn cosx')

        hist, bins = np.histogram(np.arccos(r[:,7]), bins=bins, density=True)
        rust2 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust2[0])
        labels.append('rustbca cosy')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(rf[:,7]), bins=bins, density=True)
        ftridyn2 = plt.polar(bins[:-1], hist/np.max(hist), color=rust2[0].get_color(), linestyle='--')
        plots.append(ftridyn2[0])
        labels.append('ftridyn cosy')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(r[:,8]), bins=bins, density=True)
        rust3 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust3[0])
        labels.append('rustbca cosz')
        ax = plt.gca()

        hist, bins = np.histogram(np.arccos(rf[:,8]), bins=bins, density=True)
        ftridyn3 = plt.polar(bins[:-1], hist/np.max(hist), color=rust3[0].get_color(), linestyle='--')
        plots.append(ftridyn3[0])
        labels.append('ftridyn cosz')
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

def starshot():
    beam_species = [helium]
    beam_linestyles = {
        'He': '--',
        'H': '-',
    }
    target_species = [beryllium]
    thickness = 100
    depth = 100000
    N = 1000
    N_ = 1
    theta = 0.0001
    #velocities =
    velocities = np.arange(0.025, 0.51, 0.01)
    #velocities = [0.05]

    os.system('rm rustBCA.exe')
    os.system('cargo build --release')
    os.system('mv target/release/rustBCA.exe .')

    os.system('rm *.output')

    run_sim = True
    track_recoils = True
    track_trajectories = False
    track_recoil_trajectories = True
    ffp = True
    esmode = INTERPOLATED

    plt.figure('reflection')
    plt.figure('sputtering')

    rustbca_reflection_names = []
    for beam in beam_species:
        Za = beam['Z']
        Ma = beam['m']
        Eca = beam['Ec']
        Esa = beam['Es']

        for target in target_species:
            Zb = target['Z']
            Ecb = target['Ec']
            Esb = target['Es']
            Eb = target['Eb']
            Mb = target['m']
            n = target['n']

            rustbca_reflection = []
            rustbca_yield = []

            for velocity in velocities:
                energy = Ma*AMU*C**2.*(1./np.sqrt(1. - velocity**2) - 1.)/Q
                print(f'E: {energy}')

                name = str(velocity)+'_'+str(beam['symbol'])+'_'+str(target['symbol'])

                if run_sim:
                    rustbca_yield_, rustbca_reflection_, rustbca_range_, rustbca_straggle_  = main(
                        Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy,
                        N, N_, theta, thickness, depth,
                        track_trajectories=track_trajectories, track_recoils=track_recoils,
                        track_recoil_trajectories=track_recoil_trajectories, write_files=True,
                        name=name, random=True, free_flight_path=ffp, electronic_stopping_mode=esmode)

                    rustbca_reflection.append(rustbca_reflection_)
                    rustbca_yield.append(rustbca_yield_)

                    print(f'{beam["symbol"]} {target["symbol"]} velocity: {velocity} Y: {rustbca_yield_} R: {rustbca_reflection_} depth: {rustbca_range_}')

                    do_plots(name, file_num=str(1), symmetric=False, thickness=thickness, depth=depth)
                    #breakpoint()

            rustbca_reflection_names.append(beam['symbol']+' on '+target['symbol'])

            plt.figure('reflection')
            plt.plot(velocities, rustbca_reflection, beam_linestyles[beam['symbol']])
            plt.figure('sputtering')
            plt.plot(velocities, rustbca_yield)

    plt.figure('reflection')
    plt.legend(rustbca_reflection_names)
    plt.title('Reflection Coefficients from 100cmX10um Target')
    plt.xlabel('v/c')
    plt.ylabel('R')

    plt.figure('sputtering')
    plt.legend(rustbca_reflection_names)
    plt.title('Sputtering Yields from 1cmX10um Target')
    plt.xlabel('v/c')
    plt.ylabel('Y [at/ion]')

    plt.show()

def benchmark_srim():
    beam = helium
    target = aluminum

    energies = np.logspace(1, 4, 20)*1000.
    energies = [10E6]
    theta = 0.0001
    N = 1
    N_ = 10000
    depth = 200
    thickness = 200

    Za = beam['Z']
    Ma = beam['m']
    Eca = beam['Ec']
    Esa = beam['Es']
    Zb = target['Z']
    Ecb = target['Ec']
    Esb = target['Es']
    Eb = target['Eb']
    Mb = target['m']
    n = target['n']

    track_recoils = False
    track_trajectories = False
    track_recoil_trajectories = False
    ffp = True

    os.system('rm *.output')
    os.system('rm rustBCA.exe')
    os.system('cargo build --release')
    os.system('mv target/release/rustBCA.exe .')

    ranges = []
    straggles = []
    for index, energy in enumerate(energies):
        name = beam['symbol']+'_'+target['symbol']+str(index)+'_'
        rustbca_yield_, rustbca_reflection_, rustbca_range_, rustbca_straggle_  = main(
            Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy,
            N, N_, theta, thickness, depth,
            track_trajectories=track_trajectories, track_recoils=track_recoils,
            track_recoil_trajectories=track_recoil_trajectories, write_files=True,
            name=name, random=False, free_flight_path=ffp)
        do_plots(name, thickness=thickness, depth=depth)
        ranges.append(rustbca_range_)
        straggles.append(rustbca_straggle_)
        breakpoint()

    ranges = np.array(ranges)
    plt.loglog(energies, ranges)
    plt.loglog(energies, straggles)
    geant4 = np.genfromtxt('comparison_data/geant4_as_si.dat')
    plt.loglog(geant4[:,0], geant4[:,1]*1e-3, '.')
    plt.show()

def test_rustbca():
    beam_species = [hydrogen, helium]
    target_species = [beryllium, tungsten]

    N = 1 #Number of unique ion kinds
    N_ = 10000 #Number of ions
    theta = 0.0001 #Incident angle w.r.t. surface normal
    thickness = 1 #micron
    depth = 1 #micron
    energies = np.logspace(3, 4, 10)
    tr = True
    trt = False
    tt = False
    ffp = False
    esmode = LOW_ENERGY_EQUIPARTITION

    #os.system('rm *.png')
    os.system('rm *.output')
    os.system('rm rustBCA.exe')
    os.system('cargo build --release')
    os.system('mv target/release/rustBCA.exe .')

    fig1, ax1 = plt.subplots(1, 1, num='yields', figsize=(16, 14))
    fig2, ax2 = plt.subplots(1, 1, num='ranges', figsize=(16, 14))

    name_1 = []
    name_2 = []
    rustbca_time = 0.
    ftridyn_time = 0.
    for beam in beam_species:
        Za = beam['Z']
        Ma = beam['m']
        Eca = beam['Ec']
        Esa = beam['Es']
        for target in target_species:

            #Append names to namelists for legend
            name_1.append(beam['symbol']+' on '+target['symbol']+' Y rustbca')
            name_1.append(beam['symbol']+' on '+target['symbol']+' R rustbca')
            name_1.append(beam['symbol']+' on '+target['symbol']+' Y Yamamura')
            nba
            name_1.append(beam['symbol']+' on '+target['symbol']+' R Wierzbicki-Biersack')
            name_1.append(beam['symbol']+' on '+target['symbol']+' R Thomas et al.')

            name_2.append(beam['symbol']+' on '+target['symbol']+' Range rustbca')

            #Generate empty arrays for each ion-target combo
            rustbca_yield = []
            rustbca_reflection = []
            yamamura_yield = []
            rustbca_range = []
            reflection_coefficient = []
            reflection_coefficient_thomas = []

            Zb = target['Z']
            Ecb = target['Ec']
            Esb = target['Es']
            Eb = target['Eb']
            Mb = target['m']
            n = target['n']

            for energy in energies:
                rustbca_name = str(energy)+'eV_'+str(theta)+'deg_'+beam['symbol']+'_'+target['symbol']

                start = time.time()
                rustbca_yield_, rustbca_reflection_, rustbca_range_, rustbca_straggle_ = main(Zb,
                    Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy, N, N_,
                    theta, thickness, depth, track_trajectories=tt,
                    track_recoils=tr, track_recoil_trajectories=trt,
                    write_files=True, name=rustbca_name, random=False,
                    free_flight_path=ffp, electronic_stopping_mode=esmode)
                stop = time.time()
                rustbca_time += stop - start

                rustbca_yield.append(rustbca_yield_)
                rustbca_reflection.append(rustbca_reflection_)
                rustbca_range.append(rustbca_range_)

                yamamura_yield_ = yamamura(beam, target, energy)
                yamamura_yield.append(yamamura_yield_)

                reflection_coefficient_ = wierzbicki_biersack(beam, target, energy)
                reflection_coefficient.append(reflection_coefficient_)

                reflection_coefficient_thomas_ = thomas_reflection(beam, target, energy)
                reflection_coefficient_thomas.append(reflection_coefficient_thomas_)

                do_plots(rustbca_name, file_num=str(1), symmetric=False, thickness=thickness, depth=depth)
                #breakpoint()

            handle = ax1.loglog(energies, rustbca_yield, '-*')
            ax1.loglog(energies, rustbca_reflection, ':', color=handle[0].get_color())
            ax1.loglog(energies, yamamura_yield, '-', color=handle[0].get_color())
            ax1.loglog(energies, reflection_coefficient, '-.', color=handle[0].get_color())
            ax1.loglog(energies, reflection_coefficient_thomas, '-^', color=handle[0].get_color())
            ax2.loglog(energies, rustbca_range, color=handle[0].get_color())

    ax1.legend(name_1, fontsize='small', loc='upper left', bbox_to_anchor=(0.85, 1.0))
    ax1.axis([0.1*np.min(energies), 10.0*np.max(energies), 1./N_, 10.])
    ax1.set_ylabel('Y/R [at/ion]/[ion/ion]')
    ax1.set_xlabel('E [eV]')
    fig1.savefig(rustbca_name+'yields.png')

    ax2.legend(name_2, fontsize='small', loc='upper left', bbox_to_anchor=(1.0, 0.5))
    ax2.set_xlabel('E [eV]')
    ax2.set_ylabel('Range [um]')
    ax2.axis([0.5*np.min(energies), 10.*np.max(energies), 0.1*np.min(rustbca_range), 10.*np.max(rustbca_range)])
    fig2.savefig(rustbca_name+'ranges.png')

    plt.show()

def benchmark():
    beam_species = [helium]#, argon, boron]#, hydrogn]#, helium]#, beryllium, boron, neon, silicon, argon, copper, tungsten]
    target_species = [tungsten]#, tungsten]#, tungsten]#, nickel]#, silicon, copper]#, copper]#, copper]#, boron, silicon, copper]

    N = 1
    N_ = 10000
    angles = np.array([0.001])
    thickness = 1
    depth = 1
    #energies = [1000.]
    #energies = [20.0]
    #energies = [0.1, 1., 10., 50.]
    energies = np.round(np.logspace(1., 4, 20), 1)
    #energies = [100.]
    #energies = np.array([1E7])
    #energies = [1000.]
    #angles = np.round(np.linspace(0.001, 80, 5), 1)
    #energies = [8000.]
    tr = True
    trt = False
    tt = False
    ffp = False
    esmode = LOW_ENERGY_NONLOCAL
    ck = 1.
    weak_collision_order = 3
    mfp_model = LIQUID
    do_trajectory_plot = False
    run_magic = False
    interaction_potential = TRIDYN
    scattering_integral = MAGIC

    #os.system('rm *.png')
    os.system('rm *.output')
    os.system('rm *.DAT')
    os.system('rm *.IN')
    os.system('rm *.OUT')

    os.system('rm rustBCA.exe')
    os.system('cargo build --release')
    os.system('mv target/release/rustBCA.exe .')

    for theta in angles:

        fig1, ax1 = plt.subplots(1, 1, num='yields_', figsize=(8, 6))
        fig2, ax2 = plt.subplots(1, 1, num='refl_', figsize=(8, 6))
        fig3, ax3 = plt.subplots(1, 1, num='range_', figsize=(8, 6))
        fig4, ax4 = plt.subplots(1, 1, num='perf_', figsize=(8, 6))

        name_1 = []
        name_2 = []
        name_3 = []
        rustbca_time = 0.
        ftridyn_time = 0.
        for beam in beam_species:
            Za = beam['Z']
            Ma = beam['m']
            Eca = beam['Ec']
            Esa = beam['Es']
            for target in target_species:
                Zb = target['Z']
                Ecb = target['Ec']
                Esb = target['Es']
                Eb = target['Eb']
                Mb = target['m']
                n = target['n']

                name_1.append(beam['symbol']+' on '+target['symbol']+' rustbca')
                if run_magic: name_1.append(beam['symbol']+' on '+target['symbol']+' rustbca MAGIC')

                name_2.append(beam['symbol']+' on '+target['symbol']+' rustbca')
                if run_magic: name_2.append(beam['symbol']+' on '+target['symbol']+' rustbca MAGIC')

                name_1.append(beam['symbol']+' on '+target['symbol']+' F-TRIDYN')
                name_2.append(beam['symbol']+' on '+target['symbol']+' F-TRIDYN')

                name_1.append(beam['symbol']+' on '+target['symbol']+' Yamamura')
                if Ma/Mb < 0.5: name_1.append(beam['symbol']+' on '+target['symbol']+' Bohdansky')

                #name_2.append(beam['symbol']+' on '+target['symbol']+' Wierzbicki-Biersack')
                name_2.append(beam['symbol']+' on '+target['symbol']+' Thomas et al.')

                name_3.append(beam['symbol']+' on '+target['symbol']+' Range rustbca')
                if run_magic: name_3.append(beam['symbol']+' on '+target['symbol']+' Range rustbca MAGIC')

                name_3.append(beam['symbol']+' on '+target['symbol']+' Straggle rustbca')
                if run_magic: name_3.append(beam['symbol']+' on '+target['symbol']+' Straggle rustbca MAGIC')

                name_3.append(beam['symbol']+' on '+target['symbol']+' Range F-TRIDYN')
                name_3.append(beam['symbol']+' on '+target['symbol']+' Straggle F-TRIDYN')

                rustbca_yield = []
                rustbca_reflection = []
                rustbca_yield_MAGIC = []
                rustbca_reflection_MAGIC = []

                ftridyn_yield = []
                ftridyn_reflection = []

                yamamura_yield = []
                bohdansky_yield = []

                rustbca_range = []
                rustbca_range_MAGIC = []
                ftridyn_range = []

                rustbca_straggle = []
                rustbca_straggle_MAGIC = []
                ftridyn_straggle = []

                reflection_coefficient = []
                reflection_coefficient_thomas = []

                computational_time_rustbca = []
                computational_time_rustbca_ffp = []
                computational_time_ftridyn = []

                for energy in energies:
                    rustbca_name = str(energy)+'eV_'+str(theta)+'deg_'+beam['symbol']+'_'+target['symbol']

                    start = time.time()
                    rustbca_yield_, rustbca_reflection_, rustbca_range_, rustbca_straggle_ = main(Zb,
                        Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy, N, N_,
                        theta, thickness, depth, track_trajectories=tt,
                        track_recoils=tr, track_recoil_trajectories=trt,
                        write_files=True, name=rustbca_name, random=False,
                        free_flight_path=ffp, electronic_stopping_mode=esmode,
                        weak_collision_order=weak_collision_order, ck=ck,
                        mean_free_path_model=mfp_model, scattering_integral=scattering_integral,
                        interaction_potential=interaction_potential)
                    stop = time.time()

                    rustbca_time += stop - start
                    computational_time_rustbca.append(stop - start)

                    if run_magic:
                        start = time.time()
                        os.system('rm *.output')
                        rustbca_yield_MAGIC_, rustbca_reflection_MAGIC_, rustbca_range_MAGIC_, rustbca_straggle_MAGIC_ = main(Zb,
                            Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy, N, N_,
                            theta, thickness, depth, track_trajectories=tt,
                            track_recoils=tr, track_recoil_trajectories=trt,
                            write_files=True, name=rustbca_name, random=False,
                            free_flight_path=True, electronic_stopping_mode=INTERPOLATED,
                            weak_collision_order=0, ck=ck,
                            mean_free_path_model=mfp_model, scattering_integral=MAGIC,
                            interaction_potential=TRIDYN)
                        stop = time.time()
                        computational_time_rustbca_ffp.append(stop - start)

                        rustbca_yield_MAGIC.append(rustbca_yield_MAGIC_)
                        rustbca_reflection_MAGIC.append(rustbca_reflection_MAGIC_)
                        rustbca_range_MAGIC.append(rustbca_range_MAGIC_)
                        rustbca_straggle_MAGIC.append(rustbca_straggle_MAGIC_)

                    rustbca_yield.append(rustbca_yield_)
                    rustbca_reflection.append(rustbca_reflection_)
                    rustbca_range.append(rustbca_range_)
                    rustbca_straggle.append(rustbca_straggle_)

                    ftridyn_name = beam['symbol'].ljust(2, '_') + target['symbol'].ljust(2, '_')
                    interface = g.tridyn_interface(beam['symbol'], target['symbol'])
                    #os.system('rm *.DAT')
                    #os.system('rm *.OUT')
                    start = time.time()
                    ftridyn_yield_, ftridyn_reflection_, ftridyn_range_ = interface.run_tridyn_simulations_from_iead([energy], [theta], np.array([[1.]]), number_histories=np.max([100, N_]), depth=depth*MICRON/ANGSTROM)
                    stop = time.time()
                    ftridyn_time += stop - start
                    computational_time_ftridyn.append(stop - start)
                    ftridyn_yield.append(ftridyn_yield_/np.max([100, N_]))
                    ftridyn_reflection.append(ftridyn_reflection_/np.max([100, N_]))

                    df = np.atleast_2d(np.genfromtxt(ftridyn_name+'DUMPPRJ.DAT'))

                    if np.size(df) > 0:
                        ftridyn_range.append(np.mean(df[:, 2])*1E-10)
                        ftridyn_straggle.append(np.std(df[:, 2])*1E-10)
                    else:
                        ftridyn_range.append(0.)
                        ftridyn_straggle.append(0.)

                    yamamura_yield_ = yamamura(beam, target, energy)
                    if Ma/Mb < 0.5:
                        bohdansky_yield_ = bohdansky_light_ion(beam, target, energy)
                        bohdansky_yield.append(bohdansky_yield_)
                    yamamura_yield.append(yamamura_yield_)

                    #reflection_coefficient_ = wierzbicki_biersack(beam, target, energy)
                    #reflection_coefficient.append(reflection_coefficient_)

                    reflection_coefficient_thomas_ = thomas_reflection(beam, target, energy)
                    reflection_coefficient_thomas.append(reflection_coefficient_thomas_)

                    plot_distributions(rustbca_name, ftridyn_name, beam, target,
                        file_num=1, incident_energy=energy, incident_angle=theta,
                        plot_2d_reflected_contours=False,
                        plot_reflected_energies_by_number_collisions=True,
                        plot_scattering_energy_curve=True)
                    if do_trajectory_plot: do_plots(rustbca_name, file_num=str(1), symmetric=False, thickness=thickness, depth=depth)
                    #breakpoint()
                    #os.system('rm *.output')

                handle = ax1.semilogx(energies, rustbca_yield, '-*')
                if run_magic: ax1.semilogx(energies, rustbca_yield_MAGIC, '-+')
                ax1.semilogx(energies, ftridyn_yield, '--o', color=handle[0].get_color())
                ax1.semilogx(energies, yamamura_yield, '-', color=handle[0].get_color())
                if Ma/Mb < 0.5: ax1.semilogx(energies, bohdansky_yield, '-s', color=handle[0].get_color())

                #if beam['symbol'] == 'D' and target['symbol'] == 'Be' and abs(theta) < 1:
                #    data = np.genfromtxt('sdtrimsp_d_be')
                #    data_magic = data[:20, :]
                #    data_gauss_mehler = data[20:, :]
                #    ax1.semilogx(data_magic[:, 0], data_magic[:, 1], '-X')
                #    ax1.semilogx(data_gauss_mehler[:, 0], data_gauss_mehler[:, 1], '-P')
                #    name_1.append(beam['symbol']+' on '+target['symbol']+' SD.TRIM.SP MAGIC')
                #    name_1.append(beam['symbol']+' on '+target['symbol']+' SD.TRIM.SP Gauss-Mehler')

                ax2.loglog(energies, rustbca_reflection, '-*')#, color=handle[0].get_color())
                if run_magic: ax2.loglog(energies, rustbca_reflection_MAGIC, '-+')
                ax2.loglog(energies, ftridyn_reflection, '--^')#, color=handle[0].get_color())
                #ax2.loglog(energies, reflection_coefficient, '-.')#, color=handle[0].get_color())
                ax2.plot(energies, reflection_coefficient_thomas, '-')#, color=handle[0].get_color())

                ax3.loglog(energies, rustbca_range, '-', color=handle[0].get_color())
                if run_magic: ax3.loglog(energies, rustbca_range_MAGIC, '-+', color=handle[0].get_color())

                ax3.loglog(energies, rustbca_straggle, '-*', color=handle[0].get_color())
                if run_magic: ax3.loglog(energies, rustbca_straggle_MAGIC, '-+', color=handle[0].get_color())

                ax3.loglog(energies, np.array(ftridyn_range)/MICRON, '--o', color=handle[0].get_color())
                ax3.loglog(energies, np.array(ftridyn_straggle)/MICRON, '--^', color=handle[0].get_color())

        ax1.set_ylabel('Y [at/ion]')
        ax1.set_xlabel('E [eV]')
        ax1.set_xlim(1, 2.*np.max(energies))
        #ax1.axis([0.5*np.min(energies), 2.*np.max(energies), 1./N_, 2.])

        ax2.set_ylabel('R [at/ion]')
        ax2.set_xlabel('E [eV]')
        ax2.set_xlim(0.5*np.min(energies), 2.*np.max(energies))
        #ax2.axis([0.5*np.min(energies), 2.*np.max(energies), 1./N_, 2.])

        if beam['symbol'] == 'He' and target['symbol'] == 'W' and abs(theta) < 1:
            # V. Borovikov et al. (2014)
            energies_md = np.array([4.5111375,9.609054,19.816422,29.652786,39.885143,59.930798,79.96684,99.99904])
            r_md = np.array([0.9973039,0.9975538,0.9745242,0.9083381,0.8343281,0.7411909,0.6676615,0.6019752])
            ax2.loglog(energies_md, r_md, '--1')#, color=handle[0].get_color())
            name_2.append('He on W MD*')

            energies_exp = np.array([7.2638583,14.7358055,24.606771,49.786865,99.69532,119.71021])
            r_van_gorkum = np.array([0.9817522,0.9389802,0.80220586,0.6348087,0.42156383,0.39117166])
            ax2.loglog(energies_exp, r_van_gorkum, '-x')#, color=handle[0].get_color())
            name_2.append('He on W Experiment*')

            energies_marlowe = np.array([7.2773147,9.643655,11.221856,13.570896,15.547012,17.523129,19.495398,29.704687,39.123917,49.721508,59.14458])
            r_marlowe = np.array([0.95430124,0.9269656,0.90743464,0.9153931,0.8841167,0.85284024,0.829407,0.8024558,0.78723085,0.76814204,0.74507403])
            ax2.loglog(energies_marlowe, r_marlowe, '--s')#, color=handle[0].get_color())
            name_2.append('He on W MARLOWE*')

        if target['symbol'] == 'W' and abs(theta) < 1:
            close_yarwood = np.genfromtxt('comparison_data/close_yarwood')
            close_yarwood_2 = np.genfromtxt('comparison_data/close_yarwood_2')
            if beam['symbol'] == 'He':
                data = close_yarwood_2[:15, :]
                #ax2.loglog(data[:, 0], 1. - data[:, 1], '>')
                #name_2.append('He on W Experiment+')

            if beam['symbol'] == 'Ne':
                data = close_yarwood[:13, :]
                handle = ax2.loglog(data[:, 0], 1. - data[:, 1], '<')
                name_2.append('Ne on W Experiment+')

                data = close_yarwood_2[15:, :]
                ax2.loglog(data[:, 0], 1. - data[:, 1], '>', color=handle[0].get_color())
                name_2.append('Ne on W Experiment+')

            if beam['symbol'] == 'Ar':
                data = close_yarwood[13:25, :]
                ax2.loglog(data[:, 0], 1. - data[:, 1], '>')
                name_2.append('Ar on W Experiment+')

            if beam['symbol'] == 'Kr':
                data = close_yarwood[25:, :]
                ax2.loglog(data[:, 0], 1. - data[:, 1], '>')
                name_2.append('K on W Experiment+')

        ax1.set_title('Sputtering Yields')
        ax1.legend(name_1, fontsize='small', loc='lower right')#, bbox_to_anchor=(0.85, 1.0))
        fig1.savefig(str(theta)+'_deg_yields.png')

        ax2.set_title('Reflection Coefficients')
        ax2.legend(name_2, fontsize='small', loc='lower left')#, bbox_to_anchor=(0.85, 1.0))
        ax2.set_ylim(2.E-1, 1.1)
        fig2.savefig(str(theta)+'_deg_refl.png')

        ax3.set_title('Ranges and Straggles')
        ax3.legend(name_3, fontsize='small', loc='upper left')#, bbox_to_anchor=(1.0, 0.5))
        ax3.set_xlabel('E [eV]')
        ax3.set_ylabel('Range [um]')
        #ax3.axis([0.5*np.min(energies), 10.*np.max(energies), 0.1*np.min(rustbca_range), 10.*np.max(rustbca_range)])
        fig3.savefig(str(theta)+'_deg_ranges.png')

        ax4.loglog(energies, np.array(computational_time_rustbca)/N_/N/energies)
        ax4.loglog(energies, np.array(computational_time_ftridyn)/N_/N/energies)
        if run_magic: ax4.loglog(energies, computational_time_rustbca_ffp)

        ax4.set_title('Computational Time of BCA Codes for '+beam['symbol']+' on '+target['symbol'])
        ax4.set_xlabel('Energy [eV]')
        ax4.set_ylabel('Computational Time per Ion / E0 [s/eV]')
        ax4.legend(['rustbca', 'F-TRIDYN', 'rustbca FFP'])
        fig4.savefig('performance.png')

    print(f'time rustBCA: {rustbca_time} time F-TRIDYN: {ftridyn_time}')

if __name__ == '__main__':
    benchmark()
