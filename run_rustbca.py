import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, box
from itertools import chain
import os
import toml
import generate_ftridyn_input as g
import time

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

titanium = {
    'symbol': 'Ti',
    'name': 'titanium',
    'Z': 22,
    'm': 47.867,
    'Es': 4.84,
    'Ec': 3.5,
    'Eb': 0.,
    'Q': 1.5,
    'n': 5.67e28
}

hydrogen = {
    'symbol': 'H',
    'name': 'hydrogen',
    'Z': 1,
    'm': 1.008,
    'Ec': 0.95,
    'Es': 1.5,
}
helium = {
    'symbol': 'He',
    'name': 'helium',
    'Z': 2,
    'm': 4.002602,
    'Ec': 0.1,
    'Es': 0.
}

beryllium = {
    'symbol': 'Be',
    'name': 'beryllium',
    'Z': 4,
    'm': 9.012182,
    'n': 1.235E29,
    'Es': 3.31,
    'Eb': 0.,
    'Ec': 3.,
    'Q': 2.17,
}

boron = {
    'symbol': 'B',
    'name': 'boron',
    'Z': 5,
    'm': 10.811,
    'n': 1.37E29,
    'Es': 5.76,
    'Eb': 0.,
    'Ec': 3.,
    'Q': 4.6,
}

neon = {
    'symbol': 'Ne',
    'name': 'neon',
    'Z': 10,
    'm': 20.1797,
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
    'Q': 0.78,
}

argon = {
    'symbol': 'Ar',
    'name': 'argon',
    'Z': 18,
    'm': 39.948,
    'Ec': 0.1,
    'Es': 0.
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
    'Q': 1.30,
}

tungsten = {
    'symbol': 'W',
    'name': 'tungsten',
    'Z': 74,
    'm': 183.84,
    'n': 6.306E28,
    'Es': 11.75,
    'Eb': 0.,
    'Ec': 3.,
    'Q': 1.5
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

def main(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, E0, N, N_, theta, thickness, depth, track_trajectories=False, track_recoils=False, track_recoil_trajectories=False, write_files=True, name='test_'):

    options = {
        'name': name,
        'track_trajectories': track_trajectories,
        'track_recoils': track_recoils,
        'track_recoil_trajectories': track_recoil_trajectories,
        'write_files': write_files,
        'stream_size': 256000,
        'print': True,
        'print_num': 10
    }

    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'Eb': Eb,
        'Es': Esb,
        'Ec': Ecb,
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

    cosx = np.cos(theta*np.pi/180.)
    sinx = np.sin(theta*np.pi/180.)

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
        'pos': [(-dx, 0., 0.) for _ in range(N)],
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

    os.system('./rustBCA')

    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.output', delimiter=','))
    trajectory_data = np.genfromtxt(name+'trajectory_data.output', delimiter=',').transpose().astype(int)

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

def wierzbicki_biersack(ion, target, energy_eV):
    #Wierzbicki and Biersack (1994)
    Z1 = ion['Z']
    Z2 = target['Z']
    M1 = ion['m']
    M2 = target['m']
    energy_keV = energy_eV/1E3

    #I've never seen this form of the reduced energy before
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

def yamamura(ion, target, energy_eV):
    z1 = ion['Z']
    z2 = target['Z']
    m1 = ion['m']
    m2 = target['m']
    Us = target['Es']
    Q = target['Q']

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

def do_plots(name, file_num=''):
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
    if np.size(trajectories) > 0:
        for trajectory_length in trajectory_data:

            M = trajectories[index, 0]
            Z = trajectories[index, 1]
            E = trajectories[index:(trajectory_length + index), 2]
            x = trajectories[index:(trajectory_length + index), 3]
            y = trajectories[index:(trajectory_length + index), 4]
            z = trajectories[index:(trajectory_length + index), 5]

            plt.scatter(x, y, color = colors[Z], marker='o', s=5)
            plt.plot(x, y, color = colors[Z], linewidth = 1)

            index += trajectory_length
        if np.size(sputtered) > 0:
            plt.scatter(sputtered[:,3], sputtered[:,4], s=50, color='blue', marker='*')
        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        plt.title(name+' Trajectories')
        plt.savefig(name+'trajectories_'+file_num+'.png')
        plt.close()

    if np.size(deposited) > 0:
        plt.figure(2)
        num_bins = 100
        #bins = np.linspace(minx, maxx, num_bins)
        plt.hist(deposited[:, 2], bins=num_bins)
        plt.title(name+' X Deposition')
        plt.xlabel('x [um]')
        plt.ylabel('Helium Deposition (A.U.)')
        plt.savefig(name+'x'+file_num+'.png')
        plt.close()

        plt.figure(3)
        num_bins = 100
        #bins = np.linspace(miny, maxy, num_bins)
        plt.hist(deposited[:, 3], bins=num_bins)
        plt.title(name+' Y Deposition')
        plt.xlabel('y [um]')
        plt.ylabel('Helium Deposition (A.U.)')
        plt.savefig(name+'y'+file_num+'.png')
        plt.close()

        plt.figure(4)
        num_bins_x = 100
        num_bins_y = 100
        #binx = np.linspace(0.0, maxx, num_bins_x)
        #biny = np.linspace(miny, maxy, num_bins_y)
        plt.hist2d(deposited[:, 2], deposited[:, 3], bins=(num_bins_x, num_bins_y))
        #plt.plot(*surface.exterior.xy, color='dimgray')
        #plt.plot(*energy_surface.exterior.xy, '--', color='dimgray')
        #plt.plot(*simulation_surface.exterior.xy, '--', color='dimgray')
        plt.title(name+' 2D Deposition')
        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        #plt.axis('square')
        plt.savefig(name+'2D'+file_num+'.png')
        plt.close()

def plot_distributions(name, ftridyn_name, file_num=1):
    file_num = str(file_num)
    num_bins = 200

    r = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    s = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))

    rf = np.atleast_2d(np.genfromtxt(ftridyn_name+'RFLST.DAT'))
    sf = np.atleast_2d(np.genfromtxt(ftridyn_name+'SPLST.DAT'))

    plt.figure(num='re')
    plots = []
    labels = []
    #bins = np.linspace(0.0, 1000.0, num_bins)
    if np.size(r) > 0:
        _, _, rust = plt.hist(r[:,2], histtype='step',  bins=num_bins, density=True)
        plots.append(rust[0])
        labels.append('rust')
    if np.size(rf) > 0:
        _, _, ftridyn = plt.hist(rf[:,2], histtype='step', bins=num_bins, density=True)
        plots.append(ftridyn[0])
        labels.append('ftridyn')
    if len(plots) > 0: plt.legend(plots, labels, fontsize='small', fancybox=True, shadow=True, loc='best')
    plt.title(name+' Reflected Energy Distributions', fontsize='small')
    plt.xlabel('E [eV]')
    plt.savefig(name+'ref_e.png')
    plt.close()

    plt.figure(num='se')
    plots = []
    labels = []
    #bins=np.linspace(0.0, 100.0, num_bins)
    if np.size(s) > 0:
        _, _, rust = plt.hist(s[:,2], histtype='step', bins=num_bins, density=True)
        plots.append(rust[0])
        labels.append('rust')
    if np.size(sf) > 0:
        _, _, ftridyn = plt.hist(sf[:,2], histtype='step', bins=num_bins, density=True)
        plots.append(ftridyn[0])
        labels.append('ftridyn')
    if len(plots) > 0: plt.legend(plots, labels, fontsize='small', fancybox=True, shadow=True, loc='best')
    plt.title(name+' Sputtered Energy Distributions', fontsize='small')
    plt.xlabel('E [eV]')
    plt.savefig(name+'spt_e.png')
    plt.close()

    plt.figure(num='sa')
    plots = []
    labels = []
    ax = None
    bins = np.linspace(0, np.pi, num_bins)
    #if np.size(s) > 0:
        #hist, bins = np.histogram(np.arccos(-s[:,6]), bins=bins, density=True)
        #rust1 = plt.polar(bins[:-1], hist/np.max(hist))
        #plots.append(rust1[0])
        #labels.append('rustbca cosx')
    #if np.size(sf) > 0:
        #hist, bins = np.histogram(np.arccos(sf[:,6]), bins=bins, density=True)
        #ftridyn1 = plt.polar(bins[:-1], hist/np.max(hist))
        #plots.append(ftridyn1[0])
        #labels.append('ftridyn cosx')
    if np.size(s) > 0:
        hist, bins = np.histogram(np.arccos(s[:,7]), bins=bins, density=True)
        rust2 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust2[0])
        labels.append('rustbca cosy')
        ax = plt.gca()
    if np.size(sf) > 0:
        hist, bins = np.histogram(np.arccos(sf[:,7]), bins=bins, density=True)
        ftridyn2 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(ftridyn2[0])
        labels.append('ftridyn cosy')
        ax = plt.gca()
    if np.size(s) > 0:
        hist, bins = np.histogram(np.arccos(s[:,8]), bins=bins, density=True)
        rust3 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust3[0])
        labels.append('rustbca cosz')
        ax = plt.gca()
    if np.size(sf) > 0:
        hist, bins = np.histogram(np.arccos(sf[:,8]), bins=bins, density=True)
        ftridyn3 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(ftridyn3[0])
        labels.append('ftridyn cosz')
        ax = plt.gca()
    if len(plots) > 0: plt.legend(plots, labels, fontsize='small', bbox_to_anchor=(0.0, 1.05), fancybox=True, shadow=True, loc='upper center')
    plt.title(name+' Sputtered Angular Distributions', fontsize='small')
    if ax:
        ax.set_thetamin(0.)
        ax.set_thetamax(180.)
        ax.set_yticks([0., 0.5, 1.])
        ax.set_xticks([0., np.pi/6., 2*np.pi/6., 3*np.pi/6., 4*np.pi/6., 5*np.pi/6., np.pi])
        ax.set_xticklabels(['0', '30', '60', '90', '', '', ''])
    plt.savefig(name+'spt_ang.png')
    plt.close()

    plt.figure(num='ra')
    plots = []
    labels = []
    ax = None
    bins = np.linspace(0, np.pi, num_bins)
    #if np.size(s) > 0:
        #hist, bins = np.histogram(np.arccos(-r[:,6]), bins=bins, density=True)
        #rust1 = plt.polar(bins[:-1], hist/np.max(hist))
        #plots.append(rust1[0])
        #labels.append('rustbca cosx')
    #if np.size(sf) > 0:
        #hist, bins = np.histogram(np.arccos(rf[:,6]), bins=bins, density=True)
        #ftridyn1 = plt.polar(bins[:-1], hist/np.max(hist))
        #plots.append(ftridyn1[0])
        #labels.append('ftridyn cosx')
    if np.size(r) > 0:
        hist, bins = np.histogram(np.arccos(r[:,7]), bins=bins, density=True)
        rust2 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust2[0])
        labels.append('rustbca cosy')
        ax = plt.gca()
    if np.size(rf) > 0:
        hist, bins = np.histogram(np.arccos(rf[:,7]), bins=bins, density=True)
        ftridyn2 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(ftridyn2[0])
        labels.append('ftridyn cosy')
        ax = plt.gca()
    if np.size(r) > 0:
        hist, bins = np.histogram(np.arccos(r[:,8]), bins=bins, density=True)
        rust3 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(rust3[0])
        labels.append('rustbca cosz')
        ax = plt.gca()
    if np.size(rf) > 0:
        hist, bins = np.histogram(np.arccos(rf[:,8]), bins=bins, density=True)
        ftridyn3 = plt.polar(bins[:-1], hist/np.max(hist))
        plots.append(ftridyn3[0])
        labels.append('ftridyn cosz')
        ax = plt.gca()
    if len(plots) > 0: plt.legend(plots, labels, fontsize='small', bbox_to_anchor=(0.0, 1.05), fancybox=True, shadow=True, loc='upper center')
    plt.title(name+' Reflected Angular Distributions', fontsize='small')
    if ax:
        ax.set_thetamin(0.)
        ax.set_thetamax(180.)
        ax.set_yticks([0., 0.5, 1.])
        ax.set_xticks([0., np.pi/6., 2*np.pi/6., 3*np.pi/6., 4*np.pi/6., 5*np.pi/6., np.pi])
        ax.set_xticklabels(['0', '30', '60', '90', '', '', ''])
    plt.savefig(name+'ref_ang.png')
    plt.close()

def benchmark():
        beam_species = [neon]#, hydrogn]#, helium]#, beryllium, boron, neon, silicon, argon, copper, tungsten]
        target_species = [titanium]#, copper]#, copper]#, boron, silicon, copper]

        N = 1
        N_ = 1000000
        theta = 30.
        thickness = 1000
        depth = 1000
        energies = np.round(np.logspace(1, 2, 20))
        energies = [38.]
        #energies = [15.]
        #energies = [1000.]
        tr = True
        trt = False
        tt = False

        os.system('rm *.png')
        os.system('rm *.output')

        os.system('rm ./rustBCA')
        os.system('cargo build --release')
        os.system('mv target/release/rustBCA .')

        fig1, ax1 = plt.subplots(1, 1, num='yields', figsize=(16,12))
        fig2, ax2 = plt.subplots(1, 1, num='ranges', figsize=(16,12))

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
                name_1.append(beam['symbol']+' on '+target['symbol']+' Y RustBCA')
                name_1.append(beam['symbol']+' on '+target['symbol']+' R RustBCA')
                name_1.append(beam['symbol']+' on '+target['symbol']+' Y F-TRIDYN')
                name_1.append(beam['symbol']+' on '+target['symbol']+' R F-TRIDYN')
                name_1.append(beam['symbol']+' on '+target['symbol']+' Y Yamamura')
                name_1.append(beam['symbol']+' on '+target['symbol']+' R Wierzbicki-Biersack')

                name_2.append(beam['symbol']+' on '+target['symbol']+' Range RustBCA')
                name_2.append(beam['symbol']+' on '+target['symbol']+' Range F-TRIDYN')
                rustbca_yield = []
                rustbca_reflection = []
                ftridyn_yield = []
                ftridyn_reflection = []
                yamamura_yield = []
                rustbca_range = []
                ftridyn_range = []
                reflection_coefficient = []
                Zb = target['Z']
                Ecb = target['Ec']
                Esb = target['Es']
                Eb = target['Eb']
                Mb = target['m']
                n = target['n']
                for energy in energies:

                    rustbca_name = str(energy)+'eV_'+str(theta)+'deg_'+beam['symbol']+'_'+target['symbol']
                    start = time.time()
                    rustbca_yield_, rustbca_reflection_, rustbca_range_ = main(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy, N, N_, theta, thickness, depth, track_trajectories=tt, track_recoils=tr, track_recoil_trajectories=trt, write_files=True, name=rustbca_name)
                    stop = time.time()
                    rustbca_time += stop - start
                    rustbca_yield.append(rustbca_yield_)
                    rustbca_reflection.append(rustbca_reflection_)
                    rustbca_range.append(rustbca_range_)

                    ftridyn_name = beam['symbol'].ljust(2, '_') + target['symbol'].ljust(2, '_')
                    interface = g.tridyn_interface(beam['symbol'], target['symbol'])
                    os.system('rm *.DAT')
                    os.system('rm *.OUT')
                    start = time.time()
                    ftridyn_yield_, ftridyn_reflection_, ftridyn_range_ = interface.run_tridyn_simulations_from_iead([energy], [theta], np.array([[1.]]), number_histories=np.max([100, N_]), depth=depth*MICRON/ANGSTROM)
                    stop = time.time()
                    ftridyn_time += stop - start
                    ftridyn_yield.append(ftridyn_yield_/np.max([100, N_]))
                    ftridyn_reflection.append(ftridyn_reflection_/np.max([100, N_]))
                    ftridyn_range.append(ftridyn_range_)

                    yamamura_yield_ = yamamura(beam, target, energy)
                    yamamura_yield.append(yamamura_yield_)

                    reflection_coefficient_ = wierzbicki_biersack(beam, target, energy)
                    reflection_coefficient.append(reflection_coefficient_)

                    do_plots(rustbca_name, file_num=str(1))
                    plot_distributions(rustbca_name, ftridyn_name, file_num=1)

                handle = ax1.loglog(energies, rustbca_yield, '*')
                ax1.loglog(energies, rustbca_reflection, ':', color=handle[0].get_color())

                ax1.loglog(energies, ftridyn_reflection, '.', color=handle[0].get_color())
                ax1.loglog(energies, ftridyn_yield, 'o', color=handle[0].get_color())
                ax1.loglog(energies, yamamura_yield, '-', color=handle[0].get_color())
                ax1.loglog(energies, reflection_coefficient, '-.', color=handle[0].get_color())
                #yamamura_table = yamamura_tables[beam['symbol']][target['symbol']]
                #plt.loglog(yamamura_table[:,0], yamamura_table[:,1], '-', color=handle[0].get_color())
                ax2.loglog(energies, rustbca_range, color=handle[0].get_color())
                ax2.loglog(energies, np.array(ftridyn_range)/MICRON, '--', color=handle[0].get_color())

        ax1.legend(name_1, fontsize='small', loc='best')
        ax1.axis([0.5*np.min(energies), 2.0*np.max(energies), 1./N_, 10.])
        ax1.set_ylabel('Y/R [at/ion]/[ion/ion]')
        ax1.set_xlabel('E [eV]')
        fig1.savefig(rustbca_name+'yields.png')

        ax2.legend(name_2, fontsize='small', loc='best')
        ax2.set_xlabel('E [eV]')
        ax2.set_ylabel('Range [um]')
        ax2.axis([0.5*np.min(energies), 2*np.max(energies), 0.5*np.min(rustbca_range), 2*np.max(rustbca_range)])
        fig2.savefig(rustbca_name+'ranges.png')

        print(f'time rustBCA: {rustbca_time} time F-TRIDYN: {ftridyn_time}')

if __name__ == '__main__':

    benchmark()
