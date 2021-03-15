import os
import time

import toml
import numpy as np
from shapely.geometry import Point, Polygon, box
from scipy import constants
import matplotlib.pyplot as plt
from matplotlib import rcParams, cm
import matplotlib as mpl
import matplotlib.colors as colors

from .materials import *
from .formulas import *

#from generate_ftridyn_input import *

rcParams.update({'figure.autolayout': True})

#Constants
Q = constants.physical_constants["elementary charge"][0]
PI = constants.pi
AMU = constants.physical_constants["unified atomic mass unit"][0]
ANGSTROM = constants.angstrom
MICRON = constants.micro
NM = constants.nano
CM = constants.centi
EPS0 = constants.epsilon_0
A0 = constants.physical_constants["Bohr radius"][0]
K = constants.physical_constants["atomic unit of permittivity"][0]
ME = constants.physical_constants["electron mass"][0]
SQRTPI = np.sqrt(PI)
SQRT2PI = np.sqrt(2 * PI)
C = constants.physical_constants["speed of light in vacuum"][0]

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
        colormap = cm.ScalarMappable(norm=colors.SymLogNorm( min_Z*4, vmin=min_Z, vmax=max_Z), cmap='tab20')

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
            x = [x_ for (x_, y_) in boundary]
            y = [y_ for (x_, y_) in boundary]
            x.append(x[0])
            y.append(y[0])
            plt.plot(x, y, linewidth=3, color="dimgray")

        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
        plt.title(name+' Trajectories')
        plt.axis('square')
        if show: plt.show()
        plt.savefig(name+'trajectories_.png')
        plt.close()

def do_trajectory_plot_3d(name, thickness=None, depth=None, boundary=None, plot_final_positions=True, plot_origins=True, radius=None, cube_length=None):
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

    from mayavi.mlab import points3d, plot3d, mesh

    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.output', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.output', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.output', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.output', delimiter=','))
    trajectory_data = np.atleast_1d(np.genfromtxt(name+'trajectory_data.output', delimiter=',').transpose().astype(int))

    if np.size(trajectories) > 0:
        min_Z = np.min(trajectories[:, 1])
        max_Z = np.max(trajectories[:, 1])
        colormap = cm.ScalarMappable(norm=colors.SymLogNorm( min_Z*4, vmin=min_Z, vmax=max_Z), cmap='tab20')

    index = 0
    x_max = 0
    min_length = 1
    scale_factor = 100.0
    if np.size(trajectories) > 0:
        for trajectory_length in trajectory_data:

            M = trajectories[index, 0]
            Z = trajectories[index, 1]
            E = trajectories[index:(trajectory_length + index), 2]
            x = trajectories[index:(trajectory_length + index), 3]*scale_factor
            y = trajectories[index:(trajectory_length + index), 4]*scale_factor
            z = trajectories[index:(trajectory_length + index), 5]*scale_factor

            if np.max(x) > x_max:
                x_max = np.max(x)

            if plot_origins: points3d(x[0], y[0], z[0], color=colormap.to_rgba(Z)[:3], scale_factor=1)

            if len(x) > min_length: plot3d(x, y, z, color = colormap.to_rgba(Z)[:3], tube_radius=0.5)

            index += trajectory_length

        if plot_final_positions:
            if np.size(sputtered) > 0:
                sputtered_colors = [colormap.to_rgba(Z)[:3] for Z in sputtered[:,1]]
                for x, y, z, c in zip(sputtered[:,3], sputtered[:,4], sputtered[:,5], sputtered_colors):
                    points3d(x*scale_factor, y*scale_factor, z*scale_factor, color=c, scale_factor=2)

            if np.size(reflected) > 0:
                reflected_colors = [colormap.to_rgba(Z)[:3] for Z in reflected[:,1]]
                for x, y, z, c in zip(reflected[:,3], reflected[:,4], reflected[:,5], reflected_colors):
                    points3d(x*scale_factor, y*scale_factor, z*scale_factor, color=c, scale_factor=2)

            if np.size(deposited) > 0:
                deposited_colors = [colormap.to_rgba(Z)[:3] for Z in deposited[:,1]]
                for x, y, z, c in zip(deposited[:,2], deposited[:,3], deposited[:,4], deposited_colors):
                    points3d(x*scale_factor, y*scale_factor, z*scale_factor, color=c, scale_factor=2)

        if boundary:
            x = [x_ for (x_, y_) in boundary]
            y = [y_ for (x_, y_) in boundary]
            z = [0. for (x_, y_) in boundary]
            x.append(x[0])
            y.append(y[0])
            z.append(z[0])
            plot3d(x, y, z, color=(0.1, 0.1, 0.1))

        if radius:
            [phi, theta] = np.mgrid[0:2 * np.pi:64j, 0:np.pi:64j]
            x = np.cos(phi)*np.sin(theta)
            y = np.sin(phi)*np.sin(theta)
            z = np.cos(theta)
            mesh(radius*x, radius*y, radius*z, color=(0.1,0.7,0.3), opacity=0.2)

        if cube_length:
            faces = []
            
            xmin = -cube_length/2.*scale_factor
            xmax = cube_length/2.*scale_factor
            ymin = -cube_length/2.*scale_factor
            ymax = cube_length/2.*scale_factor
            zmin = -cube_length/2.*scale_factor
            zmax = cube_length/2.*scale_factor

            x,y = np.mgrid[xmin:xmax:3j,ymin:ymax:3j]
            z = np.ones(y.shape)*zmin
            faces.append((x,y,z))

            x,y = np.mgrid[xmin:xmax:3j,ymin:ymax:3j]
            z = np.ones(y.shape)*zmax
            faces.append((x,y,z))

            x,z = np.mgrid[xmin:xmax:3j,zmin:zmax:3j]
            y = np.ones(z.shape)*ymin
            faces.append((x,y,z))

            x,z = np.mgrid[xmin:xmax:3j,zmin:zmax:3j]
            y = np.ones(z.shape)*ymax
            faces.append((x,y,z))

            y,z = np.mgrid[ymin:ymax:3j,zmin:zmax:3j]
            x = np.ones(z.shape)*xmin
            faces.append((x,y,z))

            y,z = np.mgrid[ymin:ymax:3j,zmin:zmax:3j]
            x = np.ones(z.shape)*xmax
            faces.append((x,y,z))

            for grid in faces:
                x,y,z = grid
                mesh(x, y, z, opacity=0.4, color=(0.1,0.7,0.3))


def generate_rustbca_input(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, E0, N, N_, theta,
    thickness, depth, track_trajectories=True, track_recoils=True,
    track_recoil_trajectories=True, name='test_',
    track_displacements=True, track_energy_losses=True,
    electronic_stopping_mode=LOW_ENERGY_NONLOCAL,
    weak_collision_order=3, ck=1., mean_free_path_model=LIQUID,
    interaction_potential='"KR_C"', high_energy=False, energy_barrier_thickness=(6.306E10)**(-1/3.),
    initial_particle_position = -1*ANGSTROM/MICRON, integral='"MENDENHALL_WELLER"',
    root_finder = '{"NEWTON"={max_iterations=100, tolerance=1e-3}}',
    delta_x_angstrom=5., uniformly_distributed_ions=False):

    '''
    Generates a rustbca input file. Assumes eV, amu, and microns for units.
    '''

    options = {
        'name': name,
        'track_trajectories': track_trajectories,
        'track_recoils': track_recoils,
        'track_recoil_trajectories': track_recoil_trajectories,
        'write_buffer_size': 8000,
        'weak_collision_order': weak_collision_order,
        'suppress_deep_recoils': False,
        'high_energy_free_flight_paths': high_energy,
        'num_threads': 4,
        'num_chunks': 100,
        'use_hdf5': False,
        'electronic_stopping_mode': electronic_stopping_mode,
        'mean_free_path_model': mean_free_path_model,
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
        'surface_binding_model': "AVERAGE",
        'bulk_binding_model': "AVERAGE"
    }

    dx = delta_x_angstrom*ANGSTROM/MICRON

    minx, miny, maxx, maxy = 0.0, -thickness/2., depth, thickness/2.
    surface = box(minx, miny, maxx, maxy)

    simulation_surface = surface.buffer(10*dx, cap_style=2, join_style=2)
    geometry_input = {
        'length_unit': 'MICRON',
        'triangles': [[0., depth, 0., thickness/2., -thickness/2., -thickness/2.], [0., depth, depth, thickness/2., thickness/2., -thickness/2.]],
        'densities': [np.array(n)*(MICRON)**3, np.array(n)*(MICRON)**3],
        'material_boundary_points': [[0., thickness/2.], [depth, thickness/2.], [depth, -thickness/2.], [0., -thickness/2.]],
        'simulation_boundary_points':  list(simulation_surface.exterior.coords),
        'energy_barrier_thickness': energy_barrier_thickness,
        'electronic_stopping_correction_factors': [ck, ck],
    }

    cosx = np.cos(theta*np.pi/180.)
    sinx = np.sin(theta*np.pi/180.)
    if uniformly_distributed_ions:
        positions = [(initial_particle_position, np.random.uniform(-thickness/2., thickness/2.), 0.) for _ in range(N)]
    else:
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
        'geometry_input': geometry_input,
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

def plot_energy_loss(name, N, num_bins=50, thickness=None, depth=None):
    '''
    Plots energy loss plots separated by electronic and nuclear losses from rustbca.

    Args:
        name (string): name of rustbca simulation
        N (int): number of incident ions
        num_bins (int): number of histogram bins; default 100
    '''

    #energy_loss = np.genfromtxt(f'{name}energy_loss.output', delimiter=',', dtype=np.float32)
    energy_loss = np.loadtxt(f'{name}energy_loss.output', delimiter=',', usecols=[2,3,4,5], dtype=np.float32)

    En = energy_loss[:,0]
    Ee = energy_loss[:,1]
    x = energy_loss[:,2]
    y = energy_loss[:,3]

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
    #plt.axis('square')
    plt.title(f'Nuclear Energy Loss {name}')
    Hn, bins_x, bins_y, _ = plt.hist2d(x, y, weights=En/N, bins=num_bins, norm=colors.SymLogNorm(0.1))

    plt.savefig(name+'nuclear_energy_loss.png')
    plt.close()

    plt.figure(3)
    plt.xlabel('x [um]')
    plt.ylabel('y [um]')
    plt.title(f'Electronic Energy Loss {name}')
    #plt.axis('square')
    He, bins_x, bins_y, _ = plt.hist2d(x, y, weights=Ee/N, bins=num_bins, norm=colors.SymLogNorm(0.1))
    plt.savefig(name+'electronic_energy_loss.png')
    plt.close()

    np.savetxt(f'{name}_electronic_energy_loss_2D.dat', He)
    np.savetxt(f'{name}_electronic_energy_loss_2D_bins.dat', [bins_x, bins_y])
    np.savetxt(f'{name}_nuclear_energy_loss_2D.dat', Hn)
    np.savetxt(f'{name}_nuclear_energy_loss_2D_bins.dat', [bins_x, bins_y])

    d = np.genfromtxt(f'{name}deposited.output', delimiter=',')
    plt.figure(4)
    plt.xlabel('x [um]')
    plt.ylabel('y [um]')
    plt.title(f'Deposition {name}')
    #plt.axis('square')
    plt.hist2d(d[:,2], d[:,3], bins=num_bins, norm=colors.SymLogNorm(0.1))
    plt.savefig(name+'deposited2d.png')
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
        'write_buffer_size': 8000,
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
        'surface_binding_model': "AVERAGE",
        'bulk_binding_model': "AVERAGE"
    }

    dx = 5.*ANGSTROM/MICRON

    minx, miny, maxx, maxy = 0.0, -thickness/2., depth, thickness/2.
    surface = box(minx, miny, maxx, maxy)

    simulation_surface = surface.buffer(10.*dx, cap_style=2, join_style=2)
    geometry_input = {
        'length_unit': 'MICRON',
        'energy_barrier_thickness': sum(n)**(-1./3.)/np.sqrt(2.*np.pi),
        'triangles': [[0., depth, 0., thickness/2., -thickness/2., -thickness/2.], [0., depth, depth, thickness/2., thickness/2., -thickness/2.]],
        'densities': [n, n],
        'material_boundary_points': [[0., thickness/2.], [depth, thickness/2.], [depth, -thickness/2.], [0., -thickness/2.]],
        'simulation_boundary_points':  list(simulation_surface.exterior.coords),
        'electronic_stopping_correction_factors': [1.0, 1.0],
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
        'geometry_input': geometry_input,
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

def beam_target(ions, target, energy, angle, N_=10000, N=1, run_sim=True,
    high_energy=False, integral='"GAUSS_LEGENDRE"',
    interaction_potential='"ZBL"',
    root_finder='{"NEWTON"={max_iterations=100, tolerance=1e-6}}', tag=None,
    track_trajectories = False,
    plot_trajectories = False,
    plot_distributions = False,
    track_energy_losses = False,
    track_displacements = False,
    plot_depth_distributions = False,
    track_recoils = True,
    do_plots = False,
    thickness=100,
    depth=100, ck=1.,
    weak_collision_order=3,
    electronic_stopping_mode=LOW_ENERGY_NONLOCAL,
    uniformly_distributed_ions=False,
    mean_free_path_model='"LIQUID"'):

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

    name = ions['name']+'_'+target['name']+'_'
    if tag: name += tag+'_'

    pmax = (target['n']**(-1/3)/MICRON)/np.sqrt(PI)

    generate_rustbca_input(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy, N, N_,
        angle, thickness, depth, name=name, track_trajectories=track_trajectories,
        track_recoil_trajectories=track_trajectories, track_energy_losses=track_energy_losses,
        track_displacements=track_displacements, track_recoils=track_recoils,
        electronic_stopping_mode=electronic_stopping_mode, high_energy=high_energy,
        interaction_potential=interaction_potential, integral=integral,
        root_finder=root_finder, delta_x_angstrom=2.*pmax*MICRON/ANGSTROM,
        initial_particle_position = -2.01*pmax, weak_collision_order=weak_collision_order, ck=ck,
        uniformly_distributed_ions=uniformly_distributed_ions, mean_free_path_model=mean_free_path_model)

    if run_sim: os.system(f'rustBCA.exe {name}.toml')

    if do_plots:
        if track_energy_losses: plot_energy_loss(name, N*N_, num_bins=100)

        if plot_depth_distributions and track_displacements and track_recoils: plot_all_depth_distributions(name, 38., N_*N, num_bins=200)

        if plot_distributions: plot_distributions_rustbca(name, ions, target, incident_angle=angle, incident_energy=energy)

        if track_displacements and track_recoils: plot_displacements(name, 20., num_bins=200)

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

    thickness = 10
    depth = 10

    name = ions['name']+'_'+target['name']+'_'

    track_trajectories = False
    plot_trajectories = False
    plot_distributions = False
    track_energy_losses = False
    track_displacements = False
    plot_depth_distributions = False

    integral = '"GAUSS_LEGENDRE"'
    options = ['"KR_C"', '"MOLIERE"', '"ZBL"']

    ftridyn = tridyn_interface(ions['symbol'], target['symbol'])

    legends = []
    for option_index, option in enumerate(options):
        s = []
        r = []

        sf = []
        rf = []
        for index, energy in enumerate(energies):

            generate_rustbca_input(Zb, Mb, n, Eca, Ecb, Esa, Esb, Eb, Ma, Za, energy, N, N_,
                angle, thickness, depth, name=name+str(index)+'_'+str(option_index), track_trajectories=track_trajectories,
                track_recoil_trajectories=track_trajectories, track_energy_losses=track_energy_losses,
                track_displacements=track_displacements, track_recoils=True,
                electronic_stopping_mode="LOW_ENERGY_EQUIPARTITION", integral=integral,
                interaction_potential = option)

            if run_sim: os.system(f'rustBCA.exe {name+str(index)+"_"+str(option_index)}.toml')

            sputtered = np.atleast_2d(np.genfromtxt(name+str(index)+'_'+str(option_index)+'sputtered.output', delimiter=','))
            reflected = np.atleast_2d(np.genfromtxt(name+str(index)+'_'+str(option_index)+'reflected.output', delimiter=','))

            if np.size(reflected) > 0:
                r.append(np.shape(reflected)[0]/(N_*N))
            else:
                r.append(0.0)

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

        plt.loglog(energies, s, linestyle='-', marker='*')
        legends.append(option+' Y')
        plt.loglog(energies, r, linestyle='--', marker='*')
        legends.append(option+' R')

    for energy in energies:
        sputteredf, reflectedf, _ = ftridyn.run_tridyn_simulations_from_iead([energy], [angle], [[N_]], number_histories=N)

        if np.size(reflectedf) > 0:
            rf.append(np.shape(reflectedf)[0]/(N_*N))
        else:
            rf.append(0.0)

        if np.size(sputteredf) > 0:
            sf.append(np.shape(sputteredf)[0]/(N_*N))
        else:
            sf.append(0.0)

    plt.loglog(energies, sf, linestyle='-.', marker='^')
    legends.append('F-TRIDYN Y')
    plt.loglog(energies, rf, linestyle='-.', marker='^')
    legends.append('F-TRIDYN R')

    energies_plot = np.logspace(0, 6, 10000)
    y = [yamamura(ions, target, energy) for energy in energies_plot]
    legends.append('Yamamura')
    bl = [bohdansky_light_ion(ions, target, energy) for energy in energies_plot]
    legends.append('Bohdansky Light Ion')
    bh = [bohdansky_heavy_ion(ions, target, energy) for energy in energies_plot]
    legends.append('Bohdansky Heavy Ion')
    rt = [thomas_reflection(ions, target, energy) for energy in energies_plot]
    legends.append('Thomas et al.')
    rwb = [wierzbicki_biersack(ions, target, energy) for energy in energies_plot]
    legends.append('Wierzbicki-Biersack')

    plt.loglog(energies_plot, y, alpha=0.5, linewidth=3)
    plt.loglog(energies_plot, bl, alpha=0.5, linewidth=3)
    plt.loglog(energies_plot, bh, alpha=0.5, linewidth=3)
    plt.loglog(energies_plot, rt, linestyle='--', alpha=0.5, linewidth=3)
    plt.loglog(energies_plot, rwb, linestyle='--', alpha=0.5, linewidth=3)
    plt.legend(legends)
    axis = plt.gca()
    axis.set_ylim(np.min(np.array(s)[np.array(s)>0]))
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

def metal_oxide_bilayer_iead(ions, iead, energies, angles, metal,
    metal_stoichiometry, oxygen_stoichiometry, metal_oxide_number_density,
    layer_depth, total_depth, thickness, N=0.01, tag='', run_sim=True,
    spread=0.0):

    import itertools

    name = metal['symbol']+oxygen['symbol']+str(oxygen_stoichiometry)+'_IEAD_'+tag
    energy_angle_pairs = list(itertools.product(energies, angles))
    E0 = np.array([pair[0] for pair in energy_angle_pairs])
    theta = np.array([pair[1] for pair in energy_angle_pairs])

    #skip last row of hPIC input because it's usually garbage
    N_ = np.array([int(np.floor(iead[i, j]*N)) for i in range(len(energies)) for j in range(len(angles))], dtype=int)

    #filter out zeros
    E0 = E0[N_>0]
    theta = theta[N_>0]
    N_ = N_[N_>0]

    #build geometry
    minx, miny, maxx, maxy = 0.0, -thickness/2., total_depth, thickness/2.
    material_boundary = box(minx, miny, maxx, maxy)
    dx = 5*ANGSTROM/MICRON
    simulation_boundary = material_boundary.buffer(10.*dx, cap_style=2, join_style=2)

    #x1, x2, x3, y1, y2, y3
    triangles = [
        [0., layer_depth, 0., -thickness/2., thickness/2., thickness/2.],
        [0., layer_depth, layer_depth, -thickness/2., -thickness/2., thickness/2.],
        [layer_depth, layer_depth, total_depth, -thickness/2., thickness/2., thickness/2.],
        [layer_depth, total_depth, total_depth, -thickness/2., -thickness/2., thickness/2.]
    ]

    total_stoichiometry = metal_stoichiometry + oxygen_stoichiometry
    metal_number_density = metal_oxide_number_density*metal_stoichiometry/total_stoichiometry
    oxygen_number_density = metal_oxide_number_density*oxygen_stoichiometry/total_stoichiometry

    densities = [
        [metal_number_density, oxygen_number_density],
        [metal_number_density, oxygen_number_density],
        [metal['n']*(MICRON**3), 0.],
        [metal['n']*(MICRON**3), 0.],
    ]

    E0 = E0[N_>0]
    theta = theta[N_>0]
    N_ = N_[N_>0]

    Za = ions['Z']
    Ma = ions['m']
    Esa = ions['Es']
    Eca = ions['Ec']

    options = {
        'name': name,
        'track_trajectories': True,
        'track_recoils': True,
        'track_recoil_trajectories': True,
        'write_buffer_size': 8000,
        'weak_collision_order': 3,
        'suppress_deep_recoils': False,
        'high_energy_free_flight_paths': False,
        'num_threads': 8,
        'num_chunks': 10,
        'use_hdf5': False,
        'electronic_stopping_mode': LOW_ENERGY_LOCAL,
        'mean_free_path_model': LIQUID,
        'interaction_potential': [["KR_C"]],
        'scattering_integral': [["MENDENHALL_WELLER"]],
        'track_displacements': True,
        'track_energy_losses': True,
    }

    material_parameters = {
        'energy_unit': 'EV',
        'mass_unit': 'AMU',
        'Eb': [metal['Eb'], oxygen['Eb']],
        'Es': [metal['Es'], oxygen['Es']],
        'Ec': [metal['Ec'], oxygen['Ec']],
        'Z': [metal['Z'], oxygen['Z']],
        'm': [metal['m'], oxygen['m']],
        'interaction_index': [0, 0],
        'bulk_binding_model': "AVERAGE",
        'surface_binding_model': "AVERAGE"
    }

    geometry_input = {
        'length_unit': 'MICRON',
        'energy_barrier_thickness': sum(densities[0])**(-1./3.)/np.sqrt(2.*np.pi),
        'triangles': triangles,
        'densities': densities,
        'material_boundary_points': list(material_boundary.exterior.coords),
        'simulation_boundary_points':  list(simulation_boundary.exterior.coords),
        'electronic_stopping_correction_factors': np.ones(len(densities)),
    }

    cosx = np.cos(theta*np.pi/180.)
    sinx = np.sin(theta*np.pi/180.)
    positions = [(-dx/2., np.random.uniform(-spread/2., spread/2.), 0.) for _ in range(len(N_))]

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
        'geometry_input': geometry_input,
        'options': options,
    }

    with open(f'{name}.toml', 'w') as file:
        toml.dump(input_file, file, encoder=toml.TomlNumpyEncoder())
    with open(f'{name}.toml', 'a') as file:
        file.write(r'root_finder = [[{"NEWTON"={max_iterations = 100, tolerance=1E-3}}]]')

    if run_sim: os.system(f'cargo run --release {name}.toml')

    #do_trajectory_plot(name, thickness=thickness, depth=layer_depth)
    plot_distributions_rustbca(name, ions, metal, incident_energy=np.max(energies))

    s = np.atleast_2d(np.genfromtxt(f'{name}sputtered.output', delimiter=','))
    r = np.atleast_2d(np.genfromtxt(f'{name}reflected.output', delimiter=','))
    d = np.atleast_2d(np.genfromtxt(f'{name}deposited.output', delimiter=','))

    return np.sum(N_), s, r, d

def helium_on_tungsten_oxide_layer():
    ions = helium
    iead = np.genfromtxt('cprobe_1_IEAD_sp0.dat')
    Te_eV = 1.62167190E1#, 7.69748285, 3.14783250, 9.55495830E-1, 2.90032039E-1]
    energies = np.linspace(0.0, 48.0*Te_eV, 480)
    angles = np.linspace(0.001, 89.99, 90)

    plt.figure(1)
    plt.pcolormesh(angles, energies, iead)
    plt.title('He IEAD WEST Collector Probe Sample 1')
    plt.xlabel('angles [deg]')
    plt.ylabel('energies [eV]')
    plt.savefig('iead_west_1.png')
    plt.close()

    metal = tungsten
    metal['Ec'] = 1.0
    metal_stoichiometry = 1
    oxygen_stoichiometry = 2

    porosities = [0., 0.2, 0.5]
    layer_depths = [50*ANGSTROM/MICRON, 100*ANGSTROM/MICRON, 200*ANGSTROM/MICRON]

    thickness = 10000*ANGSTROM/MICRON
    total_depth = 10000*ANGSTROM/MICRON

    num_bins = 100

    for j, layer_depth in enumerate(layer_depths):
        handles = []
        legends = []

        for i, porosity in enumerate(porosities):

            metal_oxide_number_density = (1. - porosity)*1.004E10

            N, s, r, d = metal_oxide_bilayer_iead(ions, iead, energies, angles, metal,
                metal_stoichiometry, oxygen_stoichiometry, metal_oxide_number_density,
                layer_depth, total_depth, thickness, N=0.01, tag=f'{str(i)}_{str(j)}_',
                run_sim=True, spread=5*ANGSTROM/MICRON)

            #deposition profile
            plt.figure(j)
            heights, bins, patches = plt.hist(d[:, 2], bins=num_bins, histtype='step', linewidth=2)
            handles.append(patches[0])
            legends.append(f'Porosity: {porosity*100}% R: {np.round(np.shape(r)[0]/N,2)} Y: {np.round(np.shape(s)[0]/N,2)} [at/ion]')

        plt.text(layer_depth - (layer_depths[0])/8, 0.0, 'WO2', horizontalalignment='right')
        plt.text(layer_depth + (layer_depths[0])/8, 0.0, 'W', horizontalalignment='left')
        plt.plot([layer_depth, layer_depth], [0., 2*np.max(heights)], linewidth=3, color='dimgray')
        plt.title(f'He on WO2-W, Δx = {layer_depth*MICRON/ANGSTROM} Angstrom')
        plt.xlabel('x [um]')
        plt.ylabel('Counts')
        plt.legend(handles, legends)
        plt.axis([0., 1.5*np.max(d[:,2]), 0., 1.2*np.max(heights)])
        plt.savefig(f'dep_west_1_wo2_w_{j}.png')

def main():
    '''
    Here an example usage of beam_target is shown. This code runs rustbca and produces plots.

    For helium on copper at 1 keV and 0 degrees angle of incidence,
    all the distributions and trajectories are plotted.
    '''

    energy = 2000
    angle = 0.00001
    N = 10000

    beam_target(helium, copper, energy, angle, N_=N, N=1, do_plots=True, run_sim=True,
        plot_trajectories=True, track_trajectories=False, thickness=50, depth=1000.,
        interaction_potential=r"'KR_C'", track_energy_losses=True, track_displacements=False,
        uniformly_distributed_ions=True, high_energy=False, electronic_stopping_mode=INTERPOLATED,
        weak_collision_order=0, plot_distributions=True, mean_free_path_model='LIQUID')

if __name__ == '__main__':
    main()
