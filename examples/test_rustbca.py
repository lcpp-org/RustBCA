from libRustBCA import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'/../scripts')
sys.path.append('scripts')
from materials import *
from formulas import *
import time

def main():

    #test rotation to and from RustBCA coordinates

    #nx, ny, nz is the normal vector (out of surface)
    nx = -np.sqrt(2)/2
    ny = -np.sqrt(2)/2
    nz = 0.0

    #ux, uy, uz is the particle direction (simulation coordinates)
    ux = 1.0
    uy = 0.0
    uz = 0.0

    print(f'Before rotation: ({ux}, {uy}, {uz})')
    ux, uy, uz = rotate_given_surface_normal_py(nx, ny, nz, ux, uy, uz)
    print(f'After rotation: ({ux}, {uy}, {uz})')
    ux, uy, uz = rotate_back_py(nx, ny, nz, ux, uy, uz)
    print(f'After rotation back: ({ux}, {uy}, {uz})')

    #After rotating and rotating back, effect should be where you started (minus fp error)
    assert(abs(ux - 1.0) < 1e-6)
    assert(abs(uy) < 1e-6)
    assert(abs(uz) < 1e-6)

    #test vectorized rotation to and from RustBCA coordinates

    num_rot_test = 1000000

    #nx, ny, nz is the normal vector (out of surface)
    nx = [-np.sqrt(2)/2]*num_rot_test
    ny = [-np.sqrt(2)/2]*num_rot_test
    nz = [0.0]*num_rot_test

    #ux, uy, uz is the particle direction (simulation coordinates)
    ux = [1.0]*num_rot_test
    uy = [0.0]*num_rot_test
    uz = [0.0]*num_rot_test

    start = time.time()
    ux, uy, uz = rotate_given_surface_normal_vec_py(nx, ny, nz, ux, uy, uz)
    ux, uy, uz = rotate_back_vec_py(nx, ny, nz, ux, uy, uz)
    stop = time.time()
    print(f'Time to rotate: {(stop - start)/num_rot_test} sec/vector')

    #After rotating and rotating back, effect should be where you started (minus fp error)
    assert(abs(ux[0] - 1.0) < 1e-6)
    assert(abs(uy[0]) < 1e-6)
    assert(abs(uz[0]) < 1e-6)

    #scripts/materials.py has a number of potential ions and targets
    ion = helium
    ion['Eb'] = 0.0
    target = tungsten

    #tungsten Eb = 0 is a 'worst-case' - accepted literature value is 3 eV
    target['Eb'] = 3.0

    energy = 2500.0 #eV
    angle = 0.0 #degrees
    num_samples = 1000
    #For automatic control, libRustBCA offers two helpful functions for the sputtering yield and reflection coefficients
    Y = sputtering_yield(ion, target, energy, angle, num_samples)

    print(f'Sputtering yield for {ion["symbol"]} on {target["symbol"]} at {energy} eV is {Y} at/ion. Yamamura predicts { np.round(yamamura(ion, target, energy),3)} at/ion.')

    R_N, R_E = reflection_coefficient(ion, target, energy, angle, num_samples)
    print(f'Particle reflection coefficient for {ion["symbol"]} on {target["symbol"]} at {energy} eV is {R_N}. Thomas predicts {np.round(thomas_reflection(ion, target, energy), 3)}.')
    print(f'Energy reflection coefficient for {ion["symbol"]} on {target["symbol"]} at {energy} eV is {R_E}')

    R_N, R_E = compound_reflection_coefficient(ion, [target, ion], [target['n'], 0.1*target['n']], energy, angle, num_samples)
    print(f'Particle reflection coefficient for {ion["symbol"]} on {ion["symbol"]}x{target["symbol"]} where x=0.1 at {energy} eV is {R_N}. Thomas predicts {np.round(thomas_reflection(ion, target, energy), 3)}.')
    print(f'Energy reflection coefficient for {ion["symbol"]}x{target["symbol"]} where x=0.1 at {energy} eV is {R_E}')

    vx0 = 1e5
    vy0 = 1e5
    vz0 = 0.0
    vx1, vy1, vz1 = reflect_single_ion_py(ion, target, vx0, vy0, vz0)

    print(f'(vx, vy, vz) before reflection: ({vx0}, {vy0}, {vz0})')
    print(f'(vx, vy, vz) after reflection: ({vx1}, {vy1}, {vz1})')

    #For smooth distributions and good statistics, you should use at least 10k ions
    number_ions = 10000

    #1 keV is above the He on W sputtering threshold of ~150 eV
    energies_eV = 1000.0*np.ones(number_ions)

    #Working with angles of exactly 0 is problematic due to gimbal lock
    angle = 0.0001

    #In RustBCA's 0D geometry, +x -> into the surface
    ux = np.cos(angle*np.pi/180.)*np.ones(number_ions)
    uy = np.sin(angle*np.pi/180.)*np.ones(number_ions)
    uz = np.zeros(number_ions)

    print(f'Running RustBCA for {number_ions} {ion["symbol"]} ions on {target["symbol"]} at {energies_eV[0]/1000.} keV...')
    print(f'This may take several minutes.')

    start = time.time()
    #Note that simple_bca_list_py expects number densities in 1/Angstrom^3
    output = np.array(simple_bca_list_py(energies_eV, ux, uy, uz, ion['Z'],
        ion['m'], ion['Ec'], ion['Es'], target['Z'], target['m'],
        target['Ec'], target['Es'], target['n']/10**30, target['Eb']))
    stop = time.time()
    delta_time = stop - start

    print('Simulation complete. Processing data...')

    Z = output[:, 0]
    m = output[:, 1]
    E = output[:, 2]
    x = output[:, 3]
    y = output[:, 4]
    z = output[:, 5]
    ux = output[:, 6]
    uy = output[:, 7]
    uz = output[:, 8]

    #For the python bindings, these conditionals can be used to distinguish
    #between sputtered, reflected, and implanted particles in the output list
    sputtered = output[np.logical_and(Z == target['Z'], E > 0), :]
    reflected = output[np.logical_and(Z == ion['Z'], x < 0), :]
    implanted = output[np.logical_and(Z == ion['Z'], x > 0), :]

    plt.figure(1)
    plt.title(f'Sputtered {target["symbol"]} Energy Distribution')
    plt.xlabel('E [eV]')
    plt.ylabel('f(E) [A.U.]')
    plt.hist(sputtered[:, 2], bins=100, density=True, histtype='step')

    plt.figure(2)
    plt.title(f'Implanted {ion["symbol"]} Depth Distribution')
    plt.xlabel('x [A]')
    plt.ylabel('f(x) [A.U.]')
    plt.hist(implanted[:, 3], bins=100, density=True, histtype='step')

    thomas = thomas_reflection(ion, target, energies_eV[0])
    yamamura_yield = yamamura(ion, target, energies_eV[0])

    print('Data processing complete.')
    print(f'RustBCA Y: {len(sputtered[:, 0])/number_ions} Yamamura Y: {yamamura_yield}')
    print(f'RustBCA R: {len(reflected[:, 0])/number_ions} Thomas R: {thomas}')
    print(f'Time per ion: {delta_time/number_ions*1e3} us/{ion["symbol"]}')

    #Next up is the layered target version. I'll add a 50 Angstrom layer of W-H to the top of the target.

    #1 keV is above the He on W sputtering threshold of ~150 eV
    energies_eV = 1000.0*np.ones(number_ions)

    #Working with angles of exactly 0 is problematic due to gimbal lock
    angle = 0.0001

    #In RustBCA's 0D geometry, +x -> into the surface
    ux = np.cos(angle*np.pi/180.)*np.ones(number_ions)
    uy = np.sin(angle*np.pi/180.)*np.ones(number_ions)
    uz = np.zeros(number_ions)

    print(f'Running RustBCA for {number_ions} {ion["symbol"]} ions on {target["symbol"]} with hydrogenated layer at {energies_eV[0]/1000.} keV...')
    print(f'This may take several minutes.')
    #Note the different argument order; when a breaking change is due, this will
    #be back-ported to the other bindings as well for consistency.
    output, incident, stopped = compound_bca_list_1D_py(
        ux, uy, uz, energies_eV, [ion['Z']]*number_ions,
        [ion['m']]*number_ions, [ion['Ec']]*number_ions, [ion['Es']]*number_ions, [target['Z'], 1.0], [target['m'], 1.008],
        [target['Ec'], 1.0], [target['Es'], 1.5], [target['Eb'], 0.0], [[target['n']/10**30, target['n']/10**30], [target['n']/10**30, 0.0]], [50.0, 1e6]
    )


    output = np.array(output)

    Z = output[:, 0]
    m = output[:, 1]
    E = output[:, 2]
    x = output[:, 3]
    y = output[:, 4]
    z = output[:, 5]
    ux = output[:, 6]
    uy = output[:, 7]
    uz = output[:, 8]

    plt.figure(3)
    plt.title(f'Implanted {ion["symbol"]} Depth Distribution with 50A {target["symbol"]}-H layer')
    plt.xlabel('x [A]')
    plt.ylabel('f(x) [A.U.]')
    heights, _, _ = plt.hist(x[np.logical_and(incident, stopped)], bins=100, density=True, histtype='step')
    plt.plot([50.0, 50.0], [0.0, np.max(heights)*1.1])
    plt.gca().set_ylim([0.0, np.max(heights)*1.1])

    #For smooth distributions and good statistics, you should use at least 10k ions
    number_ions = 100

    #1 keV is above the He on W sputtering threshold of ~150 eV
    energies_eV = 1000.0*np.ones(number_ions)

    #Working with angles of exactly 0 is problematic due to gimbal lock
    angle = 0.0001

    #In RustBCA's 0D geometry, +x -> into the surface
    ux = np.cos(angle*np.pi/180.)*np.ones(number_ions)
    uy = np.sin(angle*np.pi/180.)*np.ones(number_ions)
    uz = np.zeros(number_ions)

    Z1 = np.array([ion['Z']]*number_ions)
    m1 = np.array([ion['m']]*number_ions)
    Ec1 = np.array([ion['Ec']]*number_ions)
    Es1 = np.array([ion['Es']]*number_ions)

    print(f'Running tracked RustBCA for {number_ions} {ion["symbol"]} ions on {target["symbol"]} at {energies_eV[0]/1000.} keV...')
    print(f'This may take several minutes.')

    start = time.time()
    # Note that simple_bca_list_py expects number densities in 1/Angstrom^3
    output, incident, incident_index = compound_bca_list_tracked_py(
        energies_eV, ux, uy, uz, Z1, m1, Ec1, Es1,
        [target['Z'], ion['Z']], [target['m'], ion['m']],
        [target['Ec'], ion['Ec']], [target['Es'], ion['Es']], [target['n']/10**30, target['n']/10**30],
        [target['Eb'], 0.0])
    stop = time.time()
    delta_time = stop - start

    plt.figure()
    plt.plot(incident_index)
    plt.xlabel('Particle number')
    plt.ylabel('Particle index')
    plt.legend(['Incident', 'Indicies'])

    output = np.array(output)
    print('Simulation complete. Processing data...')

    plt.show()

if __name__ == '__main__':
    main()
