import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

Q = 1.602E-19
PI = 3.14159
AMU = 1.66E-27
ANGSTROM = 1E-10
MICRON = 1E-6
EPS0 = 8.85E-12
A0 = 0.52918E-10
K = 1.11265E-10
ME = 9.11E-31
SQRTPI = 1.77245385
SQRT2PI = 2.506628274631
C = 300000000.


def do_plots():
    index = 0
    name = str(index)
    reflected = np.atleast_2d(np.genfromtxt(name+'reflected.dat', delimiter=','))
    sputtered = np.atleast_2d(np.genfromtxt(name+'sputtered.dat', delimiter=','))
    deposited = np.atleast_2d(np.genfromtxt(name+'deposited.dat', delimiter=','))
    trajectories = np.atleast_2d(np.genfromtxt(name+'trajectories.dat', delimiter=','))
    trajectory_data = np.genfromtxt(name+'trajectory_data.dat', delimiter=',').transpose().astype(int)

    colors = {
        1: 'red',
        29: 'black',
        2: 'red',
    }

    linewidths = {
        1: 2,
        29: 1,
        2: 2,
    }

    species_names = {
        1: 'Hydrogen',
        2: 'Helium'
    }

    #breakpoint()

    #plt.figure(1)
    #if np.size(reflected) > 0: plt.scatter(reflected[:,1]/ANGSTROM, reflected[:,2]/ANGSTROM, s=10, color='red', marker='x')
    #if np.size(sputtered) > 0: plt.scatter(sputtered[:,1]/ANGSTROM, sputtered[:,2]/ANGSTROM, s=10, color='blue', marker='x')

    fig2, axis2 = plt.subplots()
    n = 8.491E28
    dx =  2.*n**(-1./3.)/SQRT2PI/MICRON
    thickness = 0.5
    depth = 10

    surface_x = [0.0, depth, depth, 0.0, 0.0]
    surface_y = [-thickness/2., -thickness/2, thickness/2, thickness/2, -thickness/2]

    energy_surface_x = [-dx, depth + dx, depth + dx, -dx, -dx]
    energy_surface_y = [-thickness/2 - dx, -thickness/2 - dx, thickness/2 + dx, thickness/2 + dx, -thickness/2 - dx]

    simulation_boundary_x = [-2.*dx, depth + 2.*dx, depth + 2.*dx, -2.*dx, -2.*dx]
    simulation_boundary_y = [-thickness/2 - 2.*dx, -thickness/2 - 2.*dx, thickness/2 + 2.*dx, thickness/2 + 2.*dx, -thickness/2 - 2.*dx]

    plt.plot(surface_x, surface_y, linewidth=2, color='dimgray')
    plt.plot(energy_surface_x, energy_surface_y, '--', linewidth=1, color='dimgray')
    plt.plot(simulation_boundary_x, simulation_boundary_y, '--', linewidth=1, color='dimgray')

    index = 0
    if np.size(trajectories) > 0:
        for trajectory_length in trajectory_data:

            Z = trajectories[index, 0]
            E = trajectories[index:(trajectory_length + index), 1]
            x = trajectories[index:(trajectory_length + index), 2]/MICRON
            y = trajectories[index:(trajectory_length + index), 3]/MICRON
            z = trajectories[index:(trajectory_length + index), 4]/MICRON

            index += trajectory_length

            #points = np.array([x, y]).transpose().reshape(-1, 1, 2)
            #segments = np.concatenate([points[:-1], points[1:]], axis=1)
            #lc = LineCollection(segments, linewidths = 1 + 2*E/E[0], color=colors[Z])
            #axis2.add_collection(lc)
            plt.plot(x, y, color = colors[Z], linewidth = linewidths[Z])
            plt.axis('tight')

    if np.size(sputtered) > 0: plt.scatter(sputtered[:,1]/MICRON, sputtered[:,2]/MICRON, s=50, color='blue', marker='*')
    plt.xlabel('x [um]')
    plt.ylabel('y [um]')
    plt.title('Hydrogen trajectories in copper')

    if np.size(trajectories) > 1:
        marker_color = [colors[Z] for Z in trajectories[:, 0]]
        plt.scatter(trajectories[:, 1]*1E10, trajectories[:, 2]*1E10, s=1, color=marker_color)
    plt.axis('square')

    if np.size(deposited) > 0:
        plt.figure(2)
        bins = np.linspace(-thickness/2., thickness/2., 50)
        plt.hist(deposited[:,2]/MICRON, bins=bins)
        plt.title('Helium y-Deposition at 1 MeV on 1X10 um Target')
        plt.xlabel('y [um]')
        plt.ylabel('Hydrogen Deposition (A.U.)')

        plt.figure(3)
        binx = np.linspace(0., depth, 120)
        biny = np.linspace(-thickness/2., thickness/2., 120)
        plt.hist2d(deposited[:, 1]*1E6, deposited[:, 2]*1E6, bins=((binx, biny)))
        plt.axis([0.0, 10.0, -0.5, 0.5])
        plt.title('Helium Deposition at 1 MeV on 1X10 um Target')
        plt.xlabel('x [um]')
        plt.ylabel('y [um]')
    plt.show()

def SBB(E, Za, Zb, Ma, n):
    v = np.sqrt(2.*E/Ma)
    beta = v/C

    if Zb < 13:
        I0 = (12. + 7./Zb)
    else:
        I0 = (9.76 + 58.5*Zb**(-1.19))

    I = Zb*I0*Q

    if Za < 3:
        B = 100.*Za/Zb
    else:
        B = 5

    prefactor = 8.0735880e-42 * Zb * Za**2 / beta**2
    eb = 2*ME*v*v/(I)

    return prefactor*np.log(eb + 1 + B/eb)*n

def SLS(E, Za, Zb, Ma, n):
    return 1.212*(Za**(7./6.)*Zb)/((Za**(2./3.) + Zb**(2./3.))**(3./2.))*np.sqrt(E/Ma*AMU/Q)*ANGSTROM**2*Q*n

def test_stopping():
    energies = np.logspace(3, 10, 100)*Q
    Za = 1
    Zb = 29
    Ma = 1*AMU
    n = 8.49E28
    Mb = 64*AMU

    S_BB = np.array([SBB(energy, Za, Zb, Ma, n) for energy in energies])*6.262E10/8.96
    S_LS = np.array([SLS(energy, Za, Zb, Ma, n) for energy in energies])*6.262E10/8.96

    plt.loglog(energies/Q/1E6, S_LS)
    plt.loglog(energies/Q/1E6, S_BB)
    plt.loglog(energies/Q/1E6, 1./(1./S_LS + 1./S_BB))
    plt.show()

def main():
    do_plots()

def S(E):
    e = 1.602E-19
    amu = 1.66e-27
    Za = 1
    Zb = 29

    N2 = 28
    Ma = 1*amu
    v = np.sqrt(2*E/Ma)
    hbar = 1.055e-34
    I02 = e*10.*Zb


    I01 = e*10.*Za
    me = 9.11E-31
    eps0 = 8.854E-12

    v0 = 2.19e6 #m/s
    a0 = 5.291E-11 #m

    b = (8./np.pi)**(2./3.)
    y = v/(v0*Za**(2./3.))
    a = b**2/0.60647*y**2
    Q0 = (a/2. + (b/3.)**3 + np.sqrt((a/2.)**2 + a*(b/3.)**3))**(1./3.)
    xc = -2.*(b/3.) + 1./Q0*(b/3.)**2 + Q0**(1./3.)

    c = 0.969376
    l = 5./7.

    N1 = Za*(1. - (b**2*(3.*xc + b))/(xc + b)**3.)

    L1 = c*a0/(Za**(1./3.)*b*(1. - (l/5.)*(N1/Za))) * (N1/Za)**(2./3.)

    alpha = (Zb - N2)**2*N1/(Za - N1)**2/N2

    epsilon = I02**(1./(1. + alpha))*I01**(alpha/(1. + alpha))

    qmax = 2.*Ma*v/hbar
    qmin = epsilon/hbar/v

    i_f = 1 - N1/Za
    K = 4.*np.pi*Zb*Za*Za/(me*v*v)*(e*e/4./np.pi/eps0)**2
    Q1 = ((qmax*L1)**2 + 0.27819)/((qmin*L1)**2 + 0.27819)
    Q2 = ((qmax*L1)**2 + 3.99187)/((qmin*L1)**2 + 3.99187)
    Q3 = 1./((qmax*L1)**2 + 0.27819) - 1./((qmin*L1)**2 + 0.27819)
    Q4 = 1./((qmax*L1)**2 + 3.99187) - 1./((qmin*L1)**2 + 3.99187)

    Se = K*N1*Za**2*(i_f**2*np.log(qmax/qmin) + (1. - i_f)*(0.448685*(i_f + 0.402031)*np.log(Q1) + 0.0513151*(i_f + 6.22848)*np.log(Q2)) + (1. - i_f)**2*(0.0550439*(Q3) + 0.274619*(Q4)))

if __name__ == '__main__':
    main()
