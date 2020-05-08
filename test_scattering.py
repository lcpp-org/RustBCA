import matplotlib.pyplot as plt
import numpy as np

PI = np.pi
Q = 1.602E-19
EV = Q
AMU = 1.66E-27
ANGSTROM = 1E-10
MICRON = 1E-6
NM = 1E-9
CM = 1E-2
EPS0 = 8.85E-12
A0 = 0.52918E-10
K = 1.11265E-10
ME = 9.11E-31
SQRTPI = 1.772453850906
SQRT2PI = 2.506628274631
C = 299792000.
BETHE_BLOCH_PREFACTOR = 4.*PI*(Q*Q/(4.*PI*EPS0))*(Q*Q/(4.*PI*EPS0))/ME/C/C
LINDHARD_SCHARFF_PREFACTOR = 1.212*ANGSTROM*ANGSTROM*Q
LINDHARD_REDUCED_ENERGY_PREFACTOR = 4.*PI*EPS0/Q/Q

def phi(xi):
    return 0.190945*np.exp(-0.278544*xi) + 0.473674*np.exp(-0.637174*xi) + 0.335381*np.exp(-1.919249*xi)

def dphi(xi):
    return -0.278544*0.190945*np.exp(-0.278544*xi) - 0.637174*0.473674*np.exp(-0.637174*xi) - 0.335381*1.919249*np.exp(-1.919249*xi)

def doca_function(x0, beta, reduced_energy):
    return x0 - phi(x0)/reduced_energy - beta**2/x0

def diff_doca_function(x0, beta, reduced_energy):
    return beta**2/x0**2 - dphi(x0)/reduced_energy + 1.

def f(x, beta, reduced_energy):
    return (1. - phi(x)/x/reduced_energy - beta**2/x**2)**(-0.5)

def screening_length(Za, Zb):
    return 0.8853*A0*(np.sqrt(Za) + np.sqrt(Zb))**(-2./3.)

def distance_of_closest_approach(Za, Zb, Ma, Mb, E, impact_parameter, max_iter=100, tol=1E-12):
    mu = Mb/(Ma + Mb)
    a = screening_length(Za, Zb)
    reduced_energy = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*mu/Za/Zb*E
    beta = impact_parameter/a

    x0 = 1.
    if reduced_energy > 5:
        inv_er_2 = 0.5/reduced_energy
        x0 = inv_er_2 + np.sqrt(inv_er_2**2 + beta**2)

    error = 1
    for _ in range(max_iter):
        xn = x0 - doca_function(x0, beta, reduced_energy)/diff_doca_function(x0, beta, reduced_energy)
        err = (xn - x0)**2
        x0 = xn
        if err < tol:
            break

    return xn

def magic(Za, Zb, Ma, Mb, E, impact_parameter, x0, C_):
    c = [0.190945, 0.473674, 0.335381, 0.0 ]
    d = [-0.278544, -0.637174, -1.919249, 0.0 ]
    a_ = 0.8853*A0/(np.sqrt(Za) + np.sqrt(Zb))**(2./3.)
    V0 = Za*Zb*Q*Q/4.0/PI/EPS0/a_
    E_c = E*Mb/(Ma + Mb)
    E_r = E/V0
    b = impact_parameter/a_
    SQE = np.sqrt(E_r)
    R = a_*x0
    sum = c[0]*np.exp(d[0]*x0) + c[1]*np.exp(d[1]*x0)+ c[2]*np.exp(d[2]*x0)
    V = V0*a_/R*sum
    sum = d[0]*c[0]*np.exp(d[0]*x0) + d[1]*c[1]*np.exp(d[1]*x0) + d[2]*c[2]*np.exp(d[2]*x0)
    dV = -V/R + V0/R*sum
    rho = -2.0*(E_c - V)/dV
    D = 2.0*(1.0+C_[0]/SQE)*E_r*b**((C_[1]+SQE)/(C_[2]+SQE))
    G = (C_[4]+E_r)/(C_[3]+E_r)*(np.sqrt(1.0+D*D) - D)
    delta =  D*G/(1.0+G)*(x0-b)
    ctheta2 = (b + rho/a_ + delta)/(x0 + rho/a_)
    return 2.*np.arccos((b + rho/a_ + delta)/(x0 + rho/a_))

def scattering(Za, Zb, Ma, Mb, E, impact_parameter, xn):
    mu = Mb/(Ma + Mb)
    a = screening_length(Za, Zb)
    reduced_energy = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*mu/Za/Zb*E
    beta = impact_parameter/a

    lambda_0 = (0.5 + beta*beta/xn/xn/2. - dphi(xn)/2./reduced_energy)**(-1./2.)
    alpha = 1./12.*(1. + lambda_0 + 5.*(0.4206*f(xn/0.9072, beta, reduced_energy) + 0.9072*f(xn/0.4206, beta, reduced_energy)))
    theta = PI*(1. - beta*alpha/xn)

    return theta

def main():
    Za = 1
    Zb = 29
    Ma = 1
    Mb = 63.54
    impact_parameter = 1E-10
    N = 1000
    n = 8E28

    mfp = n**(-1./3.)
    pmax = mfp/SQRTPI
    a = screening_length(Za, Zb)
    mu = Mb/(Ma + Mb)
    betas = [35, 30, 25, 20, 15, 10, 7.5, 5, 2.5, 0.]
    impact_parameters = a*np.array(betas)
    reduced_energies = np.logspace(-5, -1, N)
    energies = reduced_energies/(LINDHARD_REDUCED_ENERGY_PREFACTOR*a*mu/Za/Zb)

    plt.figure(1)
    theta = np.zeros(N)
    beta = np.zeros(N)
    xn = np.zeros(N)
    for i, impact_parameter in enumerate(impact_parameters):
        for j, energy in enumerate(energies):

            xn[j] = distance_of_closest_approach(Za, Zb, Ma, Mb, energy, impact_parameter)
            #theta[index] = scattering(Za, Zb, Ma, Mb, E, impact_parameter, xn)

        plt.semilogx(reduced_energies, xn)
    plt.legend(betas)

    plt.figure(2)

    theta = np.zeros(N)
    theta_magic_1 = np.zeros(N)
    theta_magic_2 = np.zeros(N)
    C1 = [1.0144, 0.235800, 0.126, 63950.0, 83550.0]
    C2 = [0.7887, 0.01166, 0.006913, 17.16, 10.79]
    labels = []
    pmax = mfp/np.sqrt(np.pi*2.)

    for j, reduced_energy in enumerate([1E-1, 1E-2, 1E-3, 1E-4, 1E-5]):
        for i, impact_parameter in enumerate(np.linspace(0., pmax, N)):
            energy = reduced_energy/(LINDHARD_REDUCED_ENERGY_PREFACTOR*a*mu/Za/Zb)
            xn = distance_of_closest_approach(Za, Zb, Ma, Mb, energy, impact_parameter)
            theta[i] = scattering(Za, Zb, Ma, Mb, energy, impact_parameter, xn)
            theta_magic_1[i] = magic(Za, Zb, Ma, Mb, energy, impact_parameter, xn, C1)
            theta_magic_2[i] = magic(Za, Zb, Ma, Mb, energy, impact_parameter, xn, C2)
        handle = plt.plot(np.linspace(0., pmax/mfp, N), theta*180./np.pi)
        plt.plot(np.linspace(0., pmax/mfp, N), theta_magic_1*180./np.pi, '--', color=handle[0].get_color())
        plt.plot(np.linspace(0., pmax/mfp, N), theta_magic_2*180./np.pi, ':', color=handle[0].get_color())
        labels.append('Lobatto 6th Order, E_r='+str(reduced_energy))
        labels.append('MAGIC F-TRIDYN, E_r='+str(reduced_energy))
        labels.append('MAGIC SRIM, E_r='+str(reduced_energy))
    plt.legend(labels)
    plt.xlabel('Impact Parameter in mfp')
    plt.ylabel('theta [deg]')

    plt.figure(3)
    theta = np.zeros(N)
    psi = np.zeros(N)
    pmax = mfp/np.sqrt(np.pi*2.)
    for j, reduced_energy in enumerate([1E-1, 1E-2, 1E-3, 1E-4, 1E-5]):
        for i, impact_parameter in enumerate(np.linspace(0., pmax, N)):
            energy = reduced_energy/(LINDHARD_REDUCED_ENERGY_PREFACTOR*a*mu/Za/Zb)
            xn = distance_of_closest_approach(Za, Zb, Ma, Mb, energy, impact_parameter)
            theta[i] = scattering(Za, Zb, Ma, Mb, energy, impact_parameter, xn)
            psi[i] = np.arctan2(np.sin(theta[i]), Ma/Mb + np.cos(theta[i]))*180./np.pi
        plt.plot(np.linspace(0., pmax, N)/mfp, psi)
    plt.legend(['Reduced E=1E-1', '1E-2', '1E-3', '1E-4', '1E-5'])
    plt.yticks([0, 30, 60, 90, 120, 150, 180])
    plt.xlabel('Impact Parameter in MFP')
    plt.ylabel('Lab Angle Deflection from Initial Direction [deg]')

    plt.show()

if __name__ == '__main__':
    main()
