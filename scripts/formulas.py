from scipy import constants
from scipy.interpolate import interp1d
import numpy as np

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

def thomas_reflection(ion, target, energy_eV):
    '''
    Thomas et al. (1992) semi-empirical reflection coefficient.

    https://doi.org/10.1016/0168-583X(92)95298-6

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

    https://doi.org/10.1080/10420159408221042

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

    #Thomas-Fermi reduced energy
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

def bohdansky_heavy_ion(ion, target, energy_eV):
    '''
    Bohdansky sputtering yield formula in the heavy ion (M1/M2 > 0.5) regime.
    Returns None if the target does not have a surface binding energy.

    https://doi.org/10.1063/1.327954
    https://doi.org/10.1016/0168-583X(84)90271-4

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

    if Us == 0.: return None
    alpha = 0.3*(m2/m1)**(2./3.)

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

    sputtering_yield = alpha*Sn*(1 - (Eth/energy_eV)**(2./3.))*(1. - Eth/energy_eV)**2

    return sputtering_yield

def bohdansky_light_ion(ion, target, energy_eV):
    '''
    Bohdansky sputtering yield formula in the light ion (M1/M2 < 0.5) limit.
    Returns None if the target does not have a surface binding energy.

    https://doi.org/10.1063/1.327954
    https://doi.org/10.1016/0168-583X(84)90271-4

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

    if Us == 0.: return None
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

    https://doi.org/10.1080/01422448208226913

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
