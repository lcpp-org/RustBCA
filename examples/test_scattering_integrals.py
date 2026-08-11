import numpy as np
import matplotlib.pyplot as plt
import sys, os

from libRustBCA import *
#This should allow the script to find materials and formulas from anywhere
sys.path.append(os.path.dirname(__file__)+'/../scripts')
sys.path.append('scripts')
from materials import *
from formulas import *

energies = np.logspace(0, 4, 4)
impact_parameters = np.logspace(-3, 2, 100)

ion = helium
target = boron

Za = ion['Z']
Zb = target['Z']
Ma = ion['m']
Mb = target['m']

show_plots = True

linestyles = ['-', '--', ':', '-.']

for potential in ["KR_C", "MOLIERE", "ZBL"]:
    plt.figure()
    plt.title(f'Scattering Angles for {potential}')
    for linestyle, energy in zip(linestyles, energies):
        gm = np.zeros(100)
        gl = np.zeros(100)
        mw = np.zeros(100)
        magic = np.zeros(100)
        for index, p in enumerate(impact_parameters):
            gm[index], gl[index], mw[index], magic[index] = scattering_integrals(Za, Zb, Ma, Mb, energy, p, interaction_potential=potential)

        plt.semilogx(impact_parameters, gm, label=f'Gauss-Mehler, E={np.round(energy/1000, 3)} keV', linestyle=linestyle)
        plt.semilogx(impact_parameters, gl, label=f'Gauss-Legendre, E={np.round(energy/1000, 3)} keV', linestyle=linestyle)
        plt.semilogx(impact_parameters, mw, label=f'Mendenhall-Weller, E={np.round(energy/1000, 3)} keV', linestyle=linestyle)
        plt.semilogx(impact_parameters, magic, label=f'MAGIC, E={np.round(energy/1000, 3)} keV', linestyle=linestyle)
        plt.gca().set_prop_cycle(None)

        np.testing.assert_allclose(gm, gl, atol=5e-3) # 0.5% seems reasonable? max is ~0.3%
        np.testing.assert_allclose(gm, mw, atol=5e-3)
        np.testing.assert_allclose(mw, gl, atol=5e-3)
        plt.legend()
        plt.xlabel('p [A]')
        plt.ylabel('theta [rad]')

if show_plots: plt.show()