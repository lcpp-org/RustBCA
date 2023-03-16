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

R_T_Be = np.array(
    [
    [9.67742,	0.016754618],
    [20.32258,	0.09511873],
    [30.0,	0.13153034],
    [49.677418,	0.102242745],
    [74.83871,	0.08403694],
    [100.32258,	0.031002639],
    ]
)

R_D_Be = np.array(
    [
    [10.32258,	0.11649077],
    [20.32258,	0.19168866],
    [30.32258,	0.19722955],
    [49.677418,	0.17269129],
    [75.16129,	0.11174142],
    [100.0,	0.044459105],
    ]
)

R_H_Be = np.array(
    [
    [9.67742,	0.19485489],
    [19.67742,	0.32150397],
    [30.0,	0.33575198],
    [50.322582,	0.25105542],
    [75.16129,	0.1378628],
    [100.645164,	0.044459105],
    ]
)

#Plotting the MD data points.
energies = R_H_Be[:, 0]
R = R_H_Be[:, 1]
plt.plot(energies, R, marker='o', linestyle='', label='H on Be MD+NNP')

#Plotting the MD data points.
energies = R_D_Be[:, 0]
R = R_D_Be[:, 1]
plt.plot(energies, R, marker='x', linestyle='', label='D on Be MD+NNP')

#Plotting the MD data points.
energies = R_T_Be[:, 0]
R = R_T_Be[:, 1]
plt.plot(energies, R, marker='^', linestyle='', label='T on Be MD+NNP')


for ion in [hydrogen, deuterium, tritium]:
    energies = np.logspace(1, 2, 25)
    r_rustbca = [reflection_coefficient(ion, beryllium, energy, 0.0, 10000)[0] for energy in energies]
    r_thomas = [thomas_reflection(ion, beryllium, energy) for energy in energies]
    line = plt.plot(energies, r_rustbca, label=f'{ion["symbol"]} on Be RustBCA')[0]
    plt.plot(energies, r_thomas, label=f'{ion["symbol"]} on Be Thomas', linestyle='--', color=line.get_color())

plt.xlabel('E [eV]')
plt.ylabel('R_N')
plt.title('Reflection Coefficients')
plt.legend()
plt.savefig('reflection_md_nnp_rustbca.png')

#Data digitized from https://doi.org/10.1088/1741-4326/ac592a
#MD+NNP Sputtering Yields
y_H_Be = np.array([
    [20.004726, 0.0029361406],
    [29.880848, 0.013341782],
    [49.791046, 0.024575423],
    [75.00529,	0.016207783],
    [99.68844,	0.017530797],
])

y_D_Be = np.array([
    [19.9699, 0.005002439],
    [29.601622, 0.021064475],
    [49.76866,	0.03461476],
    [74.71549,	0.030082114],
    [99.954605, 0.013559438],
])

y_T_Be = np.array([
    [20.013433, 0.0025699555],
    [29.85846,	0.01879205],
    [49.756844, 0.04147379],
    [74.707405, 0.034042966],
    [99.68098, 0.01965117],
])


plt.figure()

#Plotting the MD data points.
energies = y_H_Be[:, 0]
y = y_H_Be[:, 1]
plt.semilogy(energies, y, marker='o', linestyle='', label='H on Be MD+NNP')

#Plotting the MD data points.
energies = y_D_Be[:, 0]
y = y_D_Be[:, 1]
plt.semilogy(energies, y, marker='x', linestyle='', label='D on Be MD+NNP')

#Plotting the MD data points.
energies = y_T_Be[:, 0]
y = y_T_Be[:, 1]
plt.semilogy(energies, y, marker='^', linestyle='', label='T on Be MD+NNP')

for ion in [hydrogen, deuterium, tritium]:
    energies = np.logspace(1, 2, 25)
    y_rustbca = [sputtering_yield(ion, beryllium, energy, 0.0, 100000) for energy in energies]
    plt.semilogy(energies, y_rustbca, label=f'{ion["symbol"]} on Be RustBCA')

yamamura_yield = [yamamura(hydrogen, beryllium, energy) for energy in energies]
plt.semilogy(energies, yamamura_yield, label='Yamamura H on Be', linestyle='--')

plt.title('Sputtering Yields of Hydrogenic Species on Be')
plt.xlabel('E [eV]')
plt.ylabel('Y [at/ion]')
plt.gca().set_ylim([1e-4, 1e-1])
plt.legend(loc='lower right')

plt.savefig('sputtering_md_nnp_rustbca.png')