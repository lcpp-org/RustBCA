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

#Data digitized from Fig. 1 of https://doi.org/10.1016/0022-3115(84)90433-1
#Modeled reflection coefficients using EAM at 100 eV and below and experimentally at 100 eV and above
#Used to show that TRIM breaks down at higher incident energies
data = np.array(
    [[0.2962771, 0.18217821],
    [0.96344626, 0.4871287],
    [3.0372448, 0.8930693],
    [9.876638, 0.86534655],
    [30.184526, 0.7841584],
    [96.644066, 0.6217822], #note: this is the last EAM data point
    [102.83226, 0.4950495], #note: this is the first experimental data point
    [203.52855, 0.37623763],
    [516.34265, 0.3980198],
    [809.75903, 0.32277226],
    [1006.2257, 0.2990099],
    [2054.3218, 0.17821783],
    [4129.5522, 0.13069306],
    [6890.884, 0.0990099]]
)

#Setting the cutoff low is necessary since the incident energies are also low.
hydrogen['Ec'] = 0.001

#Plotting the EAM data points.
energies = data[:6, 0]
r_benchmark = data[:6, 1]
plt.semilogx(energies, r_benchmark, marker='o', linestyle='', label='EAM')

#Plotting the experimental data points.
energies = data[6:, 0]
r_benchmark = data[6:, 1]
plt.semilogx(energies, r_benchmark, marker='^', linestyle='', label='Exp.')

#Plotting RustBCA data points, using the ergonomic helper function reflection_coefficient().
energies = np.logspace(-1, 4, 25)
r_rustbca = [reflection_coefficient(hydrogen, nickel, energy, 0.0, 1000)[0] for energy in energies]
plt.semilogx(energies, r_rustbca, label='RustBCA', color='black')

r_thomas = [thomas_reflection(hydrogen, nickel, energy) for energy in energies]
plt.semilogx(energies, r_thomas, label='Thomas', color='red', linestyle='--')

plt.title('Reflection Coefficients of Normal Incidence H on Ni')
plt.xlabel('E [eV]')
plt.ylabel('R')
plt.gca().set_ylim([0.0, 1.0])
plt.legend(loc='upper right')

plt.show()

