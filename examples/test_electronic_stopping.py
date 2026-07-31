from libRustBCA import electronic_stopping_cross_sections
import numpy as np
import matplotlib.pyplot as plt

Ma = 1.008
Za = 1.0
Zb = 13.0
n_Al = 6.026e28

# load PSTAR data for H on Al
# PSTAR data is in MeV cm^2/g
# to convert to eV/m, multiply by 2.7 [g/cm^3] * 100 [cm/m] * 1E-6 [eV/MeV]
PSTAR = np.genfromtxt('examples/H_Al_PSTAR.dat')
energies = PSTAR[:, 0]*1e6
S_PSTAR = PSTAR[:, 1]*270.*1e6
plt.loglog(energies, S_PSTAR, label='PSTAR')

num_energies = len(energies)

lindhard_scharff = np.zeros(num_energies)
bethe_bloch = np.zeros(num_energies)
biersack_varelas = np.zeros(num_energies)

for index, energy in enumerate(energies):
    lindhard_scharff[index], bethe_bloch[index], biersack_varelas[index] = electronic_stopping_cross_sections(Za, Zb, energy, Ma, 1.4)

plt.loglog(energies, lindhard_scharff*n_Al, label='Lindhard-Scharff (ck=1.4 [Eckstein 1991 Tab. 5.1])')
plt.loglog(energies, bethe_bloch*n_Al, label='Bethe-Bloch (ck=1.4)')
plt.loglog(energies, biersack_varelas*n_Al, label='Biersack-Varelas (ck=1.4)')

plt.legend(loc='lower left')
plt.title('Nonlocal Electronic Stopping Modes in RustBCA (H -> Al)')
plt.xlabel('E [eV]')
plt.ylabel('S [eV/m]')
plt.gca().set_ylim([3e8, 1e12])
plt.savefig('electronic_stopping_cross_sections.png')

plt.figure()
index_valid = np.argmax(energies > 1e9)
residual = np.sqrt((biersack_varelas[:index_valid]*n_Al - S_PSTAR[:index_valid])**2)/S_PSTAR[:index_valid]
# Biersack-Varelas interpolation misses the peak with ~19.9% error, but has ~1% error in LS region and max 3% error in B-B region below 1e8 ev
np.testing.assert_allclose(residual, 0.0, atol=0.2)


plt.loglog(energies[:index_valid], residual)
plt.title('Biersack-Varelas residual against PSTAR (H->Al, ck=1.4 [Eckstein 1991 Tab. 5.1])')
plt.xlabel('E [eV]')
plt.ylabel('Error in BV vs PSTAR')
plt.savefig('biersack_varelas_residual.png')