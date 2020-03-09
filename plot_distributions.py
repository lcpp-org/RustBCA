import numpy as np
import matplotlib.pyplot as plt

r = np.genfromtxt('1000.0H_Breflected.output', delimiter=',')
s = np.genfromtxt('1000.0H_Bsputtered.output', delimiter=',')

rf = np.genfromtxt('H_B_RFLST.DAT')
sf = np.genfromtxt('H_B_SPLST.DAT')

plt.figure(1)
plt.hist(r[:,2], histtype='step', bins=50)
plt.hist(rf[:,2], histtype='step', bins=50)
plt.legend(['rust', 'ftridyn'])

plt.figure(2)
plt.hist(s[:,2], histtype='step', bins=50)
plt.hist(sf[:,2], histtype='step', bins=50)
plt.legend(['rust', 'ftridyn'])

plt.show()
