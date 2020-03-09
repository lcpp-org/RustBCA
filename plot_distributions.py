import numpy as np
import matplotlib.pyplot as plt

file_num = str(1)

r = np.genfromtxt('1000.0H_Breflected.output', delimiter=',')
s = np.genfromtxt('1000.0H_Bsputtered.output', delimiter=',')

rf = np.genfromtxt('H_B_RFLST.DAT')
sf = np.genfromtxt('H_B_SPLST.DAT')

plt.figure(1)
plt.hist(r[:,2], histtype='step', bins=50, density=True)
plt.hist(rf[:,2], histtype='step', bins=50, density=True)
plt.legend(['rust', 'ftridyn'])
plt.title('Reflected Energy Distributions')

plt.figure(2)
plt.hist(s[:,2], histtype='step', bins=50, density=True)
plt.hist(sf[:,2], histtype='step', bins=50, density=True)
plt.legend(['rust', 'ftridyn'])
plt.title('Sputtered Energy Distributions')

plt.figure(3)
plt.hist(-s[:,6], histtype='step', bins=50, density=True)
plt.hist(sf[:,6], histtype='step', bins=50, density=True)
plt.hist(s[:,7], histtype='step', bins=50, density=True)
plt.hist(sf[:,7], histtype='step', bins=50, density=True)
plt.hist(s[:,8], histtype='step', bins=50, density=True)
plt.hist(sf[:,8], histtype='step', bins=50, density=True)
plt.legend(['rust cosx', 'ftridyn cosx', 'rust cosy', 'ftridyn cosy', 'rust cosz', 'ftridyn cosz'])
plt.title('Sputtered Angular Distributions')

plt.figure(4)
plt.hist(-r[:,6], histtype='step', bins=50, density=True)
plt.hist(rf[:,6], histtype='step', bins=50, density=True)
plt.hist(r[:,7], histtype='step', bins=50, density=True)
plt.hist(rf[:,7], histtype='step', bins=50, density=True)
plt.hist(r[:,8], histtype='step', bins=50, density=True)
plt.hist(rf[:,8], histtype='step', bins=50, density=True)
plt.legend(['rust cosx', 'ftridyn cosx', 'rust cosy', 'ftridyn cosy', 'rust cosz', 'ftridyn cosz'])
plt.title('Reflected Angular Distributions')

plt.show()
