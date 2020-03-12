import numpy as np
import matplotlib.pyplot as plt


def plot_distributions(name, ftridyn_name):
    file_num = str(1)
    num_bins = 100

    r = np.genfromtxt(name+'reflected.output', delimiter=',')
    s = np.genfromtxt(name+'sputtered.output', delimiter=',')

    rf = np.genfromtxt(ftridyn_name+'RFLST.DAT')
    sf = np.genfromtxt(ftridyn_time+'SPLST.DAT')

    plt.figure(1)
    bins = np.linspace(0.0, 1000.0, num_bins)
    plt.hist(r[:,2], histtype='step', bins=bins, density=True)
    plt.hist(rf[:,2], histtype='step', bins=bins, density=True)
    plt.legend(['rust', 'ftridyn'])
    plt.title(name+' Reflected Energy Distributions')
    plt.xlabel('E [eV]')

    plt.figure(2)
    bins=np.linspace(0.0, 100.0, num_bins)
    plt.hist(s[:,2], histtype='step', bins=bins, density=True)
    plt.hist(sf[:,2], histtype='step', bins=bins, density=True)
    plt.legend(['rust', 'ftridyn'])
    plt.title(name+' Sputtered Energy Distributions')
    plt.xlabel('E [eV]')

    plt.figure(3)
    bins = np.linspace(0., 90., num_bins)
    plt.hist(180./np.pi*np.arccos(-s[:,6]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(sf[:,6]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(s[:,7]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(sf[:,7]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(s[:,8]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(sf[:,8]), histtype='step', bins=bins, density=True)
    plt.legend(['rust cosx', 'ftridyn cosx', 'rust cosy', 'ftridyn cosy', 'rust cosz', 'ftridyn cosz'], loc='upper left')
    plt.title(name+' Sputtered Angular Distributions')
    plt.xlabel('angle [deg]')

    plt.figure(4)
    bins = np.linspace(0., 90., num_bins)
    plt.hist(180./np.pi*np.arccos(-r[:,6]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(rf[:,6]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(r[:,7]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(rf[:,7]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(r[:,8]), histtype='step', bins=bins, density=True)
    plt.hist(180./np.pi*np.arccos(rf[:,8]), histtype='step', bins=bins, density=True)
    plt.legend(['rust cosx', 'ftridyn cosx', 'rust cosy', 'ftridyn cosy', 'rust cosz', 'ftridyn cosz'], loc='upper left')
    plt.title(name+' Reflected Angular Distributions')
    plt.xlabel('angle [deg]')
