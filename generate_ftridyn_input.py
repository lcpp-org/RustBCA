"""
Provides classes and functions to generate F-TRIDYN input files.
"""

from __future__ import print_function
import numpy as np
#import matplotlib.pyplot as plt
import os
import pickle
import matplotlib.pyplot as plt
import itertools
import matplotlib
import matplotlib.patheffects as path_effects
from scipy.interpolate import interp2d

font = {
    'family': 'sans-serif',
    'weight': 'normal',
    'size': 16
}

title_font = {
    'family': 'sans-serif',
    'weight': 'bold',
    'size': 18
}

matplotlib.rc('font', **font)

class tridyn_interface:
    """
        Experimental class to serve as an on-line interface to F-TRIDYN (run in
        TRIDYN-like mode, i.e., FD=1.0).

        Args:
            beam_species (:obj: 'str'): Chemical symbol of beam species
            target_species (:obj: 'str'): Chemical symbol of target species
    """
    def __init__(self, beam_species, target_species):
        self.beam_species = beam_species
        self.target_species = target_species
        #Create species lookup table and use it to find species from symbols given
        self.lookup_table = species_lookup_table()
        #Beam species is assumed to be only incident species
        self.beam = self.lookup_table.find_species(self.beam_species, QU=0.0, QUBEAM=1.0)
        #Target is homogeneous and formed entirely from target species
        self.target = self.lookup_table.find_species(self.target_species, QU=1.0, QUBEAM=0.0, QUMAX=1.0)
        self.species_list = [self.beam, self.target]
        self.number_species = len(self.species_list)
        #Surface file forced to be flat (FD=1.0). Width is arbitrary given periodic BCs in F-TRIDYN
        self.fractal_surface = fractal_surface(FD=1.0, width=400.0, iterations=1)
        self.surface_name = self.fractal_surface.print_fractal_surface_file()
        #Form 4-charcater name from truncating symbols and filling with underscores
        self.name = beam_species.ljust(2, '_')[:2] + target_species.ljust(2, '_')[:2]
    #end def __init__

    def run_tridyn_simulations_from_iead(self, energies, angles, iead, number_histories=1e3, depth=400.0):
        """
            Given a list of energies, angles, and an iead, run a TRIDYN-like
            F-TRIDYN simulation for every nonzero energy-angle pair in the iead.
            Returns a list of ALL particles emitted by the material in response
            to the plasma.

            Args:
                energies (:obj: 'numpy.ndarray'): list of energies in iead
                angles (:obj: 'numpy.ndarray'): list of angles in iead
                iead (:obj: 'numpy.ndarray'): ion energy-angle distribution
                number_histories (int, optional): number of BCA psuedoparticles.
                    This will be scaled by the count in the iead.
                depth (float, optional): depth of surface in F-TRIDYN simulation
                    (may need to adjust upwards from default if high energy ions
                    are present)

            Returns:
                new_particle_list_sputtered (:obj: list): list of energy, ca,
                    cb, cg, and mass in atomic units (eV for energy and u for
                    mass)
                new_particle_list_reflected (:obj: list): same as above for
                    reflected particles
        """
        #Initialize new_particle_list as empty. I'm going to use append here
        #because there's no way to estimate a priori the final length.
        new_particle_list_sputtered = []
        new_particle_list_reflected = []

        sim_number = 0
        num_deposited = 0
        num_sputtered = 0
        num_reflected = 0

        #Pull out each row of the iead, and then pull out each element from the
        #rows. Using enumerate allows us to map indices to the appropriate
        #energy and angle values.

        total_sim_number = np.sum(iead > 0)
        for energy_index, row in enumerate(iead):
            for angle_index, count in enumerate(row):
                #Skip the energy-angle pair if there are no incident particles
                if count > 0:
                    #Pull energy and angle value from approriate arrays
                    energy = energies[energy_index]
                    angle = angles[angle_index]

                    #Assign the energy and angle values to the incident species
                    self.species_list[0].E0 = energy
                    self.species_list[0].ALPHA0 = angle
                    self.species_list[0].BE = 0.
                    self.species_list[1].BE = 0.

                    #Initialize sim_params object. Each step it's overwritten.
                    simulation_parameters = sim_params(name=self.name, sim_num=sim_number,
                        IFOUT=count*number_histories//20, NH=count*number_histories,
                        IDOUT=count*number_histories//5, IQOUT=int(count)*number_histories//5,
                        NCP=self.number_species, IDREL=0, IQ0=0, IRC0=0, IRAND=np.random.randint(1, 12855897),
                        JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1e-16, INEL=1, IWC=3, IDIFF=1,
                        CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth,
                        NQX=250, DSF=100.0, IQXN=0, IQXX=250, IMCP=0,
                        surf_name=self.surface_name, species_list=self.species_list,
                        output_pka=0, output_ska=0, output_prj=0, output_splst=1,
                        output_rflst=1, output_prof=0, output_simple=1)

                    #Print simulation input file
                    simulation_parameters.print_input_file()

                    #Run F-TRIDYN using just-created input file
                    os.system('./FTridyn_Clean <'+str(self.name)+str(sim_number).zfill(4)+'.IN > /dev/null')
                    print('simulation number: '+str(sim_number)+' / '+str(total_sim_number))
                    simple_output = np.genfromtxt(str(self.name)+str(sim_number).zfill(4)+'.OUT',
                        skip_footer = 2, usecols = (2))
                    num_prj = simple_output[0]
                    num_spt = simple_output[4]
                    num_bks = simple_output[1]
                    num_dep = num_prj - num_bks

                    #Read sputtered and reflected particle lists and force 2D
                    sputtered_list = np.atleast_2d(np.genfromtxt(str(self.name)+'SPLST.DAT'))
                    reflected_list = np.atleast_2d(np.genfromtxt(str(self.name)+'RFLST.DAT'))

                    #Energy angle coordinates for new list are:
                    #Energy: 2, cos(alpha), cos(beta), cos(gamma): 6:8, mass

                    #If the list is not empty, continue
                    if np.size(sputtered_list) > 0:
                        #For each row, pull out energy angle coordinates
                        for row in sputtered_list:
                            mass = self.species_list[int(row[1]-1)].M
                            Z = int(self.species_list[int(row[1]-1)].ZZ)
                            new_particle_list_sputtered.append([row[2], row[6],
                                row[7], row[8], mass, Z])
                        #end for
                    #end if

                    if np.size(reflected_list) > 0:
                        for row in reflected_list:
                            mass = self.species_list[int(row[1]-1)].M
                            Z = int(self.species_list[int(row[1]-1)].ZZ)
                            new_particle_list_reflected.append([row[2], row[6],
                                row[7], row[8], mass, Z])
                        #end for
                    #end if

                    #TODO implement implanted/deposited projectile list
                    sim_number += 1
                    num_deposited += num_dep
                    num_sputtered += num_spt
                    num_reflected += num_bks
                #end if
            #end for
        #end for
        return num_spt, num_bks
    #end def run_tridyn_simulations_from_iead
#end class tridyn_interface

class tridyn_lookup:
    """
    Experimental class-based lookup table for F-TRIDYN results. This version is
    intended for use with hPIC and is not very flexible yet. Functionality is
    provided to generate M*N F-TRIDYN simulations over a range of M energies and
    N angles, and, given an arbitrary energy/angle pair, return all F-TRIDYN
    output parameters and files. Assumes evenly spaced energies and angles.
    """
    def __init__(self, energies, angles, beam_species, target_species):
        self.energies = energies
        self.angles = angles
        self.de = energies[1] - energies[0]
        self.da = angles[1] - angles[0]
        self.beam_species = beam_species
        self.target_species = target_species
        self.sputtering_yield = np.zeros((len(energies), len(angles)))
        self.reflection_yield = np.zeros((len(energies), len(angles)))
        self.sim_numbers = np.zeros((len(energies), len(angles)))
        self.sim_number_to_energy_angle_tuple = [None]*len(energies)*len(angles)
        self.sputtered_lists = [None]*len(energies)*len(angles)
        self.reflected_lists = [None]*len(energies)*len(angles)
        self.implanted_lists = [None]*len(energies)*len(angles)
    #end def __init__

    def fill_lookup_table(self, name, depth=400.0):
        number_histories = 1e4
        width = 400.0
        number_layers = 250

        surface = fractal_surface(FD=1.0, width=width, iterations=1)
        surface.print_fractal_surface_file()

        species_lookup_table_ = species_lookup_table()
        beam_species = species_lookup_table_.find_species(self.beam_species,
            QU=0.0, QUBEAM=1.0, QUMAX=0.0)
        target_species = species_lookup_table_.find_species(self.target_species,
            QU=1.0, QUBEAM=0.0, QUMAX=1.0)
        species_list = [beam_species, target_species]
        number_species = len(species_list)

        sim_number = 0
        for energy_index, energy in enumerate(self.energies):
            for angle_index, angle in enumerate(self.angles):
                self.sim_number_to_energy_angle_tuple[sim_number] = (energy,angle)

                species_list[0].E0 = energy
                species_list[0].ALPHA0 = angle

                simulation_parameters = sim_params(name=name, sim_num=sim_number,
                    IFOUT=number_histories//20, NH=number_histories,
                    IDOUT=number_histories//5, IQOUT=number_histories//5,
                    NCP=number_species, IDREL=0, IQ0=0, IRC0=0, IRAND=12855897,
                    JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1e-16, INEL=1, IWC=3,
                    IDIFF=1, CSBE=0, ANGMN=0, ANGMX=90, TT=depth,
                    TTDYN=1.2*depth, NQX=number_layers, DSF=100.0, IQXN=0,
                    IQXX=250, IMCP=0,surf_name=surface.name+'.surf',
                    species_list=species_list, output_pka=0, output_ska=0,
                    output_prj=1, output_splst=1, output_rflst=1, output_prof=0,
                    output_simple=1)

                simulation_parameters.print_input_file()

                os.system('./FTridyn_Clean <'+str(name)+str(sim_number).zfill(4)+'.IN')
                simple_output = np.genfromtxt(str(name)+str(sim_number).zfill(4)+'.OUT',
                    skip_footer = 2, usecols = (2))
                num_prj = simple_output[0]
                num_spt = simple_output[4]
                num_bks = simple_output[1]
                num_dep = num_prj - num_spt
                self.sputtering_yield[energy_index, angle_index] = num_spt/number_histories
                self.reflection_yield[energy_index, angle_index] = num_bks/number_histories
                self.sim_numbers[energy_index, angle_index] = sim_number
                self.sputtered_lists[sim_number] = np.genfromtxt(str(name)+'SPLST.DAT')
                self.reflected_lists[sim_number] = np.genfromtxt(str(name)+'RFLST.DAT')
                self.implanted_lists[sim_number] = np.genfromtxt(str(name)+'DUMPPRJ.DAT')

                sim_number += 1
            #end for energy_index, energy
        #end for angle_index, angle
    #end def fill_lookup_table

    def find_tridyn_results(self, energy, angle):

        energy_index = min(range(len(self.energies)), key=lambda i: abs(self.energies[i] - energy))
        angle_index = min(range(len(self.angles)), key=lambda i: abs(self.angles[i] - angle))
        sputtering_yield = self.sputtering_yield[energy_index, angle_index]
        reflection_yield = self.reflection_yield[energy_index, angle_index]
        sim_number = int(self.sim_numbers[energy_index, angle_index])
        sputtered_list = self.sputtered_lists[sim_number]
        reflected_list = self.reflected_lists[sim_number]
        implanted_list = self.implanted_lists[sim_number]

        return sim_number, sputtering_yield, reflection_yield, sputtered_list, reflected_list, implanted_list
    #end def find_tridyn_results

    def find_tridyn_particle_results_iead(self, iead):
        new_sputtered_particles = [None]*int(np.sum(iead))
        new_reflected_particles = [None]*int(np.sum(iead))
        #p2c_multiplier, energy_sputtered, c_alpha_sputtered, c_beta_sputtered, c_gamma_sputtered
        particle_index = 0
        for energy_index, row in enumerate(iead):
            for angle_index, element in enumerate(row):
                if element > 0:
                    _, Y, R, Y_list, R_list, I_list = find_tridyn_results(self.energies[energy_index], self.angles[angle_index])
                    particle_list_index = 0
                    for new_particle_index in range(element):
                        y = new_particle_index%len(Y_list)
                        r = new_particle_index%len(R_list)

                        new_sputtered_particles[particle_index] = [Y, Y_list[y, 2], Y_list[y, 6], Y_list[y, 7], Y_list[y, 8]]

                        new_reflected_particles[particle_index] = [R, R_list[r, 2], R_list[r, 6], R_list[r, 7], R_list[r, 8]]

                        particle_index += 1
                    #end for
                #end if
            #end for
        #end for
    #end def find_tridyn_results_iead

    def find_tridyn_integral_results_iead(self, iead):
        sputtering_yields = np.zeros(np.shape(iead))
        reflection_yields = np.zeros(np.shape(iead))
        print(np.shape(self.energies))
        print(np.shape(self.angles))

        for energy_index, row in enumerate(iead):
            for angle_index, element in enumerate(row):
                if element > 0:
                    tridyn_results = self.find_tridyn_results(self.energies[energy_index], self.angles[angle_index])
                    sputtering_yields[energy_index, angle_index] = tridyn_results[1]
                    reflection_yields[energy_index, angle_index] = tridyn_results[2]
                #end if
            #end for angle_index, element
        #end for energy_index, row
        return sputtering_yields, reflection_yields
    #end def find_tridyn_integral_results_iead

    def save(self, filename):
        save_file = open(filename, 'wb')
        pickle.dump(self, save_file)
    #end def save
#end class tridyn_lookup

class species:
    """
    Stores species-relevant parameters for F-TRIDYN input.

    Args:
    ZZ (int): Atomic number.
    M (float): Mass in amu.
    BE (float): Bulk binding energy in eV.
    ED (float): Displacement energy in eV.
    EF (float): Threshold (cutoff) energy in eV.
    QU (float): Proportion of species in target.
    DNS0 (float): Number density of species in inverse Angstroms cubed.
    CK (float): Electronic stopping correction factor.
    E0 (float): Incident energy in eV.
    ALPHA0 (float): Incident energy in degrees.
    QUBEAM (float): Proportion of species in beam.
    QUMAX (float): Max proportion of species allowed in target.
    SBV (float): Surface binding energy in eV.
    """
    def __init__(self, ZZ, M, BE, ED, EF, QU, DNS0, CK, E0, ALPHA0, QUBEAM, QUMAX, SBV):
        self.ZZ = ZZ					#Atomic number/charge
        self.M = M						#Mass in [u]
        self.BE = BE					#Binding energy in [eV]
        self.ED = ED					#Displacement energy in [eV]
        self.EF = EF					#Cutoff energy in [eV]
        self.QU = QU					#Concentration of species in homogeneous target
        self.DNS0 = DNS0				#Number density in [at/A^3]
        self.CK = CK					#Electronic stopping correction factor (0.0<CK<1.0)
        self.E0 = E0					#Incident energy in [eV]
        self.ALPHA0 = ALPHA0			#Incident angle in [deg] (0.0<ALPHA0<90.0)
        self.QUBEAM = QUBEAM			#Concentration in beam (0.0<QUBEAM<1.0)
        self.QUMAX = QUMAX				#Maximum allowed concentration in target (0.0<QUMAX<1.0)
        self.SBV = SBV					#Surface binding energy in [eV]
    #end def __init__

    def __repr__(self):
        #return f'species({ZZ}, {M}, {BE}, {ED}, {EF}, {QU}, {DNS0}, {CK}, {E0}, {ALPHA0}, {QUBEAM}, {QUMAX}, {SBV})'
        pass
    #end def __repr__
#end class species

class species_lookup_table:
    """
    Provides functions to lookup and generate species objects for F-TRIDYN
    input. Searches for file 'table1.txt' and stores the first section.
    """
    def __init__(self):
        cwd = os.path.dirname(os.path.realpath(__file__))
        self.table1 = open(cwd+'/'+'table1.txt','r').readlines()[11:130]
        self.species_name = []
        self.ZZ = []
        self.M = []
        self.BE = []
        self.ED = []
        self.EF = []
        self.DNS0 = []
        self.CK = []
        self.SBV = []
        for line in self.table1:
            self.species_name.append(line[0:5].strip())
            self.ZZ.append(float(line[6:10].strip()))
            self.M.append(float(line[11:22].strip()))
            self.BE.append(float(line[114:120].strip()))
            self.ED.append(float(line[53:59].strip()))
            self.EF.append(float(line[62:68].strip()))
            self.DNS0.append(float(line[32:41].strip()))
            self.SBV.append(float(line[43:49].strip()))
        #end for
    #end def __init__

    def __repr__(self):
        return 'species_lookup_table()'

    def find_species(self,name,QU=0.0,CK=1.0,E0=0.0,ALPHA0=0.0,QUBEAM=0.0,QUMAX=0.0):
        """
        Lookup species parameters from included data file and return an appropriate
        species object for F-TRIDYN input.

        Args:
        name (:obj:'str'): species name (i.e., chemical symbol).
        QU (float, optional): portion of species in target.
        CK (float, optional): Electronic stopping correction factor.
        E0 (float, optional): Incident energy of species in eV.
        ALPHA0 (float, optional): Incident angle of species in degrees.
        QUBEAM (float, optional): Portion of species in beam.
        QUMAX (float, optional): Maximum allowed proportion in target.
        """
        species_index = self.species_name.index(name)
        return species(ZZ=self.ZZ[species_index],
        M=self.M[species_index], BE=self.BE[species_index],
        ED=self.ED[species_index], EF=self.EF[species_index],
        QU=QU, DNS0=self.DNS0[species_index], CK=CK, E0=E0, ALPHA0=ALPHA0,
        QUBEAM=QUBEAM, QUMAX=QUMAX, SBV=self.SBV[species_index])
    #end def find_species
#end def species_lookup

class sim_params:
    """
    Holds simulation parameters for F-TRIDYN simulations and provides methods for
    generating and producing input files.

    Args:
    name (:obj:'str'): a 4-character name that represents the simulation (e.g., 'He_W')
    sim_num (int): Number of simulation. Useful for ordering multiple simulations with same name.
    IFOUT (int): Number of projectiles after which to write info to terminal.
    NH (int): Number of projectile histories (trajectories) to perform.
    IDOUT (int): Number of projectiles after which to produce integral output.
    IQOUT (int): Number of projectiles after which to produce profile files (.PR##).
    NCP (int): Number of speicies in simulation.
    IDREL (int): When > 0, relaxation is suppressed. When < 0, relaxation and cascades are suppressed.
    IQ0 (int): When == 0, homogeneous initial composition. When != 0, inhomogeneous initial composition from .LAY file.
    IRC0 (int): When < 0, subthreshold atoms are free. When >= 0, subthreshold atoms are bound.
    IRAND (int): Seed for PRNG.
    JSP1 (int): Suppression of recoils of species JSP1 to JSP2.
    JSP2 (int): Suppression of recoils of species JSP1 to JSP2.
    JFRP (int): Frenkel Pair generation of species 1 to JFRP./
    JNRM (int): Normalization of print output.
    FLC (float): Fluence for dynamic composition simulation.
    INEL (int): Electronic stopping model. 1: inelastic nonlocal; 2: local; 3: equipartition
    IWC (int): Order of weak collisions (i.e., how many possible weak collisions per collision event)
    IDIFF (int): Extra particle handling. 1: Reemission; 2: "Diffusion"
    CSBE (int): Surface binding energy model.
    ANGMN (float): Minimum angle to include in particle list output.
    ANGMX (float): Maximum angle to include in particle list output.
    TT (float): Depth of surface.
    TTDYN (float): Depth of dynamic surface. Defaults to 20% more than TT.
    NQX (int): Number of composition layers.
    DSF (float): Averaging depth for surface composition output.
    IQXN (int): Limiting depth intervals for profile output.
    IQXX (int): Unused.
    IMCP (int): Species for which moments are calculated.
    surf_name (:obj:'str'): Name of surface file.
    species_list (:obj:'list'): List of species in simulation.
    output_pka (int, optional): Whether to output PKA list. Defaults to 0. [0,1]
    output_ska (int, optional): Whether to output SKA list. Defaults to 0. [0,1]
    output_prj (int, optional): Whether to output implanted projectile list. Defaults to 1. [0,1]
    output_out (int, optional): Whether to output human-readable output. Defaults to 1. [0,1]
    output_srfc (int, optional): Whether to output surface composition. Defaults to 0. [0,1]
    output_spyl (int, optional): Whether to output sputtering yield. Defaults to 0. [0,1]
    output_ardn (int, optional): Unused output. Defaults to 0. [0,1]
    output_srrs (int, optional): Unused output. Defaults to 0. [0,1]
    output_reem (int, optional): Whether to output reemitted particles. Defaults to 0. [0,1]
    output_splst (int, optional): Whether to output sputtered list. Defaults to 1. [0,1]
    output_rflst (int, optional): Whether to output reflected list. Defaults to 0. [0,1]
    output_trlst (int, optional): Whether to output transmitted list. Defaults to 0. [0,1]
    output_simple (int, optional): Whether to output simple output. Defaults to 0. [0,1]
    output_fp (int, optional): Whether to output Frenkel Pair list. Defaults to 0. [0,1]
    output_prof (int, optional): Whether to output composition profiles. Defaults to 0. [0,1]

    """
    def __init__(self,name,sim_num,IFOUT,NH,IDOUT,IQOUT,NCP,IDREL,IQ0,IRC0,IRAND,
    JSP1,JSP2,JFRP,JNRM,FLC,INEL,IWC,IDIFF,CSBE,ANGMN,ANGMX,TT,TTDYN,NQX,DSF,
    IQXN,IQXX,IMCP,surf_name,species_list, output_pka=0, output_ska=0, output_prj=1,
    output_out=1,output_srfc=0,output_spyl=0,output_ardn=0,output_srrs=0,output_reem=0,
    output_splst=1,output_rflst=0,output_trlst=0,output_simple=0,output_fp=0,output_prof=0):
        self.name = name 						#4-character file name, e.g., He_W, C_C_, NaNa, etc.
        self.sim_num = str(sim_num).zfill(4) 	#4-number simulation number, e.g., 0001, 0002
        self.IFOUT = int(IFOUT) 				#Number of projectiles after which info is written to terminal
        self.NH = int(NH)						#Number of incident histories (i.e. projectiles)
        self.IDOUT = int(IDOUT) 				#Integral data ouput after every IDOUT'th history
        self.IQOUT = int(IQOUT) 				#Additional .PR## profile file output after every IQOUT'th history
        self.NCP = int(NCP) 					#Number of components (i.e. species)
        self.IDREL = int(IDREL) 				#IDREL>0:Suppression of relaxation, IDREL<0:Suppression of relaxation and cascades
        self.IQ0 = int(IQ0) 					#IQ0=0:homogeneous initial composition, IQ0!=0:inhomogeneous initial composition from .LAY
        self.IRC0 = int(IRC0) 					#IRC0<0: Subthreshold atoms free, IRC0>=0: Subthreshold atoms bound
        self.IRAND = int(IRAND) 				#Random seed
        self.JSP1 = int(JSP1) 					#Suppression of recoils of type JSP1,...,JSP2
        self.JSP2 = int(JSP2) 					#JSP1=0: JSP2 ignored, all recoils traced
        self.JFRP = int(JFRP) 					#Frenkel Pair generation for JFRP,...,NCP
        self.JNRM = int(JNRM) 					#Normalzation of print output
        self.FLC = FLC 							#Fluence [1e16/cm2]
        self.INEL = int(INEL) 					#INEL=[1,2,3]:[inelastic nonlocal, local, equipartition]
        self.IWC = int(IWC) 					#Max order of weak projectile-target collisions {1,2,3}
        self.IDIFF = int(IDIFF) 				#IDIFF=[0,1]:[reemission,"diffusion"]
        self.CSBE = CSBE 						#Surface binding energy model (see documentation)
        self.ANGMN = ANGMN 						#Minimum angle for sputtered energy/angle distributions
        self.ANGMX = ANGMX 						#Maximum angle for sputtered energy/angle distributions
        self.TT = TT 							#Depth of surface
        self.TTDYN = TTDYN 						#Depth of dynamic surface (>=TT)
        self.NQX = NQX 							#Number of composition layers within F-TRIDYN
        self.DSF = DSF 							#Averaging depth for surface composition
        self.IQXN = IQXN 						#Limiting depth intervals for .PR## profile output
        self.IQXX = IQXX						#
        self.IMCP = int(IMCP) 					#Component for which moments should be calcluated

        self.surf_name = surf_name 				#Name of surface file
        self.species_list = species_list 		#List of species objects

        self.output_pka = int(output_pka) 		#PKA output from F-TRIDYN
        self.output_ska = int(output_ska) 		#SKA output from F-TRIDYN
        self.output_prj = int(output_prj) 		#Projectile output from F-TRIDYN

        self.output_out = int(output_out) 		#Human-readable output file
        self.output_srfc = int(output_srfc) 	#
        self.output_spyl = int(output_spyl) 	#Sputtering yield calculation
        self.output_ardn = int(output_ardn) 	#
        self.output_srrs = int(output_srrs) 	#
        self.output_reem = int(output_reem) 	# Re-emission
        self.output_splst = int(output_splst) 	#Sputtered particle list
        self.output_rflst = int(output_rflst) 	#Reflected particle list
        self.output_trlst = int(output_trlst) 	#Transmitted particle list
        self.output_simple = int(output_simple) #Simple output file
        self.output_fp = int(output_fp) 		#Frenkel Pair creation site list
        self.output_prof = int(output_prof) 	#Dynamic composition profiles
        #todo: error checking, NCP=length(species_list), all parameters within bounds, etc.
    #end def __init__

    def print_input_file(self):
        """
        Generates an F-TRIDYN input file with a .IN extension.

        Input file format is as follows:
        NAMExxxx.IN
        IFOUT
        NH  IDOUT  IQOUT  NCP  IDREL  IQ0  IRC0  IRAND  JCP1  JCP2  JFRP  JNRM
        FLC  INEL  IWC  IDIFF  CSBE  ANGMN  ANGMX
        TT  TTDYN  NQX  DSF  IQXN  IQXX  IMCP
        surfname.txt output_flags ...
        ZZ(1)  M(1)  BE(1)  ED(1)  EF(1)  QU(1)  DNS0(1)  CK(1)
        E0(1)  ALPHA0(1)  QUBEAM(1)  QUMAX(1)  SBV(1)  SBV(1)
        ZZ(2)  M(2)  BE(2)  ED(2)  EF(2)  QU(2)  DNS0(2)  CK(2)
        E0(2)  ALPHA0(2)  QUBEAM(2)  QUMAX(2)  SBV(2)  SBV(2)
        ZZ(3)...
        ZZ(NCP)...
        """
        input_file = open(self.name+self.sim_num+'.IN','w+')
        print(self.name+self.sim_num+'.IN', file=input_file)
        print(self.IFOUT,file=input_file)
        print(self.NH, self.IDOUT, self.IQOUT, self.NCP, self.IDREL, self.IQ0,
        self.IRC0, self.IRAND, self.JSP1, self.JSP2, self.JFRP, self.JNRM,
        file=input_file)
        print(self.FLC, self.INEL, self.IWC, self.IDIFF, self.CSBE, self.ANGMN,
        self.ANGMX, file=input_file)
        print(self.TT, self.TTDYN, self.NQX, self.DSF, self.IQXN, self.IQXX,
        self.IMCP, file=input_file)
        print(self.surf_name,self.output_pka, self.output_ska, self.output_prj,
        self.output_out,self.output_srfc,self.output_spyl,self.output_ardn,
        self.output_srrs,self.output_reem,self.output_splst,self.output_rflst,
        self.output_trlst,self.output_simple,self.output_fp,self.output_prof,
        file=input_file)
        for i in range(self.NCP):
            current_species = self.species_list[i]

            print(current_species.ZZ, current_species.M, current_species.BE,
            current_species.ED, current_species.EF, current_species.QU,
            current_species.DNS0, current_species.CK, file=input_file)

            print(current_species.E0, current_species.ALPHA0,
            current_species.QUBEAM, current_species.QUMAX, end='',
            file=input_file)

            for j in range(self.NCP):
                print(' '+str(current_species.SBV), end='', file=input_file)
            #end for

            print(file=input_file)
        #end for
    #end def print_input_file
#end class simulation

class fractal_surface:
    """
    Stores and provides methods for generating Mandelbrot-style fractal surfaces.

    Args:
    FD (float, optional): fractal dimension of generated fractal surface.
    Defaults to 1.0, a smooth surface.
    width (float, optional): Total width of fractal surface in Angstroms.
    Defaults to 400 Angstroms.
    iterations (int, optional): Number of generations to iterate. Defaults to 1.

    """
    def __init__(self, FD=1.0, width=400.0, iterations=1):
        self.FD = FD
        self.width = width
        if FD==1.0:
            self.iterations = 1
        else:
            self.iterations = iterations
        self.fgen()
        self.name = str(1)+'p'+str(np.int((FD-0.99999)*1000.0))
    #end def __init__

    def fgen(self,shape=[0,1,0,-1,0,-1,0,1,0]):
        """
        Given a shape, produce a generated fractal using Mandelbrot's fractal
        generator method, whereby an initial shape is transformed onto each
        segment of itself over a specified number of generations to produce a
        curve with an arbitrary fractal dimension given an initial shape type.

        Args:
        shape (:obj:'list', optional) This function defines the initial
        fractal generator. Values in the list are 0 (horizontal line
        segment), 1 (upward angled line segment) and -1 (downward angled
        line segment). Sum of shape should be zero for consistency. Angles
        are calculated from FD.

        Defaults to [0,1,0,-1,0,-1,0,1,0]. This produces surfaces that look
        like the following illustration:

        y=0       y=width
        |  _      |
        |_/ \_   _|_ x=0
        |     \_/ |

        As an example of a different shape, the shape [0,1,0,-1,0] looks like
        the following illustration:

        y=0     y=width
        |       |
        |_/\_  _|_ x=0
        |    \/ |
        """
        b = np.sum(np.abs(shape))
        a = np.size(shape)-b
        beta = np.arccos((-a+(a+b)**(1.0/self.FD))/b)
        L = 1.0/(a+b*np.cos(beta))
        num_gen_seg = np.size(shape)
        self.num_gen_points = num_gen_seg+1

        self.x = np.zeros(self.num_gen_points)
        self.y = np.zeros(self.num_gen_points)
        segX = np.zeros(self.num_gen_points)
        segY = np.zeros(self.num_gen_points)
        self.gx = np.zeros(self.num_gen_points)
        self.gy = np.zeros(self.num_gen_points)

        self.x[0] = 0.0
        self.y[0] = 0.0

        for i in range(1,self.num_gen_points):
            self.x[i] = self.x[i-1] + L * np.cos(beta*shape[i-1])
            self.y[i] = self.y[i-1] + L * np.sin(beta*shape[i-1])
        #end for

        self.gx = self.x
        self.gy = self.y

        for i in range(2,self.iterations+1):
            counter = -1
            storeX = self.x
            storeY = self.y
            self.x = np.zeros(num_gen_seg**(i)+1)
            self.y = np.zeros(num_gen_seg**(i)+1)

            for j in range(np.size(storeX)-1):
                segX[0] = storeX[j]
                segY[0] = storeY[j]
                gamma = np.arctan2((storeY[j+1]-segY[0]),(storeX[j+1]-segX[0]))

                for k in range(1,num_gen_seg+1):
                    coordX = L**i*np.cos(beta*shape[k-1])
                    coordY = L**i*np.sin(beta*shape[k-1])
                    segX[k] = segX[k-1] + coordX*np.cos(gamma) - coordY*np.sin(gamma)
                    segY[k] = segY[k-1] + coordX*np.sin(gamma) + coordY*np.cos(gamma)
                #end for k

                if(j<np.size(storeX)-2):
                    segTrim = -1
                else:
                    segTrim = 0
                #end if

                for n in range(0,np.size(segX)+segTrim):
                    counter = counter + 1
                    self.x[counter] = segX[n]
                    self.y[counter] = segY[n]
                #end for n
            #end for j
        #end for i
        self.num_surf_points = np.size(self.x)
        self.x,self.y = self.width*self.x,self.width*self.y #scale the final surface
    #end def fgen

    def print_fractal_surface_file(self):
        """
        Print fractal surface file for F-TRIDYN input with .surf extension.

        File format is as follows:
        <number of generator points> <number of iterations> <number of points>
        <fractal dimension>
        <generator x points> <generator y points>
        ...
        <fractal x points> <fractal y points>
        ...

        Returns:
            filename (:obj:'str'): name of created file.
        """
        surface_file = open(self.name+'.surf','w+')

        print(self.num_gen_points, self.iterations, self.num_surf_points, file=surface_file)
        print(self.FD,file=surface_file)

        for i in range(self.num_gen_points):
            print(self.gx[i], self.gy[i], file=surface_file)
        #end for

        for i in range(self.num_surf_points):
            print(self.x[i], self.y[i], file=surface_file)
        #end for

        return self.name+'.surf'
    #end def print_fractal_surface_file
#end class fractal_surface

def He_W_xolotl(sim_number=1, number_histories=1000000,
    incident_energy=200.0, incident_angle=0.0, fractal_dimension=1.0,
    IQ0=0, number_layers=100, depth=200.0):

    #Lookup species physical parameters from table
    lookup_table = species_lookup_table()
    Helium = lookup_table.find_species('He', E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
    Tungsten = lookup_table.find_species('W', QU=1.0, QUMAX=1.0)

    #Collect species into species_list
    species_list = [Helium, Tungsten]
    number_species = np.size(species_list)

    #Create fractal and print surface file to 1p###.surf
    surface = fractal_surface(FD=fractal_dimension, width=400.0, iterations=3)
    surface.print_fractal_surface_file()

    #Define simulation parameters for F-TRIDYN
    simulation_parameters = sim_params(name='He_W', sim_num=sim_number, IFOUT=int(number_histories//20),
    NH=int(number_histories), IDOUT=int(number_histories//5), IQOUT=int(number_histories//5), NCP=number_species, IDREL=1,
    IQ0=IQ0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1E-16, INEL=1,
    IWC=3, IDIFF=1, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=number_layers,
    DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
    output_pka=0, output_ska=0, output_prj=1)

    #Print simulation parameters to <name><sim_number>.IN
    simulation_parameters.print_input_file()
#end def He_W

def Prj_Tg_xolotl(sim_number=1, number_histories=1000000,
                  incident_energy=200.0, incident_angle=0.0, fractal_dimension=1.0,
                  IQ0=0, number_layers=100, depth=200.0, projectile_name='He',target1_name='W', target2_name='',target3_name='',target4_name=''):

        #Lookup species physical parameters from table
        lookup_table = species_lookup_table()
        projectile=lookup_table.find_species(projectile_name, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
        target1=lookup_table.find_species(target1_name, QU=1.0, QUMAX=1.0)

        #Collect species into species_list
        species_list = [projectile,target1]

        #add target species if other than empty:
        #currently max of 5 species handled by F-Tridyn
        additional_targets=[target2_name, target3_name, target4_name]

        for  sp in additional_targets:
                if sp:
                        add_target=lookup_table.find_species(sp, QU=1.0, QUMAX=1.0)
                        species_list.append(add_target)


        number_species = np.size(species_list)

        #Create fractal and print surface file to 1p###.surf
        surface = fractal_surface(FD=fractal_dimension, width=400.0, iterations=3)
        surface.print_fractal_surface_file()

        #string of 4-characters required for name
        if (len(projectile_name)+len(target1_name))==4:
                parameterNameString=projectile_name+target1_name #already 4 char, e.g. HeTa
        elif (len(projectile_name)+len(target1_name))==3:
                parameterNameString=projectile_name+'_'+target1_name #e.g.: He_W
        else: #len(parameterNameString) ==2
                parameterNameString=projectile_name+'_'+target1_name+'_' #e.g., W_W -> W_W_

        #Define simulation parameters for F-TRIDYN
        simulation_parameters = sim_params(name=parameterNameString, sim_num=sim_number, IFOUT=int(number_histories//20),
                NH=int(number_histories), IDOUT=int(number_histories//5), IQOUT=int(number_histories//5), NCP=number_species, IDREL=1,
                IQ0=IQ0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1E-16, INEL=1,
                IWC=3, IDIFF=1, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=number_layers,
                DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
                output_pka=0, output_ska=0, output_prj=1)

        #Print simulation parameters to <name><sim_number>.IN
        simulation_parameters.print_input_file()
#end def Prj_Tg

def beam_and_target(name,beam_species,target_species,sim_number=1,
    number_histories=1E5,incident_energy=100.0,incident_angle=0.0,
    fractal_dimension=1.0,width=200.0,depth=200.0):
    """
    A function that generates an input file given a monoenergetic beam on a
    homogeneous target.

    Args:
    name (:obj:'str'): a 4-character name that represents the simulation (e.g., 'He_W')
    beam_species (:obj:'str'): Chemical symbol of beam species.
    target_species (:obj:'str'): Chemical symbol of target species.
    sim_number (int, optional): Number of simulation. Useful for ordering multiple simulations with same name.
    number_histories (int, optional): Number of histories (trajectories). Defaults to 1E5.
    incident_energy (float, optional): Incident beam energy. Defaults to 100.0. [eV]
    incident_angle (float, optional): Incident angle of beam. Defaults to 0.0. [degrees]
    fractal_dimension (float, optional): FD of surface. Defaults to 1.0 (smooth).
    width (float, optional): Total width of surface. [A]
    depth (float, optional): Total depth of surface. [A]

    Returns:
    simulation_parameters (:obj:'sim_params'): sim_params object.

    Prints:
    A single F-TRIDYN input file (extension: .IN) for use in F-TRIDYN simulations.

    """

    lookup_table = species_lookup_table()
    beam = lookup_table.find_species(beam_species, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0)
    target = lookup_table.find_species(target_species, QU=1.0, QUMAX=1.0)
    species_list = [beam,target]
    number_species = np.size(species_list)

    surface = fractal_surface(FD=fractal_dimension, width=width, iterations=3)
    surface.print_fractal_surface_file()

    simulation_parameters = sim_params(name=name, sim_num=sim_number, IFOUT=number_histories//20,
    NH=number_histories, IDOUT=number_histories//5, IQOUT=number_histories//5, NCP=number_species, IDREL=0,
    IQ0=0, IRC0=-1, IRAND=12855897, JSP1=0, JSP2=1, JFRP=1, JNRM=1, FLC=1e-16, INEL=1,
    IWC=3, IDIFF=0, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=500,
    DSF=100.0, IQXN=0, IQXX=250, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
    output_pka=0, output_ska=0, output_prj=1, output_rflst=1, output_splst=1)
    simulation_parameters.print_input_file()

    return simulation_parameters
#end def beam_and_target

def preferential_sputtering_boron_carbon(boron_conc, carbon_conc, beam_species='D', number_histories=10000, sim_number=0):
    target_species = ['B', 'C_a']
    flux = 1e-16
    #fluxes_normalized = [flux/np.sum(fluxes) for flux in fluxes]

    incident_energy = 2500 #3 kTe from Stangeby
    incident_angle = 0.

    lookup_table = species_lookup_table()
    beam = lookup_table.find_species(beam_species, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=1.0, QUMAX=0.0, QU=0.0)

    targets = [lookup_table.find_species(target, QU=1.0, QUMAX=1.0, E0=incident_energy, ALPHA0=incident_angle, QUBEAM=0.0) for target in target_species]

    targets[0].QU = boron_conc
    targets[1].QU = carbon_conc

    species_list = [beam, targets[0], targets[1]]
    number_species = len(species_list)

    surface = fractal_surface(FD=1.0, width=400., iterations=1)
    surface.print_fractal_surface_file()

    name = 'CONC'
    IDREL = -1 #Suppression of dynamic relaxation + cascades
    IDIFF = 1 #Turn off diffusion
    IQ0 = 0 #Start homogeneous
    IRC0 = -1 #Subtrhreshold atoms bound
    FLC = 1e-16 #Fluence*1e16
    depth = 100.
    IWC = 1 #max 1 collision
    NQX = 250

    simulation_parameters = sim_params(name=name, sim_num=sim_number, IFOUT=number_histories//20,
    NH=number_histories, IDOUT=number_histories//25, IQOUT=number_histories//25, NCP=number_species, IDREL=IDREL,
    IQ0=IQ0, IRC0=IRC0, IRAND=12855897 - int(boron_conc*100) + int(carbon_conc*100), JSP1=0, JSP2=0, JFRP=1, JNRM=1, FLC=FLC, INEL=1,
    IWC=IWC, IDIFF=IDIFF, CSBE=0, ANGMN=0, ANGMX=90, TT=depth, TTDYN=1.2*depth, NQX=NQX,
    DSF=100.0, IQXN=0, IQXX=NQX, IMCP=0, surf_name=surface.name+'.surf', species_list=species_list,
    output_pka=0, output_ska=0, output_prj=0, output_prof=1, output_srfc=0, output_splst=1)
    simulation_parameters.print_input_file()

def test_tridyn_lookup_and_interface():
    energies = np.linspace(100.0, 200.0, 4)
    angles = np.linspace(0.0, 89.0, 3)
    beam_species = 'W'
    target_species = 'B'
    #tl = tridyn_lookup(energies, angles, beam_species, target_species)
    #tl.fill_lookup_table('H_B_')
    #tl.save('tridyn_lookup.pickle')
    #N, S, R, splst, rflst, implst = tl.find_tridyn_results(200.0, 70.0)
    #print(N, S, R, splst, rflst, implst)
    #S_array, R_array = tl.find_tridyn_integral_results_iead(np.floor(np.random.uniform(0.0, 200.0, (4,3))))
    #print(np.max(S_array))
    #print(np.max(R_array))
    iead = np.floor(np.random.uniform(0.0, 200.0, (len(energies), len(angles))))
    ti = tridyn_interface('W', 'B')
    new_particle_list, new_particle_number = ti.run_tridyn_simulations_from_iead(energies, angles, iead)
#end def test_tridyn_lookup

def main():
    load_from_file = True
    run_hpic = False
    Te = 10.0
    Ti = 10.0

    psi = 86.0
    B0 = 2.0

    m_background = 2.0
    n_background = 1e19

    m_B = 10.18
    n_B1 = 1e12
    n_B2 = 2e13
    n_B3 = 4e13
    n_B4 = 9e13
    n_B5 = 2e13
    n_B = [n_B1, n_B2, n_B3, n_B4, n_B5]

    num_debye_lengths = 200
    num_gridpoints_per_debye_length = 3
    num_timesteps_per_gyroperiod = 3
    num_ion_transit_times = 10.0
    num_particles_per_cell = 100

    BC_left = 0.0
    BC_right = 0.0

    kinfo = 10
    kgrid = 0
    kpart = 10

    sim_name = 'boron_5'
    command_line_string = f'{sim_name} {num_debye_lengths} {num_gridpoints_per_debye_length} {num_timesteps_per_gyroperiod} {num_ion_transit_times} {num_particles_per_cell} {B0} {psi} {Te} {Ti} {BC_left} {BC_right} {kinfo} {kgrid} {kpart} '
    species_string = f'{m_background} 1 {n_background} {m_B} 1 {n_B1} {m_B} 2 {n_B2} {m_B} 3 {n_B3} {m_B} 4 {n_B4} {m_B} 5 {n_B5}'

    if run_hpic:
        os.system('./hpic -command_line ' + command_line_string + species_string)
        os.system('rm *PARTICLEDATA*')
        os.system('rm *PHI*')
        os.system('cp *.dat dat/.')

    #n_B1, n_B2, n_B3, n_B4, n_B5 = 1, 1, 1, 1, 1

    iead_b = n_B1*np.genfromtxt(f'dat/{sim_name}_IEAD_sp1.dat')
    iead_b += n_B2*np.genfromtxt(f'dat/{sim_name}_IEAD_sp2.dat')
    iead_b += n_B3*np.genfromtxt(f'dat/{sim_name}_IEAD_sp3.dat')
    iead_b += n_B4*np.genfromtxt(f'dat/{sim_name}_IEAD_sp4.dat')
    iead_b += n_B5*np.genfromtxt(f'dat/{sim_name}_IEAD_sp5.dat')
    #np.savetxt('iead_b.dat', iead_b)
    #iead_b_unwrapped = np.reshape(iead_b, (np.size(iead_b)))

    #iead_b[iead_b<1e16]=0
    iead_d = np.genfromtxt(f'dat/{sim_name}_IEAD_sp0.dat')
    #np.savetxt('iead_d.dat', iead_d)

    #IEAD D
    plt.figure(figsize=(8, 6), dpi=300.0)
    energies = np.linspace(0.0, 24.0*Te, 240)
    max_energy_plot = energies[199]
    dE = energies[1] - energies[0]
    energies += dE/2.0
    angles = np.linspace(0.0, 90.0, 90)
    plt.pcolormesh(angles, energies[0:200], iead_d[0:200, :], cmap='plasma')

    #plt.title('D IEAD', fontdict=title_font)
    plt.xlabel('Angle [deg]')
    plt.ylabel('Energy [eV]')
    plt.axis([0.0, 90.0, 0.0, max_energy_plot])
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        top=False,      # ticks along the bottom edge are off
        bottom=True,         # ticks along the top edge are off
        labelbottom=True,
        labeltop=False
    )
    plt.savefig('iead_d.png')
    #plt.xticks(np.arange(0.0, max_energy_plot, 50))
    #plt.yticks(np.arange(0.0, 1.1, 0.5))

    #TRIDYN overlay
    energies_tridyn = np.linspace(10.0, 300.0, 25)
    energies_tridyn = np.linspace(0.0, 290.0, 25)
    angles_tridyn = np.linspace(0.0, 89.0, 25)
    yield_tridyn = np.genfromtxt('D_B__yield.sav')
    reflection_tridyn = np.genfromtxt('D_B__reflected.sav')

    c1 = plt.contour(angles_tridyn, energies_tridyn, yield_tridyn, colors='red', levels=[0.0, 0.005, 0.01, 0.02, 0.05, 0.07], linewidths=2, zorder=2)
    c1_labels = plt.clabel(c1, fontsize=11)
    for contour_label in c1_labels:
        contour_label.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'), path_effects.Normal()])
    plt.savefig('iead_d_overlay_yield.png')

    c2 = plt.contour(angles_tridyn, energies_tridyn, reflection_tridyn, colors='white', levels=[0.4, 0.6, 0.8, 0.95], linewidths=2, zorder=2)
    c2_labels = plt.clabel(c2, fontsize=11)
    for contour_label in c2_labels:
        contour_label.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'), path_effects.Normal()])
    #c2.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'), path_effects.Normal()])

    #Flattened A, E distributions
    angle_dist = np.sum(iead_d[0:200, :], axis=0)
    energy_dist = np.sum(iead_d[0:200, :], axis=1)
    peak_angle = angles[np.argmax(angle_dist)]
    peak_energy = energies[np.argmax(energy_dist)]
    peak_value = 1.0

    #Most probably yields
    R = interp2d(angles_tridyn, energies_tridyn, reflection_tridyn)
    Y = interp2d(angles_tridyn, energies_tridyn, yield_tridyn)
    print(peak_angle, peak_energy)
    print(R(peak_angle, peak_energy))
    print(Y(peak_angle, peak_energy))

    plt.savefig('iead_d_overlay_both.png')
    plt.scatter(peak_angle, peak_energy, marker='*', color='purple', s=100, zorder=1)
    plt.savefig('iead_d_overlay_both_star.png')

    #ED
    plt.figure(figsize=(3, 6), dpi=300.0)
    plt.plot(energy_dist/np.max(energy_dist), energies[0:200], c='black')
    #plt.title('D Energies', fontdict=title_font)
    plt.xlabel('f(E)')

    plt.xticks(np.arange(0.0, 1.1, 0.5))
    plt.plot( [0.0, 1.2], [peak_energy, peak_energy], ':', color='black')
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right=False,      # ticks along the bottom edge are off
        left=False,         # ticks along the top edge are off
        labelleft=False
    )
    plt.axis([0.0, 1.20, 0.0, max_energy_plot])
    plt.text(0.0, peak_energy+1, f'max: {int(peak_energy)} eV')
    plt.savefig('ed_d.png')

    #AD
    plt.figure(figsize=(8, 3), dpi=300.0)
    plt.plot(angles, angle_dist/np.max(angle_dist), c='black')
    #plt.title('D Angles', fontdict=title_font)
    plt.ylabel('f(θ)')
    plt.yticks(np.arange(0.0, 1.1, 0.5))

    plt.plot([peak_angle, peak_angle], [0.0, 1.2], ':', color='black')
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False
    )
    plt.axis([0.0, 90.0, 0.0, 1.2])
    plt.text(45, 1.05, f'max: {int(peak_angle)} deg')
    plt.savefig('ad_d.png')

    #IEAD B
    plt.figure(figsize=(8, 6), dpi=300.0)
    energies = np.linspace(0.0, 24.0*Te, 240)
    max_energy_plot = energies[199]
    dE = energies[1] - energies[0]
    energies += dE/2.0
    angles = np.linspace(0.0, 90.0, 90)
    plt.pcolormesh(angles, energies[0:200], iead_b[0:200, :], cmap='plasma')

    #plt.title('B IEAD', fontdict=title_font)
    plt.xlabel('Angle [deg]')
    plt.ylabel('Energy [eV]')
    plt.axis([0.0, 90.0, 0.0, max_energy_plot])
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        top=False,      # ticks along the bottom edge are off
        bottom=True,         # ticks along the top edge are off
        labelbottom=True,
        labeltop=False
    )
    plt.savefig('iead_b.png')
    #plt.xticks(np.arange(0.0, max_energy_plot, 50))
    #plt.yticks(np.arange(0.0, 1.1, 0.5))

    #Flattened A, E distributions
    angle_dist = np.sum(iead_b[0:200, :], axis=0)
    energy_dist = np.sum(iead_b[0:200, :], axis=1)
    peak_angle = angles[np.argmax(angle_dist)]
    peak_energy = energies[np.argmax(energy_dist)]
    peak_value = 1.0

    #TRIDYN overlay
    energies_tridyn = np.linspace(10.0, 300.0, 25)
    angles_tridyn = np.linspace(0.0, 89.0, 25)
    yield_tridyn = np.genfromtxt('B_B__yield.sav')
    reflection_tridyn = np.genfromtxt('B_B__reflected.sav')

    R = interp2d(angles_tridyn, energies_tridyn, reflection_tridyn)
    Y = interp2d(angles_tridyn, energies_tridyn, yield_tridyn)
    print(peak_angle, peak_energy)
    print(R(peak_angle, peak_energy))
    print(Y(peak_angle, peak_energy))

    c1 = plt.contour(angles_tridyn, energies_tridyn, yield_tridyn, colors='red', levels=[0.01, 0.1, 0.2, 0.3, 0.5], linewidths=2)
    c1_labels = plt.clabel(c1, fontsize=11)
    for contour_label in c1_labels:
        contour_label.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'), path_effects.Normal()])
    plt.savefig('iead_b_overlay_yield.png')

    c2 = plt.contour(angles_tridyn, energies_tridyn, reflection_tridyn, colors='white', levels=[0.1, 0.3, 0.5, 0.9], linewidths=2)
    c2_labels = plt.clabel(c2, fontsize=11)
    for contour_label in c2_labels:
        contour_label.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'), path_effects.Normal()])
    #c2.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'), path_effects.Normal()])

    plt.savefig('iead_b_overlay_both.png')
    plt.scatter(peak_angle, peak_energy, marker='*', color='purple', s=100, zorder=1)
    plt.savefig('iead_b_overlay_both_star.png')

    #plt.show()
    #exit()

    #iead_b = np.genfromtxt('dat/iead_b.dat')
    iead_b *= 10 / np.max(iead_b)
    iead_b = iead_b.astype('int')

    #iead_d = np.genfromtxt('dat/iead_d.dat')
    iead_d *= 10 / np.max(iead_d)
    iead_d = iead_d.astype('int')

    num_incident_b = np.sum(iead_b)
    num_incident_d = np.sum(iead_d)

    bw = tridyn_interface('B', 'W')
    bb = tridyn_interface('B', 'B')
    bc = tridyn_interface('B', 'C_a')

    dw = tridyn_interface('D', 'W')
    db = tridyn_interface('D', 'B')
    dc = tridyn_interface('D', 'C_a')

    num_hist = 1e2

    if load_from_file:
        output_file = open('output.sav', 'rb')
        output = pickle.load(output_file)
        sputtered_bw, reflected_bw, sputtered_bb, reflected_bb, sputtered_bc, reflected_bc, sputtered_dw, reflected_dw, sputtered_dw, reflected_dw, sputtered_db, reflected_db, sputtered_dc, reflected_dc = output
        output_file.close()
    else:
        sputtered_bw, reflected_bw = bw.run_tridyn_simulations_from_iead(energies, angles, iead_b, number_histories=num_hist)
        sputtered_bb, reflected_bb = bb.run_tridyn_simulations_from_iead(energies, angles, iead_b, number_histories=num_hist)
        sputtered_bc, reflected_bc = bc.run_tridyn_simulations_from_iead(energies, angles, iead_d, number_histories=num_hist)

        sputtered_dw, reflected_dw = dw.run_tridyn_simulations_from_iead(energies, angles, iead_d, number_histories=num_hist)
        sputtered_db, reflected_db = db.run_tridyn_simulations_from_iead(energies, angles, iead_d, number_histories=num_hist)
        sputtered_dc, reflected_dc = dc.run_tridyn_simulations_from_iead(energies, angles, iead_d, number_histories=num_hist)

        output = [sputtered_bw, reflected_bw, sputtered_bb, reflected_bb,
            sputtered_bc, reflected_bc, sputtered_dw, reflected_dw,
            sputtered_dw, reflected_dw, sputtered_db, reflected_db,
            sputtered_dc, reflected_dc]

        output_file = open('output.sav', 'wb')
        pickle.dump(output, output_file)
        output_file.close()

    print('incident boron: ', num_incident_b)
    print('incident plasma: ', num_incident_d)

    print('rbb: ', len(reflected_bb)/num_hist/num_incident_b)
    print('rbw: ', len(reflected_bw)/num_hist/num_incident_b)
    print('rbc: ', len(reflected_bc)/num_hist/num_incident_b)

    y_b = (len(sputtered_bb) + len(sputtered_db))/num_hist/(num_incident_b + num_incident_d)
    y_w = (len(sputtered_bw) + len(sputtered_dw))/num_hist/(num_incident_b + num_incident_d)
    y_c = (len(sputtered_bc) + len(sputtered_dc))/num_hist/(num_incident_b + num_incident_d)

    y_bb = len(sputtered_bb)/num_hist/num_incident_b
    y_db = len(sputtered_db)/num_hist/num_incident_d

    y_bw = len(sputtered_bw)/num_hist/num_incident_b
    y_dw = len(sputtered_dw)/num_hist/num_incident_d

    y_bc = len(sputtered_bc)/num_hist/num_incident_b
    y_dc = len(sputtered_dc)/num_hist/num_incident_d

    print(f'sb: {y_b}, y_bb: {y_bb}, y_db: {y_db}')
    print(f'sb: {y_c}, y_bc: {y_bc}, y_dc: {y_dc}')
    print(f'sb: {y_w}, y_bw: {y_bw}, y_dw: {y_dw}')

    print('rdb: ', len(reflected_db)/num_hist/num_incident_d)
    print('rdw: ', len(reflected_dw)/num_hist/num_incident_d)
    print('rdc: ', len(reflected_dc)/num_hist/num_incident_d)

    max_energy_reflected = 200.0
    max_energy_deuterium_reflected = 100.0
    max_energy_sputtering = 100.0
    bins = (80, 90)
    #boron

    #boron wall boron reflected
    plt.figure(figsize=(8, 6), dpi=300.0)
    boron_out = reflected_bb
    boron_out = np.array(boron_out)
    #boron_out = np.loadtxt('boron_out.dat')
    energy = boron_out[:,0]
    angle = np.arccos(np.abs(boron_out[:,1]))*180.0/np.pi
    #iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_reflected),(0.0, 90.0)))
    #print(energies, angles)
    #plt.pcolormesh(angles[:-1], energies[:-1], iead_out)

    #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    #plt.clabel(c, inline=1, fontsize=9)
    energy_dist_out, energy_bins = np.histogram(energy, bins=bins[0], range=(0., max_energy_reflected), density=True)
    de = energy_bins[1] - energy_bins[0]
    energies = energy_bins + de/2.
    np.savetxt('energy_dist_b_b_reflected.txt', np.array([energies[:-1], energy_dist_out/np.max(energy_dist_out)]).transpose())
    plt.plot(energy_bins[:-1], energy_dist_out/np.max(energy_dist_out))
    plt.axis([energy_bins[0], energy_bins[-1], 0.0, 1.2])
    plt.xticks(np.arange(0., max_energy_reflected+1, 50))
    plt.yticks(np.arange(0., 1.1, 0.5))
    plt.title('Boron Wall Boron Reflected Energies', fontdict=title_font)
    #plt.xlabel('angle [deg]')
    plt.xlabel('Energy [eV]')
    plt.ylabel('f(E)')
    plt.savefig('bb_reflected_2.png', format='png')

    #boron wall sputtered
    plt.figure(figsize=(8, 6), dpi=300.0)
    boron_out = sputtered_bb + sputtered_db
    boron_out = np.array(boron_out)
    energy = boron_out[:,0]
    angle = np.arccos(np.abs(boron_out[:,1]))*180.0/np.pi
    #iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_sputtering),(0.0, 90.0)))
    #print(energies, angles)
    #plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    #plt.clabel(c, inline=1, fontsize=9)
    energy_dist_out, energy_bins = np.histogram(energy, bins=bins[0], range=(0., max_energy_sputtering), density=True)
    de = energy_bins[1] - energy_bins[0]
    energies = energy_bins + de/2.
    np.savetxt('energy_dist_b_b_sputtered.txt', np.array([energies[:-1], energy_dist_out/np.max(energy_dist_out)]).transpose())
    plt.plot(energy_bins[:-1], energy_dist_out/np.max(energy_dist_out))
    plt.axis([energy_bins[0], energy_bins[-1], 0.0, 1.2])
    plt.xticks(np.arange(0., max_energy_sputtering+1, 50))
    plt.yticks(np.arange(0., 1.1, 0.5))
    plt.title('Boron Wall Boron Sputtered Energies', fontdict=title_font)
    #plt.xlabel('angle [deg]')
    plt.xlabel('Energy [eV]')
    plt.ylabel('f(E)')
    plt.savefig('b_sputtered_2.png', format='png')

    # #boron wall deuterium reflected
    # plt.figure(figsize=(8, 6), dpi=300.0)
    # deuterium_out = reflected_db
    # deuterium_out = np.array(deuterium_out)
    # energy = deuterium_out[:,0]
    # angle = np.arccos(np.abs(deuterium_out[:,1]))*180.0/np.pi
    # iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_deuterium_reflected),(0.0, 90.0)))
    # #print(energies, angles)
    # plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    # #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    # #plt.clabel(c, inline=1, fontsize=9)
    # plt.title('Boron Wall Deuterium Reflected IEAD')
    # plt.xlabel('angle [deg]')
    # plt.ylabel('energy [eV]')
    # plt.savefig('db_reflected_2.png', format='png')

    #boron total flux out
    plt.figure(figsize=(8, 6), dpi=300.0)
    boron_total = np.array(sputtered_bb + sputtered_bb + reflected_bb)
    energy = boron_total[:, 0]
    angle = np.arccos(np.abs(boron_total[:,1]))*180.0/np.pi
    #iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_reflected),(0.0, 90.0)))
    #print(energies, angles)
    #plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    #plt.clabel(c, inline=1, fontsize=9)''
    energy_dist_out, energy_bins = np.histogram(energy, bins=bins[0], range=(0., max_energy_reflected), density=True)
    de = energy_bins[1] - energy_bins[0]
    energies = energy_bins + de/2.
    np.savetxt('energy_dist_b_b_total.txt', np.array([energies[:-1], energy_dist_out/np.max(energy_dist_out)]).transpose())
    plt.plot(energy_bins[:-1], energy_dist_out/np.max(energy_dist_out))
    plt.axis([energy_bins[0], energy_bins[-1], 0.0, 1.2])
    plt.xticks(np.arange(0., max_energy_reflected+1, 50))
    plt.yticks(np.arange(0., 1.1, 0.5))
    plt.title('Boron Wall Boron Total Emitted Energies', fontdict=title_font)
    #plt.xlabel('angle [deg]')
    plt.xlabel('Energy [eV]')
    plt.ylabel('f(E)')
    plt.savefig('b_total_2.png', format='png')

    #tungsten

    #tungsten wall boron reflected
    plt.figure(figsize=(8, 6), dpi=300.0)
    tungsten_out = reflected_bw
    tungsten_out = np.array(tungsten_out)
    #boron_out = np.loadtxt('boron_out.dat')
    energy = tungsten_out[:,0]
    angle = np.arccos(np.abs(tungsten_out[:,1]))*180.0/np.pi
    #iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_reflected),(0.0, 90.0)))
    #print(energies, angles)
    #plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    #plt.clabel(c, inline=1, fontsize=9)
    energy_dist_out, energy_bins = np.histogram(energy, bins=bins[0], range=(0., max_energy_reflected), density=True)
    de = energy_bins[1] - energy_bins[0]
    energies = energy_bins + de/2.
    np.savetxt('energy_dist_b_w_reflected.txt', np.array([energies[:-1], energy_dist_out/np.max(energy_dist_out)]).transpose())
    plt.plot(energy_bins[:-1], energy_dist_out/np.max(energy_dist_out))
    plt.axis([energy_bins[0], energy_bins[-1], 0.0, 1.2])
    plt.xticks(np.arange(0., max_energy_reflected+1, 50))
    plt.yticks(np.arange(0., 1.1, 0.5))
    plt.title('Tungsten Wall Boron Reflected Energies', fontdict=title_font)
    #plt.xlabel('angle [deg]')
    plt.xlabel('Energy [eV]')
    plt.ylabel('f(E)')
    plt.savefig('bw_reflected_2.png', format='png')

    # #tungsten wall sputtered
    # plt.figure(figsize=(8, 6), dpi=300.0)
    # tungsten_out = sputtered_bw + sputtered_dw
    # tungsten_out = np.array(tungsten_out)
    # energy = tungsten_out[:,0]
    # angle = np.arccos(np.abs(tungsten_out[:,1]))*180.0/np.pi
    # iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_sputtering),(0.0, 90.0)))
    # #print(energies, angles)
    # plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    # #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    # #plt.clabel(c, inline=1, fontsize=9)
    # plt.title('Tungsten Wall Sputtered EAD')
    # plt.xlabel('angle [deg]')
    # plt.ylabel('energy [eV]')
    # plt.savefig('w_sputtered_2.png', format='png')

    # #tungsten wall deuterium reflected
    # plt.figure(figsize=(8, 6), dpi=300.0)
    # deuterium_out = reflected_dw
    # deuterium_out = np.array(deuterium_out)
    # energy = deuterium_out[:,0]
    # angle = np.arccos(np.abs(deuterium_out[:,1]))*180.0/np.pi
    # iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_deuterium_reflected),(0.0, 90.0)))
    # #print(energies, angles)
    # plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    # #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    # #plt.clabel(c, inline=1, fontsize=9)
    # plt.title('Tungsten Wall Deuterium Reflected IEAD')
    # plt.xlabel('angle [deg]')
    # plt.ylabel('energy [eV]')
    # plt.savefig('dw_reflected_2.png', format='png')
    #
    # #tungsten total flux out
    # plt.figure(figsize=(8, 6), dpi=300.0)
    # tungsten_total = np.array(sputtered_bw + sputtered_bw + reflected_bw)
    # energy = tungsten_total[:, 0]
    # angle = np.arccos(np.abs(tungsten_total[:,1]))*180.0/np.pi
    # iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_reflected),(0.0, 90.0)))
    # #print(energies, angles)
    # plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    # #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    # #plt.clabel(c, inline=1, fontsize=9)
    # plt.title('Tungsten Wall Sputtered and Reflected EAD')
    # plt.xlabel('angle [deg]')
    # plt.ylabel('energy [eV]')
    # plt.savefig('w_total_2.png', format='png')

    #carbon

    #carbon wall boron reflected
    plt.figure(figsize=(8, 6), dpi=300.0)
    boron_out = reflected_bc
    boron_out = np.array(boron_out)
    #boron_out = np.loadtxt('boron_out.dat')
    energy = boron_out[:,0]
    angle = np.arccos(np.abs(boron_out[:,1]))*180.0/np.pi
    #iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_reflected),(0.0, 90.0)))
    #print(energies, angles)
    #plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    #plt.clabel(c, inline=1, fontsize=9)
    energy_dist_out, energy_bins = np.histogram(energy, bins=bins[0], range=(0., max_energy_reflected), density=True)
    de = energy_bins[1] - energy_bins[0]
    energies = energy_bins + de/2.
    np.savetxt('energy_dist_b_c_reflected.txt', np.array([energies[:-1], energy_dist_out/np.max(energy_dist_out)]).transpose())
    plt.plot(energy_bins[:-1], energy_dist_out/np.max(energy_dist_out))
    plt.axis([energy_bins[0], energy_bins[-1], 0.0, 1.2])
    plt.xticks(np.arange(0., max_energy_reflected+1, 50))
    plt.yticks(np.arange(0., 1.1, 0.5))
    plt.title('Carbon Wall Boron Reflected Energies', fontdict=title_font)
    #plt.xlabel('angle [deg]')
    plt.xlabel('Energy [eV]')
    plt.ylabel('f(E)')
    plt.savefig('bc_reflected_2.png', format='png')

    # #carbon wall sputtered
    # plt.figure(figsize=(8, 6), dpi=300.0)
    # boron_out = sputtered_bc + sputtered_dc
    # boron_out = np.array(boron_out)
    # energy = boron_out[:,0]
    # angle = np.arccos(np.abs(boron_out[:,1]))*180.0/np.pi
    # iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_sputtering),(0.0, 90.0)))
    # #print(energies, angles)
    # plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    # #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    # #plt.clabel(c, inline=1, fontsize=9)
    # plt.title('Carbon Wall Sputtered EAD')
    # plt.xlabel('angle [deg]')
    # plt.ylabel('energy [eV]')
    # plt.savefig('c_sputtered_2.png', format='png')

    # #carbon wall deuterium reflected
    # plt.figure(figsize=(8, 6), dpi=300.0)
    # deuterium_out = reflected_dc
    # deuterium_out = np.array(deuterium_out)
    # energy = deuterium_out[:,0]
    # angle = np.arccos(np.abs(deuterium_out[:,1]))*180.0/np.pi
    # iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_deuterium_reflected),(0.0, 90.0)))
    # #print(energies, angles)
    # plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    # #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    # #plt.clabel(c, inline=1, fontsize=9)
    # plt.title('Carbon Wall Deuterium Reflected IEAD')
    # plt.xlabel('angle [deg]')
    # plt.ylabel('energy [eV]')
    # plt.savefig('dc_reflected_2.png', format='png')

    # #carbon total flux out
    # plt.figure(figsize=(8, 6), dpi=300.0)
    # carbon_total = np.array(sputtered_bc + sputtered_bc + reflected_bc)
    # energy = carbon_total[:, 0]
    # angle = np.arccos(np.abs(carbon_total[:,1]))*180.0/np.pi
    # iead_out, energies, angles = np.histogram2d(energy, angle, bins=bins, range=((0.0, max_energy_reflected),(0.0, 90.0)))
    # #print(energies, angles)
    # plt.pcolormesh(angles[:-1], energies[:-1], iead_out)
    # #plt.contourf(angles[:-1], energies[:-1], iead_out/np.max(iead_out), levels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # #c = plt.contour(angles[:-1], energies[:-1], iead_out/np.max(iead_out), colors='black', levels=[0.2, 0.4, 0.6, 0.8, 1.0])
    # #plt.clabel(c, inline=1, fontsize=9)
    # plt.title('Carbon Wall Sputtered and Reflected EAD')
    # plt.xlabel('angle [deg]')
    # plt.ylabel('energy [eV]')
    # plt.savefig('c_total_2.png', format='png')

    os.system('rm *.IN *.OUT')
    #plt.show()

def boron_carbon():
    boron_concs = np.linspace(0.0, 1.0, 10)
    carbon_concs = [1. - x for x in boron_concs]
    number_histories = 100000

    beam_species_list = ['He', 'Ar']
    total_yield_array = np.zeros((len(beam_species_list), len(boron_concs)))
    carbon_yield_array = np.zeros((len(beam_species_list), len(boron_concs)))
    boron_yield_array = np.zeros((len(beam_species_list), len(boron_concs)))

    for beam_index, beam_species in enumerate(beam_species_list):
        for simulation_index, (boron_conc, carbon_conc) in enumerate(zip(boron_concs, carbon_concs)):
            print(f'beam: {beam_species} {simulation_index}')
            os.system('rm *.DAT')
            preferential_sputtering_boron_carbon(boron_conc, carbon_conc, beam_species=beam_species, number_histories=number_histories, sim_number=simulation_index)
            os.system(f'./FTridyn_Clean < CONC{str(simulation_index).zfill(4)}.IN')

            sputtered_list = np.atleast_2d(np.genfromtxt('CONCSPLST.DAT'))

            sputtered_id = sputtered_list[:, 1]
            number_boron_sputtered = np.sum(sputtered_id==2)
            number_carbon_sputtered = np.sum(sputtered_id==3)
            total_sputtered = number_boron_sputtered + number_carbon_sputtered
            print(f'{total_sputtered} {len(sputtered_list)}')
            total_yield_array[beam_index, simulation_index] = total_sputtered / number_histories
            carbon_yield_array[beam_index, simulation_index] = number_carbon_sputtered / number_histories
            boron_yield_array[beam_index, simulation_index] = number_boron_sputtered / number_histories

            print(f'B:C: {boron_conc / carbon_conc}')
            #print(f'Y_C: {number_carbon_sputtered/total_sputtered} Y_B: {number_boron_sputtered/total_sputtered}')
            print(f'Y_total: {total_sputtered/number_histories}')

    #plt.figure(figsize=(8, 6), dpi=300.0))
    #plt.plot(boron_concs, total_yield_array[0, :])
    #plt.plot(boron_concs, total_yield_array[1, :])
    #plt.plot(boron_concs, total_yield_array[2, :])
    #plt.legend([f'beam: {beam_species}' for beam_species in beam_species_list])
    #plt.title('Total Sputtering Yield of B-C Targets (E=80eV, alpha=86deg)')
    #plt.xlabel('Boron concentration')
    #plt.ylabel('Y (at./ion)')

    for beam_species_index, beam_species in enumerate(beam_species_list):
        figure, axis1 = plt.subplots()
        axis2 = axis1.twinx()
        plt.title(f'Beam: {beam_species}, E=2.5keV, alpha=0')
        axis1.plot(boron_concs, carbon_yield_array[beam_species_index, :], c='black')
        axis1.plot(boron_concs, boron_yield_array[beam_species_index, :], c='red')
        axis1.plot(boron_concs, total_yield_array[beam_species_index, :], c='blue')
        axis2.plot(boron_concs ,boron_yield_array[beam_species_index, :] * carbon_concs / carbon_yield_array[beam_species_index, :] / boron_concs, ':')
        axis1.legend(['Y_C', 'Y_B', 'Y_tot'], loc='lower left')
        axis2.legend(['r'], loc='lower right')
        plt.xlabel('Boron concentration')
        axis1.set_ylabel('Y [at./ion]')
        axis2.set_ylabel('Preferential Sputtering Parameter')
        axis2.set_ylim((0.0, 2.0))

    plt.show()

def boron_carbon_tungsten_aps():
    energies = np.linspace(10.0, 300.0, 25)
    angles = np.linspace(0.0, 89.0, 25)

    beam_species = ['B', 'D']
    beam_names = [ 'B', 'D']

    target_species = ['B']
    target_names = ['B']

    number_histories = 50000

    for beam_index, beam in enumerate(beam_species):
        for target_index, target in enumerate(target_species):
            print(f'beam: {beam}')
            print(f'target: {target}')

            energy_angle_sputtered = np.zeros((len(energies), len(angles)))
            energy_angle_reflected = np.zeros((len(energies), len(angles)))
            energy_angle_deposited = np.zeros((len(energies), len(angles)))

            sim_number = 0

            for energy_index, energy in enumerate(energies):

                print(f'progress: {np.round(energy_index/len(energies),2)*100}%')

                for angle_index, angle in enumerate(angles):

                    sim_number += 1

                    sim_params = beam_and_target(
                        beam_names[beam_index] + '_' + target_names[target_index] + '_',
                        beam, target, sim_number=1,
                        number_histories=number_histories,
                        incident_energy=energy, incident_angle=angle,
                        depth=1000
                        )

                    input_file = sim_params.name+str(sim_params.sim_num).zfill(4)+'.IN'

                    os.system(f'./FTridyn_Clean < {input_file} > /dev/null 2>&1')

                    sputtered = np.genfromtxt(sim_params.name+'SPLST.DAT')
                    reflected = np.genfromtxt(sim_params.name+'RFLST.DAT')

                    num_sputtered = len(sputtered)
                    num_reflected = len(reflected)

                    Y = num_sputtered / number_histories
                    R = num_reflected / number_histories
                    D = (number_histories - num_reflected) / number_histories

                    energy_angle_sputtered[energy_index, angle_index] = Y
                    energy_angle_reflected[energy_index, angle_index] = R
                    energy_angle_deposited[energy_index, angle_index] = D

            np.savetxt(sim_params.name+'_yield.sav', energy_angle_sputtered)
            np.savetxt(sim_params.name+'_reflected.sav', energy_angle_reflected)
            np.savetxt(sim_params.name+'_deposited.sav', energy_angle_deposited)


            levels = np.arange(0, 2.1, 0.25)
            fig = plt.figure(figsize=(8, 6), dpi=300.0)
            plt.contourf(energies, angles, energy_angle_sputtered.transpose(), levels=levels, antialiased=False)
            contour = c = plt.contour(energies, angles, energy_angle_sputtered.transpose(), colors='black', levels=levels)
            plt.title(f'{beam_names[beam_index]} on {target_names[target_index]} Sputtering Yield')
            plt.clabel(contour, fontsize=9, inline=True)
            plt.xlabel('Energy [eV]')
            plt.ylabel('Angle [deg]')
            plt.savefig(f'{beam_names[beam_index]}_{target_names[target_index]}_Y')

            levels = np.arange(0, 1.1, 0.2)
            fig = plt.figure(figsize=(8, 6), dpi=300.0)
            plt.contourf(energies, angles, energy_angle_reflected.transpose(), levels=levels, antialiased=False)
            contour = c = plt.contour(energies, angles, energy_angle_reflected.transpose(), colors='black', levels=levels)
            plt.title(f'{beam_names[beam_index]} on {target_names[target_index]} Reflection Coefficient')
            plt.clabel(contour, fontsize=9, inline=True)
            plt.xlabel('Energy [eV]')
            plt.ylabel('Angle [deg]')
            plt.savefig(f'{beam_names[beam_index]}_{target_names[target_index]}_R')

            levels = np.arange(0, 1.1, 0.2)
            fig = plt.figure(figsize=(8, 6), dpi=300.0)
            plt.contourf(energies, angles, energy_angle_deposited.transpose(), levels=levels, antialiased=False)
            contour = c = plt.contour(energies, angles, energy_angle_deposited.transpose(), colors='black', levels=levels)
            plt.title(f'{beam_names[beam_index]} on {target_names[target_index]} Stopped Fraction')
            plt.clabel(contour, fontsize=9, inline=True)
            plt.xlabel('Energy [eV]')
            plt.ylabel('Angle [deg]')
            plt.savefig(f'{beam_names[beam_index]}_{target_names[target_index]}_D')

            os.system('rm sim_params.name*.IN')
            os.system('rm *.OUT')
            os.system('rm *.DAT')

if __name__ == '__main__':
    #boron_carbon_tungsten_aps()
    #boron_carbon_tungsten_aps()
    main()
