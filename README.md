# rustBCA

A Binary Collision Approximation (BCA) code for ion-solid interactions, written in Rust!

BCA codes are valid for incident ion energies between approximately 10 eV through MeV (or up to GeV for light ions). By discretizing the collision cascade into a sequence of binary collisions, BCA codes can accurately model reflection, implantation, sputtering, transmission, and displacement damage.

# Installation

Dependences:
* Python 3.6+
* Numpy, Matplotlib, toml, Shapely
* rustup

Steps to install:
1. Install Python 3.6+.
2. Install rustup, the Rust toolchain (includes rustc, the compiler, and cargo, the package manager) from https://rustup.rs/
3. Run python3.# -m pip install numpy matplotlib shapely
4. Install toml from source:
  * git clone https://github.com/uiri/toml.git
  * cd toml
  * python3.# setup.py install
5. Build RustBCA
  * git clone https://github.com/lcpp-org/rustBCA
  * cd rustBCA
  * cargo build --release
 6. input.toml is the input file -- example included below
 7. cargo run --relase
 8. Output files are .output

# Example input file: input.toml
~~~~
[options]
name = "test_He_Cu" #Output file name
track_trajectories = false #Tracking of full ion trajectories
track_recoils = true #Tracking of generated recoils
track_recoil_trajectories = false #Tracking of generated recoil trajectories
write_files = true #Write files to disk
stream_size = 256000 #Size in bytes of disk-writing chunk
print = true #Print to console
print_num = 100 #Number of statements to print to console
weak_collision_order = 3 #Number of weak collisions
suppress_deep_recoils = false #Eckstein's empirical deep recoil suppression - see Eq. 7.5.3 in Eckstein (1991)
high_energy_free_flight_paths = false #Use TRIM-style free flight path - must use electronic stopping mode 0
electronic_stopping_mode = 1 # 0: High energy (Biersack-Varelas interpolation of Lindhard-Scharff and Bethe-Bloch) 1: Low energy nonlocal (Lindhard-Scharff) 2: Low energy local (Oen-Robinson - currently in testing) 3: Low energy equipartition (LS + OR - currently in testing)

[material_parameters]
energy_unit = "EV"
mass_unit = "AMU"
Eb = 0.0 #Bulk binding energy
Es = 3.52 #Surface binding energy
Ec = 3.0 #Cutoff energy
n = 8.491e+28 #Atomic density
Z = 29 #Atomic number
m = 63.546 #Mass number

[geometry]
length_unit = "MICRON"
surface = [ [ 1.0, -0.5,], [ 1.0, 0.5,], [ 0.0, 0.5,], [ 0.0, -0.5,], [ 1.0, -0.5,],]
#Actual geometry of surface
energy_surface = [ [ 1.0001815322460903, -0.5001815322460903,], [ -0.00018153224609026986, -0.5001815322460903,], [ -0.00018153224609026986, 0.5001815322460903,], [ 1.0001815322460903, 0.5001815322460903,], [ 1.0001815322460903, -0.5001815322460903,],]
#Geometry of surface boundary potential barrier
simulation_surface = [ [ 1.0003630644921806, -0.5003630644921805,], [ -0.0003630644921805397, -0.5003630644921805,], [ -0.0003630644921805397, 0.5003630644921805,], [ 1.0003630644921806, 0.5003630644921805,], [ 1.0003630644921806, -0.5003630644921805,],]
#Simulation boundary

[particle_parameters]
length_unit = "MICRON"
energy_unit = "EV"
mass_unit = "AMU"
N = [ 10000,] #Number of particles in each column
m = [ 4.002602,] #Mass of particles
Z = [ 2,] #Atomci number
E = [ 1000.0,] #Incident energy
Ec = [ 0.1,] #Cutoff energy
Es = [ 0.0,] #Surface binding energy
pos = [ [ -0.00018, 0.0, 0.0,],] #Initial position
dir = [ [ 0.9999999999984769, 1.7453292519934434e-6, 0.0,],] #Initial velocity unit vectors

 ~~~~
# Usage

Modify input.toml as wanted; run ./rustBCA with input.toml in the same directory as rustBCA

Recommended usage: use run_rustbca.py to generate input files programmatically
