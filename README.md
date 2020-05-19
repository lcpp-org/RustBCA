# rustBCA

A Binary Collision Approximation (BCA) code for ion-material interactions, written in Rust!

BCA codes are valid for incident ion energies between approximately 10 eV  through GeV. By discretizing the collision cascade into a sequence of binary collisions, BCA codes can accurately and efficiently model reflection, implantation, sputtering, transmission, and displacement damage.

# Features

* Ion-solid interactions for all combinations of incident ion and target species
* Includes homogeneous multicomponent materials
* Low energy (<25 keV/nucleon) electronic stopping modes include local (Oen-Robinson), nonlocal (Lindhard-Scharff), and equipartition forms
* Includes Biersack-Varelas interpolation to extend electronic stopping validity up to ~1 GeV/nucleon
* Includes MAGIC algorithm and more accurate Mendenhall-Weller quadrature to solve scattering integral
* Arbitrary 2D geometry
* Full trajectory tracking
* Human-readable input file using the TOML format
* Output of energies and directions of emitted particles (reflected ions and sputtered atoms)
* Output of final positions of implanted ions
* Output of full trajectory tracking for incident ions and target atoms

# Installation

Dependences:
* Python 3.6+
* Numpy, Matplotlib, toml, Shapely, scipy
* rustup

Ubuntu 18.04 LTS:
1. Make sure you have Python 3.6+ (this comes natively in Ubuntu 18.04)
2. Install curl, `apt-get install curl`
3. Install rustup, the Rust toolchain (includes rustc, the compiler, and cargo, the package manager) from https://rustup.rs/
4. Install pip for Python-3, `apt-get install python3-pip`
5. Install Python librarires, `python3.# -m pip install numpy matplotlib shapely scipy`
6. Install Python toml library from source:
- git clone https://github.com/uiri/toml.git
- cd toml
- python3.# setup.py install
7. Install Cargo, `apt install cargo`
8. Build RustBCA
- git clone https://github.com/lcpp-org/rustBCA
- cd rustBCA
- cargo build --release
9. input.toml is the input file -- example included below
10. cargo run --relase
11. Output files are .output

# Example input file: input.toml
~~~~
[options]
name = "10.0eV_0.0001deg_He_N" #Output file name
track_trajectories = false #Track ion trajectories (memory intensive)
track_recoils = false #Track recoils
track_recoil_trajectories = false #Track recoil trajectories (memory intensive)
write_files = true #Write output files to disk
stream_size = 8000  #Streaming write packet size in bytes
print = true #Print to command line
print_num = 10 #Number of times to print to command line
weak_collision_order = 0 #Number of max weak collisions
suppress_deep_recoils = false #Use Eckstein's deep recoil supression routine
high_energy_free_flight_paths = true #Use TRIM-style high-energy free-flight paths between ions
electronic_stopping_mode = 0 #0: Interpolated, Biersack-Varelas (eV-GeV) 1: Low energy nonlocal, Lindhard-Scharff (eV-25keV/amu) 2: Low energy local, Oen-Robinson (eV-25 keV/amu), 3: Equipartition between LS and OR
mean_free_path_model = 1 #0: Liquid/solid target, 1: Gaseous target
interaction_potential = 2 #0: Moliere 1: Kr-C 2: ZBL
scattering_integral = 0 #0: 6th order Lobatto Quadrature 1: MAGIC (Recommended: 0, don't use 1 unless you need to)
tolerance = 1e-6 #tolerance on distance of closest approach algorithm
max_iterations = 100 #maximum number of newton-raphson iterations for distance of closest approach calculation

[material_parameters]
energy_unit = "EV"
mass_unit = "AMU"
Eb = [0.0] #Bulk binding energy
Es = [0.0] #Surface binding energy
Ec = [1.0] #Cutoff energy
n = [5.4e+25] #Number density 
Z = [7]
m = [14]
electronic_stopping_correction_factor = 1.0 #C_k, used to compensate for Z1 oscillations in Lindhard-Scharff when known

[geometry]
length_unit = "MICRON"
surface = [ [ 10.0, -5.0,], [ 10.0, 5.0,], [ 0.0, 5.0,], [ 0.0, -5.0,], [ 10.0, -5.0,],]
energy_surface = [ [ 10.0, -5.0,], [ 0.0, -5.0,], [ 0.0, 5.0,], [ 10.0, 5.0,], [ 10.0, -5.0,],]
simulation_surface = [ [ 10.0, -5.0,], [ 0.0, -5.0,], [ 0.0, 5.0,], [ 10.0, 5.0,], [ 10.0, -5.0,],]

[particle_parameters]
length_unit = "MICRON"
energy_unit = "EV"
mass_unit = "AMU"
N = [ 10000,] #Number of particles for first column
m = [ 4.002602,]
Z = [ 2,]
E = [ 10.0,] #Incident energy
Ec = [ 0.01,] #Cutoff energy
Es = [ 0.0,] #Surface binding energy
pos = [ [ -0.0, 0.0, 0.0,],]
dir = [ [ 0.9999999999984769, 1.7453292519934434e-6, 0.0,],] #Note: because of gimbal lock, dir[0] != 0

 ~~~~
# Usage

Modify input.toml to configure simulation; run ./rustBCA with input.toml in the same directory as rustBCA

Recommended usage: use run_rustbca.py to generate input files programmatically using toml
