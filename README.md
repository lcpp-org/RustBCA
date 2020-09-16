# rustBCA

A Binary Collision Approximation (BCA) code for ion-material interactions, written in Rust!

BCA codes are valid for incident ion energies between approximately 10 eV  through GeV. By discretizing the collision cascade into a sequence of binary collisions, BCA codes can accurately and efficiently model reflection, implantation, sputtering, transmission, and displacement damage.

# Features

* Ion-solid interactions for all combinations of incident ion and target species
* 2D geometry
* Arbitrary, mesh-based inhomogeneous composition
* Low energy (<25 keV/nucleon) electronic stopping modes include local (Oen-Robinson), nonlocal (Lindhard-Scharff), and equipartition forms
* Includes Biersack-Varelas interpolation to extend electronic stopping validity up to ~1 GeV/nucleon
* Includes MAGIC algorithm and Mendenhall-Weller quadrature to solve scattering integral
* Full trajectory tracking of ions and target atoms
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
name = "10.0eV_0.0001deg_Ne_TiO2_Al_Si_" #Output file name
track_trajectories = false #Track incident ion trajectories
track_recoils = true #Track target atoms
track_recoil_trajectories = false #Track target atom trajectories
write_files = true #Write output files to disk
stream_size = 8000 #Size of streaming write in bytes
print = true #Print to std out
print_num = 10 #Number of times to print to std out
weak_collision_order = 3 #Number of weak collisions (>0)
suppress_deep_recoils = false #Eckstein's empirical recoil suppression for target atoms that cannot sputter
high_energy_free_flight_paths = false #TRIM-style free flight paths between collisions
electronic_stopping_mode = 1 #0: Biersack-Varelas Interpolation (~10eV/nucleon-100MeV/nucleon) 1: Lindhard-Scharff nonlocal (<=25keV/nucleon) 2: Oen-Robinson local (<=25keV/nucleon) 3: Equipartition betw. LS and OR (<=25keV/nucleon)
mean_free_path_model = 0 #0: Liquid (1/n^3) 1: Gaseous (exponential)
interaction_potential = 2 #0: Moliere 1: Krypton-Carbon 2: ZBL 3: Lenz-Jensen
scattering_integral = 0 #0: Mendenhall-Weller Quadrature 1: MAGIC algorithm
tolerance = 0.001 #Tolerance for Newton-Raphson distance of closest approach algorithm
max_iterations = 100 #Max iterations for Newton-Raphson distance of closest approach algorithm
#Options below refer to optional HDF5, Chebyshev Proxy Rootfinder, and Polynomial Rootfinder options
use_hdf5 = false
root_finder = 0
num_threads = 4
num_chunks = 10
cpr_n0 = 0
cpr_nmax = 0
cpr_epsilon = 0
cpr_complex = 0
cpr_truncation = 0
cpr_far_from_zero = 0
cpr_interval_limit = 0
cpr_upper_bound_const = 0
polynom_complex_threshold = 0

[material_parameters]
energy_unit = "EV" #energy unit, one of: EV, J, KEV, MEV
mass_unit = "AMU" #mass unit, one of: AMU, KG
Eb = [ 3.0, 2.58, 3.0, 0.0,] #Binding energy of each species
Es = [ 4.84, 2.58, 3.39, 4.72,] #Surface binding energy of each species
Ec = [ 3.5, 2.0, 3.0, 1.5,] #Cutoff energy of each species
n = [ 5.67e+28, 4.291e+28, 6.026e+28, 4.996e+28,] #Nominal number density of each species
Z = [ 22, 8, 13, 14,] #Atomic number of each species
m = [ 47.867, 15.9994, 26.98, 28.08553,] #Atomic mass of each species
electronic_stopping_correction_factor = 1.0 #Correction factor for electronic stopping, used for Z1 oscillations when known
energy_barrier_thickness = 6.713580575220781e-5 #Distance from surface that surface binding energy is computed (typically ~mfp)

[particle_parameters]
length_unit = "MICRON" #length unit, one of ANGSTROM, NM, MICRON, CM, M
energy_unit = "EV" #energy unit, one of: EV, J, KEV, MEV
mass_unit = "AMU" #mass unit, one of: AMU, KG
N = [ 100000,] #Number of particles in column 0
m = [ 20.1797,] #Atomic mass of each particle species
Z = [ 10,] #Atomic number of each particle species
E = [ 10.0,] #Initial energy of each particle species
Ec = [ 1.0,] #Cutoff energy of each particle species
Es = [ 0.0,] #Surface binding energy of each particle species
pos = [ [ -1.7453292519934434e-8, 0.0, 0.0,],] #Initial position (x, y, z) of each particle species
dir = [ [ 0.9999999999984769, 1.7453292519934434e-6, 0.0,],] #Initial directional cosines (ux, uy, uz) of each particle species

[mesh_2d_input]
length_unit = "MICRON"
#List of triangles that make up mesh
coordinate_sets = [ [ 0.0, 0.1, 0.0, 0.5, 0.5, -0.5,], [ 0.0, 0.1, 0.1, -0.5, 0.5, -0.5,], [ 0.1, 0.1, 0.4, -0.5, 0.5, -0.5,], [ 0.1, 0.4, 0.4, 0.5, 0.5, -0.5,], [ 0.4, 0.5, 0.4, 0.5, 0.5, -0.5,], [ 0.4, 0.5, 0.5, -0.5, 0.5, -0.5,],]
densities = [ [ 3e+28, 6e+28, 0.0, 0.0,], [ 3e+28, 6e+28, 0.0, 0.0,], [ 0.0, 0.0, 6.026e+28, 0.0,], [ 0.0, 0.0, 6.026e+28, 0.0,], [ 0.0, 0.0, 0.0, 4.996e+28,], [ 0.0, 0.0, 0.0, 4.996e+28,],]
#Boundary of mesh
boundary_points = [ [ 0.0, 0.5,], [ 0.5, 0.5,], [ 0.5, -0.5,], [ 0.0, -0.5,],]
#Simulation boundary
simulation_boundary_points = [ [ 0.6, -0.6,], [ -0.1, -0.6,], [ -0.1, 0.6,], [ 0.6, 0.6,], [ 0.6, -0.6,],]

 ~~~~
# Usage

Modify input.toml to configure simulation; run ./rustBCA with input.toml in the same directory as rustBCA

Recommended usage: use run_rustbca.py to generate input files programmatically using toml
