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
2. Install rustup, the Rust toolchain (includes rustc, the compiler, and cargo, the package manager)
3. Run python3.# -m pip install numpy matplotlib shapely
4. Install toml from source:
  * git clone https://github.com/uiri/toml.git
  * cd toml
  * python3.# setup.py install
5. Build RustBCA
  * git clone https://github.com/drobnyjt/rustBCA
  * cd rustBCA
  * cargo build --release
 6. input.toml is the input file -- example included below
 7. cargo run --relase
 8. Output files are .output

# Example input file: input.toml

'[options]
name = "test_" #Name of output file
track_trajectories = false #Track incident ion trajectories (memory intensive)
track_recoils = false #Track recoils
track_recoil_trajectories = false #Track recoil trajectories (memory intensive)
write_files = true #Write output files
[material_parameters]
energy_unit = "EV"
mass_unit = "AMU"
Eb = 1.0 #Bulk binding energy
Es = 3.52 #Surface binding energy
Ec = 3.52 #Cutoff energy
n = 8.4e+28 #Number density
Z = 29 
m = 63.54
[geometry] #Must define atomic geometry, geometry representing surface energy barrier, and simulation boundary
length_unit = "MICRON"
surface = [ [ 20.0, -0.5,], [ 20.0, 0.5,], [ 0.0, 0.5,], [ 0.0, -0.5,], [ 20.0, -0.5,],]
energy_surface = [ [ 20.000182185426144, -0.5001821854261432,], [ -0.00018218542614322072, -0.5001821854261432,], [ -0.00018218542614322072, 0.5001821854261432,], [ 20.000182185426144, 0.5001821854261432,], [ 20.000182185426144, -0.5001821854261432,],]
simulation_surface = [ [ 20.000364370852285, -0.5003643708522865,], [ -0.0003643708522864414, -0.5003643708522865,], [ -0.0003643708522864414, 0.5003643708522865,], [ 20.000364370852285, 0.5003643708522865,], [ 20.000364370852285, -0.5003643708522865,],]
[particle_parameters]
length_unit = "MICRON"
energy_unit = "EV"
mass_unit = "AMU"
N = [ 100,] #Number of particles from this column to run
m = [ 4,]
Z = [ 2,]
E = [ 5000000.0,] #Energy in energy_units
pos = [ [ -0.00020040396875754278, -0.27233199674655106, 0.0,],] #Initial position
dir = [ [ 0.999999835859687, 0.0005729577637823323, 0.0,],] #Initial direction unit vector'
 
# Usage

Modify input.toml as wanted; run ./rustBCA with input.toml in the same directory as rustBCA
