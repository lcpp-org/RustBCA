# rustBCA

A Binary Collision Approximation (BCA) code for ion-material interactions, written in Rust!

BCA codes are valid for incident ion energies between approximately 10 eV  through GeV. By discretizing the collision cascade into a sequence of binary collisions, BCA codes can accurately and efficiently model reflection, implantation, sputtering, transmission, and displacement damage.

# Features

* Ion-material interactions for all combinations of incident ion and target species
* 2D geometry
* Arbitrary, mesh-based inhomogeneous composition
* Amorphous solid/liquid and gaseous targets
* Low energy (<25 keV/nucleon) electronic stopping modes include local (Oen-Robinson), nonlocal (Lindhard-Scharff), and equipartition forms
* Includes Biersack-Varelas interpolation to extend electronic stopping validity up to ~1 GeV/nucleon
* High-energy free-flight paths to greatly speed up high-energy simulations
* Includes Kr-C, ZBL, Lenz-Jensen, Moliere, Lennard-Jones 12-6, and Lennard-Jones 6.5-6 interaction potentials
* Multiple interaction potentials can be used in a simulation - for example, the He-W interaction can be specified as a Lennard-Jones 12-6 while the W-W interaction can be specified as a Kr-C, for any number of potentials or species.
* Includes Gauss-Mehler quadrature, Gauss-Legendre quadrature, Mendenall-Weller quadrature, and MAGIC algorithm to determine scattering integral
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

# Usage

Modify input.toml to configure simulation; run ./rustBCA with input.toml in the same directory as rustBCA or run with a single input argument, the name of a .toml input file. See the rustBCA Wiki for example input files.

Recommended usage: use run_rustbca.py to generate input files programmatically using toml
