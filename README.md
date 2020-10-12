# rustBCA

A Binary Collision Approximation (BCA) code for ion-material interactions, written in Rust!

Check out the [rustBCA wiki](https://github.com/lcpp-org/RustBCA/wiki) for detailed information, installation instructions, usage instructions, examples, and more.

BCA codes are valid for incident ion energies between approximately 10 eV  through GeV. By discretizing the collision cascade into a sequence of binary collisions, BCA codes can accurately and efficiently model reflection, implantation, sputtering, transmission, and displacement damage.
# Features

* Ion-material interactions for all combinations of incident ion and target species
* 2D geometry with triangular-mesh-based inhomogeneous composition
* Amorphous solid/liquid and gaseous targets
* Low energy (<25 keV/nucleon) electronic stopping modes include local (Oen-Robinson), nonlocal (Lindhard-Scharff), and equipartition forms; Includes Biersack-Varelas interpolation to extend electronic stopping validity up to ~1 GeV/nucleon
* Includes optional Biersack-Haggmark high-energy free-flight path to greatly speed up high-energy simulations by neglecting very small angle scattering
* Includes Kr-C, ZBL, Lenz-Jensen, and Moliere screened-Coulomb potentials; Includes Lennard-Jones 12-6, Lennard-Jones 6.5-6, and Morse attractive-repulsive interaction potentials
* Includes Newton-Raphson for simple root-finding and optionally includes Adaptive Chebyshev Proxy Rootfinder with Automatic Subdivision and Polynomial root-finding algorithms through the [rcpr](https://github.com/drobnyjt/rcpr) crate for solving the distance-of-closest-approach problem
* Multiple interaction potentials can be used in a single simulation - for example, the He-W interaction can be specified as a Lennard-Jones 12-6 while the W-W interaction can be specified as a Kr-C, for any number of potentials or species
* Includes Gauss-Mehler quadrature, Gauss-Legendre quadrature, Mendenall-Weller quadrature, and the MAGIC algorithm to calculate the scattering integral
* Full trajectory tracking of ions and target atoms
* Human-readable input file using the TOML format
* User-friendly error messages pinpoint the cause of errors and suggest fixes to the user
* Output of energies and directions of emitted particles (reflected ions and sputtered atoms)
* Output of final positions of implanted ions
* Output of full trajectory tracking for incident ions and target atoms

# Installation

Without optional features, rustBCA should compile with cargo on Windows, MacOS, and Linux systems. HDF5 has been tested on Windows, but HDF5 1.10.6 must be used. [rcpr](https://github.com/drobnyjt/rcpr), the Rust adaptive Chebyshev Proxy Rootfinder with automatic subdivision and polynomial rootfinder package, has not yet been successfully compiled on Windows. However, it can be compiled on Windows Subsystem for Linux (WSL) and likely on Ubuntu for Windows or Cygwin.

Manual Dependences:
* [rustup](https://rustup.rs), the Rust toolchain (includes cargo, the Rust package manager, rustc, the Rust compiler, and more).

Automatic Dependencies:
* see [Cargo.toml](https://github.com/lcpp-org/RustBCA/blob/master/Cargo.toml) for a complete list

Optional Dependencies:
* HDF5 libraries
* `rcpr`: https://github.com/drobnyjt/rcpr the CPR and polynomial rootfinder, required for using attractive-repulsive interaction potentials such as Lennard-Jones or Morse, may require the following to be installed depending on the system:
  * gcc
  * build-essential
  * cmake
  * gfortran
  * liblapack-dev
  * libblas-dev
  * liblapacke-dev
* Python 3.6+
* Numpy, Matplotlib, toml, Shapely, scipy

Ubuntu 18.04 LTS:
1. Optional: install Python 3.6+ (this comes natively in Ubuntu 18.04)
2. Install curl, `apt-get install curl`
3. Install rustup, the Rust toolchain (includes rustc, the compiler, and cargo, the package manager) from https://rustup.rs/
4. Optional: Install pip for Python-3, `apt-get install python3-pip`
5. Optional: Install Python libraries for making input files, `python3 -m pip install numpy matplotlib shapely scipy`
6. Optional: Install Python toml library from source:
- `git clone https://github.com/uiri/toml.git`
- `cd toml`
- `python3 setup.py install`
7. Install Cargo, `apt install cargo`
8. Build RustBCA
- `git clone https://github.com/lcpp-org/rustBCA`
- `cd rustBCA`
- `cargo build --release`
9. Optional: Build RustBCA with optional dependencies, hdf5 and/or rcpr (with your choice of backed: openblas, netlib, or intel-mkl):
 - `cargo build --release --features cpr_rootfinder_netlib,hdf5_input`
 - `cargo build --release --features cpr_rootfinder_openblas,hdf5_input`
 - `cargo build --release --features cpr_rootfinder_intel_mkl,hdf5_input`
10. input.toml is the input file -- see [Usage](https://github.com/lcpp-org/RustBCA/wiki/Usage,-Input-File,-and-Output-Files) for more information
11. `cargo test` will run all required tests
12. Optional: `cargo test --features cpr_rootfinder_*` will run all required and optional tests for desired backend *
# Usage

Modify input.toml to configure simulation; run ./rustBCA with input.toml in the same directory as rustBCA or run with a single input argument, the name of a .toml input file. See [Usage](https://github.com/lcpp-org/RustBCA/wiki/Usage,-Input-File,-and-Output-Files) on the rustBCA wiki for usage instructions and see the examples on the wiki for example input files.

