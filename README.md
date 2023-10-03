# `RustBCA`

`RustBCA` is a general-purpose, high-performance code for simulating
ion-material interactions including sputtering, reflection, and implantation
using the binary collision approximation ([BCA]), written in [Rust]!
RustBCA includes a standalone version and libraries for including
ion-material interactions in simulations written in C/C++, Python,
and Fortran.

By discretizing the collision cascade into a sequence of binary collisions,
[BCA] codes can accurately and efficiently model the prompt interaction
between an energetic ion and a target material. This includes reflection,
implantation, and transmission of the incident ion, as well as sputtering
and displacement damage of the target. Generally, [BCA] codes can be
valid for incident ion energies between approximately ~1 eV/nucleon 
to <1 GeV/nucleon. Improvements to RustBCA have been shown to expand the 
regime of validity for some quantities even lower than 1 eV/nucleon.

Check out the `RustBCA` [Wiki] for detailed information, installation
instructions, use cases, examples, and more. See the RustBCA paper at the
Journal of Open Source Software by clicking the badge below:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03298/status.svg)](https://doi.org/10.21105/joss.03298)

Selected citations of RustBCA as of 5/24/23:
* [Simulation of liquid lithium divertor geometry using SOLPS-ITER](https://doi.org/10.1109/TPS.2022.3166402), JD Lore et al. (2022)
* [Characterizing W sources in the all-W wall, all-RF WEST tokamak environment](https://doi.org/10.1088/1361-6587/ac8acc), CC Klepper et al. (2022)
* [hPIC2: A hardware-accelerated, hybrid particle-in-cell code for dynamic plasma-material interactions](https://doi.org/10.1016/j.cpc.2022.108569), LT Meredith et al. (2023)
* [Global sensitivity analysis of a coupled multiphysics model to predict surface evolution in fusion plasma–surface interactions](https://doi.org/10.1016/j.commatsci.2023.112229), P. Robbe et al. (2023)
* [Modeling the effect of nitrogen recycling on the erosion and leakage of tungsten impurities from the SAS-VW divertor in DIII-D during nitrogen gas injection](https://doi.org/10.1016/j.nme.2022.101254), MS Parsons et al. (2023)

## Getting started

The easiest way to get started is with the ergonomic Python functions.
Follow these steps to install, build, and run simple RustBCA simulations
for sputtering yields and reflection coefficients:
```
git clone https://github.com/lcpp-org/rustbca
cd rustbca
python -m pip install .
python
Python 3.9.6 (tags/v3.9.6:db3ff76, Jun 28 2021, 15:26:21) [MSC v.1929 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> from libRustBCA import *; from scripts.materials import *
>>> angle = 0.0 # deg
>>> energy = 1000.0 # eV
>>> num_samples = 10000
>>> sputtering_yield(argon, tungsten, energy, angle, num_samples)
1.0398
>>> reflection_coefficient(argon, tungsten, energy, angle, num_samples)
(0.3294, 0.10230906775743769) # (reflection coefficient, energy reflection coefficient)
>>>
```

For those eager to get started with the standalone code, try running one of the examples in the
`RustBCA/examples` directory. Note that to automatically manipulate input files and reproduce
the plots located on the [Wiki], these may require some optional
[Python] packages (`matplotlib`, `numpy`, `scipy`, `shapely`, and `toml`).

### H trajectories and collision cascades in a boron nitride dust grain

First, run the example using:

```bash
cargo run --release examples/boron_nitride.toml
```

Afterwords, fire up your favourite [Python] interpreter
(e.g., [IPython]) and execute:

```python
from scripts.rustbca import *
do_trajectory_plot("boron_dust_grain_")
```

### He implantation into a layered TiO<sub>2</sub>/Al/Si target

First, run the example using:

```bash
cargo run --release examples/layered_geometry.toml
```

Afterwords, fire up your favourite [Python] interpreter
(e.g., [IPython]) and execute:

```python
import numpy as np
import matplotlib.pyplot as plt

deposited_ions = np.genfromtxt(
    "2000.0eV_0.0001deg_He_TiO2_Al_Sideposited.output",
    delimiter=",",
    names=["M", "Z", "x", "y", "z", "collisions"],
)

plt.hist(deposited_ions["x"], bins=100)

plt.show()
```

## Features

The following features are implemented in `RustBCA`:

* Ion-material interactions for all combinations of incident ion and target species.
* Infinite, homogeneous targets (Mesh0D), Layered, finite-depth inhomogeneous targets (Mesh1D), arbitrary 2D composition through a triangular mesh (Mesh2D), fast homogeneous 2D geometry (Homogeneous2D), homogeneous spherical geometry (Sphere), and homogeneous 3D triangular mesh geometry (TriMesh).
* Amorphous Solid/Liquid targets, Gaseous targets, and targets with both solid/liquid and gaseous elements
* Low energy (< 25 keV/nucleon) electronic stopping modes including:
  * local (Oen-Robinson),
  * nonlocal (Lindhard-Scharff),
  * and equipartition
* Biersack-Varelas interpolation is also included for electronic stopping up to ~1 GeV/nucleon. Note that high energy physics beyond electronic stopping are not included, and that Biersack-Varelas may not be as accurate as other methods.
* Biersack-Haggmark treatment of high-energy free-flight paths between collisions can be included to greatly speed up high-energy simulations (i.e., by neglecting very small angle scattering).
* A wide range of interaction potentials are provided, including:
  * the Kr-C, ZBL, Lenz-Jensen, and Moliere universal, screened-Coulomb potentials.
  * the Lennard-Jones 12-6 and Morse attractive-repulsive potentials.
* Solving the distance-of-closest-approach problem is achieved using:
  * the Newton-Raphson method for simple root-finding,
  * or, for attractive-repulsive potentials, an Adaptive Chebyshev Proxy Rootfinder with Automatic Subdivision algorithm and a Polynomial root-finding algorithm are provided through the [rcpr] crate.
* Multiple interaction potentials can be used in a single simulation for any number of potentials/species.
  * For example, the He-W interaction can be specified using a Lennard-Jones 12-6 potential, while the W-W interaction can be defined using a Kr-C potential.
* The scattering integral can be calculated using:
  * Gauss-Mehler quadrature,
  * Gauss-Legendre quadrature,
  * Mendenall-Weller quadrature,
  * or the MAGIC algorithm (for certain screened Coulomb potentials only).
* Input files use the [TOML] format, making them both human-readable and easily parsable.
* RustBCA generates user-friendly, context-providing error messages, which help pinpoint the cause of errors and provide suggested fixes to the user.
* The simulation results are comma-delimited (`csv` format) and include:
  * the energies and directions of emitted particles (reflected ions and sputtered atoms),
  * the final positions of implanted ions,
  * displacements,
  * full trajectory tracking for both the incident ions and target atoms,
  * and many other parameters such as position of origin of sputtered particles and energy loss along trajectories.
* Optionally, the code can produce energy-angle and implantation distributions when built with the `--features distributions` flag and disable space-intensive particle list output with `--features no_list_output`.
* Library functions for modeling ion reflection, implantation, and sputtering in C++/C, Python, and Fortran codes.

## Installation

Without optional features, `RustBCA` should compile with `cargo` alone on
Windows, MacOS, and Linux systems.

[HDF5] for particle list input has been tested on Windows, but version 1.10.6 must be used.

#### Manual Dependences

* [rustup], the [Rust] toolchain (includes `cargo`, the [Rust] package manager, `rustc`, the [Rust] compiler, and more).

#### Automatic Dependencies

* see [Cargo.toml](https://github.com/lcpp-org/RustBCA/blob/master/Cargo.toml) for a complete list.

#### Optional Dependencies

* [HDF5] libraries
* [rcpr], a CPR and polynomial rootfinder, required for using attractive-repulsive interaction potentials such as Lennard-Jones or Morse. It may require additional software (see below).
* For manipulating input files and running associated scripts, the following are required:
  * [Python] 3.6+
  * The [Python] libraries: `numpy`, `matplotlib`, `toml` (must build from source), `shapely`, and `scipy`.

### Detailed instructions for Ubuntu 18.04 LTS

1. (Optional) Install Python 3.6+ (this comes natively in Ubuntu 18.04)
2. Install `curl`:
```bash
sudo apt-get install curl
```
3. Install [rustup], the Rust toolchain (includes rustc, the compiler, and cargo, the package manager) from https://rustup.rs/ by running the following command and following on-screen instructions:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
4. (Optional) Install `pip` for [Python]:
```bash
sudo apt-get install python3-pip
```
5. (Optional) Install [Python] libraries for making input files:
```bash
python3 -m pip install numpy matplotlib shapely scipy
```
6. (Optional) Install [Python] [TOML] library from source:
```bash
git clone https://github.com/uiri/toml.git
cd toml
python3 setup.py install
```
7. (Optional) Install software for [rcpr]:
```bash
sudo apt-get install gcc gfortran build-essential cmake liblapack-dev libblas-dev liblapacke-dev
```
8. (Optional - should come with rustup) Install `cargo`:
```bash
sudo apt-get install cargo
```
9. Build `RustBCA`:
```bash
git clone https://github.com/lcpp-org/rustBCA
cd RustBCA
cargo build --release
```
10. (Optional) Build `RustBCA` with optional dependencies, `hdf5` and/or `rcpr` (with your choice of backend: `openblas`, `netlib`, or `intel-mkl`):
```bash
cargo build --release --features cpr_rootfinder_netlib,hdf5_input
cargo build --release --features cpr_rootfinder_openblas,hdf5_input
cargo build --release --features cpr_rootfinder_intel_mkl,hdf5_input
 ```
11. `input.toml` is the input file - see [Usage](https://github.com/lcpp-org/RustBCA/wiki/Usage,-Input-File,-and-Output-Files) for more information
12. Run the required tests using:
```bash
cargo test
```
13. (Optional) Run the required and optional tests for the desired backend(s):
```bash
cargo test --features cpr_rootfinder_netlib
cargo test --features cpr_rootfinder_openblas
cargo test --features cpr_rootfinder_intel_mkl
```

### Detailed instructions for Fedora 33

Most of the ingredients for building `RustBCA` and running the [Python] helper
scripts are available natively in the Fedora software repository, so the setup
is relatively painless.

The [Rust] toolchain can be aquired using:

```bash
sudo dnf install rust rust-src rust-std-static rust-analysis rust-gdb rust-lldb rustfmt
```

The (optional) [Python] packages can be obtained using:

```bash
sudo dnf install python3-numpy python3-scipy python3-matplotlib python3-toml python3-shapely
```

or, alternatively, using `pip3`.

If desired, RustBCA can be built with [rcpr] to simulate attractive-repuslive interaction potentials; rcpr requires (at least) the following:

```bash
sudo dnf install gcc gcc-gfortran cmake lapack lapack-devel blas blas-devel
```

Building `RustBCA` is straightforward, and can be done using:

```bash
git clone https://github.com/lcpp-org/rustBCA
cd RustBCA
cargo build --release
```

with all of the explicit dependencies listed in `Cargo.toml` handled
automatically during the build.

## Usage

To use `RustBCA`, modify an `input.toml` file, which is used to configure each
simulation.
To run a simulation, execute:

```bash
./RustBCA
```

with `input.toml` in the same directory as `RustBCA`.

Alternatively, `RustBCA` accepts the name of a`.toml` input file as a single
command line argument.

```bash
./RustBCA /path/to/input.toml
```

Additionally, `RustBCA` accepts an input file type (one of: `0D`, `1D`, `2D`, `TRIMESH`, `SPHERE` - see the wiki for more details):
```bash
./RustBCA 0D /path/to/input.toml
```

For further details, have a look at
[Usage](https://github.com/lcpp-org/RustBCA/wiki/Usage,-Input-File,-and-Output-Files)
on the `RustBCA` [Wiki] for usage instructions.
Also have a look at the examples on the [Wiki] for writing `.toml` input files.

[BCA]: https://en.wikipedia.org/wiki/Binary_collision_approximation
[HDF5]: https://en.wikipedia.org/wiki/Hierarchical_Data_Format
[IPython]: https://en.wikipedia.org/wiki/IPython
[Python]: https://en.wikipedia.org/wiki/Python_(programming_language)
[rcpr]: https://github.com/drobnyjt/rcpr
[rustup]: https://rustup.rs
[Rust]: https://en.wikipedia.org/wiki/Rust_(programming_language)
[TOML]: https://en.wikipedia.org/wiki/TOML
[Wiki]: https://github.com/lcpp-org/RustBCA/wiki
