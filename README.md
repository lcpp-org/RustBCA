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
valid for incident ion energies between several eV/nucleon  to 
<1 GeV/nucleon. Improvements to RustBCA have expanded the regime
of validity for some quantities, such as reflection coefficients, below 
1 eV/nucleon for some ion/target pairs.

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
* [Enabling attractive-repulsive potentials in binary-collision-approximation monte-carlo codes for ion-surface interactions](https://doi.org/10.1088/2053-1591/ad1262), J Drobny and D Curreli (2023)
* [Multi-physics modeling of tungsten collector probe samples during the WEST C4 He campaign](https://doi.org.10.1088/1741-4326/ad6c5b), A. Lasa et al. (2024)

## Getting started

The easiest way to get started is with the ergonomic Python functions.
These functions use the default RustBCA options detailed on the
[Input Files](https://github.com/lcpp-org/RustBCA/wiki/Standalone-Code:-Input-File) page, which are not universally applicable.
These examples use example material parameters located in
`scripts/materials.py` that should be verified before use.
Follow these steps to install, build, and run simple RustBCA simulations
for sputtering yields and reflection coefficients:
```
git clone https://github.com/lcpp-org/rustbca
cd rustbca
python -m pip install .
python
Python 3.9.6 (tags/v3.9.6:db3ff76, Jun 28 2021, 15:26:21) [MSC v.1929 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> from libRustBCA import *; from scripts.materials import *; import numpy as np
>>> angle = 0.0 # deg
>>> energy = 1000.0 # eV
>>> num_samples = 10000
>>> 1 < sputtering_yield(argon, tungsten, energy, angle, num_samples) < 1.1 # Y approx. 1.04
True
>>> R_N, R_E = reflection_coefficient(argon, tungsten, energy, angle, num_samples)
>>> 0.3 < R_N < 0.4 # R_N approx. 0.35 
True
>>> 0.0 < R_E < 0.2 # R_E approx 0.1
True
```

For those eager to get started with the standalone code, try running one of the examples in the
`RustBCA/examples` directory. Note that to automatically manipulate input files and reproduce
the plots located on the [Wiki], these may require some optional
[Python] packages (`matplotlib`, `numpy`, `scipy`, `shapely`, and `toml`).

### H trajectories and collision cascades in boron nitride
First, run the example using:

```shell
$ cargo run --release 0D examples/boron_nitride_0D.toml 2>/dev/null  #suppress progress bar for automatic testing
Processing 10 ions...
Initializing with 4 threads...
Finished!
```

Afterwords, fire up your favourite [Python] interpreter
(e.g., [IPython]) and execute:

```python
from scripts.rustbca import *
do_trajectory_plot("boron_nitride_")
```

### He implantation into a layered TiO<sub>2</sub>/Al/Si target

First, run the example using:

```shell
$ cargo run --release examples/layered_geometry.toml 2>/dev/null #suppress progress bar for automatic testing
Processing 10000 ions...
Initializing with 4 threads...
Finished!
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
  * the Newton-Raphson method for purely repulsive potentials,
  * or, for attractive-repulsive potentials, an Adaptive Chebyshev Proxy Rootfinder with Automatic Subdivision algorithm and a polynomial root-finding algorithm are provided through [rcpr].
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
Windows, MacOS, and Linux systems:

```
cargo build --release 
```

will add an executable at `target/release/`.

[HDF5] for particle list input has been tested on Windows, but version 1.10.6 must be used.

#### Manual Dependences

* [rustup], the [Rust] toolchain (includes `cargo`, the [Rust] package manager, `rustc`, the [Rust] compiler, and more).

#### Automatic Dependencies

* see [Cargo.toml](https://github.com/lcpp-org/RustBCA/blob/master/Cargo.toml) for a complete list of required and optional dependencies managed by `cargo`.

#### Optional Dependencies

* [HDF5] libraries
* For manipulating input files and running associated scripts, the following are suggested:
  * [Python] 3.6+
  * [Python] libraries: `numpy`, `matplotlib`, `toml`, `shapely`, and `scipy`.

### Detailed instructions for Ubuntu 18.04 LTS

1. (Optional) Install Python 3.6+ (this comes natively in Ubuntu 18.04)
2. Install `curl`:
```bash
sudo apt-get install curl
```
3. Install [rustup], the Rust toolchain (includes rustc, the compiler, and cargo, the package manager) from https://rustup.rs/ by running the following command and following on-screen instructions:
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
```
4. (Optional) Install `pip` for [Python]:
```bash
sudo apt-get install python3-pip
```
5. (Optional) Install [Python] libraries for making input files:
```bash
python3 -m pip install numpy matplotlib shapely scipy toml
```
6. (Optional - should come with rustup) Install `cargo`:
```bash
sudo apt-get install cargo
```
7. Build `RustBCA`:
```bash
git clone https://github.com/lcpp-org/RustBCA
cd RustBCA
cargo build --release
```
8. (Optional) Build `RustBCA` with optional dependencies, `hdf5` and/or `rcpr`:
```bash
cargo build --release --features cpr_rootfinder,hdf5
```
9. `input.toml` is the input file - see the [Input File](https://github.com/lcpp-org/RustBCA/wiki/Standalone-Code:-Input-File) page for more information
10. Run the required tests using:
```bash
cargo test
```
11. (Optional) Run the tests for the advanced rootfinder:
```bash
cargo test --features cpr_rootfinder
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
**Warning: RustBCA defaults to the 2D triangular mesh input mode.** For more details, see [Input Files](https://github.com/lcpp-org/RustBCA/wiki/Standalone-Code:-Input-File).
Also have a look at the examples on the [Wiki] to see some examples of RustBCA input files.

[BCA]: https://en.wikipedia.org/wiki/Binary_collision_approximation
[HDF5]: https://en.wikipedia.org/wiki/Hierarchical_Data_Format
[IPython]: https://en.wikipedia.org/wiki/IPython
[Python]: https://en.wikipedia.org/wiki/Python_(programming_language)
[rcpr]: https://github.com/drobnyjt/rcpr
[rustup]: https://rustup.rs
[Rust]: https://en.wikipedia.org/wiki/Rust_(programming_language)
[TOML]: https://en.wikipedia.org/wiki/TOML
[Wiki]: https://github.com/lcpp-org/RustBCA/wiki
