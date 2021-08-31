# `RustBCA`

`RustBCA` is a general-purpose, high-performance code for simulating
ion-material interactions using the binary collision approximation ([BCA]),
written in [Rust]!

By discretizing the collision cascade into a sequence of binary collisions,
[BCA] codes can accurately and efficiently model the prompt interaction
between an energetic ion and a target material.
This includes reflection, implantation, and transmission of the incident ion,
as well as sputtering and displacement damage of the target.
Generally, [BCA] codes are valid for incident ion energies between approximately
~1 eV/nucleon to ~1 GeV/nucleon.

Check out the `RustBCA` [Wiki] for detailed information, installation
instructions, use cases, examples, and more. See the RustBCA paper at the
Journal of Open Source Software by clicking the badge below:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03298/status.svg)](https://doi.org/10.21105/joss.03298)

## Getting started

For those eager to get started, try running one of the examples in the
`rustBCA` directory. Note that to automatically manipulate input files and reproduce the plots located on the [Wiki], these require several optional, but common,
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

The following features are implemented in `rustBCA`:

* Ion-material interactions for all combinations of incident ion and target species.
* Infinite, homogeneous targets (Mesh0D), Layered, finite-depth inhomogeneous targets (Mesh1D) and arbitrary 2D geometry and composition through a triangular mesh (Mesh2D).
* Amorphous Solid/Liquid targets, Gaseous targets, and targets with both solid/liquid and gaseous elements
* Low energy (< 25 keV/nucleon) electronic stopping modes including:
  * local (Oen-Robinson),
  * nonlocal (Lindhard-Scharff),
  * and equipartition forms.
* Biersack-Varelas interpolation is also included for electronic stopping up to ~1 GeV/nucleon.
* Optionally, the Biersack-Haggmark treatment of high-energy free-flight paths between collisions can be included to greatly speed up high-energy simulations (i.e., by neglecting very small angle scattering).
* A wide range of interaction potentials are provided, including:
  * the Kr-C, ZBL, Lenz-Jensen, and Moliere universal, screened-Coulomb potentials.
  * the Lennard-Jones 12-6, Lennard-Jones 6.5-6, and Morse attractive-repulsive potentials.
* Solving the distance-of-closest-approach problem is achieved using:
  * the Newton-Raphson method for simple root-finding,
  * or, for attractive-repulsive potentials, an Adaptive Chebyshev Proxy Rootfinder with Automatic Subdivision algorithm and a Polynomial root-finding algorithm are provided through the [rcpr] crate.
* Multiple interaction potentials can be used in a single simulation for any number of potentials/species.
  * For example, the He-W interaction can be specified using a Lennard-Jones 12-6 potential, while the W-W interaction can be defined using a Kr-C potential.
* The scattering integral can be calculated using:
  * Gauss-Mehler quadrature,
  * Gauss-Legendre quadrature,
  * Mendenall-Weller quadrature,
  * or the MAGIC algorithm.
* Input files use the [TOML] format, making them both human-readable and easily parsable.
* RustBCA generates user-friendly, context-providing error messages, which help pinpoint the cause of errors and provide suggested fixes to the user.
* The simulation results are formatted as the ubiquitous `csv` format and include:
  * the energies and directions of emitted particles (reflected ions and sputtered atoms),
  * the final positions of implanted ions,
  * and full trajectory tracking for both the incident ions and target atoms.
* Optionally, the code can produce energy-angle and implantation distributions when built with the `--features distributions` flag and disable space-intensive particle list output with `--features no_list_output`.

## Installation

Without optional features, `rustBCA` should compile with `cargo` alone on
Windows, MacOS, and Linux systems.
[HDF5] has been tested on Windows, but version 1.10.6 must be used.
[rcpr], the adaptive Chebyshev Proxy Rootfinder with automatic subdivision and
polynomial rootfinder package for [Rust], has not yet been successfully compiled
on Windows.
However, it can be compiled on the Windows Subsystem for Linux (WSL) and, likely,
on Ubuntu for Windows or Cygwin.

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
8. Install `cargo`:
```bash
sudo apt-get install cargo
```
9. Build `rustBCA`:
```bash
git clone https://github.com/lcpp-org/rustBCA
cd rustBCA
cargo build --release
```
10. (Optional) Build `rustBCA` with optional dependencies, `hdf5` and/or `rcpr` (with your choice of backend: `openblas`, `netlib`, or `intel-mkl`):
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

Most of the ingredients for building `rustBCA` and running the [Python] helper
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

If the [rcpr] is desired, it's probably also a good idea to install the following:

```bash
sudo dnf install gcc gcc-gfortran cmake lapack lapack-devel blas blas-devel
```

Building `rustBCA` is straightforward and can be done using:

```bash
git clone https://github.com/lcpp-org/rustBCA
cd RustBCA
cargo build --release
```

with all of the explicit dependencies listed in `Cargo.toml` handled
automatically during the build.

## Usage

To use `rustBCA`, modify the `input.toml` file, which is used to configure each
simulation.
To run a simulation, execute:

```bash
./rustBCA
```

with `input.toml` in the same directory as `rustBCA`.
Alternatively, `rustBCA` accepts the name of a`.toml` input file as a single
command line argument:

```bash
./rustBCA /path/to/input.toml
```

For further details, have a look at
[Usage](https://github.com/lcpp-org/RustBCA/wiki/Usage,-Input-File,-and-Output-Files)
on the `rustBCA` [Wiki] for usage instructions.
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
