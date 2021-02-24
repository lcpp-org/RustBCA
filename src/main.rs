#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

#[cfg(feature = "cpr_rootfinder_openblas")]
extern crate openblas_src;
#[cfg(feature = "cpr_rootfinder_netlib")]
extern crate netlib_src;
#[cfg(feature = "cpr_rootfinder_intel_mkl")]
extern crate intel_mkl_src;

use std::{env, fmt};
use std::mem::discriminant;

//Progress bar crate - works with rayon
use indicatif::{ProgressIterator, ProgressBar, ProgressStyle};

//Error handling crate
use anyhow::Result;
use anyhow::*;

//Serializing/Deserializing crate
use serde::*;

//Array input via hdf5
#[cfg(feature = "hdf5_input")]
use hdf5::*;

//Parallelization
use rayon::prelude::*;
use rayon::*;

//I/O
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::BufWriter;

//Math
use std::f64::consts::FRAC_2_SQRT_PI;
use std::f64::consts::PI;
use std::f64::consts::SQRT_2;

//Load internal modules
pub mod material;
pub mod particle;
pub mod tests;
pub mod interactions;
pub mod bca;
pub mod mesh;
pub mod input;
pub mod output;

//Physical constants
///Fundamental charge in Coulombs.
const Q: f64 = 1.602176634E-19;
/// One electron-volt in Joules.
const EV: f64 = Q;
/// One atomic mass unit in kilograms.
const AMU: f64 = 1.66053906660E-27;
/// One Angstrom in meters.
const ANGSTROM: f64 = 1E-10;
/// One micron in meters.
const MICRON: f64 = 1E-6;
/// One nanometer in meters.
const NM: f64 = 1E-9;
/// One centimeter in meters.
const CM: f64 = 1E-2;
/// Vacuum permitivity in Farads/meter.
const EPS0: f64 = 8.8541878128E-12;
/// Bohr radius in meters.
const A0: f64 = 5.29177210903E-11;
/// Electron mass in kilograms.
const ME: f64 = 9.1093837015E-31;
/// sqrt(pi).
const SQRTPI: f64 = 2. / FRAC_2_SQRT_PI;
/// sqrt(2 * pi).
const SQRT2PI: f64 = 2. * SQRT_2 / FRAC_2_SQRT_PI;
/// Speed of light in meters/second.
const C: f64 = 299792458.;
/// Bethe-Bloch electronic stopping prefactor, in SI units.
const BETHE_BLOCH_PREFACTOR: f64 = 4.*PI*(Q*Q/(4.*PI*EPS0))*(Q*Q/(4.*PI*EPS0))/ME/C/C;
/// Lindhard-Scharff electronic stopping prefactor, in SI units.
const LINDHARD_SCHARFF_PREFACTOR: f64 = 1.212*ANGSTROM*ANGSTROM*Q;
/// Lindhard reduced energy prefactor, in SI units.
const LINDHARD_REDUCED_ENERGY_PREFACTOR: f64 = 4.*PI*EPS0/Q/Q;

/// Mode of electronic stopping to use.
#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum ElectronicStoppingMode {
    /// Biersack-Varelas interpolated electronic stopping. Valid for ~eV/nucleon to ~GeV/nucleon.
    INTERPOLATED,
    /// Oen-Robinson Firsov-type local electronic stopping. Valid up to ~25 keV/nucleon.
    LOW_ENERGY_LOCAL,
    /// Lindhard-Scharff nonlocal electronic stopping. Valid up to ~25 keV/nucleon.
    LOW_ENERGY_NONLOCAL,
    /// Equipartition between Oen-Robinson and Lindhard-Scharff electronic stopping formulas.
    LOW_ENERGY_EQUIPARTITION,
}

impl fmt::Display for ElectronicStoppingMode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            ElectronicStoppingMode::INTERPOLATED => write!(f, "Biersack-Varelas electronic stopping"),
            ElectronicStoppingMode::LOW_ENERGY_NONLOCAL => write!(f, "Lindhard-Scharff electronic stopping"),
            ElectronicStoppingMode::LOW_ENERGY_LOCAL => write!(f, "Oen-Robinson electronic stopping"),
            ElectronicStoppingMode::LOW_ENERGY_EQUIPARTITION => write!(f, "Equipartition with Lindhard-Scharff and Oen-Robinson"),
        }
    }
}

/// Mode of surface binding energy calculation.
#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum SurfaceBindingModel {
    /// Surface binding energies will be determined individually depending only on the particle's `Es`.
    INDIVIDUAL,
    /// Surface binding energies will be a concentration-weighted average of material surface-binding energies, unless `particle.Es == 0` in which case it will be zero.
    TARGET,
    /// Surface binding energies will be the average of the particle and TARGET, unless either is zero in which case it will be zero.
    AVERAGE,
}

impl fmt::Display for SurfaceBindingModel {
    fn fmt (&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            SurfaceBindingModel::INDIVIDUAL => write!(f,
                "Individual surface binding energies."),
            SurfaceBindingModel::TARGET => write!(f,
                "Concentration-dependent linear combinaion of target binding energies."),
            SurfaceBindingModel::AVERAGE => write!(f,
                "Average between particle and concentration-dependent linear combination of target binding energies."),
        }
    }
}

/// Mean-free-path model.
#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum MeanFreePathModel {
    /// Constant mean-free-path for liquids and amorphous solids.
    LIQUID,
    /// Exponentially-distributed mean-free-paths for gases.
    GASEOUS,
}

impl fmt::Display for MeanFreePathModel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            MeanFreePathModel::LIQUID => write!(f, "Amorphous Solid/Liquid Model"),
            MeanFreePathModel::GASEOUS => write!(f, "Gaseous Model"),
        }
    }
}

/// Interatomic potentials between particles in rustbca.
#[derive(Deserialize, Clone, Copy)]
pub enum InteractionPotential {
    /// TRIDYN-style Kr-C. Equivalent to KR_C, except for the MAGIC constants.
    TRIDYN,
    /// Moliere's approximation to the Thomas-Fermi interatomic potential.
    MOLIERE,
    /// Krypton-Carbon "universal" interatomic potential.
    KR_C,
    /// Ziegler-Biersack-Littmark "unviversal" semi-empirical interatomic potential.
    ZBL,
    /// Lenz-Jensen screened Coulomb potential.
    LENZ_JENSEN,
    /// Lennard-Jones 12-6 potential, with user-defined sigma and epsilon.
    LENNARD_JONES_12_6 {sigma: f64, epsilon: f64},
    /// Lennard-Jones 6.5-6 potential, with user-defined sigma and epsilon.
    LENNARD_JONES_65_6 {sigma: f64, epsilon: f64},
    /// Morse potential, with user-defined D, alpha, and r0.
    MORSE{D: f64, alpha: f64, r0: f64},
    /// Tungsten-tungsten cubic spline potential (Following Bjorkas et al.)
    WW,
    /// Unscreened Coulombic interatomic potential between ions with charges Za and Zb.
    COULOMB{Za: f64, Zb: f64}
}

impl fmt::Display for InteractionPotential {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            InteractionPotential::TRIDYN => write!(f, "TRIDYN-style Kr-C (Different MAGIC constants)"),
            InteractionPotential::MOLIERE => write!(f, "Moliere Potential"),
            InteractionPotential::KR_C => write!(f, "Kr-C Potential"),
            InteractionPotential::ZBL => write!(f, "ZBL Potential"),
            InteractionPotential::LENZ_JENSEN => write!(f, "Lenz-Jensen Potential"),
            InteractionPotential::LENNARD_JONES_12_6{sigma, epsilon} => write!(f, "Lennard-Jones 12-6 Potential with sigma = {} A, epsilon = {} eV", sigma/ANGSTROM, epsilon/EV),
            InteractionPotential::LENNARD_JONES_65_6{sigma, epsilon} => write!(f, "Lennard-Jones 6.5-6 Potential with sigma = {} A, epsilon = {} eV", sigma/ANGSTROM, epsilon/EV),
            InteractionPotential::MORSE{D, alpha, r0} => write!(f, "Morse potential with D = {} eV, alpha = {} 1/A, and r0 = {} A", D/EV, alpha*ANGSTROM, r0/ANGSTROM),
            InteractionPotential::WW => write!(f, "W-W cubic spline interaction potential."),
            InteractionPotential::COULOMB{Za, Zb} => write!(f, "Coulombic interaction with Za = {} and Zb = {}", Za, Zb)
        }
    }
}

impl PartialEq for InteractionPotential {
    fn eq(&self, other: &Self) -> bool {
        discriminant(self) == discriminant(other)
    }
}

/// Method for solving the scattering integral.
#[derive(Deserialize, Clone, Copy)]
pub enum ScatteringIntegral {
    /// Mendenhall-Weller Gauss-Lobatto 4-point quadrature.
    MENDENHALL_WELLER,
    /// Ziegler's MAGIC algorithm.
    MAGIC,
    /// Gauss-Mehler n-point quadrature.
    GAUSS_MEHLER{n_points: usize},
    /// Gauss-Legendre 5-point quadrature.
    GAUSS_LEGENDRE
}

impl fmt::Display for ScatteringIntegral {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            ScatteringIntegral::MENDENHALL_WELLER => write!(f, "Mendenhall-Weller 4-Point Lobatto Quadrature"),
            ScatteringIntegral::MAGIC => write!(f, "MAGIC Algorithm"),
            ScatteringIntegral::GAUSS_MEHLER{n_points} => write!(f, "Gauss-Mehler {}-point Quadrature", n_points),
            ScatteringIntegral::GAUSS_LEGENDRE => write!(f, "Gauss-Legendre 5-point Quadrature"),
        }
    }
}

impl PartialEq for ScatteringIntegral {
    fn eq(&self, other: &Self) -> bool {
        discriminant(self) == discriminant(other)
    }
}

/// Root-finding algorithm.
#[derive(Deserialize, Clone, Copy)]
pub enum Rootfinder {
    /// Newton root-finder with user-defined `max_iterations` and `tolerance`.
    NEWTON{max_iterations: usize, tolerance: f64},
    CPR{n0: usize, nmax: usize, epsilon: f64, complex_threshold: f64, truncation_threshold: f64,
        far_from_zero: f64, interval_limit: f64, derivative_free: bool},
    POLYNOMIAL{complex_threshold: f64},
}

impl fmt::Display for Rootfinder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Rootfinder::NEWTON{max_iterations, tolerance} => write!(f, "Newton-Raphson Rootfinder with maximum {} iterations and toleance = {}", max_iterations, tolerance),
            Rootfinder::CPR{n0, nmax, epsilon, complex_threshold, truncation_threshold, far_from_zero, interval_limit, derivative_free} =>
                write!(f, "Chebyshev-Proxy Rootfinder with {}-polishing", match derivative_free { true => "Secant", false => "Newton"}),
            Rootfinder::POLYNOMIAL{complex_threshold} => write!(f, "Frobenius Companion Matrix Polynomial Real Rootfinder with a complex tolerance of {}", complex_threshold),
        }
    }
}
impl PartialEq for Rootfinder {
    fn eq(&self, other: &Self) -> bool {
        discriminant(self) == discriminant(other)
    }
}

/// 3D vector.
#[derive(Clone)]
pub struct Vector {
    x: f64,
    y: f64,
    z: f64,
}
impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Vector {
        Vector {
            x,
            y,
            z
        }
    }

    /// Calculates vector magnitude.
    fn magnitude(&self) -> f64 {
        return (self.x*self.x + self.y*self.y + self.z*self.z).sqrt();
    }

    /// Assigns vector values from another vector.
    fn assign(&mut self, other: &Vector) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
    }

    /// Normalizes vector components to magnitude 1.
    fn normalize(&mut self) {
        let magnitude = self.magnitude();
        self.x /= magnitude;
        self.y /= magnitude;
        self.z /= magnitude;
    }

    /// Add this vector and another and return a new vector.
    pub fn add(&self, other: &Vector) -> Vector {
        Vector::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

/// Vector4 is a trajectory-tracking object that includes x, y, z, and the current energy.
#[derive(Clone)]
pub struct Vector4 {
    E: f64,
    x: f64,
    y: f64,
    z: f64,
}

impl Vector4 {
    fn new(E: f64, x: f64, y: f64, z: f64) -> Vector4 {
        Vector4 {
            E,
            x,
            y,
            z
        }
    }
}

/// Energy loss is an output tracker that tracks the separate nuclear and electronic energy losses.
#[derive(Clone)]
pub struct EnergyLoss {
    En: f64,
    Ee: f64,
    x: f64,
    y: f64,
    z: f64,
}

impl EnergyLoss {
    fn new(Ee: f64, En: f64, x: f64, y: f64, z: f64) -> EnergyLoss {
        EnergyLoss {
            En,
            Ee,
            x,
            y,
            z
        }
    }
}

/// Rustbca's internal representation of an input file.
#[derive(Deserialize)]
pub struct Input {
    options: Options,
    material_parameters: material::MaterialParameters,
    particle_parameters: particle::ParticleParameters,
    mesh_2d_input: mesh::Mesh2DInput,
}

/// Rustbca's internal representation of the simulation-level options.
#[derive(Deserialize)]
pub struct Options {
    name: String,
    track_trajectories: bool,
    track_recoils: bool,
    track_recoil_trajectories: bool,
    write_buffer_size: usize,
    weak_collision_order: usize,
    suppress_deep_recoils: bool,
    high_energy_free_flight_paths: bool,
    electronic_stopping_mode: ElectronicStoppingMode,
    mean_free_path_model: MeanFreePathModel,
    interaction_potential: Vec<Vec<InteractionPotential>>,
    scattering_integral: Vec<Vec<ScatteringIntegral>>,
    root_finder: Vec<Vec<Rootfinder>>,
    num_threads: usize,
    num_chunks: u64,
    use_hdf5: bool,
    track_displacements: bool,
    track_energy_losses: bool,
}

pub struct OutputUnits {
    length_unit: f64,
    energy_unit: f64,
    mass_unit: f64
}

fn main() {
    //Open and process input_file
    let (particle_input_array, material, options, output_units) = input::input();

    println!("Processing {} ions...", particle_input_array.len());
    let total_count: u64 = particle_input_array.len() as u64;
    assert!(total_count/options.num_chunks > 0, "Input error: chunk size == 0 - reduce num_chunks or increase particle count.");

    //Open output files
    let mut output_list_streams = output::open_output_lists(&options);
    let mut summary_stream = output::open_output_summary(&options);
    let mut summary = output::Summary::new(total_count);

    //Initialize threads with rayon
    println!("Initializing with {} threads...", options.num_threads);
    if options.num_threads > 1 {let pool = rayon::ThreadPoolBuilder::new().num_threads(options.num_threads).build_global().unwrap();};

    //Create and configure progress bar
    let bar: ProgressBar = ProgressBar::new(options.num_chunks);
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] [{bar:40.cyan/blue}] {percent}%")
        .progress_chars("#>-"));

    //Main loop
    for (chunk_index, particle_input_chunk) in particle_input_array.chunks((total_count/options.num_chunks) as usize).progress_with(bar).enumerate() {

        let mut finished_particles: Vec<particle::Particle> = Vec::new();

        if options.num_threads > 1 {
            finished_particles.par_extend(
                particle_input_chunk.into_par_iter()
                .map(|particle_input| bca::single_ion_bca(particle::Particle::from_input(*particle_input, &options), &material, &options))
                    .flatten()
            );
        } else {
            finished_particles.extend(
                particle_input_chunk.iter()
                .map(|particle_input| bca::single_ion_bca(particle::Particle::from_input(*particle_input, &options), &material, &options))
                    .flatten()
            );
        }

        for particle in finished_particles {
            summary.add(&particle);
            output::output_lists(&mut output_list_streams, particle, &options, &output_units);
        }
        //Flush all file streams before dropping to ensure all data is written
        output::output_list_flush(&mut output_list_streams);
    }

    //Write to summary file
    writeln!(summary_stream, "incident, sputtered, reflected")
        .expect(format!("Output error: could not write to {}summary.output.", options.name).as_str());
    writeln!(summary_stream, "{}, {}, {}",
    summary.num_incident, summary.num_sputtered, summary.num_reflected)
        .expect(format!("Output error: could not write to {}summary.output.", options.name).as_str());
    summary_stream.flush().unwrap();

    println!("Finished!");
}
