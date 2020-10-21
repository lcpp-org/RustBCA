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

//extern crate openblas_src;
//Progress bar crate - works wiht rayon
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
use std::f64::consts::PI;

//Load internal modules
pub mod material;
pub mod particle;
pub mod tests;
pub mod interactions;
pub mod bca;
pub mod mesh;

//Physical constants
const Q: f64 = 1.60217646E-19;
const EV: f64 = Q;
const AMU: f64 = 1.660539E-27;
const ANGSTROM: f64 = 1E-10;
const MICRON: f64 = 1E-6;
const NM: f64 = 1E-9;
const CM: f64 = 1E-2;
const EPS0: f64 = 8.85418781E-12;
const A0: f64 = 5.29177211E-11;
const ME: f64 =  9.109383632E-31;
const SQRTPI: f64 = 1.772453850906;
const SQRT2PI: f64 = 2.506628274631;
const C: f64 = 299792000.;
const BETHE_BLOCH_PREFACTOR: f64 = 4.*PI*(Q*Q/(4.*PI*EPS0))*(Q*Q/(4.*PI*EPS0))/ME/C/C;
const LINDHARD_SCHARFF_PREFACTOR: f64 = 1.212*ANGSTROM*ANGSTROM*Q;
const LINDHARD_REDUCED_ENERGY_PREFACTOR: f64 = 4.*PI*EPS0/Q/Q;

#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum ElectronicStoppingMode {
    INTERPOLATED,
    LOW_ENERGY_LOCAL,
    LOW_ENERGY_NONLOCAL,
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

#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum SurfaceBindingModel {
    INDIVIDUAL,
    TARGET,
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


#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum MeanFreePathModel {
    LIQUID,
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

#[derive(Deserialize, Clone, Copy)]
pub enum InteractionPotential {
    TRIDYN,
    MOLIERE,
    KR_C,
    ZBL,
    LENZ_JENSEN,
    LENNARD_JONES_12_6 {sigma: f64, epsilon: f64},
    LENNARD_JONES_65_6 {sigma: f64, epsilon: f64},
    MORSE{D: f64, alpha: f64, r0: f64},
    WW,
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

#[derive(Deserialize, Clone, Copy)]
pub enum ScatteringIntegral {
    MENDENHALL_WELLER,
    MAGIC,
    GAUSS_MEHLER{n_points: usize},
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

#[derive(Deserialize, Clone, Copy)]
pub enum Rootfinder {
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
    fn magnitude(&self) -> f64 {
        return (self.x*self.x + self.y*self.y + self.z*self.z).sqrt();
    }

    fn assign(&mut self, other: &Vector) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
    }

    fn normalize(&mut self) {
        let magnitude = self.magnitude();
        self.x /= magnitude;
        self.y /= magnitude;
        self.z /= magnitude;
    }

    pub fn add(&self, other: &Vector) -> Vector {
        Vector::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

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

#[derive(Deserialize)]
pub struct Input {
    options: Options,
    material_parameters: material::MaterialParameters,
    particle_parameters: particle::ParticleParameters,
    mesh_2d_input: mesh::Mesh2DInput,
}

#[derive(Deserialize)]
pub struct Options {
    name: String,
    track_trajectories: bool,
    track_recoils: bool,
    track_recoil_trajectories: bool,
    stream_size: usize,
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

fn main() {

    let args: Vec<String> = env::args().collect();

    let input_file = match args.len() {
        1 => "input.toml".to_string(),
        2 => args[1].clone(),
        _ => panic!("Too many command line arguments. RustBCA accepts 0 (use 'input.toml') or 1 (input file name).")
    };

    //Read input file, convert to string, and open with toml
    let mut input_toml = String::new();
    let mut file = OpenOptions::new()
        .read(true)
        .write(false)
        .create(false)
        .open(&input_file)
        .expect(format!("Input errror: could not open input file {}.", &input_file).as_str());
    file.read_to_string(&mut input_toml).context("Could not convert TOML file to string.").unwrap();
    let input: Input = toml::from_str(&input_toml)
        .expect("Input error: failed to parse TOML file. Check that all floats have terminal 0s and that there are no missing/extra delimiters.");

    //Unpack toml information into structs
    let material = material::Material::new(input.material_parameters, input.mesh_2d_input);
    let options = input.options;
    let particle_parameters = input.particle_parameters;

    //Check that all material arrays are of equal length.
    assert!(material.m.len() == material.Z.len(), "Input error: material input arrays of unequal length.");
    assert!(material.m.len() == material.Eb.len(), "Input error: material input arrays of unequal length.");
    assert!(material.m.len() == material.Es.len(), "Input error: material input arrays of unequal length.");
    assert!(material.m.len() == material.interaction_index.len(), "Input error: material input arrays of unequal length.");

    //Check that incompatible options are not on simultaneously

    if options.high_energy_free_flight_paths {
        assert!(options.electronic_stopping_mode == ElectronicStoppingMode::INTERPOLATED,
            "Input error: High energy free flight paths used with low energy stoppping power. Change to INTERPOLATED.");
        assert!(options.weak_collision_order == 0,
            "Input error: Cannot use weak collisions with free flight paths. Set weak_collision_order = 0.");
    }

    if options.mean_free_path_model == MeanFreePathModel::GASEOUS {
        assert!(options.weak_collision_order == 0,
            "Input error: Cannot use weak collisions with gaseous mean free path model. Set weak_collision_order = 0.");
    }

    for i in 0..options.interaction_potential.len() {
        assert!(&options.interaction_potential.len() == &options.interaction_potential[i].len(),
            "Input error: interaction matrix not square.");

        assert!(&options.scattering_integral.len() == &options.scattering_integral[i].len(),
            "Input error: scattering intergral matrix not square.");

        assert!(&options.root_finder.len() == &options.root_finder[i].len(),
            "Input error: rootfinder matrix not square.");
    }

    for ((interaction_potentials, scattering_integrals), root_finders) in options.interaction_potential.clone().iter().zip(options.scattering_integral.clone()).zip(options.root_finder.clone()) {
        for ((interaction_potential, scattering_integral), root_finder) in interaction_potentials.iter().zip(scattering_integrals).zip(root_finders) {

            if cfg!(not(any(feature="cpr_rootfinder_openblas", feature="cpr_rootfinder_netlib", feature="cpr_rootfinder_intel_mkl",))) {
                assert!( match root_finder {
                    Rootfinder::POLYNOMIAL{..} => false,
                    Rootfinder::CPR{..} => false,
                    _ => true,
                },
                "Input error: CPR rootfinder not enabled. Build with --features cpr_rootfinder");
            }

            assert!(
                match (interaction_potential, root_finder) {
                    (InteractionPotential::LENNARD_JONES_12_6{..}, Rootfinder::CPR{..}) => true,
                    (InteractionPotential::LENNARD_JONES_12_6{..}, Rootfinder::POLYNOMIAL{..}) => true,
                    (InteractionPotential::LENNARD_JONES_12_6{..}, _) => false,
                    (InteractionPotential::LENNARD_JONES_65_6{..}, Rootfinder::CPR{..}) => true,
                    (InteractionPotential::LENNARD_JONES_65_6{..}, Rootfinder::POLYNOMIAL{..}) => true,
                    (InteractionPotential::LENNARD_JONES_65_6{..}, _) => false,
                    (InteractionPotential::MORSE{..}, Rootfinder::CPR{..}) => true,
                    (InteractionPotential::MORSE{..}, _) => false,
                    (_, Rootfinder::POLYNOMIAL{..}) => false,
                    (_, _) => true,
                },
            "Input error: cannot use {} with {}. Try switching to a different rootfinder.", interaction_potential, root_finder);
        }
    }

    //Check that particle arrays are equal length
    assert_eq!(particle_parameters.Z.len(), particle_parameters.m.len(),
        "Input error: particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.E.len(),
        "Input error: particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.pos.len(),
        "Input error: particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.dir.len(),
        "Input error: particle input arrays of unequal length.");

    //Check that interaction indices are all within interaction matrices
    assert!(material.interaction_index.iter().max().unwrap() < &options.interaction_potential.len(),
        "Input error: interaction matrix too small for material interaction indices.");
    assert!(particle_parameters.interaction_index.iter().max().unwrap() < &options.interaction_potential.len(),
        "Input error: interaction matrix too small for particle interaction indices.");


    //N is the number of distinct particles.
    let N = particle_parameters.Z.len();

    //Determine the length, energy, and mass units for particle input
    let length_unit: f64 = match particle_parameters.length_unit.as_str() {
        "MICRON" => MICRON,
        "CM" => CM,
        "ANGSTROM" => ANGSTROM,
        "NM" => NM,
        "M" => 1.,
        _ => panic!("Input error: unknown unit {} in input file. Choose one of: MICRON, CM, ANGSTROM, NM, M",
            particle_parameters.length_unit.as_str())
    };

    let energy_unit: f64 = match particle_parameters.energy_unit.as_str() {
        "EV" => EV,
        "J"  => 1.,
        "KEV" => EV*1E3,
        "MEV" => EV*1E6,
        _ => panic!("Input error: unknown unit {} in input file. Choose one of: EV, J, KEV, MEV",
            particle_parameters.energy_unit.as_str())
    };

    let mass_unit: f64 = match particle_parameters.mass_unit.as_str() {
        "AMU" => AMU,
        "KG" => 1.0,
        _ => panic!("Input error: unknown unit {} in input file. Choose one of: AMU, KG",
            particle_parameters.mass_unit.as_str())
    };

    //HDF5
    #[cfg(feature = "hdf5_input")]
    let particle_input_array: Vec<particle::ParticleInput> = {
        if options.use_hdf5 {
            let particle_input_filename = particle_parameters.particle_input_filename.as_str();
            let _e = hdf5::silence_errors();
            let particle_input_file = hdf5::File::open(particle_input_filename)
                .context("Input error: cannot open HDF5 file.")
                .unwrap();
            let particle_input = particle_input_file.dataset("particles")
                .context("Input error: cannot read from HDF5 file.")
                .unwrap();
            particle_input.read_raw::<particle::ParticleInput>().unwrap()

        } else {
            let mut particle_input: Vec<particle::ParticleInput> = Vec::new();

            for particle_index in 0..N {
                let N_ = particle_parameters.N[particle_index];
                let m = particle_parameters.m[particle_index];
                let Z = particle_parameters.Z[particle_index];
                let E = particle_parameters.E[particle_index];
                let Ec = particle_parameters.Ec[particle_index];
                let Es = particle_parameters.Es[particle_index];
                let interaction_index = particle_parameters.interaction_index[particle_index];
                let (x, y, z) = particle_parameters.pos[particle_index];
                let (cosx, cosy, cosz) = particle_parameters.dir[particle_index];
                assert!(cosx < 1.,
                    "Input error: particle x-direction cannot be exactly equal to 1 to avoid numerical gimbal lock.");
                for sub_particle_index in 0..N_ {
                    //Add new particle to particle vector
                    particle_input.push(
                        particle::ParticleInput{
                            m: m*mass_unit,
                            Z: Z,
                            E: E*energy_unit,
                            Ec: Ec*energy_unit,
                            Es: Es*energy_unit,
                            x: x*length_unit,
                            y: y*length_unit,
                            z: z*length_unit,
                            ux: cosx,
                            uy: cosy,
                            uz: cosz,
                            interaction_index: interaction_index
                        }
                    );
                }
            }
            particle_input
        }
    };

    #[cfg(not(feature = "hdf5_input"))]
    let particle_input_array: Vec<particle::ParticleInput> = {
        if options.use_hdf5 {
            panic!("HDF5 particle input not enabled. Enable with: cargo build --features hdf5_input")
        } else {
            let mut particle_input: Vec<particle::ParticleInput> = Vec::new();

            for particle_index in 0..N {
                let N_ = particle_parameters.N[particle_index];
                let m = particle_parameters.m[particle_index];
                let Z = particle_parameters.Z[particle_index];
                let E = particle_parameters.E[particle_index];
                let Ec = particle_parameters.Ec[particle_index];
                let Es = particle_parameters.Es[particle_index];
                let interaction_index = particle_parameters.interaction_index[particle_index];
                let (x, y, z) = particle_parameters.pos[particle_index];
                let (cosx, cosy, cosz) = particle_parameters.dir[particle_index];
                assert!(cosx < 1.,
                    "Input error: particle x-direction cannot be exactly equal to 1 to avoid numerical gimbal lock.");
                for sub_particle_index in 0..N_ {

                    //Add new particle to particle vector
                    particle_input.push(
                        particle::ParticleInput{
                            m: m*mass_unit,
                            Z: Z,
                            E: E*energy_unit,
                            Ec: Ec*energy_unit,
                            Es: Es*energy_unit,
                            x: x*length_unit,
                            y: y*length_unit,
                            z: z*length_unit,
                            ux: cosx,
                            uy: cosy,
                            uz: cosz,
                            interaction_index
                        }
                    );
                }
            }
            particle_input
        }
    };

    //Open output files for streaming output
    let reflected_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "reflected.output"))
        .context("Could not open output file.")
        .unwrap();
    let mut reflected_file_stream = BufWriter::with_capacity(options.stream_size, reflected_file);

    let sputtered_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "sputtered.output"))
        .context("Could not open output file.")
        .unwrap();
    let mut sputtered_file_stream = BufWriter::with_capacity(options.stream_size, sputtered_file);

    let deposited_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "deposited.output"))
        .context("Could not open output file.")
        .unwrap();
    let mut deposited_file_stream = BufWriter::with_capacity(options.stream_size, deposited_file);

    let trajectory_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "trajectories.output"))
        .context("Could not open output file.")
        .unwrap();
    let mut trajectory_file_stream = BufWriter::with_capacity(options.stream_size, trajectory_file);

    let trajectory_data = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "trajectory_data.output"))
        .context("Could not open output file.")
        .unwrap();
    let mut trajectory_data_stream = BufWriter::with_capacity(options.stream_size, trajectory_data);

    let displacements_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "displacements.output"))
        .context("Could not open output file.")
        .unwrap();
    let mut displacements_file_stream = BufWriter::with_capacity(options.stream_size, displacements_file);

    let energy_loss_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "energy_loss.output"))
        .context("Could not open output file.")
        .unwrap();
    let mut energy_loss_file_stream = BufWriter::with_capacity(options.stream_size, energy_loss_file);

    println!("Processing {} ions...", particle_input_array.len());

    let total_count: u64 = particle_input_array.len() as u64;

    //Initialize threads with rayon
    println!("Initializing with {} threads...", options.num_threads);
    if options.num_threads > 1 {let pool = rayon::ThreadPoolBuilder::new().num_threads(options.num_threads).build_global().unwrap();};

    let bar: ProgressBar = ProgressBar::new(options.num_chunks);

    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] [{bar:40.cyan/blue}] {percent}%")
        .progress_chars("#>-"));

    assert!(total_count/options.num_chunks > 0, "Input error: chunk size == 0 - reduce num_chunks or increase particle count.");

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
            //
            if !particle.incident & options.track_displacements {
                    writeln!(
                        displacements_file_stream, "{},{},{},{},{},{}",
                        particle.m/mass_unit, particle.Z, particle.energy_origin/energy_unit,
                        particle.pos_origin.x/length_unit, particle.pos_origin.y/length_unit, particle.pos_origin.z/length_unit,
                    ).expect(format!("Output error: could not write to {}displacements.output.", options.name).as_str());
            }

            if particle.incident & particle.left {
                writeln!(
                    reflected_file_stream, "{},{},{},{},{},{},{},{},{},{}",
                    particle.m/mass_unit, particle.Z, particle.E/energy_unit,
                    particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit,
                    particle.dir.x, particle.dir.y, particle.dir.z,
                    particle.number_collision_events
                ).expect(format!("Output error: could not write to {}reflected.output.", options.name).as_str());
            }

            //Incident particle, stopped in material: deposited
            if particle.incident & particle.stopped {
                writeln!(
                    deposited_file_stream, "{},{},{},{},{},{}",
                    particle.m/mass_unit, particle.Z,
                    particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit,
                    particle.number_collision_events
                ).expect(format!("Output error: could not write to {}deposited.output.", options.name).as_str());
            }

            //Not an incident particle, left material: sputtered
            if !particle.incident & particle.left {
                writeln!(
                    sputtered_file_stream, "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                    particle.m/mass_unit, particle.Z, particle.E/energy_unit,
                    particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit,
                    particle.dir.x, particle.dir.y, particle.dir.z,
                    particle.number_collision_events,
                    particle.pos_origin.x/length_unit, particle.pos_origin.y/length_unit, particle.pos_origin.z/length_unit
                ).expect(format!("Output error: could not write to {}sputtered.output.", options.name).as_str());
            }

            //Trajectory output
            if particle.track_trajectories {
                writeln!(trajectory_data_stream, "{}", particle.trajectory.len())
                    .expect(format!("Output error: could not write to {}trajectory_data.output.", options.name).as_str());

                for pos in particle.trajectory {
                    writeln!(
                        trajectory_file_stream, "{},{},{},{},{},{}",
                        particle.m/mass_unit, particle.Z, pos.E/energy_unit,
                        pos.x/length_unit, pos.y/length_unit, pos.z/length_unit,
                    ).expect(format!("Output error: could not write to {}trajectories.output.", options.name).as_str());
                }
            }

            if particle.incident & options.track_energy_losses {
                for energy_loss in particle.energies {
                    writeln!(
                        energy_loss_file_stream, "{},{},{},{},{},{},{}",
                        particle.m/mass_unit, particle.Z,
                        energy_loss.En/energy_unit, energy_loss.Ee/energy_unit,
                        energy_loss.x/length_unit, energy_loss.y/length_unit, energy_loss.z/length_unit,
                    ).expect(format!("Output error: could not write to {}energy_loss.output.", options.name).as_str());
                }
            }
        }
    }

    //Flush all file streams before dropping to ensure all data is written
    reflected_file_stream.flush().unwrap();
    deposited_file_stream.flush().unwrap();
    sputtered_file_stream.flush().unwrap();
    trajectory_data_stream.flush().unwrap();
    trajectory_file_stream.flush().unwrap();

    println!("Finished!");
}
