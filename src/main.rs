#![allow(unused_variables)]
#![allow(non_snake_case)]

use itertools::Itertools;

//extern crate openblas_src;
//Progress bar crate - works wiht rayon
use indicatif::{ParallelProgressIterator, ProgressIterator, ProgressBar, ProgressStyle};

//Error handling crate
//use anyhow::Result;
use anyhow::Result;
use anyhow::*;

//Geometry crate
use geo::algorithm::contains::Contains;
//use geo::algorithm::convexhull::ConvexHull;
use geo::algorithm::closest_point::ClosestPoint;
use geo::{point, Polygon, LineString, Closest, Point};// MultiPoint, Point};

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
//const K: f64 = 1.11265E-10;
const ME: f64 =  9.109383632E-31;
const SQRTPI: f64 = 1.772453850906;
const SQRT2PI: f64 = 2.506628274631;
const C: f64 = 299792000.;
const BETHE_BLOCH_PREFACTOR: f64 = 4.*PI*(Q*Q/(4.*PI*EPS0))*(Q*Q/(4.*PI*EPS0))/ME/C/C;
const LINDHARD_SCHARFF_PREFACTOR: f64 = 1.212*ANGSTROM*ANGSTROM*Q;
const LINDHARD_REDUCED_ENERGY_PREFACTOR: f64 = 4.*PI*EPS0/Q/Q;

//Electronic stopping models
const INTERPOLATED: i32 = 0;
const LOW_ENERGY_NONLOCAL: i32 = 1;
const LOW_ENERGY_LOCAL: i32 = 2;
const LOW_ENERGY_EQUIPARTITION: i32 = 3;

//Mean free path models
const LIQUID: i32 = 0;
const GASEOUS: i32 = 1;

//Interaction potentials
const MOLIERE: i32 = 0;
const KR_C: i32 = 1;
const ZBL: i32 = 2;
const LENZ_JENSEN: i32 = 3;
const LENNARD_JONES_12_6: i32 = 4;
const TRIDYN: i32 = -1;
const BUCKINGHAM: i32 = 5;

//Scattering integral forms
const MENDENHALL_WELLER: i32 = 0;
const MAGIC: i32 = 1;
const GAUSS_MEHLER: i32 = 2;
const GAUSS_LEGENDRE: i32 = 3;

//Distance of closest approach determination
const NEWTON: i32 = 0;
const CPR: i32 = 1;
const POLYNOMIAL: i32 = 2;

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
    write_files: bool,
    stream_size: usize,
    print: bool,
    print_num: usize,
    weak_collision_order: usize,
    suppress_deep_recoils: bool,
    high_energy_free_flight_paths: bool,
    electronic_stopping_mode: i32,
    mean_free_path_model: i32,
    interaction_potential: i32,
    scattering_integral: i32,
    tolerance: f64,
    max_iterations: usize,
    num_threads: usize,
    num_chunks: u64,
    use_hdf5: bool,
    root_finder: i32,
    cpr_n0: usize,
    cpr_nmax: usize,
    cpr_epsilon: f64,
    cpr_complex: f64,
    cpr_truncation: f64,
    cpr_far_from_zero: f64,
    cpr_interval_limit: f64,
    cpr_upper_bound_const: f64,
    polynom_complex_threshold: f64
}

fn main() {

    //Read input file, convert to string, and open with toml
    let mut input_toml = String::new();
    let mut file = OpenOptions::new()
        .read(true)
        .write(false)
        .create(false)
        .open("input.toml")
        .expect("Input errror: could not open input file, input.toml.");
    file.read_to_string(&mut input_toml).unwrap();
    let input: Input = toml::from_str(&input_toml).unwrap();

    //Unpack toml information into structs
    let material = material::Material::new(input.material_parameters, input.mesh_2d_input);
    assert!(material.n.len() == material.m.len(), "Input error: material input arrays of unequal length.");
    assert!(material.n.len() == material.Z.len(), "Input error: material input arrays of unequal length.");
    assert!(material.n.len() == material.Eb.len(), "Input error: material input arrays of unequal length.");
    assert!(material.n.len() == material.Es.len(), "Input error: material input arrays of unequal length.");

    let options = input.options;
    let particle_parameters = input.particle_parameters;

    //Check that incompatible options are not on simultaneously
    if options.high_energy_free_flight_paths {
        assert!(options.electronic_stopping_mode == INTERPOLATED,
            "Input error: High energy free flight paths used with low energy stoppping power.");
    }

    if options.electronic_stopping_mode == INTERPOLATED {
        assert!(options.weak_collision_order == 0,
            "Input error: Cannot use weak collisions with free flight paths. Set weak_collision_order = 0.");
    }

    if options.mean_free_path_model == GASEOUS {
        assert!(options.weak_collision_order == 0,
            "Input error: Cannot use weak collisions with gaseous mean free path model. Set weak_collision_order = 0.");
    }

    if options.interaction_potential == LENNARD_JONES_12_6 {
        assert!((options.scattering_integral == GAUSS_MEHLER) | (options.scattering_integral == GAUSS_LEGENDRE),
        "Input error: Cannot use scattering integral {} with interaction potential {}. Use scattering integral 2: Gauss-Mehler or 3: Gauss-Legendre.",
        options.scattering_integral, options.interaction_potential);

        assert!((options.root_finder == CPR) | (options.root_finder == POLYNOMIAL),
        "Input error: Cannot use Newton root-finder with attractive-repulsive potentials. Use 1: Chebyshev-Proxy Rootfinder. or 2: Polynomial (if interatomic potential has only power-of-r terms.)");
    }

    if (options.scattering_integral == MENDENHALL_WELLER) | (options.scattering_integral == MAGIC) {
        assert!(match options.interaction_potential {
            MOLIERE | ZBL | KR_C | LENZ_JENSEN | TRIDYN => true,
            _ => false
        }, "Input error: Mendenhall-Weller quadrature and Magic formula can only be used with screened Coulomb potentials. Use scattering integral 2: Gauss-Mehler or 3: Gauss-Legendre.")
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
            let particle_input_file = hdf5::File::open(particle_input_filename).expect("Input error: cannot open HDF5 file.");
            let particle_input = particle_input_file.dataset("particles").unwrap();
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
                let (x, y, z) = particle_parameters.pos[particle_index];
                let (cosx, cosy, cosz) = particle_parameters.dir[particle_index];
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
                            uz: cosz
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
                let (x, y, z) = particle_parameters.pos[particle_index];
                let (cosx, cosy, cosz) = particle_parameters.dir[particle_index];
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
                            uz: cosz
                        }
                    );
                }
            }
            particle_input
        }
    };

    //Open output files for streaming output
    let reflected_file = OpenOptions::new()
        .write(true).
        create(true).
        open(format!("{}{}", options.name, "reflected.output")).
        unwrap();
    let mut reflected_file_stream = BufWriter::with_capacity(options.stream_size, reflected_file);

    let sputtered_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "sputtered.output"))
        .unwrap();
    let mut sputtered_file_stream = BufWriter::with_capacity(options.stream_size, sputtered_file);

    let deposited_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "deposited.output"))
        .unwrap();
    let mut deposited_file_stream = BufWriter::with_capacity(options.stream_size, deposited_file);

    let trajectory_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "trajectories.output"))
        .unwrap();
    let mut trajectory_file_stream = BufWriter::with_capacity(options.stream_size, trajectory_file);

    let trajectory_data = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "trajectory_data.output"))
        .unwrap();
    let mut trajectory_data_stream = BufWriter::with_capacity(options.stream_size, trajectory_data);

    println!("Processing {} ions...", particle_input_array.len());

    //Main loop
    let total_count: u64 = particle_input_array.len() as u64;

    //Initialize threads with rayon
    println!("Initializing with {} threads...", options.num_threads);
    if options.num_threads > 1 {let pool = rayon::ThreadPoolBuilder::new().num_threads(options.num_threads).build_global().unwrap();};

    let bar: ProgressBar = ProgressBar::new(options.num_chunks);

    bar.set_style(ProgressStyle::default_bar()
        .template("{spinnter:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {percent}%")
        .progress_chars("#>-"));

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
                particle_input_chunk.into_iter()
                .map(|particle_input| bca::single_ion_bca(particle::Particle::from_input(*particle_input, &options), &material, &options))
                    .flatten()
            );
        }

        for particle in finished_particles {
            if particle.incident & particle.left {
                writeln!(
                    reflected_file_stream, "{},{},{},{},{},{},{},{},{},{}",
                    particle.m/mass_unit, particle.Z, particle.E/energy_unit,
                    particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit,
                    particle.dir.x, particle.dir.y, particle.dir.z,
                    particle.number_collision_events
                ).expect("Output error: could not write to reflected.output.");
            }

            //Incident particle, stopped in material: deposited
            if particle.incident & particle.stopped {
                writeln!(
                    deposited_file_stream, "{},{},{},{},{},{}",
                    particle.m/mass_unit, particle.Z,
                    particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit,
                    particle.number_collision_events
                ).expect("Output error: could not write to deposited.output.");
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
                ).expect("Output error: could not write to sputtered.output.");
            }

            //Trajectory output
            if particle.track_trajectories {
                writeln!(trajectory_data_stream, "{}", particle.trajectory.len())
                    .expect("Output error: could not write trajectory length data.");

                for pos in particle.trajectory {
                    writeln!(
                        trajectory_file_stream, "{},{},{},{},{},{}",
                        particle.m/mass_unit, particle.Z, pos.E/energy_unit,
                        pos.x/length_unit, pos.y/length_unit, pos.z/length_unit,
                    ).expect("Output error: could not write to trajectories.output.");
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
