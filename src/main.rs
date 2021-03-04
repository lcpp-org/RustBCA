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
use indicatif::{ProgressBar, ProgressStyle};

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

//itertools
use itertools::izip;

//Math
use std::f64::consts::FRAC_2_SQRT_PI;
use std::f64::consts::PI;
use std::f64::consts::SQRT_2;

//rng
//use rand::{Rng, thread_rng};

//Load internal modules
pub mod material;
pub mod particle;
pub mod tests;
pub mod interactions;
pub mod bca;
pub mod geometry;
pub mod input;
pub mod output;
pub mod enums;
pub mod consts;
pub mod structs;

pub use crate::enums::*;
pub use crate::consts::*;
pub use crate::structs::*;
pub use crate::input::{Input2D, Input1D, Input0D, Options, InputFile, GeometryInput};
pub use crate::output::{OutputUnits};
pub use crate::geometry::{Geometry, GeometryElement, Mesh0D, Mesh1D, Mesh2D};

fn physics_loop<T: Geometry + Sync>(particle_input_array: Vec<particle::ParticleInput>, material: material::Material<T>, options: Options, output_units: OutputUnits) {

        println!("Processing {} ions...", particle_input_array.len());

        let total_count: u64 = particle_input_array.len() as u64;
        assert!(total_count/options.num_chunks > 0, "Input error: chunk size == 0 - reduce num_chunks or increase particle count.");

        #[cfg(not(feature = "no_list_output"))]
        let mut output_list_streams = output::open_output_lists(&options);

        let mut summary = output::SummaryPerSpecies::new(&options);

        #[cfg(feature = "distributions")]
        let mut distributions = output::Distributions::new(&options);

        //Initialize threads with rayon
        println!("Initializing with {} threads...", options.num_threads);
        if options.num_threads > 1 {let pool = rayon::ThreadPoolBuilder::new().num_threads(options.num_threads).build_global().unwrap();};

        //Create and configure progress bar
        let bar: ProgressBar = ProgressBar::new(total_count);
        bar.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}][{bar:40.cyan/blue}][{eta_precise}] {percent}%")
            .progress_chars("#>-"));

        //Main loop
        for (chunk_index, particle_input_chunk) in particle_input_array.chunks((total_count/options.num_chunks) as usize).enumerate() {

            let mut finished_particles: Vec<particle::Particle> = Vec::new();

            if options.num_threads > 1 {
                // BCA loop is implemented as parallelized extension of a per-chunk initially empty
                // finished particle array via map from particle -> finished particles via BCA
                finished_particles.par_extend(
                    particle_input_chunk.into_par_iter()
                    .map(|particle_input| {
                        bar.tick();
                        bar.inc(1);
                        bca::single_ion_bca(particle::Particle::from_input(*particle_input, &options), &material, &options)
                    }).flatten()
                );
            } else {
                finished_particles.extend(
                    particle_input_chunk.iter()
                    .map(|particle_input| {
                        bar.tick();
                        bar.inc(1);
                        bca::single_ion_bca(particle::Particle::from_input(*particle_input, &options), &material, &options)
                    }).flatten()
                );
            }

            // Process this chunk of finished particles for output
            for particle in finished_particles {

                summary.update(&particle);

                #[cfg(feature = "distributions")]
                distributions.update(&particle, &output_units, &options, total_count as usize);

                #[cfg(not(feature = "no_list_output"))]
                output::output_lists(&mut output_list_streams, particle, &options, &output_units);

            }
            //Flush all file streams before dropping to ensure all data is written
            #[cfg(not(feature = "no_list_output"))]
            output::output_list_flush(&mut output_list_streams);
        }

        summary.print(&options, &output_units);

        //Write distributions to file
        #[cfg(feature = "distributions")]
        distributions.print(&options);

        println!("Finished!");
}

fn main() {

    let args: Vec<String> = env::args().collect();

    let (input_file, geometry_type) = match args.len() {
        1 => ("input.toml".to_string(), GeometryType::MESH2D),
        2 => (args[1].clone(), GeometryType::MESH2D),
        3 => (args[2].clone(), match args[1].as_str() { "0D" => GeometryType::MESH0D, "1D" => GeometryType::MESH1D, "2D" => GeometryType::MESH2D, _ => panic!("Unimplemented geometry {}.", args[1].clone()) }),
        _ => panic!("Too many command line arguments. RustBCA accepts 0 (use 'input.toml') 1 (<input file name>) or 2 (<geometry type> <input file name>)"),
    };

     match geometry_type {
        GeometryType::MESH0D => {
            let (particle_input_array, material, options, output_units) = input::input::<geometry::Mesh0D>(input_file);
            physics_loop::<Mesh0D>(particle_input_array, material, options, output_units);
        },
        GeometryType::MESH1D => {
            let (particle_input_array, material, options, output_units) = input::input::<geometry::Mesh1D>(input_file);
            physics_loop::<Mesh1D>(particle_input_array, material, options, output_units);
        },
        GeometryType::MESH2D => {
            let (particle_input_array, material, options, output_units) = input::input::<geometry::Mesh2D>(input_file);
            physics_loop::<Mesh2D>(particle_input_array, material, options, output_units);
        },
    }
}
