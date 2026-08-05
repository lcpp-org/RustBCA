#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

use std::{env, fmt};
use std::mem::discriminant;

//Progress bar crate - works with rayon
use indicatif::{ProgressBar, ProgressStyle};

//Error handling crate
use anyhow::{Result, Context, anyhow};

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

//RNG
use rand::{SeedableRng, rngs::ChaCha8Rng};

//Load internal modules
pub mod material;
pub mod particle;
pub mod interactions;
pub mod bca;
pub mod geometry;
pub mod input;
pub mod output;
pub mod enums;
pub mod consts;
pub mod structs;
pub mod sphere;
pub mod physics;
pub mod math;

#[cfg(feature = "parry3d")]
pub mod parry;

pub use crate::enums::*;
pub use crate::consts::*;
pub use crate::structs::*;
pub use crate::input::{Input2D, InputHomogeneous2D, Input1D, Input0D, Options, InputFile, GeometryInput};
pub use crate::output::{OutputUnits};
pub use crate::geometry::{Geometry, GeometryElement, Mesh0D, Mesh1D, Mesh2D, HomogeneousMesh2D};
pub use crate::sphere::{Sphere, SphereInput, InputSphere};
pub use crate::physics::{physics_loop};
pub use crate::math::duff_orthonormal_basis;

#[cfg(feature = "parry3d")]
pub use crate::parry::{ParryBall, ParryBallInput, InputParryBall, ParryTriMesh, ParryTriMeshInput, InputParryTriMesh};


macro_rules! main_loop {
    ($geometry_type:ident, $input_file:expr) => {
        {
            let (particle_input_array, material, options, output_units) = input::input::<$geometry_type>($input_file);
            //Initialize threads with rayon
            println!("Processing {} ions...", particle_input_array.len());
            println!("Initializing with {} threads...", options.num_threads);
            let _ = rayon::ThreadPoolBuilder::new().num_threads(options.num_threads).build_global();
            physics_loop::<$geometry_type>(particle_input_array, material, options, output_units);
        }
    }
}

fn main() {

    let args: Vec<String> = env::args().collect();

    let (input_file, geometry_type) = match args.len() {
        1 => ("input.toml".to_string(), GeometryType::MESH2D),
        2 => (args[1].clone(), GeometryType::MESH2D),
        3 => (args[2].clone(), match args[1].as_str() {
            "0D" => GeometryType::MESH0D,
            "1D" => GeometryType::MESH1D,
            "2D" => GeometryType::MESH2D,
            "SPHERE" => GeometryType::SPHERE,
            #[cfg(feature = "parry3d")]
            "BALL" => GeometryType::BALL,
            #[cfg(feature = "parry3d")]
            "TRIMESH" => GeometryType::TRIMESH,
            "HOMOGENEOUS2D" => GeometryType::HOMOGENEOUS2D,
            _ => panic!("Unimplemented geometry {}.", args[1].clone())
        }),
        _ => panic!("Too many command line arguments. RustBCA accepts 0 (use 'input.toml') 1 (<input file name>) or 2 (<geometry type> <input file name>)"),
    };

     match geometry_type {
        GeometryType::MESH0D => main_loop!(Mesh0D, input_file),
        GeometryType::MESH1D => main_loop!(Mesh1D, input_file),
        GeometryType::MESH2D => main_loop!(Mesh2D, input_file),
        GeometryType::SPHERE => main_loop!(Sphere, input_file),
        #[cfg(feature = "parry3d")]
        GeometryType::BALL => main_loop!(ParryBall, input_file),
        #[cfg(feature = "parry3d")]
        GeometryType::TRIMESH => main_loop!(ParryTriMesh, input_file),
        GeometryType::HOMOGENEOUS2D => main_loop!(HomogeneousMesh2D, input_file),
    }
}
