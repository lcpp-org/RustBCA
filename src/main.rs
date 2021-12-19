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
pub mod sphere;
pub mod physics;

#[cfg(feature = "parry3d")]
pub mod parry;

pub use crate::enums::*;
pub use crate::consts::*;
pub use crate::structs::*;
pub use crate::input::{Input2D, Input1D, Input0D, Options, InputFile, GeometryInput};
pub use crate::output::{OutputUnits};
pub use crate::geometry::{Geometry, GeometryElement, Mesh0D, Mesh1D, Mesh2D};
pub use crate::sphere::{Sphere, SphereInput, InputSphere};
pub use crate::physics::{physics_loop};

#[cfg(feature = "parry3d")]
pub use crate::parry::{ParryBall, ParryBallInput, InputParryBall, ParryTriMesh, ParryTriMeshInput, InputParryTriMesh};

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
            _ => panic!("Unimplemented geometry {}.", args[1].clone())
        }),
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
        GeometryType::SPHERE => {
            let (particle_input_array, material, options, output_units) = input::input::<Sphere>(input_file);
            physics_loop::<Sphere>(particle_input_array, material, options, output_units);
        },
        #[cfg(feature = "parry3d")]
        GeometryType::BALL => {
            let (particle_input_array, material, options, output_units) = input::input::<ParryBall>(input_file);
            physics_loop::<ParryBall>(particle_input_array, material, options, output_units);
        }
        #[cfg(feature = "parry3d")]
        GeometryType::TRIMESH => {
            let (particle_input_array, material, options, output_units) = input::input::<ParryTriMesh>(input_file);
            physics_loop::<ParryTriMesh>(particle_input_array, material, options, output_units);
        }
    }
}
