#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

#[cfg(feature = "cpr_rootfinder_openblas")]
extern crate openblas_src;
#[cfg(feature = "cpr_rootfinder_netlib")]
extern crate netlib_src;
#[cfg(feature = "cpr_rootfinder_intel_mkl")]
extern crate intel_mkl_src;

use std::{fmt};
use std::mem::discriminant;

//Error handling crate
use anyhow::Result;
use anyhow::*;

//Serializing/Deserializing crate
use serde::*;

//Array input via hdf5
#[cfg(feature = "hdf5_input")]
use hdf5::*;

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

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

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

#[cfg(feature = "parry3d")]
pub mod parry;

pub use crate::enums::*;
pub use crate::consts::*;
pub use crate::structs::*;
pub use crate::input::{Input2D, Input1D, Input0D, Options, InputFile, GeometryInput};
pub use crate::output::{OutputUnits};
pub use crate::geometry::{Geometry, GeometryElement, Mesh0D, Mesh1D, Mesh2D};
pub use crate::sphere::{Sphere, SphereInput, InputSphere};

#[cfg(feature = "parry3d")]
pub use crate::parry::{ParryBall, ParryBallInput, InputParryBall, ParryTriMesh, ParryTriMeshInput, InputParryTriMesh};

#[pymodule]
pub fn pybca(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(simple_bca_py, m)?)?;
    m.add_function(wrap_pyfunction!(simple_bca_list_py, m)?)?;
    Ok(())
}

#[repr(C)]
pub struct OutputBCA {
    len: usize,
    pub particles: *mut [f64; 9],
}

#[no_mangle]
pub extern "C" fn simple_bca_c(x: f64, y: f64, z: f64, ux: f64, uy: f64, uz: f64, E1: f64, Z1: f64, m1: f64, Ec1: f64, Es1: f64, Z2: f64, m2: f64, Ec2: f64, Es2: f64, n2: f64, Eb2: f64) -> OutputBCA {
    let mut output = simple_bca(x, y, z, ux, uy, uz, E1, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2);
    let len = output.len();
    let particles = output.as_mut_ptr();

    std::mem::forget(output);
    OutputBCA {
        len,
        particles
    }
}

#[pyfunction]
pub fn simple_bca_py(x: f64, y: f64, z: f64, ux: f64, uy: f64, uz: f64, E1: f64, Z1: f64, m1: f64, Ec1: f64, Es1: f64, Z2: f64, m2: f64, Ec2: f64, Es2: f64, n2: f64, Eb2: f64) -> Vec<[f64; 9]> {
    simple_bca(x, y, z, ux, uy, uz, E1, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2)
}

#[pyfunction]
pub fn simple_bca_list_py(energies: Vec<f64>, usx: Vec<f64>, usy: Vec<f64>, usz: Vec<f64>, Z1: f64, m1: f64, Ec1: f64, Es1: f64, Z2: f64, m2: f64, Ec2: f64, Es2: f64, n2: f64, Eb2: f64) -> Vec<[f64; 9]> {

    assert_eq!(energies.len(), usx.len());
    assert_eq!(energies.len(), usy.len());
    assert_eq!(energies.len(), usz.len());

    let x = -2.*(n2*10E30).powf(-1./3.);
    let y = 0.0;
    let z = 0.0;

    let mut total_output = vec![];
    for (((E1, ux), uy), uz) in energies.iter().zip(usx).zip(usy).zip(usz) {
        let output = simple_bca(x, y, z, ux, uy, uz, *E1, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2);
        for particle in output {
            total_output.push(particle);
        }
    }
    total_output
}

pub fn simple_bca(x: f64, y: f64, z: f64, ux: f64, uy: f64, uz: f64, E1: f64, Z1: f64, m1: f64, Ec1: f64, Es1: f64, Z2: f64, m2: f64, Ec2: f64, Es2: f64, n2: f64, Eb2: f64) -> Vec<[f64; 9]> {

    assert!(E1 > 0.0, "Error: Incident energy cannot be less than or equal to 0.");
    assert!(Ec1 > 0.0, "Error: Cutoff energy cannot be less than or equal to 0.");
    assert!(Ec2 > 0.0, "Error: Cutoff energy cannot be less than or equal to 0.");

    #[cfg(feature = "distributions")]
    let options = Options {
        name: "test".to_string(),
        track_trajectories: false,
        track_recoils: true,
        track_recoil_trajectories: false,
        write_buffer_size: 8000,
        weak_collision_order: 3,
        suppress_deep_recoils: false,
        high_energy_free_flight_paths: false,
        electronic_stopping_mode: ElectronicStoppingMode::INTERPOLATED,
        mean_free_path_model: MeanFreePathModel::LIQUID,
        interaction_potential: vec![vec![InteractionPotential::KR_C]],
        scattering_integral: vec![vec![ScatteringIntegral::MENDENHALL_WELLER]],
        num_threads: 1,
        num_chunks: 1,
        use_hdf5: false,
        root_finder: vec![vec![Rootfinder::DEFAULTNEWTON]],
        track_displacements: false,
        track_energy_losses: false,
        energy_min: 0.0,
        energy_max: 10.0,
        energy_num: 11,
        angle_min: 0.0,
        angle_max: 90.0,
        angle_num: 11,
        x_min: 0.0,
        y_min: -10.0,
        z_min: -10.0,
        x_max: 10.0,
        y_max: 10.0,
        z_max: 10.0,
        x_num: 11,
        y_num: 11,
        z_num: 11,
    };

    #[cfg(not(feature = "distributions"))]
    let options = Options {
        name: "test".to_string(),
        track_trajectories: false,
        track_recoils: true,
        track_recoil_trajectories: false,
        write_buffer_size: 8000,
        weak_collision_order: 3,
        suppress_deep_recoils: false,
        high_energy_free_flight_paths: false,
        electronic_stopping_mode: ElectronicStoppingMode::INTERPOLATED,
        mean_free_path_model: MeanFreePathModel::LIQUID,
        interaction_potential: vec![vec![InteractionPotential::KR_C]],
        scattering_integral: vec![vec![ScatteringIntegral::MENDENHALL_WELLER]],
        num_threads: 1,
        num_chunks: 1,
        use_hdf5: false,
        root_finder: vec![vec![Rootfinder::DEFAULTNEWTON]],
        track_displacements: false,
        track_energy_losses: false,
    };

    let p = particle::Particle {
        m: m1*AMU,
        Z: Z1,
        E: E1*EV,
        Ec: Ec1*EV,
        Es: Es1*EV,
        pos: Vector::new(x, y, z),
        dir: Vector::new(ux, uy, uz),
        pos_origin: Vector::new(x, y, z),
        pos_old: Vector::new(x, y, z),
        dir_old: Vector::new(ux, uy, uz),
        energy_origin: E1*EV,
        asymptotic_deflection: 0.0,
        stopped: false,
        left: false,
        incident: true,
        first_step: true,
        trajectory: vec![],
        energies: vec![],
        track_trajectories: false,
        number_collision_events: 0,
        backreflected: false,
        interaction_index : 0
    };

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![Eb2],
        Es: vec![Es2],
        Ec: vec![Ec2],
        Z: vec![Z2],
        m: vec![m2],
        interaction_index: vec![0],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: vec![n2],
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let output = bca::single_ion_bca(p, &m, &options);

    output.iter().filter(|particle| (particle.incident) | (particle.left)).map(|particle|
        [
            particle.Z,
            particle.m/AMU,
            particle.E/EV,
            particle.pos.x/ANGSTROM,
            particle.pos.y/ANGSTROM,
            particle.pos.z/ANGSTROM,
            particle.dir.x,
            particle.dir.y,
            particle.dir.z
        ]
    ).collect()
}
