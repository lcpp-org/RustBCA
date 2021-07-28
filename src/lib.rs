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

//standard slice
use std::slice;

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
    pub len: usize,
    pub particles: *mut [f64; 9],
}

#[derive(Debug)]
#[repr(C)]
pub struct InputSimpleBCA {
    pub len: usize,
    /// vx, vy, vz
    pub velocities: *mut [f64; 3],
    pub Z1: f64,
    pub m1: f64,
    pub Ec1: f64,
    pub Es1: f64,
    pub Z2: f64,
    pub m2: f64,
    pub n2: f64,
    pub Ec2: f64,
    pub Es2: f64,
    pub Eb2: f64,
}

#[derive(Debug)]
#[repr(C)]
pub struct InputCompoundBCA {
    pub len: usize,
    /// vx, vy, vz
    pub velocities: *mut [f64; 3],
    pub Z1: f64,
    pub m1: f64,
    pub Ec1: f64,
    pub Es1: f64,
    pub num_species_target: usize,
    pub Z2: *mut f64,
    pub m2: *mut f64,
    pub n2: *mut f64,
    pub Ec2: *mut f64,
    pub Es2: *mut f64,
    pub Eb2: *mut f64,
}

#[no_mangle]
pub extern "C" fn simple_bca_list_c(input: InputSimpleBCA) -> OutputBCA {

    let x = -2.*(input.n2*10E30).powf(-1./3.);
    let y = 0.0;
    let z = 0.0;

    let mut total_output = vec![];

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

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![input.Eb2],
        Es: vec![input.Es2],
        Ec: vec![input.Ec2],
        Z: vec![input.Z2],
        m: vec![input.m2],
        interaction_index: vec![0],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: vec![input.n2],
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let velocities = unsafe { slice::from_raw_parts(input.velocities, input.len) };

    for velocity in velocities {

        let vx = velocity[0];
        let vy = velocity[1];
        let vz = velocity[2];

        let v = (vx*vx + vy*vy + vz*vz).sqrt();

        let E1 = 0.5*input.m1*AMU*v*v;

        let ux = vx/v;
        let uy = vy/v;
        let uz = vz/v;

        let p = particle::Particle {
            m: input.m1*AMU,
            Z: input.Z1,
            E: E1,
            Ec: input.Ec1*EV,
            Es: input.Es1*EV,
            pos: Vector::new(x, y, z),
            dir: Vector::new(ux, uy, uz),
            pos_origin: Vector::new(x, y, z),
            pos_old: Vector::new(x, y, z),
            dir_old: Vector::new(ux, uy, uz),
            energy_origin: E1,
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


        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {

            if (particle.left) | (particle.incident) {
                total_output.push(
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
                );
            }
        }
    }

    let len = total_output.len();
    let particles = total_output.as_mut_ptr();

    std::mem::forget(total_output);
    OutputBCA {
        len,
        particles
    }
}

#[no_mangle]
pub extern "C" fn compound_bca_list_c(input: InputCompoundBCA) -> OutputBCA {

    let mut total_output = vec![];

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

    let Z2 = unsafe { slice::from_raw_parts(input.Z2, input.num_species_target).to_vec() };
    let m2 = unsafe { slice::from_raw_parts(input.m2, input.num_species_target).to_vec() };
    let n2 = unsafe { slice::from_raw_parts(input.n2, input.num_species_target).to_vec() };
    let Ec2 = unsafe { slice::from_raw_parts(input.Ec2, input.num_species_target).to_vec() };
    let Es2 = unsafe { slice::from_raw_parts(input.Es2, input.num_species_target).to_vec() };
    let Eb2 = unsafe { slice::from_raw_parts(input.Eb2, input.num_species_target).to_vec() };

    let Z2 = vec![1.0, 1.0];
    let m2 = vec![1.0, 1.0];
    let n2 = vec![0.06, 0.05];
    let Ec2 = vec![1.0, 1.0];
    let Es2 = vec![1.0, 1.0];
    let Eb2 = vec![1.0, 1.0];

    let x = -2.*(n2.iter().sum::<f64>()*10E30).powf(-1./3.);
    let y = 0.0;
    let z = 0.0;

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: Eb2,
        Es: Es2,
        Ec: Ec2,
        Z: Z2,
        m: m2,
        interaction_index: vec![0; input.num_species_target],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: n2,
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let velocities = unsafe { slice::from_raw_parts(input.velocities, input.len) };

    for velocity in velocities {

        let vx = velocity[0];
        let vy = velocity[1];
        let vz = velocity[2];

        let v = (vx*vx + vy*vy + vz*vz).sqrt();

        let E1 = 0.5*input.m1*AMU*v*v;

        let ux = vx/v;
        let uy = vy/v;
        let uz = vz/v;

        let p = particle::Particle {
            m: input.m1*AMU,
            Z: input.Z1,
            E: E1,
            Ec: input.Ec1*EV,
            Es: input.Es1*EV,
            pos: Vector::new(x, y, z),
            dir: Vector::new(ux, uy, uz),
            pos_origin: Vector::new(x, y, z),
            pos_old: Vector::new(x, y, z),
            dir_old: Vector::new(ux, uy, uz),
            energy_origin: E1,
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


        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {

            if (particle.left) | (particle.incident) {
                total_output.push(
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
                );
            }
        }
    }

    let len = total_output.len();
    let particles = total_output.as_mut_ptr();

    std::mem::forget(total_output);
    OutputBCA {
        len,
        particles
    }
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
    assert!(Ec1 > 0.0, "Error: Cutoff energy Ec1 cannot be less than or equal to 0.");
    assert!(Ec2 > 0.0, "Error: Cutoff energy Ec2 cannot be less than or equal to 0.");

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