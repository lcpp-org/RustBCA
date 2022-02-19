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
use anyhow::{Result, Context, anyhow};
//use anyhow::*;

//Serializing/Deserializing crate
use serde::*;

//Array input via hdf5
#[cfg(feature = "hdf5_input")]
use hdf5::*;

//I/O
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::BufWriter;

//C integer
use std::os::raw::c_int;

//standard slice
use std::slice;

//itertools
use itertools::izip;

//Math
use std::f64::consts::FRAC_2_SQRT_PI;
use std::f64::consts::PI;
use std::f64::consts::SQRT_2;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
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

#[derive(Debug)]
#[repr(C)]
pub struct InputTaggedBCA {
    pub len: usize,
    /// x y z
    pub positions: *mut [f64; 3],
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
    pub tags: *mut i32,
    pub weights: *mut f64,
}

#[derive(Debug)]
#[repr(C)]
pub struct OutputTaggedBCA {
    pub len: usize,
    pub particles: *mut [f64; 9],
    pub weights: *mut f64,
    pub tags: *mut i32,
    pub incident: *mut bool,
}

#[no_mangle]
pub extern "C" fn compound_tagged_bca_list_c(input: InputTaggedBCA) -> OutputTaggedBCA {

    let mut total_output = vec![];
    let mut output_tags = vec![];
    let mut output_weights = vec![];
    let mut output_incident = vec![];

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
    let positions = unsafe { slice::from_raw_parts(input.positions, input.num_species_target).to_vec() };
    let tags = unsafe { slice::from_raw_parts(input.tags, input.num_species_target).to_vec() };
    let weights = unsafe { slice::from_raw_parts(input.weights, input.num_species_target).to_vec() };

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

    let mut index: usize = 0;
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
            interaction_index : 0,
            weight: weights[index],
            tag: tags[index],
            tracked_vector: Vector::new(positions[index][0], positions[index][1], positions[index][2]),
        };

        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {

            if (particle.left) | (particle.incident) {
                total_output.push(
                    [
                        particle.Z,
                        particle.m/AMU,
                        particle.E/EV,

                        particle.tracked_vector.x/ANGSTROM,
                        particle.tracked_vector.y/ANGSTROM,
                        particle.tracked_vector.z/ANGSTROM,

                        particle.dir.x,
                        particle.dir.y,
                        particle.dir.z
                    ]
                );
                output_tags.push(particle.tag);
                output_weights.push(particle.weight);
                output_incident.push(particle.incident);
            }
        }
        index += 1;
    }

    let len = total_output.len();
    let particles = total_output.as_mut_ptr();
    let tags_ptr = output_tags.as_mut_ptr();
    let weights_ptr = output_weights.as_mut_ptr();
    let incident_ptr = output_incident.as_mut_ptr();

    std::mem::forget(total_output);

    OutputTaggedBCA {
        len,
        particles,
        tags: tags_ptr,
        weights: weights_ptr,
        incident: incident_ptr,
    }
}


#[no_mangle]
#[cfg(not(feature = "distributions"))]
pub extern "C" fn reflect_single_ion_c(num_species_target: &mut c_int, ux: &mut f64, uy: &mut f64, uz: &mut f64, E1: &mut f64, Z1: &mut f64, m1: &mut f64, Ec1: &mut f64, Es1: &mut f64, Z2: *mut f64, m2: *mut f64, Ec2: *mut f64, Es2: *mut f64, Eb2: *mut f64, n2: *mut f64) {

    assert!(E1 > &mut 0.0);

    let options = Options {
        name: "test".to_string(),
        track_trajectories: false,
        track_recoil_trajectories: false,
        track_recoils: false,
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

    let Z2 = unsafe { slice::from_raw_parts(Z2, *num_species_target as usize).to_vec() };
    let m2 = unsafe { slice::from_raw_parts(m2, *num_species_target as usize).to_vec() };
    let n2 = unsafe { slice::from_raw_parts(n2, *num_species_target as usize).to_vec() };
    let Ec2 = unsafe { slice::from_raw_parts(Ec2, *num_species_target as usize).to_vec() };
    let Es2 = unsafe { slice::from_raw_parts(Es2, *num_species_target as usize).to_vec() };
    let Eb2 = unsafe { slice::from_raw_parts(Eb2, *num_species_target as usize).to_vec() };

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
        interaction_index: vec![0; *num_species_target as usize],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: n2,
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let p = particle::Particle {
        m: *m1*AMU,
        Z: *Z1,
        E: *E1*EV,
        Ec: *Ec1*EV,
        Es: *Es1*EV,
        pos: Vector::new(x, y, z),
        dir: Vector::new(*ux, *uy, *uz),
        pos_origin: Vector::new(x, y, z),
        pos_old: Vector::new(x, y, z),
        dir_old: Vector::new(*ux, *uy, *uz),
        energy_origin: *E1*EV,
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
        interaction_index : 0,
        weight: 1.0,
        tag: 0,
        tracked_vector: Vector::new(0.0, 0.0, 0.0),
    };

    let output = bca::single_ion_bca(p, &m, &options);

    *ux = output[0].dir.x;
    *uy = output[0].dir.y;
    *uz = output[0].dir.z;
    if output[0].pos.x >= 0.0 {
        *E1 = 0.0
    } else {
        *E1 = output[0].E/EV;
    }
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
            interaction_index : 0,
            weight: 1.0,
            tag: 0,
            tracked_vector: Vector::new(0.0, 0.0, 0.0),
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
            interaction_index : 0,
            weight: 1.0,
            tag: 0,
            tracked_vector: Vector::new(0.0, 0.0, 0.0),
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
pub extern "C" fn compound_bca_list_fortran(num_incident_ions: &mut c_int, track_recoils: &mut bool,
    ux: *mut f64, uy: *mut f64, uz: *mut f64, E1: *mut f64,
    Z1: *mut f64, m1: *mut f64, Ec1: *mut f64, Es1: *mut f64,
    num_species_target: &mut c_int,
    Z2: *mut f64, m2: *mut f64, Ec2: *mut f64, Es2: *mut f64, Eb2: *mut f64, n2: *mut f64,
    num_emitted_particles: &mut c_int
    ) -> *const [f64; 6] {

    //println!("{} {}", num_incident_ions, num_species_target);

    let mut total_output = vec![];

    #[cfg(feature = "distributions")]
    let options = Options {
        name: "test".to_string(),
        track_trajectories: false,
        track_recoils: *track_recoils,
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
        track_recoils: *track_recoils,
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

    let ux = unsafe { slice::from_raw_parts(ux, *num_incident_ions as usize).to_vec() };
    let uy = unsafe { slice::from_raw_parts(uy, *num_incident_ions as usize).to_vec() };
    let uz = unsafe { slice::from_raw_parts(uz, *num_incident_ions as usize).to_vec() };
    let Z1 = unsafe { slice::from_raw_parts(Z1, *num_incident_ions as usize).to_vec() };
    let m1 = unsafe { slice::from_raw_parts(m1, *num_incident_ions as usize).to_vec() };
    let E1 = unsafe { slice::from_raw_parts(E1, *num_incident_ions as usize).to_vec() };
    let Ec1 = unsafe { slice::from_raw_parts(Ec1, *num_incident_ions as usize).to_vec() };
    let Es1 = unsafe { slice::from_raw_parts(Es1, *num_incident_ions as usize).to_vec() };

    //println!("ux: {} uy: {} uz: {} Z1: {} m1: {} E1: {} Ec1: {} Es1: {}", ux[0], uy[0], uz[0], Z1[0], m1[0], E1[0], Ec1[0], Es1[0]);

    let Z2 = unsafe { slice::from_raw_parts(Z2, *num_species_target as usize).to_vec() };
    let m2 = unsafe { slice::from_raw_parts(m2, *num_species_target as usize).to_vec() };
    let n2 = unsafe { slice::from_raw_parts(n2, *num_species_target as usize).to_vec() };
    let Ec2 = unsafe { slice::from_raw_parts(Ec2, *num_species_target as usize).to_vec() };
    let Es2 = unsafe { slice::from_raw_parts(Es2, *num_species_target as usize).to_vec() };
    let Eb2 = unsafe { slice::from_raw_parts(Eb2, *num_species_target as usize).to_vec() };

    //println!("Z2: {} m2: {} n2: {} Ec2: {} Es2: {} Eb2: {}", Z2[0], m2[0], n2[0], Ec2[0], Es2[0], Eb2[0]);

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
        interaction_index: vec![0; *num_species_target as usize],
        surface_binding_model: SurfaceBindingModel::INDIVIDUAL,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: n2,
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    for (((((((E1_, ux_), uy_), uz_), Z1_), Ec1_), Es1_), m1_) in E1.iter().zip(ux).zip(uy).zip(uz).zip(Z1).zip(Ec1).zip(Es1).zip(m1) {

        let p = particle::Particle {
            m: m1_*AMU,
            Z: Z1_,
            E: *E1_*EV,
            Ec: Ec1_*EV,
            Es: Es1_*EV,
            pos: Vector::new(x, y, z),
            dir: Vector::new(ux_, uy_, uz_),
            pos_origin: Vector::new(x, y, z),
            pos_old: Vector::new(x, y, z),
            dir_old: Vector::new(ux_, uy_, uz_),
            energy_origin: *E1_*EV,
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
            interaction_index : 0,
            weight: 1.0,
            tag: 0,
            tracked_vector: Vector::new(0.0, 0.0, 0.0)
        };

        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {

            if (particle.left) | (particle.incident) {
                total_output.push(
                    [
                        particle.Z,
                        particle.m/AMU,
                        particle.E/EV,
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

    *num_emitted_particles = len as c_int;
    particles
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

#[cfg(feature = "python")]
/// simple_bca_py( x, y, z, ux, uy, uz, energy, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2)
/// --
///
/// This function runs a 0D Binary Collision Approximation simulation for the given single incident ion and material.
/// Args:
///    x (f64): initial ion position x. Material target is x>0
///    y (f64): initial ion position y.
///    z (f64): initial ion position z.
///    ux (f64): initial ion direction x. ux != 0.0 to avoid gimbal lock
///    uy (f64): initial ion direction y.
///    uz (f64): initial ion direction z.
///    energy (f64): initial ion energy in eV.
///    Z1 (f64): initial ion atomic number.
///    m1 (f64): initial ion mass in amu.
///    Ec1 (f64): ion cutoff energy in eV. If ion energy < Ec1, it stops in the material.
///    Es1 (f64): ion surface binding energy. Assumed planar.
///    Z2 (f64): target material atomic number.
///    m2 (f64): target material mass in amu.
///    Ec2 (f64): target material cutoff energy in eV. If recoil energy < Ec2, it stops in the material.
///    Es2 (f64): target atom surface binding energy. Assumed planar.
///    n2 (f64): target material atomic number density in inverse cubic Angstroms.
///    Eb2 (f64): target material bulk binding energy in eV.
/// Returns:
///    output (NX9 list of f64): each row in the list represents an output particle (implanted,
///    sputtered, or reflected). Each row consists of:
///      [Z, m (amu), E (eV), x, y, z, (angstrom), ux, uy, uz]
#[pyfunction]
pub fn simple_bca_py(x: f64, y: f64, z: f64, ux: f64, uy: f64, uz: f64, E1: f64, Z1: f64, m1: f64, Ec1: f64, Es1: f64, Z2: f64, m2: f64, Ec2: f64, Es2: f64, n2: f64, Eb2: f64) -> Vec<[f64; 9]> {
    simple_bca(x, y, z, ux, uy, uz, E1, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2)
}

#[cfg(feature = "python")]
/// simple_bca_list_py( energies, ux, uy, uz, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2)
/// --
///
/// This function runs a 0D Binary Collision Approximation simulation for the given incident ions and material.
/// Args:
///    energy (list(f64)): initial energies in eV.
///    ux (list(f64)): initial ion directions x. ux != 0.0 to avoid gimbal lock
///    uy (list(f64)): initial ion directions y.
///    uz (list(f64)): initial ion directions z.
///    Z1 (f64): initial ion atomic number.
///    m1 (f64): initial ion mass in amu.
///    Ec1 (f64): ion cutoff energy in eV. If ion energy < Ec1, it stops in the material.
///    Es1 (f64): ion surface binding energy. Assumed planar.
///    Z2 (f64): target material atomic number.
///    m2 (f64): target material mass in amu.
///    Ec2 (f64): target material cutoff energy in eV. If recoil energy < Ec2, it stops in the material.
///    Es2 (f64): target atom surface binding energy. Assumed planar.
///    n2 (f64): target material atomic number density in inverse cubic Angstroms.
///    Eb2 (f64): target material bulk binding energy in eV.
/// Returns:
///    output (NX9 list of f64): each row in the list represents an output particle (implanted,
///    sputtered, or reflected). Each row consists of:
///      [Z, m (amu), E (eV), x, y, z, (angstrom), ux, uy, uz]
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
        interaction_index : 0,
        weight: 1.0,
        tag: 0,
        tracked_vector: Vector::new(0.0, 0.0, 0.0),
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
