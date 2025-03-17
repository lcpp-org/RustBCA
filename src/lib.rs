#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

use std::{fmt};
use std::mem::discriminant;

use std::alloc::{dealloc, Layout};
use std::mem::align_of;

//Parallelization - currently only used in python library functions
#[cfg(feature = "python")]
use rayon::prelude::*;
#[cfg(feature = "python")]
use rayon::*;

//Error handling crate
use anyhow::{Result, Context, anyhow};

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
//Mutex for multithreading in ergonomic Python library functions
#[cfg(feature = "python")]
use std::sync::Mutex;

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
#[cfg(feature = "python")]
use pyo3::types::*;

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
pub use crate::input::{Input2D, InputHomogeneous2D, Input1D, Input0D, Options, InputFile, GeometryInput};
pub use crate::output::{OutputUnits};
pub use crate::geometry::{Geometry, GeometryElement, Mesh0D, Mesh1D, Mesh2D};
pub use crate::sphere::{Sphere, SphereInput, InputSphere};

#[cfg(feature = "parry3d")]
pub use crate::parry::{ParryBall, ParryBallInput, InputParryBall, ParryTriMesh, ParryTriMeshInput, InputParryTriMesh};
#[cfg(feature = "parry3d")]
pub use parry3d_f64::na::{Point3, Vector3, Matrix3, Rotation3};

#[cfg(feature = "python")]
#[pymodule]
pub fn libRustBCA(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(simple_bca_py, m)?)?;
    m.add_function(wrap_pyfunction!(simple_bca_list_py, m)?)?;
    m.add_function(wrap_pyfunction!(compound_bca_list_py, m)?)?;
    m.add_function(wrap_pyfunction!(compound_bca_list_1D_py, m)?)?;
    m.add_function(wrap_pyfunction!(sputtering_yield, m)?)?;
    m.add_function(wrap_pyfunction!(reflection_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(compound_reflection_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(reflect_single_ion_py, m)?)?;
    #[cfg(feature = "parry3d")]
    m.add_function(wrap_pyfunction!(rotate_given_surface_normal_py, m)?)?;
    #[cfg(feature = "parry3d")]
    m.add_function(wrap_pyfunction!(rotate_given_surface_normal_vec_py, m)?)?;
    #[cfg(feature = "parry3d")]
    m.add_function(wrap_pyfunction!(rotate_back_py, m)?)?;
    #[cfg(feature = "parry3d")]
    m.add_function(wrap_pyfunction!(rotate_back_vec_py, m)?)?;
    Ok(())
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

#[repr(C)]
pub struct OutputBCA {
    pub len: usize,
    pub particles: *mut [f64; 9],
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
pub extern "C" fn drop_output_tagged_bca(output: OutputTaggedBCA) {
    let length = output.len;

    if length > 0 {

        let particles_layout = Layout::from_size_align(length, align_of::<[f64; 9]>()).unwrap();
        let weights_layout = Layout::from_size_align(length, align_of::<f64>()).unwrap();
        let tags_layout = Layout::from_size_align(length, align_of::<i32>()).unwrap();
        let incident_layout = Layout::from_size_align(length, align_of::<bool>()).unwrap();

        unsafe {
            dealloc(output.particles as *mut u8, particles_layout);
            dealloc(output.weights as *mut u8, weights_layout);
            dealloc(output.tags as *mut u8, tags_layout);
            dealloc(output.incident as *mut u8, incident_layout);
        };
    }
}

#[no_mangle]
pub extern "C" fn drop_output_bca(output: OutputBCA) {
    let length = output.len;

    if length > 0 {
        let particles_layout = Layout::from_size_align(length, align_of::<[f64; 9]>()).unwrap();

        unsafe {
            dealloc(output.particles as *mut u8, particles_layout);
        };
    }
}

#[no_mangle]
pub extern "C" fn compound_tagged_bca_list_c(input: InputTaggedBCA) -> OutputTaggedBCA {

    let mut total_output = vec![];
    let mut output_tags = vec![];
    let mut output_weights = vec![];
    let mut output_incident = vec![];

    let options = Options::default_options(true);

    let Z2 = unsafe { slice::from_raw_parts(input.Z2, input.num_species_target).to_vec() };
    let m2 = unsafe { slice::from_raw_parts(input.m2, input.num_species_target).to_vec() };
    let n2 = unsafe { slice::from_raw_parts(input.n2, input.num_species_target).to_vec() };
    let Ec2 = unsafe { slice::from_raw_parts(input.Ec2, input.num_species_target).to_vec() };
    let Es2 = unsafe { slice::from_raw_parts(input.Es2, input.num_species_target).to_vec() };
    let Eb2 = unsafe { slice::from_raw_parts(input.Eb2, input.num_species_target).to_vec() };
    let positions = unsafe { slice::from_raw_parts(input.positions, input.len).to_vec() };
    let tags = unsafe { slice::from_raw_parts(input.tags, input.len).to_vec() };
    let weights = unsafe { slice::from_raw_parts(input.weights, input.len).to_vec() };

    let x = -2.*(n2.iter().sum::<f64>()*10E30).powf(-1./3.);
    let y = 0.0;
    let z = 0.0;

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: Eb2,
        Es: Es2,
        Ec: Ec2,
        Ed: vec![0.0; input.num_species_target],
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
            Ed: 0.0,
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
    std::mem::forget(output_tags);
    std::mem::forget(output_weights);
    std::mem::forget(output_incident);

    OutputTaggedBCA {
        len,
        particles,
        tags: tags_ptr,
        weights: weights_ptr,
        incident: incident_ptr,
    }
}

#[no_mangle]
pub extern "C" fn reflect_single_ion_c(num_species_target: &mut c_int, ux: &mut f64, uy: &mut f64, uz: &mut f64, E1: &mut f64, Z1: &mut f64, m1: &mut f64, Ec1: &mut f64, Es1: &mut f64, Z2: *mut f64, m2: *mut f64, Ec2: *mut f64, Es2: *mut f64, Eb2: *mut f64, n2: *mut f64) {

    assert!(E1 > &mut 0.0);

    let options = Options::default_options(false);

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
        Ed: vec![0.0; *num_species_target as usize],
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
        Ed: 0.0,
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

    let options = Options::default_options(true);

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![input.Eb2],
        Es: vec![input.Es2],
        Ec: vec![input.Ec2],
        Ed: vec![0.0; input.len],
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
            Ed: 0.0,
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

    let options = Options::default_options(true);

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
        Ed: vec![0.0; input.num_species_target],
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
            Ed: 0.0,
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

    let options = Options::default_options(*track_recoils);

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
        Ed: vec![0.0; *num_species_target as usize],
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
            Ed: 0.0,
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
///compound_tagged_bca_list_py(ux, uy,  uz, energy, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2)
/// runs a BCA simulation for a list of particles and outputs a list of sputtered, reflected, and implanted particles.
/// Args:
///    energies (list(f64)): initial ion energies in eV.
///    ux (list(f64)): initial ion directions x. ux != 0.0 to avoid gimbal lock
///    uy (list(f64)): initial ion directions y.
///    uz (list(f64)): initial ion directions z.
///    Z1 (list(f64)): initial ion atomic numbers.
///    m1 (list(f64)): initial ion masses in amu.
///    Ec1 (list(f64)): ion cutoff energies in eV. If ion energy < Ec1, it stops in the material.
///    Es1 (list(f64)): ion surface binding energies. Assumed planar.
///    Z2 (list(f64)): target material species atomic numbers.
///    m2 (list(f64)): target material species masses in amu.
///    Ec2 (list(f64)): target material species cutoff energies in eV. If recoil energy < Ec2, it stops in the material.
///    Es2 (list(f64)): target species surface binding energies. Assumed planar.
///    n2 (list(f64)): target material species atomic number densities in inverse cubic Angstroms.
///    Eb2 (list(f64)): target material species bulk binding energies in eV.
/// Returns:
///    output (NX9 list of f64): each row in the list represents an output particle (implanted,
///    sputtered, or reflected). Each row consists of:
///      [Z, m (amu), E (eV), x, y, z, (angstrom), ux, uy, uz]
///    incident (list(bool)): whether each row of output was an incident ion or originated in the target
#[pyfunction]
pub fn compound_bca_list_py(energies: Vec<f64>, ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>, Z1: Vec<f64>, m1: Vec<f64>, Ec1: Vec<f64>, Es1: Vec<f64>, Z2: Vec<f64>, m2: Vec<f64>, Ec2: Vec<f64>, Es2: Vec<f64>, n2: Vec<f64>, Eb2: Vec<f64>) -> (Vec<[f64; 9]>, Vec<bool>) {
    let mut total_output = vec![];
    let mut incident = vec![];
    let num_species_target = Z2.len();
    let num_incident_ions = energies.len();

    assert_eq!(ux.len(), num_incident_ions, "Input error: list of x-directions is not the same length as list of incident energies.");
    assert_eq!(uy.len(), num_incident_ions, "Input error: list of y-directions is not the same length as list of incident energies.");
    assert_eq!(uz.len(), num_incident_ions, "Input error: list of z-directions is not the same length as list of incident energies.");
    assert_eq!(Z1.len(), num_incident_ions, "Input error: list of incident atomic numbers is not the same length as list of incident energies.");
    assert_eq!(m1.len(), num_incident_ions, "Input error: list of incident atomic masses is not the same length as list of incident energies.");
    assert_eq!(Es1.len(), num_incident_ions, "Input error: list of incident surface binding energies is not the same length as list of incident energies.");
    assert_eq!(Ec1.len(), num_incident_ions, "Input error: list of incident cutoff energies is not the same length as list of incident energies.");

    assert_eq!(m2.len(), num_species_target, "Input error: list of target atomic masses is not the same length as atomic numbers.");
    assert_eq!(Ec2.len(), num_species_target, "Input error: list of target cutoff energies is not the same length as atomic numbers.");
    assert_eq!(Es2.len(), num_species_target, "Input error: list of target surface binding energies is not the same length as atomic numbers.");
    assert_eq!(Eb2.len(), num_species_target, "Input error: list of target bulk binding energies is not the same length as atomic numbers.");
    assert_eq!(n2.len(), num_species_target, "Input error: list of target number densities is not the same length as atomic numbers.");

    let options = Options::default_options(true);

    let x = -2.*(n2.iter().sum::<f64>()*10E30).powf(-1./3.);
    let y = 0.0;
    let z = 0.0;

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: Eb2,
        Es: Es2,
        Ec: Ec2,
        Ed: vec![0.0; num_species_target],
        Z: Z2,
        m: m2,
        interaction_index: vec![0; num_species_target],
        surface_binding_model: SurfaceBindingModel::INDIVIDUAL,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: n2,
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let mut index: usize = 0;
    for (energy, ux_, uy_, uz_, Z1_, Ec1_, Es1_, m1_) in izip!(energies, ux, uy, uz, Z1, Ec1, Es1, m1) {

        let mut energy_out;

        let p = particle::Particle::default_incident(
            m1_,
            Z1_,
            energy,
            Ec1_,
            Es1_,
            x,
            ux_,
            uy_,
            uz_
        );

        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {
            if (particle.left) | (particle.incident) {

                incident.push(particle.incident);

                if particle.stopped {
                    energy_out = 0.
                } else {
                    energy_out = particle.E/EV
                }
                total_output.push(
                    [
                        particle.Z,
                        particle.m/AMU,
                        energy_out,
                        particle.pos.x/ANGSTROM,
                        particle.pos.y/ANGSTROM,
                        particle.pos.z/ANGSTROM,
                        particle.dir.x,
                        particle.dir.y,
                        particle.dir.z,
                    ]
                );
            }
        }
        index += 1;
    }
    (total_output, incident)
}

#[cfg(feature = "python")]
///reflect_single_ion_py(ion, target, vx, vy, vz)
///Performs a single BCA ion trajectory in target material with specified incident velocity.
///Args:
///    ion (dict): dictionary that defines ion parameters; examples can be found in scripts/materials.py.
///    target (dict): dictionary that defines target parameterrs; examples can be found in scripts/materials.py.
///    vx, vy, vz (float): initial x, y, and z velocity in m/s.
///Returns:
///    vx, vy, vz (float): final x, y, and z velocity in m/s. When ion implants in material, vx, vy, and vz will all be zero.
#[pyfunction]
pub fn reflect_single_ion_py(ion: &PyDict, target: &PyDict, vx: f64, vy: f64, vz: f64) -> (f64, f64, f64){

    let Z1 = unpack(ion.get_item("Z").expect("Cannot get ion Z from dictionary. Ensure ion['Z'] exists."));
    let m1 = unpack(ion.get_item("m").expect("Cannot get ion mass from dictionary. Ensure ion['m'] exists."));
    let Ec1 = unpack(ion.get_item("Ec").expect("Cannot get ion cutoff energy from dictionary. Ensure ion['Ec'] exists."));
    let Es1 = unpack(ion.get_item("Es").expect("Cannot get ion surface binding energy from dictionary. Ensure ion['Es'] exists."));

    let Z2 = unpack(target.get_item("Z").expect("Cannot get target Z from dictionary. Ensure target['Z'] exists."));
    let m2 = unpack(target.get_item("m").expect("Cannot get target mass from dictionary. Ensure target['m'] exists."));
    let Ec2 = unpack(target.get_item("Ec").expect("Cannot get target cutoff energy from dictionary. Ensure target['Ec'] exists."));
    let Es2 = unpack(target.get_item("Es").expect("Cannot get target surface binding energy from dictionary. Ensure target['Es'] exists."));
    let Eb2 = unpack(target.get_item("Eb").expect("Cannot get target bulk binding energy from dictionary. Ensure target['Eb'] exists."));
    let n2 = unpack(target.get_item("n").expect("Cannot get target density from dictionary. Ensure target['n'] exists."));

    assert!(vx > 0.0, "Input error: vx must be greater than zero for incident particles to hit surface at x=0.");

    let options = Options::default_options(false);

    let velocity2 = vx.powf(2.) + vy.powf(2.) + vz.powf(2.); //m/s
    let energy_eV = 0.5*m1*AMU*velocity2/EV; //EV

    let ux = vx/velocity2.sqrt();
    let uy = vy/velocity2.sqrt();
    let uz = vz/velocity2.sqrt();

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![Eb2],
        Es: vec![Es2],
        Ec: vec![Ec2],
        Ed: vec![0.0],
        Z: vec![Z2],
        m: vec![m2],
        interaction_index: vec![0],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "M".to_string(),
        densities: vec![n2],
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let x = -m.geometry.energy_barrier_thickness;
    let y = 0.0;
    let z = 0.0;

    let p = particle::Particle::default_incident(
        m1,
        Z1,
        energy_eV,
        Ec1,
        Es1,
        x,
        ux,
        uy,
        uz
    );

    let output = bca::single_ion_bca(p, &m, &options);

    let reflected_energy = output[0].E; //Joules

    let reflected_velocity = (2.*reflected_energy/(m1*AMU)).sqrt(); //m/s

    let vx2 = output[0].dir.x*reflected_velocity;
    let vy2 = output[0].dir.y*reflected_velocity;
    let vz2 = output[0].dir.z*reflected_velocity;

    if output[0].E > 0.0 && output[0].dir.x < 0.0 && output[0].left && output[0].incident {
        (vx2, vy2, vz2)
    } else {
        (0.0, 0.0, 0.0)
    }
}

#[cfg(feature = "python")]
///compound_bca_list_1D_py(ux, uy,  uz, energies, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, Eb2 n2, dx)
/// runs a BCA simulation for a list of particles and outputs a list of sputtered, reflected, and implanted particles.
/// Args:
///    ux (list(f64)): initial ion directions x. ux != 0.0 to avoid gimbal lock
///    uy (list(f64)): initial ion directions y.
///    uz (list(f64)): initial ion directions z.
///    energies (list(f64)): initial ion energies in eV.
///    Z1 (list(f64)): initial ion atomic numbers.
///    m1 (list(f64)): initial ion masses in amu.
///    Ec1 (list(f64)): ion cutoff energies in eV. If ion energy < Ec1, it stops in the material.
///    Es1 (list(f64)): ion surface binding energies. Assumed planar.
///    Z2 (list(f64)): target material species atomic numbers.
///    m2 (list(f64)): target material species masses in amu.
///    Ec2 (list(f64)): target material species cutoff energies in eV. If recoil energy < Ec2, it stops in the material.
///    Es2 (list(f64)): target species surface binding energies. Assumed planar.
///    Eb2 (list(f64)): target material species bulk binding energies in eV.
///    n2 (list(list(f64))): target material species atomic number densities in inverse cubic Angstroms.
///    dx (list(f64)): target material layer thicknesses starting at surface.
/// Returns:
///    output (NX9 list of f64): each row in the list represents an output particle (implanted,
///    sputtered, or reflected). Each row consists of:
///      [Z, m (amu), E (eV), x, y, z, (angstrom), ux, uy, uz]
///    incident (list(bool)): whether each row of output was an incident ion or originated in the target
/// stopped (list(bool)): whether each row of output is associated with a particle that stopped in the target
#[pyfunction]
pub fn compound_bca_list_1D_py(ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>, energies: Vec<f64>, Z1: Vec<f64>, m1: Vec<f64>, Ec1: Vec<f64>, Es1: Vec<f64>, Z2: Vec<f64>, m2: Vec<f64>, Ec2: Vec<f64>, Es2: Vec<f64>, Eb2: Vec<f64>, n2: Vec<Vec<f64>>,  dx: Vec<f64>) -> (Vec<[f64; 9]>, Vec<bool>, Vec<bool>) {
    let mut total_output = vec![];
    let mut incident = vec![];
    let mut stopped = vec![];
    let num_layers_target = n2.len();
    let num_species = Z2.len();
    let num_incident_ions = energies.len();

    assert_eq!(ux.len(), num_incident_ions, "Input error: list of x-directions is not the same length as list of incident energies.");
    assert_eq!(uy.len(), num_incident_ions, "Input error: list of y-directions is not the same length as list of incident energies.");
    assert_eq!(uz.len(), num_incident_ions, "Input error: list of z-directions is not the same length as list of incident energies.");
    assert_eq!(Z1.len(), num_incident_ions, "Input error: list of incident atomic numbers is not the same length as list of incident energies.");
    assert_eq!(m1.len(), num_incident_ions, "Input error: list of incident atomic masses is not the same length as list of incident energies.");
    assert_eq!(Es1.len(), num_incident_ions, "Input error: list of incident surface binding energies is not the same length as list of incident energies.");
    assert_eq!(Ec1.len(), num_incident_ions, "Input error: list of incident cutoff energies is not the same length as list of incident energies.");

    assert_eq!(m2.len(), num_species, "Input error: list of target atomic masses is not the same length as atomic numbers.");
    assert_eq!(Ec2.len(), num_species, "Input error: list of target cutoff energies is not the same length as atomic numbers.");
    assert_eq!(Es2.len(), num_species, "Input error: list of target surface binding energies is not the same length as atomic numbers.");
    assert_eq!(Eb2.len(), num_species, "Input error: list of target bulk binding energies is not the same length as atomic numbers.");
    assert_eq!(n2[0].len(), num_species, "Input error: first layer list of target number densities is not the same length as atomic numbers.");

    assert_eq!(n2[0].len(), num_species, "Input error: first layer species list of target number densities is not the same length as atomic numbers.");
    assert_eq!(dx.len(), num_layers_target, "Input error: number of layer thicknesses not the same as number of layers in atomic densities list.");

    let options = Options::default_options(true);
    let y = 0.0;
    let z = 0.0;

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: Eb2,
        Es: Es2,
        Ec: Ec2,
        Ed: vec![0.0; num_species],
        Z: Z2,
        m: m2,
        interaction_index: vec![0; num_species],
        surface_binding_model: SurfaceBindingModel::INDIVIDUAL,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh1DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: n2,
        layer_thicknesses: dx,
        electronic_stopping_correction_factors: vec![1.0; num_layers_target]
    };

    let m = material::Material::<Mesh1D>::new(&material_parameters, &geometry_input);

    let x = -m.geometry.top_energy_barrier_thickness/2.;

    let mut index: usize = 0;

    for (energy, ux_, uy_, uz_, Z1_, Ec1_, Es1_, m1_) in izip!(energies, ux, uy, uz, Z1, Ec1, Es1, m1) {();

        let mut energy_out;

        let p = particle::Particle::default_incident(
            m1_,
            Z1_,
            energy,
            Ec1_,
            Es1_,
            x,
            ux_,
            uy_,
            uz_
        );

        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {
            if (particle.left) | (particle.incident) {

                incident.push(particle.incident);
                stopped.push(particle.stopped);

                if particle.stopped {
                    energy_out = 0.
                } else {
                    energy_out = particle.E/EV
                }
                total_output.push(
                    [
                        particle.Z,
                        particle.m/AMU,
                        energy_out,
                        particle.pos.x/ANGSTROM,
                        particle.pos.y/ANGSTROM,
                        particle.pos.z/ANGSTROM,
                        particle.dir.x,
                        particle.dir.y,
                        particle.dir.z,
                    ]
                );
            }
        }
        index += 1;
    }
    (total_output, incident, stopped)
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

    let options = Options::default_options(true);

    let p = particle::Particle {
        m: m1*AMU,
        Z: Z1,
        E: E1*EV,
        Ec: Ec1*EV,
        Es: Es1*EV,
        Ed: 0.0,
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
        Ed: vec![0.0],
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

pub fn simple_compound_bca(x: f64, y: f64, z: f64, ux: f64, uy: f64, uz: f64, E1: f64, Z1: f64, m1: f64, Ec1: f64, Es1: f64, Z2: Vec<f64>, m2: Vec<f64>, Ec2: Vec<f64>, Es2: Vec<f64>, n2: Vec<f64>, Eb2: Vec<f64>) -> Vec<[f64; 9]> {

    assert!(E1 > 0.0, "Error: Incident energy cannot be less than or equal to 0.");
    assert!(Ec1 > 0.0, "Error: Cutoff energy Ec1 cannot be less than or equal to 0.");
    //assert!(Ec2 > 0.0, "Error: Cutoff energy Ec2 cannot be less than or equal to 0.");

    let options = Options::default_options(true);

    let p = particle::Particle {
        m: m1*AMU,
        Z: Z1,
        E: E1*EV,
        Ec: Ec1*EV,
        Es: Es1*EV,
        Ed: 0.0,
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
        Eb: Eb2,
        Es: Es2,
        Ec: Ec2,
        Ed: vec![0.0; Z2.len()],
        Z: Z2,
        m: m2,
        interaction_index: vec![0],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: n2,
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

#[cfg(feature = "parry3d")]
#[no_mangle]
pub extern "C" fn rotate_given_surface_normal(nx: f64, ny: f64, nz: f64, ux: &mut f64, uy: &mut f64, uz: &mut f64) {
    const DELTA: f64 = 1e-9;
    let RUSTBCA_DIRECTION: Vector3::<f64> = Vector3::<f64>::new(1.0, 0.0, 0.0);

    let into_surface = Vector3::new(-nx, -ny, -nz);
    let direction = Vector3::new(*ux, *uy, *uz);

    //Rotation to local RustBCA coordinates from global
    //Here's how this works: a rotation matrix is found that maps the rustbca
    //into-the-surface vector (1.0, 0.0, 0.0) onto the local into-the-surface vector (negative normal w.r.t. ray origin).
    //That rotation is then applied to the particle direction, and can be undone later.
    //Algorithm is from here:
    //https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/180436#180436
    let v: Vector3<f64> = into_surface.cross(&RUSTBCA_DIRECTION);
    let c = into_surface.dot(&RUSTBCA_DIRECTION);
    let vx = Matrix3::<f64>::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0);

    let rotation_matrix = if (c + 1.0).abs() > 1e-6 {
        Matrix3::identity() + vx + vx*vx/(1. + c)
    } else {
        //If c == -1.0, the correct rotation should simply be a 180 degree rotation
        //around a non-x axis; y is chosen arbitrarily
        Rotation3::from_axis_angle(&Vector3::y_axis(), PI).into()
    };

    let incident = rotation_matrix*direction;

    // ux must not be exactly 1.0 to avoid gimbal lock in RustBCA
    // simple_bca does not normalize direction before proceeding, must be done manually
    assert!(
        incident.x + DELTA > 0.0, "Error: RustBCA initial direction out of surface. Please check surface normals and incident direction. c={} n = ({}, {}, {}) u = ({}, {}, {}), unew = ({}, {}, {})",
        (c + 1.0).abs(), nx, ny, nz, ux, uy, uz, incident.x, incident.y, incident.z
    );

    *ux = incident.x + DELTA;
    *uy = incident.y - DELTA;
    *uz = incident.z;
    let mag = (ux.powf(2.) + uy.powf(2.) + uz.powf(2.)).sqrt();

    *ux /= mag;
    *uy /= mag;
    *uz /= mag;
}



#[cfg(all(feature = "python", feature = "parry3d"))]
#[pyfunction]
/// rotate_given_surface_normal_py(nx, ny, nz, ux, uy, uz)
/// --
///
/// This function takes a particle direction and a normal vector and rotates from simulation to RustBCA coordinates.
/// Args:
///     nx (f64): surface normal in global frame x-component.
///     ny (f64): surface normal in global frame y-component.
///     nz (f64): surface normal in global frame z-component.
///     ux (f64): particle direction in global frame x-component.
///     uy (f64): particle direction in global frame normal y-component.
///     uz (f64): particle direction in global frame normal z-component.
/// Returns:
///    direction (f64, f64, f64): direction vector of particle in RustBCA coordinates.
pub fn rotate_given_surface_normal_py(nx: f64, ny: f64, nz: f64, ux: f64, uy: f64, uz: f64) -> (f64, f64, f64) {
    let mut ux = ux;
    let mut uy = uy;
    let mut uz = uz;
    rotate_given_surface_normal(nx, ny, nz, &mut ux, &mut uy, &mut uz);
    (ux, uy, uz)
}

#[cfg(all(feature = "python", feature = "parry3d"))]
#[pyfunction]
/// rotate_given_surface_normal_vec_py(nx, ny, nz, ux, uy, uz)
/// --
///
/// This function takes a particle direction and a normal vector and rotates from simulation to RustBCA coordinates.
/// Args:
///     nx (list(f64)): surface normal in global frame x-component.
///     ny (list(f64)): surface normal in global frame y-component.
///     nz (list(f64)): surface normal in global frame z-component.
///     ux (list(f64)): particle direction in global frame x-component.
///     uy (list(f64)): particle direction in global frame normal y-component.
///     uz (list(f64)): particle direction in global frame normal z-component.
/// Returns:
///    direction (list(f64), list(f64), list(f64)): direction vector of particle in RustBCA coordinates.
///    Note: non-incident particles will be returned with ux, uy, uz = (0, 0, 0)
pub fn rotate_given_surface_normal_vec_py(nx: Vec<f64>, ny: Vec<f64>, nz: Vec<f64>, ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {

    let length = nx.len();

    let mut ux_new = Vec::with_capacity(length);
    let mut uy_new = Vec::with_capacity(length);
    let mut uz_new = Vec::with_capacity(length);

    (0..length).into_iter().for_each(|index| {

        let mut ux_ = ux[index];
        let mut uy_ = uy[index];
        let mut uz_ = uz[index];

        if (nx[index]*ux_ + ny[index]*uy_ + nz[index]*uz_) <= 0. {
            rotate_given_surface_normal(nx[index], ny[index], nz[index], &mut ux_, &mut uy_, &mut uz_);
            ux_new.push(ux_);
            uy_new.push(uy_);
            uz_new.push(uz_);
        } else {
            ux_new.push(0.);
            uy_new.push(0.);
            uz_new.push(0.);
        }

    });

    (ux_new, uy_new, uz_new)
}

#[cfg(feature = "parry3d")]
#[no_mangle]
pub extern "C" fn rotate_back(nx: f64, ny: f64, nz: f64, ux: &mut f64, uy: &mut f64, uz: &mut f64) {
    let RUSTBCA_DIRECTION: Vector3::<f64> = Vector3::<f64>::new(1.0, 0.0, 0.0);

    let into_surface = Vector3::new(-nx, -ny, -nz);
    let direction = Vector3::new(*ux, *uy, *uz);

    //Rotation to local RustBCA coordinates from global
    //Here's how this works: a rotation matrix is found that maps the rustbca
    //into-the-surface vector (1.0, 0.0, 0.0) onto the local into-the-surface vector (negative normal w.r.t. ray origin).
    //That rotation is then applied to the particle direction, and can be undone later.
    //Algorithm is from here:
    //https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/180436#180436
    let v: Vector3<f64> = into_surface.cross(&RUSTBCA_DIRECTION);
    let c = into_surface.dot(&RUSTBCA_DIRECTION);
    let vx = Matrix3::<f64>::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0);
    let rotation_matrix = if c != -1.0 {
        Matrix3::identity() + vx + vx*vx/(1. + c)
    } else {
        //If c == -1.0, the correct rotation should simply be a 180 degree rotation
        //around a non-x axis; y is chosen arbitrarily
        Rotation3::from_axis_angle(&Vector3::y_axis(), PI).into()
    };

    let u = rotation_matrix.transpose()*direction;

    *ux = u.x;
    *uy = u.y;
    *uz = u.z;
}

#[cfg(all(feature = "python", feature = "parry3d"))]
#[pyfunction]
/// rotate_back_py(nx, ny, nz, ux, uy, uz)
/// --
///
/// This function takes a particle direction and a normal vector and rotates from RustBCA to simulation coordinates.
/// Args:
///     nx (f64): surface normal in global frame x-component.
///     ny (f64): surface normal in global frame y-component.
///     nz (f64): surface normal in global frame z-component.
///     ux (f64): particle direction in RustBCA frame x-component.
///     uy (f64): particle direction in RustBCA frame normal y-component.
///     uz (f64): particle direction in RustBCA frame normal z-component.
/// Returns:
///    direction (f64, f64, f64): direction vector of particle in global coordinates.
pub fn rotate_back_py(nx: f64, ny: f64, nz: f64, ux: f64, uy: f64, uz: f64) -> (f64, f64, f64) {
    let mut ux = ux;
    let mut uy = uy;
    let mut uz = uz;
    rotate_back(nx, ny, nz, &mut ux, &mut uy, &mut uz);
    (ux, uy, uz)
}

#[cfg(all(feature = "python", feature = "parry3d"))]
#[pyfunction]
/// rotate_back_vec_py(nx, ny, nz, ux, uy, uz)
/// --
///
/// This function takes a RustBCA particle direction and a normal vector and rotates back from RustBCA to simulation coordinates.
/// Args:
///     nx (list(f64)): surface normal in global frame x-component.
///     ny (list(f64)): surface normal in global frame y-component.
///     nz (list(f64)): surface normal in global frame z-component.
///     ux (list(f64)): particle direction in global frame x-component.
///     uy (list(f64)): particle direction in global frame normal y-component.
///     uz (list(f64)): particle direction in global frame normal z-component.
/// Returns:
///    direction (list(f64), list(f64), list(f64)): direction vector of particle in simulation coordinates.
pub fn rotate_back_vec_py(nx: Vec<f64>, ny: Vec<f64>, nz: Vec<f64>, ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {

    let length = nx.len();

    let mut ux_new = Vec::with_capacity(length);
    let mut uy_new = Vec::with_capacity(length);
    let mut uz_new = Vec::with_capacity(length);

    (0..length).into_iter().for_each(|index| {

        let mut ux_ = ux[index];
        let mut uy_ = uy[index];
        let mut uz_ = uz[index];
        rotate_back(nx[index], ny[index], nz[index], &mut ux_, &mut uy_, &mut uz_);
        ux_new.push(ux_);
        uy_new.push(uy_);
        uz_new.push(uz_);
    });

    (ux_new, uy_new, uz_new)
}

#[cfg(feature = "python")]
/// A helper function to unpack a python float from a python any.
fn unpack(python_float: &PyAny) -> f64 {
    python_float.downcast::<PyFloat>().expect("Error unpacking Python float to f64. Check values.").value()
}

#[cfg(feature = "python")]
#[pyfunction]
/// sputteirng_yield(ion, target, energy, angle, num_samples)
/// A routine the calculates the sputtering yield in atoms per ion of energetic ions incident upon materials using RustBCA.
/// Args:
///     ion: a dictionary with the keys Z (atomic number), m (atomic mass in AMU), Ec (cutoff energy in eV), Es (surface binding energy in eV)
///     target: a dictionary with the keys Z, m, Ec, Es, Eb (bulk binding energy in eV), n2 (number density in 1/m3)
///     energy: the incident energy of the ion in eV
///     angle: incident angle of the ion in degrees from surface normal
///     num_samples: number of ion trajectories to run; precision will go as 1/sqrt(N)
pub fn sputtering_yield(ion: &PyDict, target: &PyDict, energy: f64, angle: f64, num_samples: usize) -> f64 {

    assert!(angle.abs() <= 90.0, "Incident angle w.r.t. surface normal, {}, cannot exceed 90 degrees.", angle);

    const DELTA: f64 = 1e-6;

    let Z1 = unpack(ion.get_item("Z").expect("Cannot get ion Z from dictionary. Ensure ion['Z'] exists."));
    let m1 = unpack(ion.get_item("m").expect("Cannot get ion mass from dictionary. Ensure ion['m'] exists."));
    let Ec1 = unpack(ion.get_item("Ec").expect("Cannot get ion cutoff energy from dictionary. Ensure ion['Ec'] exists."));
    let Es1 = unpack(ion.get_item("Es").expect("Cannot get ion surface binding energy from dictionary. Ensure ion['Es'] exists."));

    let Z2 = unpack(target.get_item("Z").expect("Cannot get target Z from dictionary. Ensure target['Z'] exists."));
    let m2 = unpack(target.get_item("m").expect("Cannot get target mass from dictionary. Ensure target['m'] exists."));
    let Ec2 = unpack(target.get_item("Ec").expect("Cannot get target cutoff energy from dictionary. Ensure target['Ec'] exists."));
    let Es2 = unpack(target.get_item("Es").expect("Cannot get target surface binding energy from dictionary. Ensure target['Es'] exists."));
    let Eb2 = unpack(target.get_item("Eb").expect("Cannot get target bulk binding energy from dictionary. Ensure target['Eb'] exists."));
    let n2 = unpack(target.get_item("n").expect("Cannot get target density from dictionary. Ensure target['n'] exists."));

    let options = Options::default_options(true);

    let y = 0.0;
    let z = 0.0;

    let ux = (angle/180.0*PI).cos() + DELTA;
    let uy = (angle/180.0*PI).sin() - DELTA;
    let uz = 0.0;

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![Eb2],
        Es: vec![Es2],
        Ec: vec![Ec2],
        Ed: vec![0.0],
        Z: vec![Z2],
        m: vec![m2],
        interaction_index: vec![0],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "M".to_string(),
        densities: vec![n2],
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let x = -m.geometry.energy_barrier_thickness;

    let num_sputtered = Mutex::new(0);

    (0..num_samples as u64).into_par_iter().for_each( |index| {

        let p = particle::Particle::default_incident(
            m1,
            Z1,
            energy,
            Ec1,
            Es1,
            x,
            ux,
            uy,
            uz
        );

        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {
            if particle.E > 0.0 && particle.dir.x < 0.0 && particle.left && (!particle.incident) {
                let mut num_sputtered = num_sputtered.lock().unwrap();
                *num_sputtered += 1;
            }
        }
    });
    let num_sputtered = *num_sputtered.lock().unwrap();
    num_sputtered as f64 / num_samples as f64
}

#[cfg(feature = "python")]
#[pyfunction]
/// reflection_coefficient(ion, target, energy, angle, num_samples)
/// A routine the calculates the reflection coefficient of energetic ions incident upon materials using RustBCA.
/// Args:
///     ion: a dictionary with the keys Z (atomic number), m (atomic mass in AMU), Ec (cutoff energy in eV), Es (surface binding energy in eV)
///     target: a dictionary with the keys Z, m, Ec, Es, Eb (bulk binding energy in eV), n2 (number density in 1/m3)
///     energy: the incident energy of the ion in eV
///     angle: incident angle of the ion in degrees from surface normal
///     num_samples: number of ion trajectories to run; precision will go as 1/sqrt(N)
/// Returns:
///     R_N (f64): reflection coefficient (number of particles reflected / number of incident particles)
///     R_E (f64): energy reflection coefficient (sum of reflected particle energies / total incident energy)
pub fn reflection_coefficient(ion: &PyDict, target: &PyDict, energy: f64, angle: f64, num_samples: usize) -> (f64, f64) {

    const DELTA: f64 = 1e-6;

    assert!(angle.abs() <= 90.0, "Incident angle w.r.t. surface normal, {}, cannot exceed 90 degrees.", angle);

    let Z1 = unpack(ion.get_item("Z").expect("Cannot get ion Z from dictionary. Ensure ion['Z'] exists."));
    let m1 = unpack(ion.get_item("m").expect("Cannot get ion mass from dictionary. Ensure ion['m'] exists."));
    let Ec1 = unpack(ion.get_item("Ec").expect("Cannot get ion cutoff energy from dictionary. Ensure ion['Ec'] exists."));
    let Es1 = unpack(ion.get_item("Es").expect("Cannot get ion surface binding energy from dictionary. Ensure ion['Es'] exists."));

    let Z2 = unpack(target.get_item("Z").expect("Cannot get target Z from dictionary. Ensure target['Z'] exists."));
    let m2 = unpack(target.get_item("m").expect("Cannot get target mass from dictionary. Ensure target['m'] exists."));
    let Ec2 = unpack(target.get_item("Ec").expect("Cannot get target cutoff energy from dictionary. Ensure target['Ec'] exists."));
    let Es2 = unpack(target.get_item("Es").expect("Cannot get target surface binding energy from dictionary. Ensure target['Es'] exists."));
    let Eb2 = unpack(target.get_item("Eb").expect("Cannot get target bulk binding energy from dictionary. Ensure target['Eb'] exists."));
    let n2 = unpack(target.get_item("n").expect("Cannot get target density from dictionary. Ensure target['n'] exists."));

    let options = Options::default_options(false);

    let y = 0.0;
    let z = 0.0;

    let ux = (angle/180.0*PI).cos() + DELTA;
    let uy = (angle/180.0*PI).sin() - DELTA;
    let uz = 0.0;

    let mut direction = Vector::new(ux, uy, uz);
    direction.normalize();

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![Eb2],
        Es: vec![Es2],
        Ec: vec![Ec2],
        Ed: vec![0.0],
        Z: vec![Z2],
        m: vec![m2],
        interaction_index: vec![0],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "M".to_string(),
        densities: vec![n2],
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let x = -m.geometry.energy_barrier_thickness;

    let num_reflected = Mutex::new(0);
    let energy_reflected = Mutex::new(0.0);

    (0..num_samples as u64).into_par_iter().for_each( |index| {

        let p = particle::Particle::default_incident(
            m1,
            Z1,
            energy,
            Ec1,
            Es1,
            x,
            ux,
            uy,
            uz
        );

        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {
            if particle.E > 0.0 && particle.dir.x < 0.0 && particle.left && particle.incident {
                let mut num_reflected = num_reflected.lock().unwrap();
                *num_reflected += 1;
                let mut energy_reflected = energy_reflected.lock().unwrap();
                *energy_reflected += particle.E;
            }
        }
    });
    let num_reflected = *num_reflected.lock().unwrap();
    let energy_reflected = *energy_reflected.lock().unwrap();

    (num_reflected as f64 / num_samples as f64, energy_reflected / EV / energy / num_samples as f64)
}

#[cfg(feature = "python")]
#[pyfunction]
/// compound_reflection_coefficient(ion, target_species, target_number_densities, energy, angle, num_samples)
/// A routine the calculates the reflection coefficient of energetic ions incident upon materials using RustBCA.
/// Args:
///     ion: a dictionary with the keys Z (atomic number), m (atomic mass in AMU), Ec (cutoff energy in eV), Es (surface binding energy in eV)
///     target_species: a list of dictionaries with the keys Z, m, Ec, Es, Eb (bulk binding energy in eV), n2 (number density in 1/m3)
///     target_number_densities (list(f64)): number density of each target species in the compound
///     energy: the incident energy of the ion in eV
///     angle: incident angle of the ion in degrees from surface normal
///     num_samples: number of ion trajectories to run; precision will go as 1/sqrt(N)
/// Returns:
///     R_N (f64): reflection coefficient (number of particles reflected / number of incident particles)
///     R_E (f64): energy reflection coefficient (sum of reflected particle energies / total incident energy)
pub fn compound_reflection_coefficient(ion: &PyDict, targets: Vec<&PyDict>, target_number_densities: Vec<f64>, energy: f64, angle: f64, num_samples: usize) -> (f64, f64) {

    const DELTA: f64 = 1e-6;

    assert!(angle.abs() <= 90.0, "Incident angle w.r.t. surface normal, {}, cannot exceed 90 degrees.", angle);

    let Z1 = unpack(ion.get_item("Z").expect("Cannot get ion Z from dictionary. Ensure ion['Z'] exists."));
    let m1 = unpack(ion.get_item("m").expect("Cannot get ion mass from dictionary. Ensure ion['m'] exists."));
    let Ec1 = unpack(ion.get_item("Ec").expect("Cannot get ion cutoff energy from dictionary. Ensure ion['Ec'] exists."));
    let Es1 = unpack(ion.get_item("Es").expect("Cannot get ion surface binding energy from dictionary. Ensure ion['Es'] exists."));

    let Z2: Vec<f64> = targets.iter().enumerate().map( |(index, item)| unpack(item.get_item("Z").unwrap_or_else(|| panic!("Cannot get target Z from dictionary at index {}. Ensure target['Z'] exists.", index)))).collect();
    let m2: Vec<f64> = targets.iter().enumerate().map( |(index, item)| unpack(item.get_item("m").unwrap_or_else(|| panic!("Cannot get target m from dictionary at index {}. Ensure target['m'] exists.", index)))).collect();
    let Ec2: Vec<f64> = targets.iter().enumerate().map( |(index, item)| unpack(item.get_item("Ec").unwrap_or_else(|| panic!("Cannot get target Ec from dictionary at index {}. Ensure target['Ec'] exists.", index)))).collect();
    let Es2: Vec<f64> = targets.iter().enumerate().map( |(index, item)| unpack(item.get_item("Es").unwrap_or_else(|| panic!("Cannot get target Es from dictionary at index {}. Ensure target['Es'] exists.", index)))).collect();
    let Eb2: Vec<f64> = targets.iter().enumerate().map( |(index, item)| unpack(item.get_item("Eb").unwrap_or_else(|| panic!("Cannot get target Eb from dictionary at index {}. Ensure target['Eb'] exists.", index)))).collect();

    let number_target_species = Z2.len();

    let options = Options::default_options(false);

    let y = 0.0;
    let z = 0.0;

    let ux = (angle/180.0*PI).cos() + DELTA;
    let uy = (angle/180.0*PI).sin() - DELTA;
    let uz = 0.0;

    let mut direction = Vector::new(ux, uy, uz);
    direction.normalize();

    let material_parameters = material::MaterialParameters {
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: Eb2,
        Es: Es2,
        Ec: Ec2,
        Ed: vec![0.0; number_target_species],
        Z: Z2,
        m: m2,
        interaction_index: vec![0; number_target_species],
        surface_binding_model: SurfaceBindingModel::AVERAGE,
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let geometry_input = geometry::Mesh0DInput {
        length_unit: "M".to_string(),
        densities: target_number_densities,
        electronic_stopping_correction_factor: 1.0
    };

    let m = material::Material::<Mesh0D>::new(&material_parameters, &geometry_input);

    let x = -m.geometry.energy_barrier_thickness;

    let num_reflected = Mutex::new(0);
    let energy_reflected = Mutex::new(0.0);

    (0..num_samples as u64).into_par_iter().for_each( |index| {

        let p = particle::Particle::default_incident(
            m1,
            Z1,
            energy,
            Ec1,
            Es1,
            x,
            ux,
            uy,
            uz
        );

        let output = bca::single_ion_bca(p, &m, &options);

        for particle in output {
            if particle.E > 0.0 && particle.dir.x < 0.0 && particle.left && particle.incident {
                let mut num_reflected = num_reflected.lock().unwrap();
                *num_reflected += 1;
                let mut energy_reflected = energy_reflected.lock().unwrap();
                *energy_reflected += particle.E;
            }
        }
    });
    let num_reflected = *num_reflected.lock().unwrap();
    let energy_reflected = *energy_reflected.lock().unwrap();

    (num_reflected as f64 / num_samples as f64, energy_reflected / EV / energy / num_samples as f64)
}
