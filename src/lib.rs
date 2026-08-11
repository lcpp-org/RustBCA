#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

use std::{fmt, env};
use std::mem::discriminant;

use std::alloc::{dealloc, Layout};
use std::mem::align_of;

//Parallelization - currently only used in python library functions
//#[cfg(feature = "python")]
//use rayon::ThreadPoolBuilder;
#[cfg(feature = "python")]
use rayon::iter::{IndexedParallelIterator, ParallelExtend, IntoParallelIterator, ParallelIterator};

//Error handling crate
use anyhow::{Result, Context, anyhow};

//Serializing/Deserializing crate
use serde::*;

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
use itertools::{izip};

//RNG
use rand::{SeedableRng, rngs::ChaCha8Rng};

//Math
use std::f64::consts::FRAC_2_SQRT_PI;
use std::f64::consts::PI;
use std::f64::consts::SQRT_2;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::*;
#[cfg(feature = "python")]
use pythonize::*;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyValueError, PyRuntimeError};

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
pub mod math;
pub mod physics;

pub use crate::enums::*;
pub use crate::consts::*;
pub use crate::structs::*;
pub use crate::input::{Input2D, InputHomogeneous2D, Input1D, Input0D, Options, InputFile, GeometryInput};
pub use crate::output::{OutputUnits};
pub use crate::geometry::{Geometry, GeometryElement, Mesh0D, Mesh1D, Mesh2D};
pub use crate::sphere::{Sphere, SphereInput, InputSphere};
pub use crate::math::*;
pub use crate::material::*;
pub use crate::physics::*;

#[cfg(feature = "parry3d")]
pub mod parry;

#[cfg(feature = "parry3d")]
pub use crate::parry::{ParryBall, ParryBallInput, InputParryBall, ParryTriMesh, ParryTriMeshInput, InputParryTriMesh};
#[cfg(feature = "parry3d")]
pub use parry3d_f64::na::{Point3, Vector3, Matrix3, Rotation3};

#[cfg(feature = "python")]
#[pymodule]
mod libRustBCA {

    #[pymodule_export]
    use super::simple_bca_py;

    #[pymodule_export]
    use super::simple_bca_list_py;

    #[pymodule_export]
    use super::compound_bca_list_py;

    #[pymodule_export]
    use super::compound_bca_list_1D_py;

    #[pymodule_export]
    use super::compound_bca_list_tracked_py;

    #[pymodule_export]
    use super::reflect_single_ion_py;

    #[pymodule_export]  
    use super::reflection_coefficient;

    #[pymodule_export]
    use super::compound_reflection_coefficient;

    #[pymodule_export]
    use super::sputtering_yield;

    #[cfg(feature = "parry3d")]
    #[pymodule_export]
    use super::rotate_given_surface_normal_py;

    #[cfg(feature = "parry3d")]
    #[pymodule_export]
    use super::rotate_back_py;

    #[cfg(feature = "parry3d")]
    #[pymodule_export]
    use super::rotate_back_vec_py;

    #[cfg(feature = "parry3d")]
    #[pymodule_export]
    use super::rotate_given_surface_normal_vec_py;

    #[pymodule_export]
    use super::electronic_stopping_cross_sections;

    #[pymodule_export]
    use super::scattering_integrals;

    #[pymodule_export]
    use super::rustbca_py;

    #[pymodule_export]
    use super::rustbca_local_py;
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

#[unsafe(no_mangle)]
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

#[unsafe(no_mangle)]
pub extern "C" fn drop_output_bca(output: OutputBCA) {
    let length = output.len;

    if length > 0 {
        let particles_layout = Layout::from_size_align(length, align_of::<[f64; 9]>()).unwrap();

        unsafe {
            dealloc(output.particles as *mut u8, particles_layout);
        };
    }
}

#[unsafe(no_mangle)]
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

    let x = -2.*(n2.iter().sum::<f64>()*1E30).powf(-1./3.);
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

        let mut rng = ChaCha8Rng::seed_from_u64(index as u64);
        let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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

#[unsafe(no_mangle)]
pub unsafe extern "C" fn reflect_single_ion_c(num_species_target: &mut c_int, ux: &mut f64, uy: &mut f64, uz: &mut f64, E1: &mut f64, Z1: &mut f64, m1: &mut f64, Ec1: &mut f64, Es1: &mut f64, Z2: *mut f64, m2: *mut f64, Ec2: *mut f64, Es2: *mut f64, Eb2: *mut f64, n2: *mut f64) {

    assert!(E1 > &mut 0.0);

    let options = Options::default_options(false);

    let Z2 = unsafe { slice::from_raw_parts(Z2, *num_species_target as usize).to_vec() };
    let m2 = unsafe { slice::from_raw_parts(m2, *num_species_target as usize).to_vec() };
    let n2 = unsafe { slice::from_raw_parts(n2, *num_species_target as usize).to_vec() };
    let Ec2 = unsafe { slice::from_raw_parts(Ec2, *num_species_target as usize).to_vec() };
    let Es2 = unsafe { slice::from_raw_parts(Es2, *num_species_target as usize).to_vec() };
    let Eb2 = unsafe { slice::from_raw_parts(Eb2, *num_species_target as usize).to_vec() };

    let x = -2.*(n2.iter().sum::<f64>()*1E30).powf(-1./3.);
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

    let mut rng = ChaCha8Rng::from_rng(&mut rand::rng());
    let output = bca::single_ion_bca(p, &m, &options, &mut rng);

    *ux = output[0].dir.x;
    *uy = output[0].dir.y;
    *uz = output[0].dir.z;
    if output[0].pos.x >= 0.0 {
        *E1 = 0.0
    } else {
        *E1 = output[0].E/EV;
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn simple_bca_list_c(input: InputSimpleBCA) -> OutputBCA {

    let x = -2.*(input.n2*1E30).powf(-1./3.);
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

    let seed: u64 = get_seed().unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(seed);
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

        let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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

#[unsafe(no_mangle)]
pub extern "C" fn compound_bca_list_c(input: InputCompoundBCA) -> OutputBCA {

    let mut total_output = vec![];

    let options = Options::default_options(true);

    let Z2 = unsafe { slice::from_raw_parts(input.Z2, input.num_species_target).to_vec() };
    let m2 = unsafe { slice::from_raw_parts(input.m2, input.num_species_target).to_vec() };
    let n2 = unsafe { slice::from_raw_parts(input.n2, input.num_species_target).to_vec() };
    let Ec2 = unsafe { slice::from_raw_parts(input.Ec2, input.num_species_target).to_vec() };
    let Es2 = unsafe { slice::from_raw_parts(input.Es2, input.num_species_target).to_vec() };
    let Eb2 = unsafe { slice::from_raw_parts(input.Eb2, input.num_species_target).to_vec() };

    let x = -2.*(n2.iter().sum::<f64>()*1E30).powf(-1./3.);
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

    let seed: u64 = get_seed().unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(seed);
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

        let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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

#[unsafe(no_mangle)]
pub unsafe extern "C" fn compound_bca_list_fortran(num_incident_ions: &mut c_int, track_recoils: &mut bool,
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

    let x = -2.*(n2.iter().sum::<f64>()*1E30).powf(-1./3.);
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

    let seed: u64 = get_seed().unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(seed);
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

        
        let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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

#[unsafe(no_mangle)]
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

#[cfg(feature="python")]
#[pyfunction]
#[pyo3(signature = (Za, Zb, E, Ma, ck=1.0, ci=1.0))]
///electronic_stopping_cross_sections(Za, Zb, E, Ma, ck)
/// uses RustBCA internal functions to calculate electronic stopping power cross-sections
/// Args:
///     Za (f64): atomic number of ion
///     Zb (f64): atomic number of target
///     E (f64): ion energy in eV
///     Ma (f64): ion mass in AMU
///     ck (f64): LS correction factor
///     ci (f64): BV custom interp. weight
/// Returns:
///     (Lindhard-Scharff [eV m^2], Bethe-Bloch [eV m^2], Biersack-Varelas [eV m^2], Biersack-Varelas with custom interp. weight [eV m^2])
pub fn electronic_stopping_cross_sections<'py>(Za: f64, Zb: f64, E: f64, Ma: f64, ck: f64, ci: f64) -> (f64, f64, f64, f64) {

    let S_low = lindhard_scharff_stopping_power_cross_section(Za, Zb, E*EV, Ma*AMU);
    let S_high = bethe_bloch_stopping_power_cross_section(Za, Zb, E*EV, Ma*AMU);

    (S_low*ck/EV, S_high/EV, 1./(1./(S_high) + 1./(S_low*ck))/EV, (S_high.powf(-ci) + (S_low*ck).powf(-ci)).powf(-1./ci)/EV)
}

#[cfg(feature = "python")]
///compound_bca_list_py(ux, uy,  uz, energy, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2)
/// runs a BCA simulation for a list of particles and outputs a list of sputtered, reflected, and implanted particles.
/// Args:
///    energies (list(f64)): initial ion energies in eV.
///    ux (list(f64)): initial ion directions x.
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
pub fn compound_bca_list_py<'py>(energies: Vec<f64>, ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>, Z1: Vec<f64>, m1: Vec<f64>, Ec1: Vec<f64>, Es1: Vec<f64>, Z2: Vec<f64>, m2: Vec<f64>, Ec2: Vec<f64>, Es2: Vec<f64>, n2: Vec<f64>, Eb2: Vec<f64>) -> (Vec<[f64; 9]>, Vec<bool>) {
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

    let x = -2.*(n2.iter().sum::<f64>()*1E30).powf(-1./3.);
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

    let seed: u64 = get_seed().unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(seed);

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

        
        let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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
    }
    (total_output, incident)
}

#[cfg(feature = "python")]
///compound_bca_list_tracked_py(ux, uy,  uz, energy, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2)
/// runs a BCA simulation for a list of particles and outputs a list of sputtered, reflected, and implanted particles.
/// Args:
///    energies (list(f64)): initial ion energies in eV.
///    ux (list(f64)): initial ion directions x.
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
///    incident_index (list(usize)): index of incident particle that caused this particle to be emitted
#[pyfunction]
pub fn compound_bca_list_tracked_py<'py>(energies: Vec<f64>, ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>, Z1: Vec<f64>, m1: Vec<f64>, Ec1: Vec<f64>, Es1: Vec<f64>, Z2: Vec<f64>, m2: Vec<f64>, Ec2: Vec<f64>, Es2: Vec<f64>, n2: Vec<f64>, Eb2: Vec<f64>) -> (Vec<[f64; 9]>, Vec<bool>, Vec<usize>) {
    let mut total_output = vec![];
    let mut incident = vec![];
    let mut incident_index = vec![];
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
    //options.high_energy_free_flight_paths = true;

    let x = -2.*(n2.iter().sum::<f64>()*1E30).powf(-1./3.);
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

    let mut finished_particles: Vec<particle::Particle> = Vec::new();

    let seed: u64 = get_seed().unwrap();

    let incident_particles: Vec<particle::Particle> = izip!(energies, ux, uy, uz, Z1, Ec1, Es1, m1)
        .enumerate()
        .map(|(index, (energy, ux_, uy_, uz_, Z1_, Ec1_, Es1_, m1_))| {
            let mut p = particle::Particle::default_incident(
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
            p.tag = index as i32;
            p
        }).collect();

        finished_particles.par_extend(
            incident_particles.into_par_iter()
            .enumerate()
            .map_init(
                || ChaCha8Rng::seed_from_u64(seed),
                | rng, (particle_index, incident_particle)| {
                    rng.set_stream(particle_index as u64);
                    bca::single_ion_bca(incident_particle, &m, &options, rng)
                } 
            )
            .flatten()
        );

        for particle in finished_particles {
            if (particle.left) | (particle.incident) {
                incident.push(particle.incident);
                incident_index.push(particle.tag as usize);
                let energy_out;
                if particle.stopped {
                    energy_out = 0.;
                } else {
                    energy_out = particle.E/EV;
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
                )
            }
        }

    (total_output, incident, incident_index)
}

#[cfg(feature = "python")]
///reflect_single_ion_py(ion, target, vx, vy, vz)
///Performs a single BCA ion trajectory in target material with specified incident veloci
///    ion (dict): dictionary that defines ion parameters; examples can be found in scripts/materials.py.
///    target (dict): dictionary that defines target parameterrs; examples can be found in scripts/materials.py.
///    vx, vy, vz (float): initial x, y, and z velocity in m/s.
///Returns:
///    vx, vy, vz (float): final x, y, and z velocity in m/s. When ion implants in material, vx, vy, and vz will all be zero.
#[pyfunction]
pub fn reflect_single_ion_py<'py>(ion: &Bound<'py, PyDict>, target: &Bound<'py, PyDict>, vx: f64, vy: f64, vz: f64) -> (f64, f64, f64){
    
    let Z1: f64 = ion.get_item("Z").unwrap().expect("Error: Cannot get key 'Z' from ion dict.").extract().unwrap();
    let m1: f64 = ion.get_item("m").unwrap().expect("Error: Cannot get key 'm' from ion dict.").extract().unwrap();
    let Es1: f64 = ion.get_item("Es").unwrap().expect("Error: Cannot get key 'Es' from ion dict.").extract().unwrap();
    let Ec1: f64 = ion.get_item("Ec").unwrap().expect("Error: Cannot get key 'Ec' from ion dict.").extract().unwrap();

    let Z2: f64 = target.get_item("Z").unwrap().expect("Error: Cannot get key 'Z' from target dict.").extract().unwrap();
    let m2: f64 = target.get_item("m").unwrap().expect("Error: Cannot get key 'm' from target dict.").extract().unwrap();
    let Es2: f64 = target.get_item("Es").unwrap().expect("Error: Cannot get key 'Es' from target dict.").extract().unwrap();
    let Ec2: f64 = target.get_item("Ec").unwrap().expect("Error: Cannot get key 'Ec' from target dict.").extract().unwrap();
    let Eb2: f64 = target.get_item("Eb").unwrap().expect("Error: Cannot get key 'Eb' from target dict.").extract().unwrap();
    let n2: f64 = target.get_item("n").unwrap().expect("Error: Cannot get key 'n' from target dict.").extract().unwrap();

    assert!(vx > 0.0, "Input error: vx must be greater than zero for incident particles to hit surface at x=0.");

    let options = Options::default_options(false);

    let velocity2 = vx*vx + vy*vy + vz*vz; //m^2/s^2
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

    let mut rng = ChaCha8Rng::from_rng(&mut rand::rng());
    let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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
///    ux (list(f64)): initial ion directions x.
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
pub fn compound_bca_list_1D_py<'py>(ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>, energies: Vec<f64>, Z1: Vec<f64>, m1: Vec<f64>, Ec1: Vec<f64>, Es1: Vec<f64>, Z2: Vec<f64>, m2: Vec<f64>, Ec2: Vec<f64>, Es2: Vec<f64>, Eb2: Vec<f64>, n2: Vec<Vec<f64>>,  dx: Vec<f64>) -> PyResult<(Vec<[f64; 9]>, Vec<bool>, Vec<bool>)> {
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

    let seed: u64 = get_seed().map_err(|error| PyValueError::new_err(""))?;

    let mut rng = ChaCha8Rng::seed_from_u64(seed);
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

        let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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
    }
    Ok((total_output, incident, stopped))
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
///    ux (f64): initial ion direction x.
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
pub fn simple_bca_py<'py>(x: f64, y: f64, z: f64, ux: f64, uy: f64, uz: f64, E1: f64, Z1: f64, m1: f64, Ec1: f64, Es1: f64, Z2: f64, m2: f64, Ec2: f64, Es2: f64, n2: f64, Eb2: f64) -> PyResult<Vec<[f64; 9]>> {
    Ok(simple_bca(x, y, z, ux, uy, uz, E1, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2))
}

#[cfg(feature = "python")]
/// simple_bca_list_py( energies, ux, uy, uz, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2)
/// --
///
/// This function runs a 0D Binary Collision Approximation simulation for the given incident ions and material.
/// Args:
///    energy (list(f64)): initial energies in eV.
///    ux (list(f64)): initial ion directions x.
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
pub fn simple_bca_list_py<'py>(energies: Vec<f64>, usx: Vec<f64>, usy: Vec<f64>, usz: Vec<f64>, Z1: f64, m1: f64, Ec1: f64, Es1: f64, Z2: f64, m2: f64, Ec2: f64, Es2: f64, n2: f64, Eb2: f64) -> PyResult<Vec<[f64; 9]>> {

    assert_eq!(energies.len(), usx.len());
    assert_eq!(energies.len(), usy.len());
    assert_eq!(energies.len(), usz.len());

    let x = -2.*(n2*1E30).powf(-1./3.);
    let y = 0.0;
    let z = 0.0;

    let mut total_output = vec![];
    for (((E1, ux), uy), uz) in energies.iter().zip(usx).zip(usy).zip(usz) {
        let output = simple_bca(x, y, z, ux, uy, uz, *E1, Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, n2, Eb2);
        for particle in output {
            total_output.push(particle);
        }
    }
    Ok(total_output)
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

    let mut rng = ChaCha8Rng::from_rng(&mut rand::rng());
    let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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

    let mut rng = ChaCha8Rng::from_rng(&mut rand::rng());
    let output = bca::single_ion_bca(p, &m, &options, &mut rng);

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
#[unsafe(no_mangle)]
pub extern "C" fn rotate_given_surface_normal(nx: f64, ny: f64, nz: f64, ux: &mut f64, uy: &mut f64, uz: &mut f64) {

    let direction = Vector3::new(*ux, *uy, *uz);
    let n = Vector3::new(nx, ny, nz);

    let (b1, b2) = duff_orthonormal_basis(Vector::new(-nx, -ny, -nz));
    let e1 = Vector3::new(b1.x, b1.y, b1.z);
    let e2 = Vector3::new(b2.x, b2.y, b2.z);
    let rotation_matrix_duff = Matrix3::from_columns(&[-n, e1, e2]).transpose();
    // Duff et al. provide a robust algorithm that constructs an orthonormal basis from n
    // That basis is used to construct an R such that R ex = -n, R ey = e1, R ez = e2.
    // The transpose of this matrix gives the matrix we want, R^T n = ex.
    // That is, R maps the global normal vector onto the RustBCA normal vector.
    // And thus R maps a global particle velocity into the RustBCA frame.

    let incident = rotation_matrix_duff*direction;

    *ux = incident.x;
    *uy = incident.y;
    *uz = incident.z;
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
pub fn rotate_given_surface_normal_py<'py>(nx: f64, ny: f64, nz: f64, ux: f64, uy: f64, uz: f64) -> PyResult<(f64, f64, f64)> {
    let mut ux = ux;
    let mut uy = uy;
    let mut uz = uz;
    rotate_given_surface_normal(nx, ny, nz, &mut ux, &mut uy, &mut uz);
    Ok((ux, uy, uz))
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
pub fn rotate_given_surface_normal_vec_py<'py>(nx: Vec<f64>, ny: Vec<f64>, nz: Vec<f64>, ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>)> {

    let length = nx.len();

    let mut ux_new = Vec::with_capacity(length);
    let mut uy_new = Vec::with_capacity(length);
    let mut uz_new = Vec::with_capacity(length);

    (0..length).into_iter().for_each(|index| {

        let mut ux_ = ux[index];
        let mut uy_ = uy[index];
        let mut uz_ = uz[index];

        rotate_given_surface_normal(nx[index], ny[index], nz[index], &mut ux_, &mut uy_, &mut uz_);
        ux_new.push(ux_);
        uy_new.push(uy_);
        uz_new.push(uz_);

    });

    Ok((ux_new, uy_new, uz_new))
}

#[cfg(feature = "parry3d")]
#[unsafe(no_mangle)]
pub extern "C" fn rotate_back(nx: f64, ny: f64, nz: f64, ux: &mut f64, uy: &mut f64, uz: &mut f64) {

    let direction = Vector3::new(*ux, *uy, *uz);
    let n = Vector3::new(nx, ny, nz);

    let (b1, b2) = duff_orthonormal_basis(Vector::new(-nx, -ny, -nz));
    let e1 = Vector3::new(b1.x, b1.y, b1.z);
    let e2 = Vector3::new(b2.x, b2.y, b2.z);
    let rotation_matrix_duff = Matrix3::from_columns(&[-n, e1, e2]);
    // Duff et al. provide a robust algorithm that constructs an orthonormal basis from n
    // That basis is used to construct an R such that R ex = -n, R ey = e1, R ez = e2.
    // This is the transpose of the matrix in rotate_given_surface_normal.
    // Since, for rotation matrices, R^T = R^-1, this is the inverse transform.

    let incident = rotation_matrix_duff*direction;

    *ux = incident.x;
    *uy = incident.y;
    *uz = incident.z;
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
pub fn rotate_back_py<'py>(nx: f64, ny: f64, nz: f64, ux: f64, uy: f64, uz: f64) -> PyResult<(f64, f64, f64)> {
    let mut ux = ux;
    let mut uy = uy;
    let mut uz = uz;
    rotate_back(nx, ny, nz, &mut ux, &mut uy, &mut uz);
    Ok((ux, uy, uz))
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
pub fn rotate_back_vec_py<'py>(nx: Vec<f64>, ny: Vec<f64>, nz: Vec<f64>, ux: Vec<f64>, uy: Vec<f64>, uz: Vec<f64>) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>)> {

    let (ux_new, (uy_new, uz_new)) = (nx, ny, nz, ux, uy, uz).into_par_iter().map(|(nx_, ny_, nz_, ux_, uy_, uz_)| {

        let mut ux_mut = ux_;
        let mut uy_mut = uy_;
        let mut uz_mut = uz_;
        rotate_back(nx_, ny_, nz_, &mut ux_mut, &mut uy_mut, &mut uz_mut);

        (ux_mut, (uy_mut, uz_mut))
    }).unzip();

    Ok((ux_new, uy_new, uz_new))
}

#[cfg(feature = "python")]
#[pyfunction]
/// sputtering_yield(ion, target, energy, angle, num_samples)
/// A routine the calculates the sputtering yield in atoms per ion of energetic ions incident upon materials using RustBCA.
/// Args:
///     ion: a dictionary with the keys Z (atomic number), m (atomic mass in AMU), Ec (cutoff energy in eV), Es (surface binding energy in eV)
///     target: a dictionary with the keys Z, m, Ec, Es, Eb (bulk binding energy in eV), n2 (number density in 1/m3)
///     energy: the incident energy of the ion in eV
///     angle: incident angle of the ion in degrees from surface normal
///     num_samples: number of ion trajectories to run; precision will go as 1/sqrt(N)
pub fn sputtering_yield<'py>(ion: &Bound<'py, PyDict>, target: &Bound<'py, PyDict>, energy: f64, angle: f64, num_samples: usize) -> PyResult<f64> {

    assert!(angle.abs() <= 90.0, "Incident angle w.r.t. surface normal, {}, cannot exceed 90 degrees.", angle);

    let Z1: f64 = ion.get_item("Z")?.expect("Error: Cannot get key 'Z' from ion dict.").extract()?;
    let m1: f64 = ion.get_item("m")?.expect("Error: Cannot get key 'm' from ion dict.").extract()?;
    let Es1: f64 = ion.get_item("Es")?.expect("Error: Cannot get key 'Es' from ion dict.").extract()?;
    let Ec1: f64 = ion.get_item("Ec")?.expect("Error: Cannot get key 'Ec' from ion dict.").extract()?;

    let Z2: f64 = target.get_item("Z")?.expect("Error: Cannot get key 'Z' from target dict.").extract()?;
    let m2: f64 = target.get_item("m")?.expect("Error: Cannot get key 'm' from target dict.").extract()?;
    let Es2: f64 = target.get_item("Es")?.expect("Error: Cannot get key 'Es' from target dict.").extract()?;
    let Ec2: f64 = target.get_item("Ec")?.expect("Error: Cannot get key 'Ec' from target dict.").extract()?;
    let Eb2: f64 = target.get_item("Eb")?.expect("Error: Cannot get key 'Eb' from target dict.").extract()?;
    let n2: f64 = target.get_item("n")?.expect("Error: Cannot get key 'n' from target dict.").extract()?;

    let options = Options::default_options(true);

    let y = 0.0;
    let z = 0.0;

    let ux = (angle/180.0*PI).cos();
    let uy = (angle/180.0*PI).sin();
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

    let seed: u64 = get_seed().map_err(|error| PyValueError::new_err(""))?;

    (0..num_samples as u64).into_par_iter()
    .for_each_init(
        || ChaCha8Rng::seed_from_u64(seed), |rng, index| {

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
        
        rng.set_stream(index);
        let output = bca::single_ion_bca(p, &m, &options, rng);

        for particle in output {
            if particle.E > 0.0 && particle.dir.x < 0.0 && particle.left && (!particle.incident) {
                let mut num_sputtered = num_sputtered.lock().unwrap();
                *num_sputtered += 1;
            }
        }
    });
    let num_sputtered = *num_sputtered.lock().unwrap();
    Ok(num_sputtered as f64 / num_samples as f64)
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
pub fn reflection_coefficient<'py>(ion: &Bound<'py, PyDict>, target: &Bound<'py, PyDict>, energy: f64, angle: f64, num_samples: usize) -> PyResult<(f64, f64)> {

    assert!(angle.abs() <= 90.0, "Incident angle w.r.t. surface normal, {}, cannot exceed 90 degrees.", angle);

    let Z1: f64 = ion.get_item("Z")?.expect("Error: Cannot get key 'Z' from ion dict.").extract()?;
    let m1: f64 = ion.get_item("m")?.expect("Error: Cannot get key 'm' from ion dict.").extract()?;
    let Es1: f64 = ion.get_item("Es")?.expect("Error: Cannot get key 'Es' from ion dict.").extract()?;
    let Ec1: f64 = ion.get_item("Ec")?.expect("Error: Cannot get key 'Ec' from ion dict.").extract()?;

    let Z2: f64 = target.get_item("Z")?.expect("Error: Cannot get key 'Z' from target dict.").extract()?;
    let m2: f64 = target.get_item("m")?.expect("Error: Cannot get key 'm' from target dict.").extract()?;
    let Es2: f64 = target.get_item("Es")?.expect("Error: Cannot get key 'Es' from target dict.").extract()?;
    let Ec2: f64 = target.get_item("Ec")?.expect("Error: Cannot get key 'Ec' from target dict.").extract()?;
    let Eb2: f64 = target.get_item("Eb")?.expect("Error: Cannot get key 'Eb' from target dict.").extract()?;
    let n2: f64 = target.get_item("n")?.expect("Error: Cannot get key 'n' from target dict.").extract()?;

    let options = Options::default_options(false);

    let y = 0.0;
    let z = 0.0;

    let ux = (angle/180.0*PI).cos();
    let uy = (angle/180.0*PI).sin();
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
    let residue = Mutex::new(0.0);

    let seed: u64 = get_seed().map_err(|error| PyValueError::new_err(""))?;

    (0..num_samples as u64).into_par_iter()
    .for_each_init(
        || ChaCha8Rng::seed_from_u64(seed), |rng, index| {
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
        
        rng.set_stream(index);
        let output = bca::single_ion_bca(p, &m, &options, rng);

        for particle in output {
            if particle.E > 0.0 && particle.dir.x < 0.0 && particle.left && particle.incident {
                let mut num_reflected = num_reflected.lock().unwrap();
                *num_reflected += 1;

                let mut energy_reflected = energy_reflected.lock().unwrap();

                let residue_part;

                // Use Moller-Knuth TwoSum to preserve deterministic fp reduce
                (*energy_reflected, residue_part) = moller_knuth_two_sum(*energy_reflected, particle.E);

                let mut residue = residue.lock().unwrap();
                *residue = *residue + residue_part;
            }
        };
    });
    if let (Ok(num_reflected), Ok(energy_reflected), Ok(residue)) = (num_reflected.lock(), energy_reflected.lock(), residue.lock()) {
        return Ok((*num_reflected as f64 / num_samples as f64, (*energy_reflected + *residue) / EV / energy / num_samples as f64))
    } else {
        return Err(PyValueError::new_err("Check input values."))
    }

    
}

fn get_seed() -> Result<u64> {
    match env::var("LIBRUSTBCA_SEED") {
        Ok(seed) if seed == "-1" => Ok(rand::random()),
        Ok(seed) => Ok(seed.parse::<u64>()?),
        Err(env::VarError::NotPresent) => Ok(0_u64),
        Err(env::VarError::NotUnicode(_)) => Err(anyhow!("LIBRUSTBCA_SEED not unicode."))
    }
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
pub fn compound_reflection_coefficient<'py>(ion: &Bound<'py, PyDict>, targets: Vec<Bound<'py, PyDict>>, target_number_densities: Vec<f64>, energy: f64, angle: f64, num_samples: usize) -> PyResult<(f64, f64)> {

    assert!(angle.abs() <= 90.0, "Incident angle w.r.t. surface normal, {}, cannot exceed 90 degrees.", angle);

    let Z1: f64 = ion.get_item("Z")?.expect("Error: Cannot get key 'Z' from ion dict.").extract()?;
    let m1: f64 = ion.get_item("m")?.expect("Error: Cannot get key 'm1' from ion dict.").extract()?;
    let Es1: f64 = ion.get_item("Es")?.expect("Error: Cannot get key 'Es' from ion dict.").extract()?;
    let Ec1: f64 = ion.get_item("Ec")?.expect("Error: Cannot get key 'Ec' from ion dict.").extract()?;

    let Z2: Vec<f64> = targets.iter()
        .enumerate()
        .map(|(index, target)| target.get_item("Z").unwrap()
        .unwrap_or_else(|| panic!(
            "Error: cannot get key 'Z' from target dict at index {}.", index
        ))
        .extract().unwrap()).collect::<Vec<f64>>();
    let m2: Vec<f64> = targets.iter()
        .enumerate()
        .map(|(index, target)| target.get_item("m").unwrap()
        .unwrap_or_else(|| panic!(
            "Error: cannot get key 'm' from target dict at index {}.", index
        ))
        .extract().unwrap()).collect::<Vec<f64>>();
    let Es2: Vec<f64> = targets.iter()
        .enumerate()
        .map(|(index, target)| target.get_item("Es").unwrap()
        .unwrap_or_else(|| panic!(
            "Error: cannot get key 'Es' from target dict at index {}.", index
        ))
        .extract().unwrap()).collect::<Vec<f64>>();
    let Ec2: Vec<f64> = targets.iter()
        .enumerate()
        .map(|(index, target)| target.get_item("Ec").unwrap()
        .unwrap_or_else(|| panic!(
            "Error: cannot get key 'Ec' from target dict at index {}.", index
        ))
        .extract().unwrap()).collect::<Vec<f64>>();
    let Eb2: Vec<f64> = targets.iter()
        .enumerate()
        .map(|(index, target)| target.get_item("Eb").unwrap()
        .unwrap_or_else(|| panic!(
            "Error: cannot get key 'Eb' from target dict at index {}.", index
        ))
        .extract().unwrap()).collect::<Vec<f64>>();

    let number_target_species = Z2.len();

    let options = Options::default_options(false);

    let y = 0.0;
    let z = 0.0;

    let ux = (angle/180.0*PI).cos();
    let uy = (angle/180.0*PI).sin();
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
    let residue = Mutex::new(0.0);

    let seed: u64 = get_seed().map_err(|error| PyValueError::new_err(""))?;

    (0..num_samples as u64).into_par_iter()
    .for_each_init(
        || ChaCha8Rng::seed_from_u64(seed), |rng, index| {

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
        
        rng.set_stream(index);
        let output = bca::single_ion_bca(p, &m, &options, rng);

        for particle in output {
            if particle.E > 0.0 && particle.dir.x < 0.0 && particle.left && particle.incident {
                let mut num_reflected = num_reflected.lock().unwrap();
                *num_reflected += 1;
                let mut energy_reflected = energy_reflected.lock().unwrap();

                let residue_part;

                // Use Moller-Knuth TwoSum to preserve deterministic fp reduce
                (*energy_reflected, residue_part) = moller_knuth_two_sum(*energy_reflected, particle.E);

                let mut residue = residue.lock().unwrap();
                *residue = *residue + residue_part;
            }
        }
    });
    let num_reflected = *num_reflected.lock().unwrap();
    let energy_reflected = *energy_reflected.lock().unwrap();
    let residue = *residue.lock().unwrap();

    Ok((num_reflected as f64 / num_samples as f64, (energy_reflected + residue) / EV / energy / num_samples as f64))
}

/// Moller-Knuth TwoSum Floating-Point Adder with Residual (FPAR)
/// This function allows one to use the identity: 
/// Given two floating point numbers a, b;
/// And the sum s = IEEE754RoundToNearest(a + b);
/// And the residual from floating point error r = (a + b) - s;
/// The following is invariant: s + r = a + b
/// citation: Accurate Parallel Floating-Point Accumulation
/// E. Kadric et al., IEEE Transactions on Computers 65 11
/// doi: 10.1109/TC.2016.2532874
#[cfg(feature = "python")]
fn moller_knuth_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let b_prime = s - a;
    let a_prime = s - b_prime;
    let delta_b = b - b_prime;
    let delta_a = a - a_prime;
    let r = delta_a + delta_b;
    (s, r)
}

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (Za, Zb, Ma, Mb, E0, p, n_gl_points=100, interaction_potential="KR_C"))]
fn scattering_integrals(Za: f64, Zb: f64, Ma: f64, Mb: f64, E0: f64, p: f64, n_gl_points: usize, interaction_potential: &str) -> PyResult<(f64, f64, f64, f64)> {
    let E0 = E0*EV;
    let p = p*ANGSTROM;

    let potential = match interaction_potential {
        "KR_C" => InteractionPotential::KR_C,
        "LENZ_JENSEN" => InteractionPotential::LENZ_JENSEN,
        "MOLIERE" => InteractionPotential::MOLIERE,
        "ZBL" => InteractionPotential::ZBL,
        _ => return Err(PyValueError::new_err(format!("Unimplemented interaction potential {}; try 'KR_C'", interaction_potential)))
    };

    let screening_length = interactions::screening_length(Za, Zb, potential);

    let x0_newton = bca::newton_rootfinder(Za, Zb, Ma, Mb, E0, p, potential, 1000, 1E-12).map_err(
        |error| PyRuntimeError::new_err(format!("Rootfinder failed to find distance of closest approach; check input values."))
    )?;

    //Compute center of mass deflection angle with each algorithm
    let theta_gm = bca::gauss_mehler(Za, Zb, Ma, Mb, E0, p, x0_newton, screening_length, potential, n_gl_points);
    let theta_gl = bca::gauss_legendre(Za, Zb, Ma, Mb, E0, p, x0_newton, screening_length, potential);
    let theta_mw = bca::mendenhall_weller(Za, Zb, Ma, Mb, E0, p, x0_newton, screening_length, potential);
    let theta_magic = bca::magic(Za, Zb, Ma, Mb, E0, p, x0_newton, screening_length, potential);

    Ok((theta_gm, theta_gl, theta_mw, theta_magic))
}
#[cfg(feature = "python")]
macro_rules! geometry_typed_loops {
    ($geometry_type:ty, $input:expr, $python:expr) => {
        {
            let input: <$geometry_type as geometry::Geometry>::InputFileFormat = depythonize(&$input).unwrap();
            let (particle_input_array, material, options, output_units) = input::process_input_file(input);
            let pool = rayon::ThreadPoolBuilder::new().num_threads(options.num_threads).build().unwrap();
            pool.install( ||
                physics::physics_loop::<$geometry_type>(particle_input_array, material, options, output_units)
            );
            Ok(())
        }
    }
}

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature=(input, geometry_mode="1D"))]
fn rustbca_py<'py>(python: Python<'py>, input: &Bound<'py, PyDict>, geometry_mode: &str) -> PyResult<()> {
    match geometry_mode {
        "0D" => geometry_typed_loops!(Mesh0D, input, python),
        "1D" => geometry_typed_loops!(Mesh1D, input, python),
        "2D" => geometry_typed_loops!(Mesh2D, input, python),
        "HOMOGENEOUS2D" => geometry_typed_loops!(Mesh2D, input, python),
        "SPHERE" => geometry_typed_loops!(Sphere, input, python),
        #[cfg(feature="parry3d")]
        "BALL" => geometry_typed_loops!(ParryBall, input, python),
        #[cfg(feature="parry3d")]
        "TRIMESH" => geometry_typed_loops!(ParryTriMesh, input, python),
       _ => Err(PyValueError::new_err(format!("Input Error: Unimplemented geometry mode {}; try '1D'", geometry_mode)))
    }
}

/*
Notes on macros - this is the first I have written, so I'm taking notes here as I go.
macro_rules! makes a macro - here, the macro is called geometry_types_silent_loops
macros pattern match an argument and replace it with anything you want
I want it to take a tuple of a string (e.g., "1D") and a type (e.g., Mesh1D)
and plop those into corresponding match arms.
The first line tells the macro to expect an argument with that pattern.
arguments are $<name>:<designator>. Designators:
block
expr is used for expressions
ident is used for variable/function names
item
literal is used for literal constants
pat (pattern)
path
stmt (statement)
tt (token tree)
ty (type)
vis (visibility qualifier)
*/
#[cfg(feature = "python")]
macro_rules! geometry_typed_silent_loops {
    ($geometry_type:ty, $input:expr, $python:expr) => {
        {
            let input: <$geometry_type as geometry::Geometry>::InputFileFormat = depythonize(&$input).unwrap();
            let (particle_input_array, material, options, output_units) = input::process_input_file(input);
            let pool = rayon::ThreadPoolBuilder::new().num_threads(options.num_threads).build().unwrap();
            let finished_particles = pool.install( ||
                physics::silent_physics_loop::<$geometry_type>(particle_input_array, material, options, output_units.clone())
            );
            let finished_particles_container = physics::process_finished_particles_to_arrays(finished_particles, output_units);
            Ok(pythonize($python, &finished_particles_container)?)
        }
    }
}

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature=(input, geometry_mode="1D"))]
fn rustbca_local_py<'py>(python: Python<'py>, input: &Bound<'py, PyDict>, geometry_mode: &str) -> PyResult<Bound<'py, PyAny>> {

    match geometry_mode {
        "0D" => geometry_typed_silent_loops!(Mesh0D, input, python),
        "1D" => geometry_typed_silent_loops!(Mesh1D, input, python),
        "2D" => geometry_typed_silent_loops!(Mesh2D, input, python),
        "HOMOGENEOUS2D" => geometry_typed_silent_loops!(Mesh2D, input, python),
        "SPHERE" => geometry_typed_silent_loops!(Sphere, input, python),
        #[cfg(feature="parry3d")]
        "BALL" => geometry_typed_silent_loops!(ParryBall, input, python),
        #[cfg(feature="parry3d")]
        "TRIMESH" => geometry_typed_silent_loops!(ParryTriMesh, input, python),
        _ => Err(PyValueError::new_err(format!("Input Error: Unimplemented geometry mode {}; try '1D'", geometry_mode)))
    }
}