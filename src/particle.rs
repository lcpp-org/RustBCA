use super::*;

/// Rustbca's internal representation of the particle_parameters input.

fn default_vec_zero() -> Vec<usize> {
    vec![0]
}

#[cfg(feature = "hdf5_input")]
#[derive(Deserialize, Clone)]
pub struct ParticleParameters {
    pub particle_input_filename: String,
    pub length_unit: String,
    pub energy_unit: String,
    pub mass_unit: String,
    pub N: Vec<usize>,
    pub m: Vec<f64>,
    pub Z: Vec<f64>,
    pub E: Vec<Distributions>,
    pub Ec: Vec<f64>,
    pub Es: Vec<f64>,
    pub pos: Vec<(Distributions, Distributions, Distributions)>,
    pub dir: Vec<(Distributions, Distributions, Distributions)>,
    #[serde(default = "default_vec_zero")]
    pub interaction_index: Vec<usize>,
}

#[cfg(not(feature = "hdf5_input"))]
#[derive(Deserialize, Clone)]
pub struct ParticleParameters {
    pub length_unit: String,
    pub energy_unit: String,
    pub mass_unit: String,
    pub N: Vec<usize>,
    pub m: Vec<f64>,
    pub Z: Vec<f64>,
    pub E: Vec<Distributions>,
    pub Ec: Vec<f64>,
    pub Es: Vec<f64>,
    pub pos: Vec<(Distributions, Distributions, Distributions)>,
    pub dir: Vec<(Distributions, Distributions, Distributions)>,
    #[serde(default = "default_vec_zero")]
    pub interaction_index: Vec<usize>,
}

/// HDF5 version of particle input.
#[derive(Clone, PartialEq, Debug, Copy)]
#[cfg_attr(feature = "hdf5_input", derive(hdf5::H5Type))]
#[repr(C)]
pub struct ParticleInput {
    pub m: f64,
    pub Z: f64,
    pub E: f64,
    pub Ec: f64,
    pub Es: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub ux: f64,
    pub uy: f64,
    pub uz: f64,
    pub interaction_index: usize,
    pub weight: f64,
    pub tag: i32,
}

/// Particle object. Particles in rustbca include incident ions and material atoms.
#[derive(Clone)]
pub struct Particle {
    pub m: f64,
    pub Z: f64,
    pub E: f64,
    pub Ec: f64,
    pub Es: f64,
    pub pos: Vector,
    pub dir: Vector,
    pub pos_old: Vector,
    pub dir_old: Vector,
    pub pos_origin: Vector,
    pub energy_origin: f64,
    pub asymptotic_deflection: f64,
    pub stopped: bool,
    pub left: bool,
    pub incident: bool,
    pub first_step: bool,
    pub trajectory: Vec<Vector4>,
    pub energies: Vec<EnergyLoss>,
    pub track_trajectories: bool,
    pub number_collision_events: usize,
    pub backreflected: bool,
    pub interaction_index: usize,
    pub weight: f64,
    pub tag: i32,
    pub tracked_vector: Vector,
}
impl Particle {
    /// Construct a particle object from input.
    pub fn from_input(input: ParticleInput, options: &Options) -> Particle {
        let dirx = input.ux;
        let diry = input.uy;
        let dirz = input.uz;

        let dir_mag = (dirx*dirx + diry*diry + dirz*dirz).sqrt();

        assert!((dirx/dir_mag).abs() < 1.0 - f64::EPSILON, "Input error: incident direction cannot round to exactly (1, 0, 0) due to gimbal lock. Use a non-zero y-component.");
        assert!(input.E > 0., "Input error: incident energy {}; must be greater than zero.", input.E/EV);

        Particle {
            m: input.m,
            Z: input.Z,
            E: input.E,
            Ec: input.Ec,
            Es: input.Es,
            pos: Vector::new(input.x, input.y, input.z),
            dir: Vector::new(dirx/dir_mag, diry/dir_mag, dirz/dir_mag),
            pos_old: Vector::new(input.x, input.y, input.z),
            dir_old: Vector::new(dirx/dir_mag, diry/dir_mag, dirz/dir_mag),
            pos_origin: Vector::new(input.x, input.y, input.z),
            energy_origin: input.E,
            asymptotic_deflection: 0.,
            stopped: false,
            left: false,
            incident: true,
            first_step: true,
            trajectory: vec![Vector4::new(input.E, input.x, input.y, input.z)],
            energies: vec![],
            track_trajectories: options.track_trajectories,
            number_collision_events: 0,
            backreflected: false,
            interaction_index: input.interaction_index,
            weight: 1.0,
            tag: 0,
            tracked_vector: Vector::new(0.0, 0.0, 0.0),
        }
    }

    /// Particle constructor from raw inputs.
    pub fn new(m: f64, Z: f64, E: f64, Ec: f64, Es: f64, x: f64, y: f64, z: f64, dirx: f64, diry: f64, dirz: f64, incident: bool, track_trajectories: bool, interaction_index: usize) -> Particle {
        let dir_mag = (dirx*dirx + diry*diry + dirz*dirz).sqrt();

        Particle {
            m,
            Z,
            E,
            Ec,
            Es,
            pos: Vector::new(x, y, z),
            dir: Vector::new(dirx/dir_mag, diry/dir_mag, dirz/dir_mag),
            pos_old: Vector::new(x, y, z),
            dir_old: Vector::new(dirx/dir_mag, diry/dir_mag, dirz/dir_mag),
            pos_origin: Vector::new(x, y, z),
            energy_origin: E,
            asymptotic_deflection: 0.,
            stopped: false,
            left: false,
            incident,
            first_step: incident,
            trajectory: vec![],
            energies: vec![EnergyLoss::new(0., 0., x, y, z)],
            track_trajectories,
            number_collision_events: 0,
            backreflected: false,
            interaction_index,
            weight: 1.0,
            tag: 0,
            tracked_vector: Vector::new(0.0, 0.0, 0.0),
        }
    }

    /// If `track_trajectories`, add the current (E, x, y, z) to the trajectory.
    pub fn add_trajectory(&mut self) {
        if self.track_trajectories {
            self.trajectory.push(Vector4 {E: self.E, x: self.pos.x, y: self.pos.y, z: self.pos.z});
        }
    }

    /// If `track_energy_losses`, add the most recent electronic and nuclear energy loss terms and (x, y, z) to the energy loss tracker.
    pub fn energy_loss(&mut self, options: &Options, En: f64, Ee: f64) {
        if self.incident & options.track_energy_losses {
            self.energies.push(EnergyLoss {Ee, En, x: self.pos.x, y: self.pos.y, z: self.pos.z});
        }
    }

    /// Get the current momentum.
    pub fn get_momentum(&mut self) -> Vector {
        let speed = (2.*self.E/self.m).sqrt();
        Vector::new(
            self.m*speed*self.dir.x,
            self.m*speed*self.dir.y,
            self.m*speed*self.dir.z,
        )
    }
}

/// Rotate a particle by a deflection psi at an azimuthal angle phi.
pub fn rotate_particle(particle_1: &mut particle::Particle, psi: f64, phi: f64) {
    let cosx: f64 = particle_1.dir.x;
    let cosy: f64 = particle_1.dir.y;
    let cosz: f64 = particle_1.dir.z;
    let cphi: f64 = phi.cos();
    let sphi: f64 = phi.sin();
    let sa = (1. - cosx*cosx).sqrt();

    //Particle direction update formulas from original TRIDYN paper, see Moeller and Eckstein 1988
    let cpsi: f64 = psi.cos();
    let spsi: f64 = psi.sin();
    let cosx_new: f64 = cpsi*cosx + spsi*cphi*sa;
    let cosy_new: f64 = cpsi*cosy - spsi/sa*(cphi*cosx*cosy - sphi*cosz);
    let cosz_new: f64 = cpsi*cosz - spsi/sa*(cphi*cosx*cosz + sphi*cosy);

    let dir_new = Vector {x: cosx_new, y: cosy_new, z: cosz_new};

    particle_1.dir.assign(&dir_new);
    particle_1.dir.normalize();
}

/// Push particle in space according to previous direction and return the distance traveled.
pub fn particle_advance(particle_1: &mut particle::Particle, mfp: f64, asymptotic_deflection: f64) -> f64 {

    if particle_1.E > particle_1.Ec {
        particle_1.add_trajectory();
    }

    //Update previous position
    particle_1.pos_old.x = particle_1.pos.x;
    particle_1.pos_old.y = particle_1.pos.y;
    particle_1.pos_old.z = particle_1.pos.z;

    //In order to keep average denisty constant, must add back previous asymptotic deflection
    let distance_traveled = mfp + particle_1.asymptotic_deflection - asymptotic_deflection;
    //let distance_traveled = mfp - asymptotic_deflection;

    //dir has been updated, so use previous direction to advance in space
    particle_1.pos.x += particle_1.dir_old.x*distance_traveled;
    particle_1.pos.y += particle_1.dir_old.y*distance_traveled;
    particle_1.pos.z += particle_1.dir_old.z*distance_traveled;
    particle_1.asymptotic_deflection = asymptotic_deflection;

    //Update previous direction
    particle_1.dir_old.x = particle_1.dir.x;
    particle_1.dir_old.y = particle_1.dir.y;
    particle_1.dir_old.z = particle_1.dir.z;

    return distance_traveled;
}

pub fn surface_refraction(particle: &mut Particle, normal: Vector, Es: f64) {
    let E = particle.E;

    let costheta = particle.dir.dot(&normal);

    let a = (E/(E + Es)).sqrt();
    let b = -(E).sqrt()*costheta;
    let c = (E*costheta.powi(2) + Es).sqrt();

    let u1x = (E/(E + Es)).sqrt()*particle.dir.x + ((-(E).sqrt()*costheta + (E*costheta.powi(2) + Es).sqrt())/(E + Es).sqrt())*normal.x;
    let u1y = (E/(E + Es)).sqrt()*particle.dir.y + ((-(E).sqrt()*costheta + (E*costheta.powi(2) + Es).sqrt())/(E + Es).sqrt())*normal.y;
    let u1z = (E/(E + Es)).sqrt()*particle.dir.z + ((-(E).sqrt()*costheta + (E*costheta.powi(2) + Es).sqrt())/(E + Es).sqrt())*normal.z;
    particle.dir.x = u1x;
    particle.dir.y = u1y;
    particle.dir.z = u1z;
    particle.E += Es;
}

/// Calcualte the refraction angle based on the surface binding energy of the material.
pub fn refraction_angle(costheta: f64, energy_old: f64, energy_new: f64) -> f64 {
    let costheta = if costheta.abs() > 1. {costheta.signum()} else {costheta};
    let sintheta0 = (1. - costheta*costheta).sqrt();
    let sintheta1 = sintheta0*(energy_old/energy_new).sqrt();
    let delta_theta = sintheta1.asin() - sintheta0.asin();
    assert!(!delta_theta.is_nan(), "Numerical error: refraction returned NaN.");
    let sign = -costheta.signum();
    return sign*delta_theta;
}
