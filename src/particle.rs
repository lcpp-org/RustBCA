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
    pub Ed: f64,
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
    pub trajectory: Vec<TrajectoryElement>,
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
        assert!(input.E > 0., "Input error: incident energy {}; must be greater than zero.", input.E/EV);

        Particle {
            m: input.m,
            Z: input.Z,
            E: input.E,
            Ec: input.Ec,
            Es: input.Es,
            Ed: 0.0,
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
            trajectory: vec![TrajectoryElement::new(input.E, input.x, input.y, input.z)],
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
    pub fn new(m: f64, Z: f64, E: f64, Ec: f64, Es: f64, Ed: f64, x: f64, y: f64, z: f64, dirx: f64, diry: f64, dirz: f64, incident: bool, track_trajectories: bool, interaction_index: usize) -> Particle {
        let dir_mag = (dirx*dirx + diry*diry + dirz*dirz).sqrt();

        Particle {
            m,
            Z,
            E,
            Ec,
            Es,
            Ed,
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

    /// Particle constructor from simplified input.
    pub fn default_incident(m_amu: f64, Z: f64, E_eV: f64, Ec_eV: f64, Es_eV: f64, x: f64, dirx: f64, diry: f64, dirz: f64) -> Particle {
        let dir_mag = (dirx*dirx + diry*diry + dirz*dirz).sqrt();

        let y = 0.;
        let z = 0.;

        assert!(E_eV > 0., "Input error: incident energy {}; must be greater than zero.", E_eV);

        Particle {
            m: m_amu*AMU,
            Z: Z,
            E: E_eV*EV,
            Ec: Ec_eV*EV,
            Es: Es_eV*EV,
            Ed: 0.0,
            pos: Vector::new(x, y, z),
            dir: Vector::new(dirx/dir_mag, diry/dir_mag, dirz/dir_mag),
            pos_old: Vector::new(x, y, z),
            dir_old: Vector::new(dirx/dir_mag, diry/dir_mag, dirz/dir_mag),
            pos_origin: Vector::new(x, y, z),
            energy_origin: E_eV,
            asymptotic_deflection: 0.,
            stopped: false,
            left: false,
            incident: true,
            first_step: true,
            trajectory: vec![],
            energies: vec![EnergyLoss::new(0., 0., x, y, z)],
            track_trajectories: false,
            number_collision_events: 0,
            backreflected: false,
            interaction_index: 0,
            weight: 1.0,
            tag: 0,
            tracked_vector: Vector::new(0.0, 0.0, 0.0),
        }
    }

    /// If `track_trajectories`, add the current (E, x, y, z) to the trajectory.
    pub fn update_trajectory_tracker(&mut self) {
        if self.track_trajectories {
            self.trajectory.push(TrajectoryElement {E: self.E, x: self.pos.x, y: self.pos.y, z: self.pos.z});
        }
    }

    /// If `track_energy_losses`, add the most recent electronic and nuclear energy loss terms and (x, y, z) to the energy loss tracker.
    pub fn update_energy_loss_tracker(&mut self, options: &Options, En: f64, Ee: f64) {
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

    /// Rotate a particle by a deflection psi at an azimuthal angle phi.
    pub fn rotate(&mut self, psi: f64, phi: f64) {
        let cosx: f64 = self.dir.x;
        let cosy: f64 = self.dir.y;
        let cosz: f64 = self.dir.z;
        let cosphi: f64 = (phi + PI).cos();
        let sinphi: f64 = (phi + PI).sin();

        let cpsi: f64 = psi.cos();
        let spsi: f64 = psi.sin();

        // To resolve the singularity, a different set of rotations is used when cosx == -1
        // Because of this, the recoil location is not consistent between the two formulas
        // Since phi is sampled uniformly from (0, 2pi), this does not matter
        // However, if a crystalline structure is ever added, this needs to be considered
        let cosx_new = if cosx > -1. {
            cpsi*cosx - spsi*(cosz*sinphi + cosy*cosphi)
        } else {
            cpsi*cosx - spsi*((1. + cosz - cosx*cosx)*cosphi - cosx*cosy*sinphi)/(1. + cosz)
        };

        let cosy_new = if cosx > -1. {
            cpsi*cosy + spsi*((1. + cosx - cosy*cosy)*cosphi - cosy*cosz*sinphi)/(1. + cosx)
        } else {
            cpsi*cosy + spsi*((1. + cosz - cosy*cosy)*sinphi - cosx*cosy*cosphi)/(1. + cosz)
        };

        let cosz_new = if cosx > -1. {
            cpsi*cosz + spsi*((1. + cosx - cosz*cosz)*sinphi - cosy*cosz*cosphi)/(1. + cosx)
        } else {
            cpsi*cosz + spsi*(cosx*cosphi + cosy*sinphi)
        };

        let dir_new = Vector {x: cosx_new, y: cosy_new, z: cosz_new};

        self.dir.assign(&dir_new);
        self.dir.normalize();
    }

    /// Push particle in space according to previous direction and return the distance traveled.
    pub fn advance(&mut self, mfp: f64, asymptotic_deflection: f64) -> f64 {

        if self.E > self.Ec {
            self.update_trajectory_tracker();
        }

        //Update previous position
        self.pos_old.x = self.pos.x;
        self.pos_old.y = self.pos.y;
        self.pos_old.z = self.pos.z;

        //In order to keep average denisty constant, must add back previous asymptotic deflection
        let distance_traveled = mfp + self.asymptotic_deflection - asymptotic_deflection;

        //dir has been updated, so use previous direction to advance in space
        self.pos.x += self.dir_old.x*distance_traveled;
        self.pos.y += self.dir_old.y*distance_traveled;
        self.pos.z += self.dir_old.z*distance_traveled;
        self.asymptotic_deflection = asymptotic_deflection;

        //Update previous direction
        self.dir_old.x = self.dir.x;
        self.dir_old.y = self.dir.y;
        self.dir_old.z = self.dir.z;

        return distance_traveled;
    }
}

pub fn surface_refraction(particle: &mut Particle, normal: Vector, Es: f64) {
    let E = particle.E;

    let costheta = particle.dir.dot(&normal);

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
