use super::*;

#[derive(Deserialize)]
pub struct ParticleParameters {
    //pub particle_input_filename: String,
    pub length_unit: String,
    pub energy_unit: String,
    pub mass_unit: String,
    pub N: Vec<usize>,
    pub m: Vec<f64>,
    pub Z: Vec<f64>,
    pub E: Vec<f64>,
    pub Ec: Vec<f64>,
    pub Es: Vec<f64>,
    pub pos: Vec<(f64, f64, f64)>,
    pub dir: Vec<(f64, f64, f64)>,
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
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
    pub uz: f64
}

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
    pub asympototic_deflection: f64,
    pub stopped: bool,
    pub left: bool,
    pub incident: bool,
    pub first_step: bool,
    pub trajectory: Vec<Vector4>,
    pub track_trajectories: bool,
    pub number_collision_events: usize,
    pub backreflected: bool,
}
impl Particle {
    pub fn new(m: f64, Z: f64, E: f64, Ec: f64, Es: f64, x: f64, y: f64, z: f64, dirx: f64, diry: f64, dirz: f64, incident: bool, track_trajectories: bool) -> Particle {
        let dir_mag = (dirx*dirx + diry*diry + dirz*dirz).sqrt();

        Particle {
            m: m,
            Z: Z,
            E: E,
            Ec: Ec,
            Es: Es,
            pos: Vector::new(x, y, z),
            dir: Vector::new(dirx/dir_mag, diry/dir_mag, dirz/dir_mag),
            pos_old: Vector::new(x, y, z),
            dir_old: Vector::new(dirx/dir_mag, diry/dir_mag, dirz/dir_mag),
            pos_origin: Vector::new(x, y, z),
            asympototic_deflection: 0.,
            stopped: false,
            left: false,
            incident: incident,
            first_step: incident,
            trajectory: vec![Vector4::new(E, x, y, z)],
            track_trajectories: track_trajectories,
            number_collision_events: 0,
            backreflected: false
        }
    }
    pub fn add_trajectory(&mut self) {
        if self.track_trajectories {
            self.trajectory.push(Vector4 {E: self.E, x: self.pos.x, y: self.pos.y, z: self.pos.z});
        }
    }

    pub fn get_momentum(&mut self) -> Vector {
        let speed = (2.*self.E/self.m).sqrt();
        Vector::new(
            self.m*speed*self.dir.x,
            self.m*speed*self.dir.y,
            self.m*speed*self.dir.z,
        )
    }
}

pub fn rotate_particle(particle_1: &mut particle::Particle, psi: f64, phi: f64) {
    let cosx: f64 = particle_1.dir.x;
    let cosy: f64 = particle_1.dir.y;
    let cosz: f64 = particle_1.dir.z;
    let cphi: f64 = phi.cos();
    let sphi: f64 = phi.sin();
    let sa = (1. - cosx*cosx).sqrt();

    //Particle direction update formulas from TRIDYN, see Moeller and Eckstein 1988
    let cpsi: f64 = psi.cos();
    let spsi: f64 = psi.sin();
    let cosx_new: f64 = cpsi*cosx + spsi*cphi*sa;
    let cosy_new: f64 = cpsi*cosy - spsi/sa*(cphi*cosx*cosy - sphi*cosz);
    let cosz_new: f64 = cpsi*cosz - spsi/sa*(cphi*cosx*cosz + sphi*cosy);

    let dir_new = Vector {x: cosx_new, y: cosy_new, z: cosz_new};

    particle_1.dir.assign(&dir_new);
    particle_1.dir.normalize();
}

pub fn particle_advance(particle_1: &mut particle::Particle, mfp: f64, asympototic_deflection: f64) -> f64 {

    if particle_1.E > particle_1.Ec {
        particle_1.add_trajectory();
    }

    //Update previous position
    particle_1.pos_old.x = particle_1.pos.x;
    particle_1.pos_old.y = particle_1.pos.y;
    particle_1.pos_old.z = particle_1.pos.z;

    //In order to keep average denisty constant, must add back previous asymptotic deflection
    let distance_traveled = mfp + particle_1.asympototic_deflection - asympototic_deflection;

    particle_1.pos.x += particle_1.dir_old.x*distance_traveled;
    particle_1.pos.y += particle_1.dir_old.y*distance_traveled;
    particle_1.pos.z += particle_1.dir_old.z*distance_traveled;
    particle_1.asympototic_deflection = asympototic_deflection;

    //Update previous direction
    particle_1.dir_old.x = particle_1.dir.x;
    particle_1.dir_old.y = particle_1.dir.y;
    particle_1.dir_old.z = particle_1.dir.z;

    return distance_traveled;
}

pub fn refraction_angle(costheta: f64, energy_old: f64, energy_new: f64) -> f64 {
    //println!("energy_old: {} energy_new: {} costheta: {}", energy_old/EV, energy_new/EV, costheta);
    let costheta = if costheta.abs() > 1. {costheta.signum()} else {costheta};
    let sintheta0 = (1. - costheta*costheta).sqrt();
    let sintheta1 = sintheta0*(energy_old/energy_new).sqrt();
    let delta_theta = sintheta1.asin() - sintheta0.asin();
    assert!(!delta_theta.is_nan(), "Numerical error: refraction returned NaN.");
    let sign = -costheta.signum();
    return sign*delta_theta;
}
