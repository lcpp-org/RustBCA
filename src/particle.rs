use super::*;

#[derive(Deserialize)]
pub struct ParticleParameters {
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

pub struct Particle {
    pub m: f64,
    pub Z: f64,
    pub E: f64,
    pub Ec: f64,
    pub Es: f64,
    pub pos: Vector,
    pub dir: Vector,
    pub pos_old: Vector,
    pub pos_origin: Vector,
    pub asympototic_deflection: f64,
    pub stopped: bool,
    pub left: bool,
    pub incident: bool,
    pub first_step: bool,
    pub trajectory: Vec<Vector4>,
    pub track_trajectories: bool,
    pub number_collision_events: usize
}
impl Particle {
    pub fn new(m: f64, Z: f64, E: f64, Ec: f64, Es: f64, x: f64, y: f64, z: f64, dirx: f64, diry: f64, dirz: f64, incident: bool, track_trajectories: bool) -> Particle {
        Particle {
            m: m,
            Z: Z,
            E: E,
            Ec: Ec,
            Es: Es,
            pos: Vector::new(x, y, z),
            dir: Vector::new(dirx, diry, dirz),
            pos_old: Vector::new(x, y, z),
            pos_origin: Vector::new(x, y, z),
            asympototic_deflection: 0.,
            stopped: false,
            left: false,
            incident: incident,
            first_step: incident,
            trajectory: Vec::new(),
            track_trajectories: track_trajectories,
            number_collision_events: 0
        }
    }
    pub fn add_trajectory(&mut self) {
        if self.track_trajectories {
            self.trajectory.push(Vector4 {E: self.E, x: self.pos.x, y: self.pos.y, z: self.pos.z});
        }
    }
}

pub fn rotate_particle(particle_1: &mut particle::Particle, psi: f64, phi: f64) {
    let ca: f64 = particle_1.dir.x;
    let cb: f64 = particle_1.dir.y;
    let cg: f64 = particle_1.dir.z;
    let cphi: f64 = phi.cos();
    let sphi: f64 = phi.sin();
    let sa = (1. - ca*ca).sqrt();

    //Particle direction update formulas from TRIDYN, see Moeller and Eckstein 1988
    let cpsi: f64 = psi.cos();
    let spsi: f64 = psi.sin();
    let ca_new: f64 = cpsi*ca + spsi*cphi*sa;
    let cb_new: f64 = cpsi*cb - spsi/sa*(cphi*ca*cb - sphi*cg);
    let cg_new: f64 = cpsi*cg - spsi/sa*(cphi*ca*cg + sphi*cb);

    let dir_new = Vector {x: ca_new, y: cb_new, z: cg_new};
    particle_1.dir.assign(&dir_new);
    particle_1.dir.normalize();
}

pub fn particle_advance(particle_1: &mut particle::Particle, mfp: f64, asympototic_deflection: f64) -> f64 {

    //Update previous position
    particle_1.pos_old.x = particle_1.pos.x;
    particle_1.pos_old.y = particle_1.pos.y;
    particle_1.pos_old.z = particle_1.pos.z;

    //In order to keep average denisty constant, must add back previous asymptotic deflection
    let distance_traveled = mfp + particle_1.asympototic_deflection - asympototic_deflection;

    particle_1.pos.x += particle_1.dir.x*distance_traveled;
    particle_1.pos.y += particle_1.dir.y*distance_traveled;
    particle_1.pos.z += particle_1.dir.z*distance_traveled;
    particle_1.asympototic_deflection = asympototic_deflection;

    return distance_traveled;
}
