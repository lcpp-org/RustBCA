extern crate geo;
extern crate serde;
extern crate toml;

pub use crate::geo::*;
//use geo::prelude::*;
use geo::algorithm::contains::Contains;
use geo::algorithm::closest_point::ClosestPoint;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::Read;
use std::f64::consts::PI;
use serde::*;
use std::io::BufWriter;
use std::process;
use std::thread;

//Physical constants
const Q: f64 = 1.60217646E-19;
const EV: f64 = Q;
const AMU: f64 = 1.660539E-27;
const ANGSTROM: f64 = 1E-10;
const MICRON: f64 = 1E-6;
const NM: f64 = 1E-9;
const CM: f64 = 1E-2;
const EPS0: f64 = 8.85418781E-12;
const A0: f64 = 5.29177211E-11;
//const K: f64 = 1.11265E-10;
const ME: f64 =  9.109383632E-31;
const SQRTPI: f64 = 1.772453850906;
const SQRT2PI: f64 = 2.506628274631;
const C: f64 = 299792000.;
const BETHE_BLOCH_PREFACTOR: f64 = 4.*PI*(Q*Q/(4.*PI*EPS0))*(Q*Q/(4.*PI*EPS0))/ME/C/C;
const LINDHARD_SCHARFF_PREFACTOR: f64 = 1.212*ANGSTROM*ANGSTROM*Q;
const LINDHARD_REDUCED_ENERGY_PREFACTOR: f64 = 4.*PI*EPS0/Q/Q;

const INTERPOLATED: i32 = 0;
const LOW_ENERGY_NONLOCAL: i32 = 1;
const LOW_ENERGY_LOCAL: i32 = 2;
const LOW_ENEGY_EQUIPARTITION: i32 = 3;

const LIQUID: i32 = 0;
const GASEOUS: i32 = 1;

const MOLIERE: i32 = 0;
const KR_C: i32 = 1;
const ZBL: i32 = 2;

#[derive(Deserialize)]
pub struct Input {
    options: Options,
    material_parameters: MaterialParameters,
    geometry: Geometry,
    particle_parameters: ParticleParameters
}

#[derive(Deserialize)]
pub struct Options {
    name: String,
    track_trajectories: bool,
    track_recoils: bool,
    track_recoil_trajectories: bool,
    write_files: bool,
    stream_size: usize,
    print: bool,
    print_num: usize,
    weak_collision_order: usize,
    suppress_deep_recoils: bool,
    high_energy_free_flight_paths: bool,
    electronic_stopping_mode: i32,
    mean_free_path_model: i32,
    interaction_potential: i32,
}

#[derive(Deserialize)]
pub struct MaterialParameters {
    energy_unit: String,
    mass_unit: String,
    Eb: f64,
    Es: f64,
    Ec: f64,
    n: f64,
    Z: f64,
    m: f64,
    electronic_stopping_correction_factor: f64
}

#[derive(Deserialize)]
pub struct ParticleParameters {
    length_unit: String,
    energy_unit: String,
    mass_unit: String,
    N: Vec<usize>,
    m: Vec<f64>,
    Z: Vec<f64>,
    E: Vec<f64>,
    Ec: Vec<f64>,
    Es: Vec<f64>,
    pos: Vec<(f64, f64, f64)>,
    dir: Vec<(f64, f64, f64)>,

}

#[derive(Deserialize)]
pub struct Geometry {
    length_unit: String,
    surface: Vec<(f64, f64)>,
    energy_surface: Vec<(f64, f64)>,
    simulation_surface: Vec<(f64, f64)>
}

pub struct Particle {
    m: f64,
    Z: f64,
    E: f64,
    Ec: f64,
    Es: f64,
    pos: Vector,
    dir: Vector,
    pos_old: Vector,
    pos_origin: Vector,
    asympototic_deflection: f64,
    stopped: bool,
    left: bool,
    incident: bool,
    first_step: bool,
    trajectory: Vec<Vector4>,
    track_trajectories: bool,
    number_collision_events: usize
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
    fn add_trajectory(&mut self) {
        if self.track_trajectories {
            self.trajectory.push(Vector4 {E: self.E, x: self.pos.x, y: self.pos.y, z: self.pos.z});
        }
    }
}

pub struct Vector {
    x: f64,
    y: f64,
    z: f64,
}
impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Vector {
        Vector {
            x: x,
            y: y,
            z: z
        }
    }
    fn magnitude(&self) -> f64 {
        return (self.x*self.x + self.y*self.y + self.z*self.z).sqrt();
    }
    fn assign(&mut self, other: &Vector) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
    }

    fn dot(&self, other: &Vector) -> f64 {
        return self.x*other.x + self.y*other.y + self.z*other.z;
    }

    fn normalize(&mut self) {
        let magnitude = self.magnitude();
        self.x /= magnitude;
        self.y /= magnitude;
        self.z /= magnitude;
    }
}

pub struct Vector4 {
    E: f64,
    x: f64,
    y: f64,
    z: f64,
}
impl Vector4 {
    fn new(E: f64, x: f64, y: f64, z: f64) -> Vector4 {
        Vector4 {
            E: E,
            x: x,
            y: y,
            z: z
        }
    }
}

pub struct Material {
    n: f64,
    m: f64,
    Z: f64,
    Eb: f64,
    Es: f64,
    Ec: f64,
    surface: Polygon<f64>,
    energy_surface: Polygon<f64>,
    simulation_surface: Polygon<f64>,
    electronic_stopping_correction_factor: f64
}
impl Material {
    pub fn new(material_parameters: MaterialParameters, geometry: Geometry) -> Material {

        //Multiply all coordinates by value of geometry unit.
        let length_unit: f64 = match geometry.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => panic!("Incorrect unit {} in input file. Choose one of: MICRON, CM, ANGSTROM, NM, M",
                geometry.length_unit.as_str())
        };
        let energy_unit: f64 = match material_parameters.energy_unit.as_str() {
            "EV" => EV,
            "J"  => 1.,
            "KEV" => EV*1E3,
            "MEV" => EV*1E6,
            _ => panic!("Incorrect unit {} in input file. Choose one of: EV, J, KEV, MEV",
                material_parameters.energy_unit.as_str())
        };
        let mass_unit: f64 = match material_parameters.mass_unit.as_str() {
            "AMU" => AMU,
            "KG" => 1.,
            _ => panic!("Incorrect unit {} in input file. Choose one of: AMU, KG",
                material_parameters.mass_unit.as_str())
        };

        let mut unit_coords = geometry.surface.clone();
        for pair in &mut unit_coords {
            pair.0 *= length_unit;
            pair.1 *= length_unit;
        }
        let mut unit_e_coords = geometry.energy_surface.clone();
        for pair in &mut unit_e_coords {
            pair.0 *= length_unit;
            pair.1 *= length_unit;
        }
        let mut unit_b_coords = geometry.simulation_surface.clone();
        for pair in &mut unit_b_coords {
            pair.0 *= length_unit;
            pair.1 *= length_unit;
        }

        Material {
            n: material_parameters.n,
            m: material_parameters.m*mass_unit,
            Z: material_parameters.Z,
            Eb: material_parameters.Eb*energy_unit,
            Es: material_parameters.Es*energy_unit,
            Ec: material_parameters.Ec*energy_unit,
            electronic_stopping_correction_factor: material_parameters.electronic_stopping_correction_factor,
            surface: Polygon::new(LineString::from(unit_coords), vec![]),
            energy_surface: Polygon::new(LineString::from(unit_e_coords), vec![]),
            simulation_surface: Polygon::new(LineString::from(unit_b_coords), vec![])
        }
    }

    fn inside(&self, x: f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.surface.contains(&p);
    }

    fn inside_1D(&self, x: f64) -> bool {
        return x > 0.;
    }

    fn inside_energy_barrier_1D(&self, x: f64) -> bool {
        let dx = -2.*self.n.powf(-1./3.)/SQRT2PI;
        return x > dx;
    }

    fn inside_simulation_boundary_1D(&self, x: f64) -> bool {
        let dx = -2.*self.n.powf(-1./3.)/SQRT2PI;
        return x > 2.*dx;
    }

    fn mfp(&self, x: f64, y: f64) -> f64 {
        return self.n.powf(-1./3.);
    }

    fn number_density(&self, x: f64, y: f64) -> f64 {
        return self.n;
    }

    fn inside_energy_barrier(&self, x: f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.energy_surface.contains(&p);
    }

    fn inside_simulation_boundary(&self, x:f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.simulation_surface.contains(&p);
    }

    fn closest_point(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.surface.closest_point(&p)
    }

    fn closest_point_on_energy_barrier(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.energy_surface.closest_point(&p)
    }

    fn closest_point_on_simulation_surface(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.simulation_surface.closest_point(&p)
    }

    fn Z_eff(&self, x: f64, y: f64) -> f64 {
        return self.Z;
    }

    fn m_eff(&self, x: f64, y: f64) -> f64 {
        return self.m;
    }

    fn Eb(&self, x: f64, y: f64) -> f64 {
        return self.Eb;
    }

    fn choose(&self, x: f64, y: f64) -> (f64, f64, f64, f64) {
        return (self.Z, self.m, self.Ec, self.Es);
    }

    fn electronic_stopping_power(&self, particle_1: &Particle, electronic_stopping_mode: i32) -> f64 {

        let n = self.number_density(particle_1.pos.x, particle_1.pos.y);
        let E = particle_1.E;
        let Ma = particle_1.m;
        let Za = particle_1.Z;
        let Zb = self.Z_eff(particle_1.pos.x, particle_1.pos.y);

        let beta = (1. - (1. + E/Ma/C.powf(2.)).powf(-2.)).sqrt();
        let v = beta*C;

        //See Biersack and Haggmark - empirical relation for mean ionization potential I
        let mut I0 = 0.;
        if Zb < 13. {
            I0 = 12. + 7./Zb
        } else {
            I0 = 9.76 + 58.5*Zb.powf(-1.19);
        }
        let I = Zb*I0*Q;

        //See Biersack and Haggmark - this looks like an empirical shell correction
        let mut B = 0.;
        if Zb < 3. {
            B = 100.*Za/Zb;
        } else {
            B = 5.;
        }

    //Bethe stopping modified by Biersack and Varelas
        let prefactor = BETHE_BLOCH_PREFACTOR*Zb*Za*Za/beta/beta;
        let eb = 2.*ME*v*v/I;
        let S_high = prefactor*(eb + 1. + B/eb).ln();

        //Lindhard-Scharff electronic stopping
        let S_low = LINDHARD_SCHARFF_PREFACTOR*(Za.powf(7./6.)*Zb)/(Za.powf(2./3.) + Zb.powf(2./3.)).powf(3./2.)*(E/Q/Ma*AMU).sqrt();

        let stopping_power = match electronic_stopping_mode {
            //Biersack-Varelas Interpolation
            INTERPOLATED => 1./(1./S_high + 1./S_low),
            //Oen-Robinson
            LOW_ENERGY_LOCAL => S_low,
            //Lindhard-Scharff
            LOW_ENERGY_NONLOCAL => S_low,
            //Lindhard-Scharff and Oen-Robinson, using Lindhard Equipartition
            LOW_ENERGY_EQUIPARTITION => S_low
        };

        return stopping_power;
    }
}

fn phi(xi: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => 0.35*(-0.3*xi).exp() + 0.55*(-1.2*xi).exp() + 0.10*(-6.0*xi).exp(),
        KR_C => 0.190945*(-0.278544*xi).exp() + 0.473674*(-0.637174*xi).exp() + 0.335381*(-1.919249*xi).exp(),
        ZBL => 0.02817*(-0.20162*xi).exp() + 0.28022*(-0.40290*xi).exp() + 0.50986*(-0.94229*xi).exp() + 0.18175*(-3.1998*xi).exp(),
        _ => panic!("Unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
    }
}

fn dphi(xi: f64, interaction_potential: i32) -> f64 {
    match interaction_potential {
        MOLIERE => -0.35*0.3*(-0.3*xi).exp() + -0.55*1.2*(-1.2*xi).exp() + -0.10*6.0*(-6.0*xi).exp(),
        KR_C => -0.278544*0.190945*(-0.278544*xi).exp() - 0.637174*0.473674*(-0.637174*xi).exp() - 0.335381*1.919249*(-1.919249*xi).exp(),
        ZBL => -0.20162*0.02817*(-0.20162*xi).exp() + -0.40290*0.28022*(-0.40290*xi).exp() + -0.94229*0.50986*(-0.94229*xi).exp() + -3.1998*0.18175*(-3.1998*xi).exp(),
        _ => panic!("Unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
    }
}

fn screening_length(Za: f64, Zb: f64) -> f64 {
    //Lindhard/Firsov screening length
    return 0.8853*A0*(Za.sqrt() + Zb.sqrt()).powf(-2./3.);
}

fn doca_function(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //Transcendental function to determine distance of closest approach
    return x0 - phi(x0, interaction_potential)/reduced_energy - beta*beta/x0;
}

fn diff_doca_function(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //First differential of distance of closest approach function
    return beta*beta/x0/x0 - dphi(x0, interaction_potential)/reduced_energy + 1.
}

fn f(x: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //Function for scattering integral - see Mendenhall and Weller, 1991 & 2005
    return (1. - phi(x, interaction_potential)/x/reduced_energy - beta*beta/x/x).powf(-0.5);
}

fn determine_mfp_phi_impact_parameter(particle_1: &mut Particle, material: &Material,
    collision_order: usize, high_energy_free_flight_paths: bool, mean_free_path_model: i32) -> (Vec<f64>, Vec<f64>, f64) {

    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let z = particle_1.pos.z;

    let mut mfp = material.mfp(x, y);

    //azimuthal angle randomly selected (0..2pi)
    let mut phis_azimuthal = Vec::with_capacity(collision_order + 1);

    //Each weak collision gets its own aziumuthal angle in annuli around collision point
    for k in 0..collision_order + 1 {
        phis_azimuthal.push(2.*PI*rand::random::<f64>());
    }

    if high_energy_free_flight_paths {

        let Ma: f64 = particle_1.m;
        let Mb: f64  = material.m_eff(x, y);
        let Za: f64  = particle_1.Z;
        let Zb: f64  = material.Z_eff(x, y);
        let n: f64  = material.number_density(x, y);
        let E: f64  = particle_1.E;
        let Ec: f64 = particle_1.Ec;
        let a: f64 = screening_length(Za, Zb);
        let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E;

        //Minimum energy transfer for scattering event set to cutoff energy
        let E_min = Ec*(Ma + Mb).powf(2.)/4./Ma/Mb;
        let reduced_energy_min: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E_min;

        //Free flight path formulation here described in SRIM textbook chapter 7, and Eckstein 1991 7.5.2
        let ep = (reduced_energy*reduced_energy_min).sqrt();
        let mut pmax = a/(ep + ep.sqrt() + 0.125*ep.powf(0.1));
        let mut ffp = 1./(n*pmax*pmax*PI);
        let delta_energy_electronic = material.electronic_stopping_power(particle_1, INTERPOLATED)*n*ffp*material.electronic_stopping_correction_factor;

        //If losing too much energy, scale free-flight-path down
        //5 percent limit set in original TRIM paper, Biersack and Haggmark 1980
        if delta_energy_electronic > 0.05*E {
            ffp = 0.05*E/delta_energy_electronic*ffp;
            pmax = (1./(n*PI*ffp)).sqrt()
        }

        //If free-flight-path less than the interatomic spacing, revert to solid model
        //Mentioned in Eckstein 1991, Ziegler, Biersack, and Ziegler 2008 (SRIM textbook 7-8)
        if ffp < mfp {

            ffp = mfp;
            //Cylindrical geometry
            pmax = mfp/SQRTPI;
            let mut impact_parameter = Vec::with_capacity(1);
            let random_number = rand::random::<f64>();
            let p = pmax*random_number.sqrt();
            impact_parameter.push(p);

            //Atomically rough surface - scatter initial collisions using mfp near interface
            if particle_1.first_step {
                ffp = mfp*rand::random::<f64>();
                particle_1.first_step = false;
            }

            if mean_free_path_model == GASEOUS {
                ffp *= -rand::random::<f64>().ln();
            }

            return (phis_azimuthal, impact_parameter, ffp);

        } else {

            //Impact parameter chosen as sqrt(-ln(R))*pmax, as in Biersack and Haggmark 1980,
            //And Mendenhall Weller 2005
            //And SRIM textbook chapter 7
            //And Eckstein 1991
            let mut impact_parameter = Vec::with_capacity(1);
            let random_number = rand::random::<f64>();
            let p = pmax*(-random_number.ln()).sqrt();
            impact_parameter.push(p);

            //Atomically rough surface - scatter initial collisions using mfp near interface
            if particle_1.first_step {
                ffp = mfp*rand::random::<f64>();
                particle_1.first_step = false;
            }

            if mean_free_path_model == GASEOUS {
                ffp *= -rand::random::<f64>().ln();
            }

            return (phis_azimuthal, impact_parameter, ffp);
        }

    } else {

        //If not using free flight paths, use weak collision model
        let pmax = mfp/SQRTPI;

        //Cylindrical geometry
        let mut impact_parameters = Vec::with_capacity(collision_order + 1);
        for k in 0..(collision_order + 1) {
            let random_number = rand::random::<f64>();
            let p = pmax*(random_number + k as f64).sqrt();
            impact_parameters.push(p)
        }

        //Atomically rough surface - scatter initial collisions
        if particle_1.first_step {
            mfp *= rand::random::<f64>();
            particle_1.first_step = false;
        }

        if mean_free_path_model == GASEOUS {
            mfp *= -rand::random::<f64>().ln();
        }

        return (phis_azimuthal, impact_parameters, mfp);

    }
}

fn choose_collision_partner(particle_1: &Particle, material: &Material, phi_azimuthal: f64, impact_parameter: f64, mfp: f64) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let z = particle_1.pos.z;

    //Determine cosines and sines
    let sphi: f64 = phi_azimuthal.sin();
    let ca: f64 = particle_1.dir.x;
    let cb: f64 = particle_1.dir.y;
    let cg: f64 = particle_1.dir.z;
    let sa: f64 = (1. - ca*ca).sqrt();
    let cphi: f64 = phi_azimuthal.cos();

    //Find recoil location
    let x_recoil: f64 = x + mfp*ca - impact_parameter*cphi*sa;
    let y_recoil: f64 = y + mfp*cb - impact_parameter*(sphi*cg - cphi*cb*ca)/sa;
    let z_recoil: f64 = z + mfp*cg + impact_parameter*(sphi*cb - cphi*ca*cg)/sa;

    //Choose recoil Z, M
    let (Z_recoil, M_recoil, Ec_recoil, Es_recoil) = material.choose(x, y);

    return (Z_recoil, M_recoil, Ec_recoil, Es_recoil, x_recoil, y_recoil, z_recoil, ca, cb, cg);
}

fn calculate_binary_collision(particle_1: &Particle, particle_2: &Particle, impact_parameter: f64, max_iter: usize, tol: f64, interaction_potential: i32) -> (f64, f64, f64, f64, f64, f64) {
    let Za: f64 = particle_1.Z;
    let Zb: f64 = particle_2.Z;
    let Ma: f64 = particle_1.m;
    let Mb: f64 = particle_2.m;
    let E0: f64 = particle_1.E;
    let mu: f64 = Mb/(Ma + Mb);

    //Lindhard screening length and reduced energy
    let a: f64 = screening_length(Za, Zb);
    //let reduced_energy: f64 = 4.*PI*EPS0*a*Mb*E0/(Ma + Mb)/Za/Zb/Q/Q;
    let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
    let beta: f64 = impact_parameter/a;

    //Guess for large reduced energy from Mendenhall and Weller 1991
    //For small energies, use pure Newton-Raphson with arbitrary guess of 1
    let mut x0 = 1.;
    let mut xn: f64;
    if reduced_energy > 5. {
        let inv_er_2 = 0.5/reduced_energy;
        x0 = inv_er_2 + (inv_er_2*inv_er_2 + beta*beta).sqrt();
    }

    //Newton-Raphson to determine distance of closest approach
    let mut err: f64;
    for _ in 0..max_iter {
        xn = x0 - doca_function(x0, beta, reduced_energy, interaction_potential)/diff_doca_function(x0, beta, reduced_energy, interaction_potential);
        err = (xn - x0).powf(2.);
        x0 = xn;
        if err < tol {
            break;
        }
    }

    //Scattering integral quadrature from Mendenhall and Weller 2005
    let lambda_0 = (0.5 + beta*beta/x0/x0/2. - dphi(x0, interaction_potential)/2./reduced_energy).powf(-1./2.);
    let alpha = 1./12.*(1. + lambda_0 + 5.*(0.4206*f(x0/0.9072, beta, reduced_energy, interaction_potential) + 0.9072*f(x0/0.4206, beta, reduced_energy, interaction_potential)));
    let theta = PI*(1. - beta*alpha/x0);

    //See Eckstein 1991 for details on center of mass and lab frame angles
    let asympototic_deflection = x0*a*(theta/2.).sin();
    let psi = (theta.sin().atan2(Ma/Mb + theta.cos())).abs();
    let psi_recoil = (theta.sin().atan2(1. - theta.cos())).abs();
    let recoil_energy = 4.*(Ma*Mb)/(Ma + Mb).powf(2.)*E0*(theta/2.).sin().powf(2.);

    return (theta, psi, psi_recoil, recoil_energy, asympototic_deflection, x0);
}

fn rotate_particle(particle_1: &mut Particle, psi: f64, phi: f64) {
    let ca: f64 = particle_1.dir.x;
    let cb: f64 = particle_1.dir.y;
    let cg: f64 = particle_1.dir.z;
    let cphi: f64 = phi.cos();
    let sphi: f64 = phi.sin();
    let sa = (1. - ca*ca).sqrt();

    //Particle direction update from TRIDYN, see Moeller and Eckstein 1988
    let cpsi: f64 = psi.cos();
    let spsi: f64 = psi.sin();
    let ca_new: f64 = cpsi*ca + spsi*cphi*sa;
    let cb_new: f64 = cpsi*cb - spsi/sa*(cphi*ca*cb - sphi*cg);
    let cg_new: f64 = cpsi*cg - spsi/sa*(cphi*ca*cg + sphi*cb);

    let dir_new = Vector {x: ca_new, y: cb_new, z: cg_new};
    particle_1.dir.assign(&dir_new);
    particle_1.dir.normalize();
}

fn update_particle_energy(particle_1: &mut Particle, material: &Material, distance_traveled: f64, recoil_energy: f64, xi: f64, electronic_stopping_mode: i32) {

    //If particle energy  drops below zero before electronic stopping calcualtion, it produces NaNs
    particle_1.E = particle_1.E - recoil_energy;
    if particle_1.E < 0. {
        particle_1.E = 0.;
    }

    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let ck = material.electronic_stopping_correction_factor;

    if material.inside_energy_barrier(x, y) {

        let electronic_stopping_power = material.electronic_stopping_power(particle_1, electronic_stopping_mode);
        let n = material.number_density(x, y);

        let delta_energy = match electronic_stopping_mode {
            INTERPOLATED => electronic_stopping_power*n*distance_traveled*ck,
            LOW_ENERGY_NONLOCAL => electronic_stopping_power*n*distance_traveled*ck,
            LOW_ENERGY_LOCAL => {
                let Za: f64  = particle_1.Z;
                let Zb: f64 = material.Z_eff(x, y);

                //Oen-Robinson local electronic stopping power
                let a = screening_length(Za, Zb);
                //d1 is the first interior constant of the screening function
                //let d1 = 0.3; //Moliere
                let d1 = 0.278544; //Kr-C

                d1*d1/2./PI*electronic_stopping_power*(-d1*xi).exp()/a/a*ck
                },
            LOW_ENERGY_EQUIPARTITION => {
                let Za: f64  = particle_1.Z;
                let Zb: f64 = material.Z_eff(particle_1.pos.x, particle_1.pos.y);

                //Oen-Robinson local electronic stopping power
                let a = screening_length(Za, Zb);
                //d1 is the first interior constant of the screening function
                let d1 = 0.278544; //Moliere
                let delta_energy_local = d1*d1/2./PI*electronic_stopping_power*(-d1*xi).exp()/a/a;
                //println!("S_local: {}, S_low: {}, xi: {}", S_local/Q, S_low/Q, xi*a/ANGSTROM);
                let delta_energy_nonlocal = electronic_stopping_power*n*distance_traveled;

                (0.5*delta_energy_local + 0.5*delta_energy_nonlocal)*ck
            },
        };

        particle_1.E += -delta_energy;
    }

    //Make sure particle energy doesn't become negative again
    if particle_1.E < 0. {
        particle_1.E = 0.;
    }
}

fn particle_advance(particle_1: &mut Particle, mfp: f64, asympototic_deflection: f64) -> f64 {

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

fn boundary_condition_2D_planar(particle_1: &mut Particle, material: &Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let cosx = particle_1.dir.x;
    let cosy = particle_1.dir.y;
    let cosz = particle_1.dir.z;
    let E = particle_1.E;
    let Es = particle_1.Es;
    let Ec = particle_1.Ec;

    if !material.inside_energy_barrier(x, y) {
        if let Closest::SinglePoint(p2) = material.closest_point(x, y) {
            let dx = p2.x() - x;
            let dy = p2.y() - y;
            let mag = (dx*dx + dy*dy).sqrt();

            let costheta = dx*cosx/mag + dy*cosy/mag;
            let leaving_energy = E*costheta*costheta;

            if (costheta < 0.) {
                if leaving_energy > Es {

                    particle_1.left = true;
                    particle_1.E += -Es;

                    //Surface refraction via Snell's law
                    let sintheta0 = (1. - costheta*costheta).sqrt();
                    let sintheta1 = sintheta0*(E/(E - Es)).sqrt();
                    let delta_theta = sintheta1.asin() - sintheta0.asin();

                    rotate_particle(particle_1, delta_theta, 0.);
                    particle_1.add_trajectory();

                } else {
                    //Specular reflection at local surface normal
                    particle_1.dir.x = -2.*(costheta)*dx/mag + cosx;
                    particle_1.dir.y = -2.*(costheta)*dy/mag + cosy;

                    particle_1.add_trajectory();
                }
            }
        } else {
            panic!("Surface boundary algorithm encountered an error. Check geometry.");
        }
    }

    if (E < Ec) & !particle_1.left {
        particle_1.stopped = true;
        particle_1.add_trajectory();
    }
}

fn boundary_condition_1D_planar(particle_1: &mut Particle, material: &Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let cosx = particle_1.dir.x;
    let cosy = particle_1.dir.y;
    let cosz = particle_1.dir.z;
    let E = particle_1.E;
    let xc = -2.*material.number_density(x, y).powf(-1./3.)/SQRT2PI;
    let Es = particle_1.Es;
    let Ec = particle_1.Ec;

    if (x < xc) & (cosx < 0.) {
        let leaving_energy = E*cosx*cosx;

        if leaving_energy > Es {
            particle_1.left = true;
            particle_1.E += -Es;

            let cosx_new = ((E*cosx*cosx - Es)/(E - Es)).sqrt();
            let sinx = (1. - cosx*cosx).sqrt();
            let sinx_new = (1. - cosx_new*cosx_new).sqrt();
            particle_1.dir.x = -cosx_new;
            particle_1.dir.y *= sinx_new/sinx;
            particle_1.dir.z *= sinx_new/sinx;
            particle_1.dir.normalize();
            particle_1.add_trajectory();

        } else {
            particle_1.dir.x = cosx.abs();
        }
    }

    if (E < Ec) & !particle_1.left {
        particle_1.stopped = true;
        particle_1.add_trajectory();
    }
}

fn bca_from_file() {
    //Read input file, convert to string, and open with toml
    let mut input_toml = String::new();
    let mut file = OpenOptions::new()
        .read(true)
        .write(false)
        .create(false)
        .open("input.toml")
        .expect("Could not open input file.");
    file.read_to_string(&mut input_toml).unwrap();
    let input: Input = toml::from_str(&input_toml).unwrap();

    //Unpack toml information into structs
    let material = Material::new(input.material_parameters, input.geometry);
    let options = input.options;
    let particle_parameters = input.particle_parameters;

    //Check that incompatible options are not on simultaneously
    assert!(options.high_energy_free_flight_paths == (options.electronic_stopping_mode == INTERPOLATED),
        "High energy free flight paths used with low energy stoppping power.");

    if options.electronic_stopping_mode == INTERPOLATED {
        assert!(options.weak_collision_order == 0,
            "Cannot use weak collision loop with free flight paths.");
        //assert!(options.mean_free_path_model == LIQUID,
        //    "Gaseous model not currently implemented for high energy free flight paths.");
    }
    if options.mean_free_path_model == GASEOUS {
        assert!(options.weak_collision_order == 0,
            "Cannot use weak collisions with gaseous mean free path model.");
    }

    //Check that particle arrays are equal length
    assert_eq!(particle_parameters.Z.len(), particle_parameters.m.len(),
        "Particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.E.len(),
        "Particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.pos.len(),
        "Particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.dir.len(),
        "Particle input arrays of unequal length.");
    let N = particle_parameters.Z.len();

    //Determine the length, energy, and mass units for particle input
    let length_unit: f64 = match particle_parameters.length_unit.as_str() {
        "MICRON" => MICRON,
        "CM" => CM,
        "ANGSTROM" => ANGSTROM,
        "NM" => NM,
        "M" => 1.,
        _ => panic!("Unknown unit {} in input file. Choose one of: MICRON, CM, ANGSTROM, NM, M",
            particle_parameters.length_unit.as_str())
    };
    let energy_unit: f64 = match particle_parameters.energy_unit.as_str() {
        "EV" => EV,
        "J"  => 1.,
        "KEV" => EV*1E3,
        "MEV" => EV*1E6,
        _ => panic!("Unknown unit {} in input file. Choose one of: EV, J, KEV, MEV", particle_parameters.energy_unit.as_str())
    };
    let mass_unit: f64 = match particle_parameters.mass_unit.as_str() {
        "AMU" => AMU,
        "KG" => 1.0,
        _ => panic!("Unknown unit {} in input file. Choose one of: AMU, KG", particle_parameters.mass_unit.as_str())
    };

    //Estimate maximum number of recoils produced per ion
    let mut max_energy: f64 = 0.;
    let mut total_particles: usize = 0;
    for particle_index in 0..N {
        let E = particle_parameters.E[particle_index];
        let N_ = particle_parameters.N[particle_index];
        if E > max_energy {
            max_energy = E*energy_unit;
        }
        total_particles += N_;
    }

    //Create particle vector from input file
    let estimated_num_particles: usize = match options.track_recoils {
        true => total_particles + ((max_energy/material.Ec).ceil() as usize),
        false => total_particles,
    };

    let mut particles: Vec<Particle> = Vec::with_capacity(estimated_num_particles);
    for particle_index in 0..N {

        let N_ = particle_parameters.N[particle_index];
        let m = particle_parameters.m[particle_index];
        let Z = particle_parameters.Z[particle_index];
        let E = particle_parameters.E[particle_index];
        let Ec = particle_parameters.Ec[particle_index];
        let Es = particle_parameters.Es[particle_index];
        let (x, y, z) = particle_parameters.pos[particle_index];
        let (cosx, cosy, cosz) = particle_parameters.dir[particle_index];
        for sub_particle_index in 0..N_ {

            //Surface refraction
            let Es = particle_parameters.Es[particle_index];
            let E_new = E*energy_unit + Es*energy_unit;
            let cosx_new = ((E*cosx*cosx + Es)/(E + Es)).sqrt();
            let sinx = (1. - cosx*cosx).sqrt();
            let sinx_new = (1. - cosx_new*cosx_new).sqrt();
            let cosy_new = cosy*sinx_new/sinx;
            let cosz_new = cosz*sinx_new/sinx;

            //Add new particle to particle vector
            particles.push(Particle::new(
                m*mass_unit, Z, E_new, Ec*energy_unit, Es*energy_unit,
                x*length_unit, y*length_unit, z*length_unit,
                cosx_new, cosy_new, cosz_new, true, options.track_trajectories
            ));
        }
    }

    let mut num_sputtered: usize = 0;
    let mut num_deposited: usize = 0;
    let mut num_reflected: usize = 0;
    let mut energy_sputtered: f64 = 0.;
    let mut energy_reflected: f64 = 0.;
    let mut range: f64 = 0.;

    //Open output files for streaming output
    let mut reflected_file = OpenOptions::new()
        .write(true).
        create(true).
        open(format!("{}{}", options.name, "reflected.output")).
        unwrap();
    let mut reflected_file_stream = BufWriter::with_capacity(options.stream_size, reflected_file);

    let mut sputtered_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "sputtered.output"))
        .unwrap();
    let mut sputtered_file_stream = BufWriter::with_capacity(options.stream_size, sputtered_file);

    let mut deposited_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "deposited.output"))
        .unwrap();
    let mut deposited_file_stream = BufWriter::with_capacity(options.stream_size, deposited_file);

    let mut trajectory_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "trajectories.output"))
        .unwrap();
    let mut trajectory_file_stream = BufWriter::with_capacity(options.stream_size, trajectory_file);

    let mut trajectory_data = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "trajectory_data.output"))
        .unwrap();
    let mut trajectory_data_stream = BufWriter::with_capacity(options.stream_size, trajectory_data);

    //Main loop
    let mut particle_index: usize = particles.len();
    'particle_loop: while particle_index > 0 {
        //Remove particle from top of vector
        let mut particle_1 = particles.pop().unwrap();

        //if particle_1.track_trajectories {
        //    particle_1.add_trajectory();
        //}

        //Print to stdout
        if options.print & particle_1.incident & (particle_index % (total_particles / options.print_num) == 0){
            println!("Incident Ion {} of {}", particle_index, total_particles);
        }

        'trajectory_loop: while !particle_1.stopped & !particle_1.left {

            //BCA loop
            //Choose recoil partner location and species
            let (phis_azimuthal, impact_parameters, mfp) = determine_mfp_phi_impact_parameter(
                    &mut particle_1, &material, options.weak_collision_order,
                    options.high_energy_free_flight_paths, options.mean_free_path_model);

            //let mut total_deflection_angle = 0.;
            let mut total_energy_loss = 0.;
            let mut total_asymptotic_deflection = 0.;
            let mut distance_of_closest_approach = 0.;

            if particle_1.track_trajectories {
                particle_1.add_trajectory();
            }

            'collision_loop: for k in 0..options.weak_collision_order + 1 {

                let (Z_recoil, M_recoil, Ec_recoil, Es_recoil, xr, yr, zr, cxr, cyr, czr) = choose_collision_partner(&mut particle_1, &material, phis_azimuthal[k], impact_parameters[k], mfp);
                //If recoil location is inside, proceed with binary collision loop
                if material.inside(xr, yr) {

                    //Generate new particle at recoil position
                    let mut particle_2 = Particle::new(
                        M_recoil, Z_recoil, 0., Ec_recoil, Es_recoil,
                        xr, yr, zr,
                        cxr, cyr, czr,
                        false, options.track_recoil_trajectories
                    );

                    //Determine scattering angle from binary collision
                    let (theta, psi, psi_recoil, recoil_energy, asymptotic_deflection, xi) = calculate_binary_collision(&particle_1, &particle_2, impact_parameters[k], 100, 1E-3, options.interaction_potential);

                    //Only use 0th order collision for local stopping
                    if k == 0 {
                        distance_of_closest_approach = xi;
                    }

                    //Energy transfer to recoil
                    particle_2.E = recoil_energy - material.Eb;

                    //Accumulate asymptotic deflections for primary particle
                    total_energy_loss += recoil_energy;

                    //total_deflection_angle += psi;
                    //total_asymptotic_deflection += asymptotic_deflection;
                    total_asymptotic_deflection = 0.;

                    //Rotate particle 1, 2 by lab frame scattering angles
                    rotate_particle(&mut particle_1, psi, phis_azimuthal[k]);
                    rotate_particle(&mut particle_2, -psi_recoil, phis_azimuthal[k]);
                    if psi > 0. {
                        particle_1.number_collision_events += 1;
                    }

                    //Deep recoil suppression
                    //See Eckstein 1991 7.5.3 for recoil suppression function
                    if options.track_recoils & options.suppress_deep_recoils {
                        let E = particle_1.E;
                        let Za = particle_1.Z;
                        let Zb = particle_2.Z;

                        let Ma = particle_1.m;
                        let Mb = particle_2.m;

                        let n = material.number_density(xr, yr);
                        let a: f64 = screening_length(Za, Zb);
                        //let reduced_energy: f64 = 4.*PI*EPS0*a*Mb*E/(Ma + Mb)/Za/Zb/Q/Q;
                        let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E;
                        let estimated_range_of_recoils = (reduced_energy.powf(0.3) + 0.1).powf(3.)/n/a/a;

                        if let Closest::SinglePoint(p2) = material.closest_point(xr, yr) {
                            let dx = p2.x() - xr;
                            let dy = p2.y() - yr;
                            let distance_to_surface = (dx*dx + dy*dy).sqrt();

                            if (distance_to_surface < estimated_range_of_recoils) & (recoil_energy > particle_2.Ec) {
                                particle_2.add_trajectory();
                                particles.push(particle_2);
                            }
                        }
                    //If transferred energy > cutoff energy, add recoil to particle vector
                    } else if options.track_recoils & (recoil_energy > particle_2.Ec) {
                        particles.push(particle_2);
                    }
                }
            }

            //Advance particle in space and track total distance traveled
            let distance_traveled = particle_advance(&mut particle_1, mfp, total_asymptotic_deflection);

            if particle_1.track_trajectories {
                particle_1.add_trajectory();
            }

            //Subtract total energy from all simultaneous collisions and electronic stopping
            update_particle_energy(&mut particle_1, &material, distance_traveled, total_energy_loss, distance_of_closest_approach, options.electronic_stopping_mode);

            //Check boundary conditions on leaving and stopping
            boundary_condition_2D_planar(&mut particle_1, &material);

            //Set particle index to topmost particle
            particle_index = particles.len();
        }

        //Once particle finishes, begin data output
        //Stream current particle output to files

        //Incident particle, left material: reflected
        if particle_1.incident & particle_1.left {
            writeln!(
                reflected_file_stream, "{},{},{},{},{},{},{},{},{},{}",
                particle_1.m/mass_unit, particle_1.Z, particle_1.E/energy_unit,
                particle_1.pos.x/length_unit, particle_1.pos.y/length_unit, particle_1.pos.z/length_unit,
                particle_1.dir.x, particle_1.dir.y, particle_1.dir.z,
                particle_1.number_collision_events
            ).expect("Could not write to reflected.output.");

        }

        //Incident particle, stopped in material: deposited
        if particle_1.incident & particle_1.stopped {
            writeln!(
                deposited_file_stream, "{},{},{},{},{},{}",
                particle_1.m/mass_unit, particle_1.Z,
                particle_1.pos.x/length_unit, particle_1.pos.y/length_unit, particle_1.pos.z/length_unit,
                particle_1.number_collision_events
            ).expect("Could not write to deposited.output.");
        }

        //Not an incident particle, left material: sputtered
        if !particle_1.incident & particle_1.left {
            writeln!(
                sputtered_file_stream, "{},{},{},{},{},{},{},{},{},{},{},{},{}",
                particle_1.m/mass_unit, particle_1.Z, particle_1.E/energy_unit,
                particle_1.pos.x/length_unit, particle_1.pos.y/length_unit, particle_1.pos.z/length_unit,
                particle_1.dir.x, particle_1.dir.y, particle_1.dir.z,
                particle_1.number_collision_events,
                particle_1.pos_origin.x, particle_1.pos_origin.y, particle_1.pos_origin.z
            ).expect("Could not write to sputtered.output.");
        }

        //Trajectory output
        if particle_1.track_trajectories {
            writeln!(trajectory_data_stream, "{}", particle_1.trajectory.len())
                .expect("Could not write trajectory length data.");

            for pos in particle_1.trajectory {
                writeln!(
                    trajectory_file_stream, "{},{},{},{},{},{}",
                    particle_1.m/mass_unit, particle_1.Z, pos.E/energy_unit,
                    pos.x/length_unit, pos.y/length_unit, pos.z/length_unit,
            ).expect("Could not write to trajectories.output.");
            }
        }
    }
    //Flush all file streams before dropping to ensure all data is written
    reflected_file_stream.flush().unwrap();
    deposited_file_stream.flush().unwrap();
    sputtered_file_stream.flush().unwrap();
    trajectory_data_stream.flush().unwrap();
    trajectory_file_stream.flush().unwrap();
}

fn main() {
    bca_from_file();
}
