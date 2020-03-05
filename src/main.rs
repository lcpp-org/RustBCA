extern crate ndarray;
extern crate geo;
extern crate csv;
extern crate serde;
extern crate text_io;
extern crate scan_fmt;
extern crate toml;

pub use crate::geo::*;
use ndarray::prelude::*;
use geo::algorithm::contains::Contains;
use geo::algorithm::closest_point::ClosestPoint;
use geo::algorithm::bounding_rect::BoundingRect;
use geo::extremes::ExtremePoints;
//use rand::Rng;
//use std::error::Error;
//use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::{Read, BufReader};
use std::f64::consts::PI;
use text_io::{read, scan};
use toml::{Value};
use serde::*;

const Q: f64 = 1.602E-19;
const EV: f64 = Q;
const AMU: f64 = 1.66E-27;
const ANGSTROM: f64 = 1E-10;
const MICRON: f64 = 1E-6;
const NM: f64 = 1E-9;
const CM: f64 = 1E-2;
const EPS0: f64 = 8.85E-12;
const A0: f64 = 0.52918E-10;
const K: f64 = 1.11265E-10;
const ME: f64 = 9.11E-31;
const SQRTPI: f64 = 1.772453850906;
const SQRT2PI: f64 = 2.506628274631;
const C: f64 = 299792000.;
const BETHE_BLOCH_PREFACTOR: f64 = 4.*PI*(Q*Q/(4.*PI*EPS0))*(Q*Q/(4.*PI*EPS0))/ME/C/C;
const LINDHARD_SCHARFF_PREFACTOR: f64 = 1.212*ANGSTROM*ANGSTROM*Q;

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
    write_files: bool
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
    m: f64
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
    pos: Vector,
    dir: Vector,
    pos_old: Vector,
    pos_origin: Vector,
    t: f64,
    stopped: bool,
    left: bool,
    incident: bool,
    first_step: bool,
    trajectory: Vec<Vector4>,
    track_trajectories: bool,
    backreflected: bool
}
impl Particle {
    pub fn new(m: f64, Z: f64, E: f64, x: f64, y: f64, z: f64, dirx: f64, diry: f64, dirz: f64, incident: bool, track_trajectories: bool) -> Particle {
        Particle {
            m: m,
            Z: Z,
            E: E,
            pos: Vector::new(x, y, z),
            dir: Vector::new(dirx, diry, dirz),
            pos_old: Vector::new(x, y, z),
            pos_origin: Vector::new(x, y, z),
            t: 0.,
            stopped: false,
            left: false,
            incident: incident,
            first_step: incident,
            trajectory: vec![Vector4::new(E, x, y, z)],
            track_trajectories: track_trajectories,
            backreflected: false
        }
    }
    fn add_trajectory(&mut self) {
        self.trajectory.push(Vector4 {E: self.E, x: self.pos.x, y: self.pos.y, z: self.pos.z});
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
}
impl Material {
    pub fn new(material_parameters: MaterialParameters, geometry: Geometry) -> Material {

        //Multiply all coordinates by value of geometry unit.
        let length_unit: f64 = match geometry.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            _ => panic!("Incorrect unit {} in input file.", geometry.length_unit.as_str())
        };
        let energy_unit: f64 = match material_parameters.energy_unit.as_str() {
            "EV" => EV,
            "J"  => 1.,
            "KEV" => EV*1E3,
            "MEV" => EV*1E6,
            _ => panic!("Incorrect unit {} in input file.", material_parameters.energy_unit.as_str())
        };
        let mass_unit: f64 = match material_parameters.mass_unit.as_str() {
            "AMU" => AMU,
            "KG" => 1.0,
            _ => panic!("Incorrect unit {} in input file.", material_parameters.mass_unit.as_str())
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
            surface: Polygon::new(LineString::from(unit_coords), vec![]),
            energy_surface: Polygon::new(LineString::from(unit_e_coords), vec![]),
            simulation_surface: Polygon::new(LineString::from(unit_b_coords), vec![])
        }
    }

    fn inside(&self, x: f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.surface.contains(&p);
    }
    fn mfp(&self, x: f64, y: f64) -> f64 {
        return self.n.powf(-1_f64/3_f64);
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
}

fn phi(xi: f64) -> f64 {
    //Kr-C potential
    return 0.190945*(-0.278544*xi).exp() + 0.473674*(-0.637174*xi).exp() + 0.335381*(-1.919249*xi).exp();
}

fn dphi(xi: f64) -> f64 {
    //First differential of Kr-C potential
    return -0.278544*0.190945*(-0.278544*xi).exp() - 0.637174*0.473674*(-0.637174*xi).exp() - 0.335381*1.919249*(-1.919249*xi).exp();
}

fn screening_length(Za: f64, Zb: f64) -> f64 {
    //Lindhard screening length
    return 0.8853*A0/(Za.sqrt() + Zb.sqrt()).powf(2./3.);
}

fn doca_function(x0: f64, beta: f64, reduced_energy: f64) -> f64 {
    //Transcendental function to determine distance of closest approach
    return x0 - phi(x0)/reduced_energy - beta*beta/x0;
}

fn diff_doca_function(x0: f64, beta: f64, reduced_energy: f64) -> f64 {
    //First differential of distance of closest approach function
    return beta*beta/x0/x0 - dphi(x0)/reduced_energy + 1.
}

fn f(x: f64, beta: f64, reduced_energy: f64) -> f64 {
    //Function for scattering integral - see Mendenhall and Weller, 1991 & 2005
    return (1_f64 - phi(x)/x/reduced_energy - beta*beta/x/x).powf(-0.5);
}

fn rotate_around_axis_by_pi(vector: &Vector, axis: &Vector) -> Vector {
    //Rotate a vector around an axis by pi (180 degrees) - used for reflection at surface
    let axis_mag = axis.magnitude();
    let ux = axis.x/axis_mag;
    let uy = axis.y/axis_mag;
    let uz = axis.z/axis_mag;

    let ndarray_vector = array![vector.x, vector.y, vector.z];

    let R = array![
        [2.*ux*ux - 1., 2.*ux*uy, 2.*ux*uz],
        [2.*uy*ux, 2.*uy*uy - 1., 2.*uy*uz],
        [2.*uz*ux, 2.*uz*uy, 2.*uz*uz - 1.]
    ];

    let temp_result =  R.dot(&ndarray_vector);
    return Vector {x: temp_result[0], y: temp_result[1], z: temp_result[2]};
}

fn binary_collision(particle_1: &Particle, particle_2: &Particle,  material: &Material, impact_parameter: f64, tol: f64, max_iter: i32) -> (f64, f64, f64, f64) {
    if !material.inside(particle_2.pos.x, particle_2.pos.y) {
        //If the chosen collision partner is not within the target, skip collision algorithm
        return (0., 0., 0., 0.);
    } else {
        let Za: f64 = particle_1.Z;
        let Zb: f64 = particle_2.Z;
        let Ma: f64 = particle_1.m;
        let Mb: f64 = particle_2.m;
        let E0: f64 = particle_1.E;
        let mu: f64 = Mb/(Ma + Mb);

        let a: f64 = screening_length(Za, Zb);
        let reduced_energy: f64 = K*a*mu/(Za*Zb*Q*Q)*E0;
        let beta: f64 = impact_parameter/a;

        //Guess for large reduced energy from Mendenhall and Weller 1991
        //For small energies, use pure Newton-Raphson with arbitrary guess of 1
        let mut x0 = 1_f64;
        let mut xn: f64;
        if reduced_energy > 5_f64 {
            let inv_er_2 = 0.5_f64/reduced_energy;
            x0 = inv_er_2 + (inv_er_2.powf(2_f64) + beta.powf(2_f64)).sqrt();
        }

        //Newton-Raphson to determine distance of closest approach
        let mut err: f64;
        for _ in 0..max_iter {
            xn = x0 - doca_function(x0, beta, reduced_energy)/diff_doca_function(x0, beta, reduced_energy);
            err = (xn - x0).abs()/xn;
            x0 = xn;
            if err < tol {
                break;
            }
        }

        //Scattering integral quadrature from Mendenhall and Weller 2005
        let lambda_0 = (0.5 + beta*beta/x0/x0/2. - dphi(x0)/2_f64/reduced_energy).powf(-0.5);
        let alpha = 1_f64/12_f64*(1_f64 + lambda_0 + 5_f64*(0.4206_f64*f(x0/0.9072_f64, beta, reduced_energy) + 0.9072_f64*f(x0/0.4206, beta, reduced_energy)));
        let theta = PI*(1_f64 - beta*alpha/x0);

        //See Eckstein 1991 for details on center of mass and lab frame angles here
        let t = x0*a*(theta/2.).sin();
        let psi = (theta.sin()).atan2(Ma/Mb + theta.cos());
        let T = 4.*(Ma*Mb)/(Ma + Mb).powf(2.)*E0*((theta/2.).sin()).powf(2.);

        return (theta, psi, T, t);
    }
}

fn update_coordinates(particle_1: &mut Particle, particle_2: &mut Particle, material: &Material, phi_azimuthal: f64, theta: f64, psi: f64, T: f64, t: f64) {
    let mut mfp = material.mfp(particle_1.pos.x, particle_1.pos.y);

    //Store previous position - this is used for boundary checking
    particle_1.pos_old.assign(&particle_1.pos);

    // "atomically rough surface" - see Moeller and Eckstein 1988
    //if particle_1.first_step {
    //    mfp *= rand::random::<f64>();
    //    particle_1.first_step = false;
    //}

    //path length - must subtract next asymptotic deflection and add previous to obtain correct trajectory
    let ffp: f64 = mfp - t + particle_1.t;

    particle_1.pos.x += ffp*particle_1.dir.x;
    particle_1.pos.y += ffp*particle_1.dir.y;
    particle_1.pos.z += ffp*particle_1.dir.z;

    particle_1.t = t;

    let ca: f64 = particle_1.dir.x;
    let cb: f64 = particle_1.dir.y;
    let cg: f64 = particle_1.dir.z;
    let cphi: f64 = phi_azimuthal.cos();
    let sphi: f64 = phi_azimuthal.sin();
    let sa = (1_f64 - ca.powf(2_f64)).sqrt();

    //Update particle 1 direction
    //Particle direction update from TRIDYN, see Moeller and Eckstein 1988
    {
        let cpsi: f64 = psi.cos();
        let spsi: f64 = psi.sin();
        let ca_new: f64 = cpsi*ca + spsi*cphi*sa;
        let cb_new: f64 = cpsi*cb - spsi/sa*(cphi*ca*cb - sphi*cg);
        let cg_new: f64 = cpsi*cg - spsi/sa*(cphi*ca*cg + sphi*cb);

        let dir_new = Vector {x: ca_new, y: cb_new, z: cg_new};
        particle_1.dir.assign(&dir_new);
        particle_1.dir.normalize();
    }

    //Update particle 2 direction
    {
        let psi_b: f64 = (-(theta.sin())).atan2(1. - theta.cos());
        //println!("{}", psi_b);
        let cpsi_b: f64 = psi_b.cos();
        let spsi_b: f64 = psi_b.sin();
        let ca_new: f64 = cpsi_b*ca + spsi_b*cphi*sa;
        let cb_new: f64 = cpsi_b*cb - spsi_b/sa*(cphi*ca*cb - sphi*cg);
        let cg_new: f64 = cpsi_b*cg - spsi_b/sa*(cphi*ca*cg + sphi*cb);

        let dir_new = Vector {x: ca_new, y: cb_new, z: cg_new};
        particle_2.dir.assign(&dir_new);
        particle_2.dir.normalize();
    }

    //Update particle energies
    {
        let Za = particle_1.Z;
        let Ma = particle_1.m;
        let Zb = material.Z_eff(particle_1.pos.x, particle_1.pos.y);
        let E = particle_1.E;

        let mut stopping_power = 0.;

        //Electronic stopping only inside material
        if material.inside(particle_1.pos_old.x, particle_1.pos_old.y) {
            //Bethe-Bloch stopping, as modified by Biersack and Haggmark (1980) to fit experiment
            let n = material.number_density(particle_1.pos.x, particle_1.pos.y);
            let v = (2.*E/Ma).sqrt();
            let beta = v/C;

            let mut I0; //Mean excitation potential approximation form Biersack and Haggmark
            if Zb < 13. {
                I0 = 12. + 7./Zb;
            } else {
                I0 = 9.76 + 58.5*Zb.powf(-1.19);
            }
            let I = Zb*I0*Q;

            let mut B; //Effective empirical shell correction from Biersack and Haggmark
            if Zb < 3. {
                B = 100.*Za/Zb;
            } else {
                B = 5.;
            }

            let prefactor = BETHE_BLOCH_PREFACTOR*Zb*Za*Za/beta/beta;
            //let prefactor = 8.0735880E-42*Zb*Za*Za/beta/beta;
            let eb = 2.*ME*v*v/I;
            let S_BB = prefactor*(eb + 1. + B/eb).ln()*n; //Bethe-Bloch stopping, as modified by Biersack and Haggmark (1980) to fit experiment

            //Lindhard-Scharff stopping for low energies
            let S_LS = LINDHARD_SCHARFF_PREFACTOR*(Za.powf(7./6.)*Zb)/(Za.powf(2./3.) + Zb.powf(2./3.)).powf(3./2.)*(E/Ma*AMU/Q).sqrt()*n;

            //Biersack-Varelas stopping interpolation scheme
            stopping_power = 1./(1./S_BB + 1./S_LS);

        }
        let Enl = ffp*stopping_power;

        //Ensure that energies do not go negative
        particle_1.E = E - T - Enl;
        if particle_1.E < 0_f64 {
            particle_1.E = 0_f64;
        }

        //Ensure that energies do not go negative
        particle_2.E = T - material.Eb;
        if particle_2.E < 0_f64 {
            particle_2.E = 0_f64;
        }
    }
}

fn pick_collision_partner(particle_1: &Particle, material: &Material) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    //Liquid model of amorphous material - see Eckstein 1991
    let mfp: f64 = material.mfp(particle_1.pos.x, particle_1.pos.y);

    let pmax : f64 = mfp/SQRTPI;
    let impact_parameter: f64 = pmax*(rand::random::<f64>()).sqrt();
    let phi_azimuthal: f64 = 2.*PI*rand::random::<f64>();

    let sphi: f64 = phi_azimuthal.sin();
    let ca: f64 = particle_1.dir.x;
    let cb: f64 = particle_1.dir.y;
    let cg: f64 = particle_1.dir.z;
    let sa: f64 = (1_f64 - ca*ca).sqrt();
    let cphi: f64 = phi_azimuthal.cos();

    //Collision partner at locus of next collision, displaced by chosen impact parameter and rotated by phi azimuthal
    let x_recoil: f64 = particle_1.pos.x + mfp*ca - impact_parameter*cphi*sa;
    let y_recoil: f64 = particle_1.pos.y + mfp*cb - impact_parameter*(sphi*cg - cphi*cb*ca)/sa;
    let z_recoil: f64 = particle_1.pos.z + mfp*cg + impact_parameter*(sphi*cb - cphi*ca*cg)/sa;

    return (impact_parameter, phi_azimuthal, x_recoil, y_recoil, z_recoil, ca, cb, cg);
}

fn surface_boundary_condition(particle_1: &mut Particle, material: &Material) {
    let leaving: bool = !material.inside_energy_barrier(particle_1.pos.x, particle_1.pos.y) & material.inside_energy_barrier(particle_1.pos_old.x, particle_1.pos_old.y);
    let entering: bool = material.inside_energy_barrier(particle_1.pos.x, particle_1.pos.y) & !material.inside_energy_barrier(particle_1.pos_old.x, particle_1.pos_old.y);

    if leaving | entering {
        let cosx = particle_1.dir.x;
        let cosy = particle_1.dir.y;
        let cosz = particle_1.dir.z;

        if let Closest::SinglePoint(p2) = material.closest_point(particle_1.pos.x, particle_1.pos.y) {
            let dx = particle_1.pos.x - p2.x();
            let dy = particle_1.pos.y - p2.y();
            let mag = (dx*dx + dy*dy).sqrt();

            let E = particle_1.E;
            let M = particle_1.m;
            let Es = material.Es;

            let costheta = dx*cosx/mag + dy*cosy/mag;

            //println!("dx: {} dy: {} costheta: {}", dx, dy, costheta);

            if leaving & (E*costheta*costheta < Es) {
                //reflect back onto surface
                println!("backreflected");

                let dot_product = dx/mag*cosx + dy/mag*cosy;
                particle_1.dir.x = -2.*dot_product*dx/mag + cosx;
                particle_1.dir.y = -2.*dot_product*dy/mag + cosy;

                particle_1.backreflected = true;

            } else if !particle_1.backreflected {
                //Refract through surface
                //println!("{} {} {} {} {}", particle_1.E/Q, particle_1.dir.x, particle_1.dir.y, particle_1.dir.z, particle_1.dir.magnitude());
                let v = (2.*E/M).sqrt();
                let delta_v = costheta.signum()*(2./M).sqrt()*((E + costheta.signum()*Es).sqrt() - E.sqrt());

                let vx = v*cosx;
                let vy = v*cosy;
                let vz = v*cosz;

                let vx_new = (cosx).signum()*(vx*vx + ((dx*cosx).signum()*2.*v*delta_v + delta_v*delta_v)*dx*dx/mag/mag).sqrt();
                let vy_new = (cosy).signum()*(vy*vy + ((dy*cosy).signum()*2.*v*delta_v + delta_v*delta_v)*dy*dy/mag/mag).sqrt();
                let vz_new = vz;
                let v_new = (vx_new*vx_new + vy_new*vy_new + vz_new*vz_new).sqrt();

                particle_1.dir.x = vx_new/v_new;
                particle_1.dir.y = vy_new/v_new;
                particle_1.dir.z = vz_new/v_new;
                particle_1.E = v_new*v_new/2.*M;
                //println!("{} {} {} {} {}", particle_1.E/Q, particle_1.dir.x, particle_1.dir.y, particle_1.dir.z, particle_1.dir.magnitude());
            } else {
                //println!("backreflection flag reset");
                particle_1.backreflected = false;
            }
        }
    }
}

fn bca_input() {

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

    //Check that particle arrays are equal length
    assert_eq!(particle_parameters.Z.len(), particle_parameters.m.len());
    assert_eq!(particle_parameters.Z.len(), particle_parameters.E.len());
    assert_eq!(particle_parameters.Z.len(), particle_parameters.pos.len());
    assert_eq!(particle_parameters.Z.len(), particle_parameters.dir.len());
    let N = particle_parameters.Z.len();

    //Determine the length, energy, and mass units for particle input
    let length_unit: f64 = match particle_parameters.length_unit.as_str() {
        "MICRON" => MICRON,
        "CM" => CM,
        "ANGSTROM" => ANGSTROM,
        "NM" => NM,
        _ => panic!("Incorrect unit {} in input file.", particle_parameters.length_unit.as_str())
    };
    let energy_unit: f64 = match particle_parameters.energy_unit.as_str() {
        "EV" => EV,
        "J"  => 1.,
        "KEV" => EV*1E3,
        "MEV" => EV*1E6,
        _ => panic!("Incorrect unit {} in input file.", particle_parameters.energy_unit.as_str())
    };
    let mass_unit: f64 = match particle_parameters.mass_unit.as_str() {
        "AMU" => AMU,
        "KG" => 1.0,
        _ => panic!("Incorrect unit {} in input file.", particle_parameters.mass_unit.as_str())
    };

    //Create particle array
    let mut particles: Vec<Particle> = Vec::new();
    for particle_index in 0..N {
        let N_ = particle_parameters.N[particle_index];
        let m = particle_parameters.m[particle_index];
        let Z = particle_parameters.Z[particle_index];
        let E = particle_parameters.E[particle_index];
        let (x, y, z) = particle_parameters.pos[particle_index];
        let (dirx, diry, dirz) = particle_parameters.dir[particle_index];
        for sub_particle_index in 0..N_ {
            particles.push(Particle::new(m*mass_unit, Z, E*energy_unit, x*length_unit, y*length_unit, z*length_unit, dirx, diry, dirz, true, options.track_trajectories));
        }
    }

    let mut num_sputtered: usize = 0;
    let mut num_deposited: usize = 0;
    let mut num_reflected: usize = 0;
    let mut energy_sputtered: f64 = 0.;
    let mut energy_reflected: f64 = 0.;
    let mut range: f64 = 0.;

    //Main BCA loop
    let mut particle_index: usize = 0;
    while particle_index < particles.len() {
        if particle_index % 100 == 0 {
            println!("particle {} of {}", particle_index, particles.len());
        }
        while !particles[particle_index].stopped & !particles[particle_index].left {
            //println!("{} {} {} {}", particles[particle_index].E, particles[particle_index].pos.x, particles[particle_index].pos.y, particles[particle_index].dir.x);
            //Check simulation boundary conditions, add to energy fluxes out
            if !material.inside_simulation_boundary(particles[particle_index].pos.x, particles[particle_index].pos.y) {
                particles[particle_index].left = true;
                if particles[particle_index].incident {
                    num_reflected += 1;
                    energy_reflected += particles[particle_index].E;
                } else {
                    num_sputtered += 1;
                    energy_sputtered += particles[particle_index].E;
                }
                //Skip BCA loop
                continue;
            }

            //Check stopping condition, if incident, contribute to average range
            if particles[particle_index].E < material.Ec {
                particles[particle_index].stopped = true;
                if particles[particle_index].incident {
                    num_deposited += 1;
                    range += particles[particle_index].pos.x;
                }
                //Skip BCA loop
                continue;
            }

            //BCA loop
            //Choose collision partner
            let (impact_parameter, phi_azimuthal, xr, yr, zr, dirx, diry, dirz) = pick_collision_partner(&particles[particle_index], &material);
            let mut particle_2 = Particle::new(material.m, material.Z, 0., xr, yr, zr, dirx, diry, dirz, false, options.track_recoil_trajectories);

            //Calculate scattering and update coordinates
            let (theta, psi, T, t) = binary_collision(&particles[particle_index], &particle_2, &material, impact_parameter, 1E-3, 100);
            update_coordinates(&mut particles[particle_index], &mut particle_2, &material, phi_azimuthal, theta, psi, T, t);

            //Reflection and refraction from surface energy barrier
            surface_boundary_condition(&mut particles[particle_index], &material);

            if particles[particle_index].track_trajectories{
                particles[particle_index].add_trajectory();
            }

            //Generate recoils if transferred kinetic energy larger than cutoff energy
            if (T > material.Ec) & options.track_recoils {
                particles.push(particle_2);
            }
        }
        particle_index += 1
    }

    //Write output files
    if options.write_files {
        //println!("Writing output files...");
        let mut reflected_file = OpenOptions::new()
            .write(true).
            create(true).
            open(format!("{}{}", options.name, "reflected.output")).
            unwrap();
        let mut sputtered_file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(format!("{}{}", options.name, "sputtered.output"))
            .unwrap();
        let mut deposited_file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(format!("{}{}", options.name, "deposited.output"))
            .unwrap();
        let mut trajectory_file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(format!("{}{}", options.name, "trajectories.output"))
            .unwrap();
        let mut trajectory_data = OpenOptions::new()
            .write(true)
            .create(true)
            .open(format!("{}{}", options.name, "trajectory_data.output"))
            .unwrap();

        for particle in particles {
            if particle.incident & particle.left {
                writeln!(reflected_file, "{},{},{},{},{},{},{},{},{}", particle.m/mass_unit, particle.Z, particle.E/energy_unit, particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit, particle.dir.x, particle.dir.y, particle.dir.z)
                    .expect("Could not write to reflected.output.");
            }

            if particle.incident & particle.stopped {
                writeln!(deposited_file, "{},{},{},{},{}", particle.m/mass_unit, particle.Z, particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit)
                    .expect("Could not write to deposited.output.");
            }

            if !particle.incident & particle.left {
                writeln!(sputtered_file, "{},{},{},{},{},{},{},{},{}", particle.m/mass_unit, particle.Z, particle.E/energy_unit, particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit, particle.dir.x, particle.dir.y, particle.dir.z)
                    .expect("Could not write to sputtered.output.");
            }

            if particle.track_trajectories {
                writeln!(trajectory_data, "{}", particle.trajectory.len())
                    .expect("Could not write trajectory length data.");
                for pos in particle.trajectory {
                    writeln!(trajectory_file, "{},{},{},{},{},{}", particle.m/mass_unit, particle.Z, pos.E/energy_unit, pos.x/length_unit, pos.y/length_unit, pos.z/length_unit)
                        .expect("Could not write to trajectories.output.");
                }
            }
        }
    }
}

fn main() {
    bca_input();
}
