extern crate geo;
extern crate serde;
extern crate toml;

pub use crate::geo::*;
use geo::algorithm::contains::Contains;
use geo::algorithm::closest_point::ClosestPoint;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::Read;
use std::f64::consts::PI;
use serde::*;
use std::io::BufWriter;

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
    write_files: bool,
    stream_size: usize,
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
            trajectory: Vec::new(),
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

    fn Es(&self, Z: f64) -> f64 {
        if Z == self.Z {
            return self.Es;
        }
        else {
            return 0.;
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
    return (1. - phi(x)/x/reduced_energy - beta*beta/x/x).powf(-0.5);
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
        let mut x0 = 1.;
        let mut xn: f64;
        if reduced_energy > 5. {
            let inv_er_2 = 0.5/reduced_energy;
            x0 = inv_er_2 + (inv_er_2.powf(2.) + beta.powf(2.)).sqrt();
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
        let lambda_0 = (0.5 + beta*beta/x0/x0/2. - dphi(x0)/2./reduced_energy).powf(-0.5);
        let alpha = 1./12.*(1. + lambda_0 + 5.*(0.4206*f(x0/0.9072, beta, reduced_energy) + 0.9072*f(x0/0.4206, beta, reduced_energy)));
        let theta = PI*(1. - beta*alpha/x0);

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
    let sa = (1. - ca.powf(2.)).sqrt();

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
        let psi_b: f64 = (-theta.sin()).atan2(1. - theta.cos());
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

    //Momentum test
    //let test_E = particle_1.E - T;
    //let v0 = (2.*particle_1.E/particle_1.m).sqrt();
    //let v1 = (2.*test_E/particle_1.m).sqrt();
    //let v2 = (2.*T/particle_2.m).sqrt();

    //let px0 = particle_1.m*v0*ca;
    //let px1 = particle_1.m*v1*particle_1.dir.x;
    //let px2 = particle_2.m*v2*particle_2.dir.x;

    //let py0 = particle_1.m*v0*cb;
    //let py1 = particle_1.m*v1*particle_1.dir.y;
    //let py2 = particle_2.m*v2*particle_2.dir.y;

    //let pz0 = particle_1.m*v0*cg;
    //let pz1 = particle_1.m*v1*particle_1.dir.z;
    //let pz2 = particle_2.m*v2*particle_2.dir.z;

    //println!("x: {} {} y: {} {} z: {} {}", px0, px1 + px2, py0, py1 + py2, pz0, pz1 + pz2);

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
        if particle_1.E < 0. {
            particle_1.E = 0.;
        }

        //Ensure that energies do not go negative
        particle_2.E = T - material.Eb;
        if particle_2.E < 0. {
            particle_2.E = 0.;
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
    let sa: f64 = (1. - ca*ca).sqrt();
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
            let Es = material.Es(particle_1.Z);

            let costheta = dx*cosx/mag + dy*cosy/mag;

            if leaving & (E*costheta*costheta < Es) {
                //reflect back onto surface
                //println!("backreflected");
                particle_1.dir.x = -2.*costheta*dx/mag + cosx;
                particle_1.dir.y = -2.*costheta*dy/mag + cosy;

                particle_1.backreflected = true;

            } else if !particle_1.backreflected {

                let sx = (cosx*dx).signum();
                let sy = (cosy*dy).signum();
                let sign = costheta.signum();
                let sign_x = (E*cosx - Es*dx/mag).signum();
                let sign_y = (E*cosy - Es*dy/mag).signum();
                let sign_z = (cosz).signum();

                particle_1.dir.x = sign_x*((E*cosx*cosx + Es*dx*dx/mag/mag)/(E + sign*Es)).sqrt();
                particle_1.dir.y = sign_y*((E*cosy*cosy + Es*dy*dy/mag/mag)/(E + sign*Es)).sqrt();
                particle_1.dir.z = sign_z*(E*cosz*cosz/(E + sign*Es)).sqrt();
                particle_1.E = E + sign*Es;

                particle_1.dir.normalize();

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

    //Estimate maximum number of recoils produced
    let mut total_energy: f64 = 0.;
    let mut total_particles: usize = 0;
    for particle_index in 0..N {
        let E = particle_parameters.E[particle_index];
        let N_ = particle_parameters.N[particle_index];
        total_energy += E*(N_ as f64)*energy_unit;
    }
    let estimated_num_particles: usize = total_particles + ((total_energy/material.Ec).ceil() as usize);

    //Create particle array
    let mut particles: Vec<Particle> = Vec::with_capacity(estimated_num_particles);
    for particle_index in 0..N {
        let N_ = particle_parameters.N[particle_index];
        let m = particle_parameters.m[particle_index];
        let Z = particle_parameters.Z[particle_index];
        let E = particle_parameters.E[particle_index];
        let (x, y, z) = particle_parameters.pos[particle_index];
        let (dirx, diry, dirz) = particle_parameters.dir[particle_index];
        for sub_particle_index in 0..N_ {
            particles.push(Particle::new(
                m*mass_unit, Z, E*energy_unit,
                x*length_unit, y*length_unit, z*length_unit,
                dirx, diry, dirz, true, options.track_trajectories
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
    let mut reflected_file_stream = BufWriter::with_capacity(16000, reflected_file);

    let mut sputtered_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "sputtered.output"))
        .unwrap();
    let mut sputtered_file_stream = BufWriter::with_capacity(16000, sputtered_file);

    let mut deposited_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "deposited.output"))
        .unwrap();
    let mut deposited_file_stream = BufWriter::with_capacity(16000, deposited_file);

    let mut trajectory_file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "trajectories.output"))
        .unwrap();
    let mut trajectory_file_stream = BufWriter::with_capacity(16000, trajectory_file);

    let mut trajectory_data = OpenOptions::new()
        .write(true)
        .create(true)
        .open(format!("{}{}", options.name, "trajectory_data.output"))
        .unwrap();
    let mut trajectory_data_stream = BufWriter::with_capacity(16000, trajectory_data);

    //Main BCA loop
    let mut particle_index: usize = particles.len();

    'particle_loop: while particle_index > 0 {
        if particle_index % 1000 == 0{
            //println!("particle {} of {}", particle_index, particles.len());
        }

        let mut particle_1 = particles.pop().unwrap();

        'trajectory_loop: while !particle_1.stopped & !particle_1.left {

            if particle_1.track_trajectories {
                particle_1.add_trajectory();
            }

            //Check simulation boundary conditions, add to energy fluxes out
            if !material.inside_simulation_boundary(particle_1.pos.x, particle_1.pos.y) {
                particle_1.left = true;
                if particle_1.incident {
                    num_reflected += 1;
                    energy_reflected += particle_1.E;
                } else {
                    num_sputtered += 1;
                    energy_sputtered += particle_1.E;
                }
                //Skip BCA loop
                break 'trajectory_loop;
            }

            //Check stopping condition, if incident, contribute to average range
            if particle_1.E < material.Ec {
                particle_1.stopped = true;
                if particle_1.incident {
                    num_deposited += 1;
                    range += particle_1.pos.x;
                }
                //Skip BCA loop
                break 'trajectory_loop;
            }

            //BCA loop
            //Choose collision partner
            let (impact_parameter, phi_azimuthal, xr, yr, zr, dirx, diry, dirz) = pick_collision_partner(&particle_1, &material);
            let mut particle_2 = Particle::new(
                material.m, material.Z, 0.,
                xr, yr, zr,
                dirx, diry, dirz,
                false, options.track_recoil_trajectories
            );

            //Calculate scattering and update coordinates
            let (theta, psi, T, t) = binary_collision(&particle_1, &particle_2, &material, impact_parameter, 1E-6, 100);
            update_coordinates(&mut particle_1, &mut particle_2, &material, phi_azimuthal, theta, psi, T, t);

            //Reflection and refraction from surface energy barrier
            surface_boundary_condition(&mut particle_1, &material);

            //Generate recoils if transferred kinetic energy larger than cutoff energy
            if (T > material.Ec) & options.track_recoils {
                particle_2.add_trajectory();
                particles.push(particle_2);
            }
            //End BCA loop
        }

        //Stream current particle output to files
        if particle_1.incident & particle_1.left {
            writeln!(reflected_file_stream, "{},{},{},{},{},{},{},{},{}", particle_1.m/mass_unit, particle_1.Z, particle_1.E/energy_unit, particle_1.pos.x/length_unit, particle_1.pos.y/length_unit, particle_1.pos.z/length_unit, particle_1.dir.x, particle_1.dir.y, particle_1.dir.z)
                .expect("Could not write to reflected.output.");

        }
        if particle_1.incident & particle_1.stopped {
            writeln!(deposited_file_stream, "{},{},{},{},{}", particle_1.m/mass_unit, particle_1.Z, particle_1.pos.x/length_unit, particle_1.pos.y/length_unit, particle_1.pos.z/length_unit)
                .expect("Could not write to deposited.output.");
        }
        if !particle_1.incident & particle_1.left {
            writeln!(sputtered_file_stream, "{},{},{},{},{},{},{},{},{}", particle_1.m/mass_unit, particle_1.Z, particle_1.E/energy_unit, particle_1.pos.x/length_unit, particle_1.pos.y/length_unit, particle_1.pos.z/length_unit, particle_1.dir.x, particle_1.dir.y, particle_1.dir.z)
                .expect("Could not write to sputtered.output.");
        }
        if particle_1.track_trajectories {
            writeln!(trajectory_data_stream, "{}", particle_1.trajectory.len())
                .expect("Could not write trajectory length data.");

            for pos in particle_1.trajectory {
                writeln!(trajectory_file_stream, "{},{},{},{},{},{}", particle_1.m/mass_unit, particle_1.Z, pos.E/energy_unit, pos.x/length_unit, pos.y/length_unit, pos.z/length_unit)
                    .expect("Could not write to trajectories.output.");
            }
        }
        //Set particle index to topmost particle
        particle_index = particles.len();
    }
    //Flush all file streams before closing to ensure all data is written
    reflected_file_stream.flush().unwrap();
    deposited_file_stream.flush().unwrap();
    sputtered_file_stream.flush().unwrap();
    trajectory_data_stream.flush().unwrap();
    trajectory_file_stream.flush().unwrap();
}

fn main() {
    bca_input();
}
