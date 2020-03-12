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

use rand::SeedableRng;

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
    print: bool,
    print_num: usize,
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
    backreflected: bool
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
            backreflected: false
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
            "KG" => 1.,
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

    fn Eb(&self, x: f64, y: f64) -> f64 {
        return self.Eb;
    }

    fn choose(&self, x: f64, y: f64) -> (f64, f64, f64, f64) {
        return (self.Z, self.m, self.Ec, self.Es);
    }

    fn electronic_stopping_power(&self, particle_1: &Particle) -> f64 {
        let n = self.number_density(particle_1.pos.x, particle_1.pos.y);
        let E = particle_1.E;
        let Ma = particle_1.m;
        let Za = particle_1.Z;
        let Zb = self.Z_eff(particle_1.pos.x, particle_1.pos.y);

        let beta = (1. - (1. + E/Ma/C/C).powf(-2.)).sqrt();
        //let v = beta/C;

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
        let eb = 2.*ME*beta*beta/C/C/I;
        let S_BB = prefactor*(eb + 1. + B/eb).ln()*n; //Bethe-Bloch stopping, as modified by Biersack and Haggmark (1980) to fit experiment

        //Lindhard-Scharff stopping for low energies
        let S_LS = LINDHARD_SCHARFF_PREFACTOR*(Za.powf(7./6.)*Zb)/(Za.powf(2./3.) + Zb.powf(2./3.)).powf(3./2.)*(E/Ma*AMU/Q).sqrt()*n;

        //Biersack-Varelas stopping interpolation scheme
        let stopping_power = 1./(1./S_BB + 1./S_LS);

        return stopping_power;
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
    return 0.8853*A0*(Za.sqrt() + Zb.sqrt()).powf(-2./3.);
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

fn choose_collision_partner(particle_1: &mut Particle, material: &Material) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let z = particle_1.pos.z;

    //Determine mfp
    let mut mfp = material.mfp(x, y);
    let pmax : f64 = mfp/SQRTPI;//Cylindrical

    //Determine phi and impact parameter
    let mut impact_parameter: f64 = pmax*(rand::random::<f64>()).sqrt();
    let mut phi_azimuthal: f64 = 2.*PI*rand::random::<f64>();

    //Atomically rough surface
    if particle_1.first_step {
        mfp = 2.*material.n.powf(-1./3.)/SQRT2PI + mfp*rand::random::<f64>();
        particle_1.first_step = false;
    }

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

    return (mfp, impact_parameter, phi_azimuthal, Z_recoil, M_recoil, Ec_recoil, Es_recoil, x_recoil, y_recoil, z_recoil, ca, cb, cg);
}

fn calculate_binary_collision(particle_1: &Particle, particle_2: &Particle, impact_parameter: f64, max_iter: usize, tol: f64) -> (f64, f64, f64, f64, f64) {
    let Za: f64 = particle_1.Z;
    let Zb: f64 = particle_2.Z;
    let Ma: f64 = particle_1.m;
    let Mb: f64 = particle_2.m;
    let E0: f64 = particle_1.E;
    let mu: f64 = Mb/(Ma + Mb);

    //Lindhard screening length and reduced energy
    let a: f64 = screening_length(Za, Zb);
    let reduced_energy: f64 = 4.*PI*EPS0*a*Mb*E0/(Ma + Mb)/Za/Zb/Q/Q;
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


    //MAGIC algorithm
    let C_ = vec![ 1.0144, 0.235800, 0.126, 63950.0, 83550.0 ];
    let c = vec![ 0.190945, 0.473674, 0.335381, 0.0 ];
    let d = vec![ -0.278544, -0.637174, -1.919249, 0.0 ];
    let a_ = 0.8853*A0/((Za).sqrt() + (Zb).sqrt()).powf(2./3.);
    let V0 = Za*Zb*Q*Q/4.0/PI/EPS0/a_;
    let V00 = V0*a_;
    let E_c = E0*Mb/(Ma + Mb);
    let E_r = E0/V0;
    let b = impact_parameter/a_;
    let SQE = E_r.sqrt();
    let R = a_*x0;
    let sum = c[0]*(d[0]*x0).exp() + c[1]*(d[1]*x0).exp() + c[2]*(d[2]*x0).exp();
    let V = V00/R*sum;
    let sum = d[0]*c[0]*(d[0]*x0).exp() + d[1]*c[1]*(d[1]*x0).exp() + d[2]*c[2]*(d[2]*x0).exp();
    let dV = -V/R + V0/R*sum;
    let rho   = -2.0*(E_c - V)/dV;
    //let D = 2.0*(1.0+C_[0]/SQE)*E_r* pow(b,((C_[1]+SQE)/(C_[2]+SQE)));
    let D = 2.0*(1.0+C_[0]/SQE)*E_r*b.powf((C_[1]+SQE)/(C_[2]+SQE));
    let G = (C_[4]+E_r)/(C_[3]+E_r)*((1.0+D*D).sqrt()-D);
    let delta =  D*G/(1.0+G)*(x0-b);
    let ctheta2 = (b + rho/a_ + delta)/(x0 + rho/a_);
    let theta = 2. * ( (b + rho/a_ + delta)/(x0 + rho/a_) ).acos(); // add a 1e-3 exception?


    //Scattering integral quadrature from Mendenhall and Weller 2005
    let lambda_0 = (0.5 + beta*beta/x0/x0/2. - dphi(x0)/2./reduced_energy).powf(-1./2.);
    let alpha = 1./12.*(1. + lambda_0 + 5.*(0.4206*f(x0/0.9072, beta, reduced_energy) + 0.9072*f(x0/0.4206, beta, reduced_energy)));
    let theta_ = PI*(1. - beta*alpha/x0);
    //println!("theta M-W: {}", theta_mw);
    //println!("theta MAGIC: {}", theta);
    //println!("{} {}", theta, theta_mw);

    //See Eckstein 1991 for details on center of mass and lab frame angles
    let asympototic_deflection = x0*a*(theta/2.).sin();
    let psi = theta.sin().atan2(Ma/Mb + theta.cos());
    let psi_recoil = theta.sin().atan2(1. - theta.cos());
    let recoil_energy = 4.*(Ma*Mb)/(Ma + Mb).powf(2.)*E0*(theta/2.).sin().powf(2.);

    return (theta, psi, psi_recoil, recoil_energy, asympototic_deflection);
}

fn rotate_particle(particle_1: &mut Particle, psi: f64, phi: f64) {
    let ca: f64 = particle_1.dir.x;
    let cb: f64 = particle_1.dir.y;
    let cg: f64 = particle_1.dir.z;
    let cphi: f64 = phi.cos();
    let sphi: f64 = phi.sin();
    let sa = (1. - ca.powf(2.)).sqrt();

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

fn update_particle_energy(particle_1: &mut Particle, material: &Material, path_length: f64, recoil_energy: f64) {
    particle_1.E = particle_1.E - recoil_energy;
    if material.inside_1D(particle_1.pos.x) {
        let electronic_stopping_power = material.electronic_stopping_power(particle_1);
        particle_1.E += -electronic_stopping_power*path_length;
    }
}

fn particle_advance(particle_1: &mut Particle, mfp: f64, asympototic_deflection: f64) -> f64 {
    particle_1.pos_old.x = particle_1.pos.x;
    particle_1.pos_old.y = particle_1.pos.y;
    particle_1.pos_old.z = particle_1.pos.z;

    let path_length = mfp + particle_1.asympototic_deflection - asympototic_deflection;

    particle_1.pos.x += particle_1.dir.x*path_length;
    particle_1.pos.y += particle_1.dir.y*path_length;
    particle_1.pos.z += particle_1.dir.z*path_length;

    particle_1.asympototic_deflection = asympototic_deflection;

    return path_length;
}

fn boundary_condition(particle_1: &mut Particle, material: &Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let x0 = particle_1.pos_old.x;
    let y0 = particle_1.pos_old.y;
    let cosx = particle_1.dir.x;
    let cosy = particle_1.dir.y;
    let cosz = particle_1.dir.z;
    let E = particle_1.E;
    let xc = -2.*material.number_density(x, y).powf(-1./3.)/SQRT2PI;
    let Es = particle_1.Es;
    let Ec = particle_1.Ec;

    if (x < xc) & (cosx < 0.) {
        let leaving_energy = E*cosx.powf(2.);

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

        } else {
            particle_1.backreflected = true;
            particle_1.dir.x = cosx.abs();
        }
    }

    if (E < Ec) & !particle_1.left {
        particle_1.stopped = true;
    }
}

fn bca_file_input() {
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
    let mut max_energy: f64 = 0.;
    let mut total_particles: usize = 0;
    for particle_index in 0..N {
        let E = particle_parameters.E[particle_index];
        let Ec = particle_parameters.Ec[particle_index];
        let N_ = particle_parameters.N[particle_index];
        if E > max_energy {
            max_energy = E*energy_unit;
        }
        total_particles += N_;
    }

    //Create particle array
    let estimated_num_particles: usize = total_particles + ((max_energy/material.Ec).ceil() as usize);
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

            particles.push(Particle::new(
                m*mass_unit, Z, E_new, Ec*energy_unit, Es*energy_unit,
                x*length_unit, y*length_unit, z*length_unit,
                cosx_new, cosy_new, cosz_new, true, options.track_trajectories
            ));
            //println!("{} {} {}", particles[particles.len() - 1].E, particles[particles.len() - 1].m, particles[particles.len() - 1].dir.x);
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
        let mut particle_1 = particles.pop().unwrap();

        if options.print & particle_1.incident & (particle_index % (total_particles / options.print_num) == 0){
            println!("Particle {}", particle_index);
        }

        'trajectory_loop: while !particle_1.stopped & !particle_1.left {

            if particle_1.track_trajectories {
                particle_1.add_trajectory();
            }

            //BCA loop
            //Choose recoil partner location and species
            let (mfp, impact_parameter, phi_azimuthal, Z_recoil, M_recoil, Ec_recoil, Es_recoil, xr, yr, zr, cxr, cyr, czr) = choose_collision_partner(&mut particle_1, &material);

            //If recoil location is inside, proceed with binary collision loop
            if material.inside_1D(xr) {

                particle_1.backreflected = false;

                let mut particle_2 = Particle::new(
                    M_recoil, Z_recoil, 0., Ec_recoil, Es_recoil,
                    xr, yr, zr,
                    cxr, cyr, czr,
                    false, options.track_recoil_trajectories
                );

                //Determine scattering angle from binary collision
                let (theta, psi, psi_recoil, recoil_energy, asympototic_deflection) = calculate_binary_collision(&particle_1, &particle_2, impact_parameter, 100, 1E-6);
                //Energy transfer to recoil
                particle_2.E = recoil_energy - material.Eb;

                //Rotate particle 1 by lab frame scattering angle
                rotate_particle(&mut particle_1, psi, phi_azimuthal);
                //Rotate particle 2 by lab frame scattering angle
                rotate_particle(&mut particle_2, psi_recoil, PI - phi_azimuthal);

                //If transferred energy > cutoff energy, add recoil to particle vector
                if (recoil_energy > material.Ec) & options.track_recoils {
                    particle_2.add_trajectory();
                    particles.push(particle_2);
                }

                //Push particle, find actual path length from mfp and asymptotic deflection
                let path_length = particle_advance(&mut particle_1, mfp, asympototic_deflection);

                //Update moving particle energy
                update_particle_energy(&mut particle_1, &material, path_length, recoil_energy);

            } else {
                //Push particle, find actual path length from mfp and asymptotic deflection
                let path_length = particle_advance(&mut particle_1, mfp, 0.);

                //Update moving particle energy
                update_particle_energy(&mut particle_1, &material, mfp, 0.);
            }
            //Check boundary conditions on leaving and stopping
            boundary_condition(&mut particle_1, &material);
        }

        //Once particle finishes, begin data output
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
    //Flush all file streams before dropping streams to ensure all data is written
    reflected_file_stream.flush().unwrap();
    deposited_file_stream.flush().unwrap();
    sputtered_file_stream.flush().unwrap();
    trajectory_data_stream.flush().unwrap();
    trajectory_file_stream.flush().unwrap();
}


fn main() {
    bca_file_input();
}
