extern crate ndarray;
extern crate geo;

pub use crate::geo::*;
use ndarray::prelude::*;
use geo::algorithm::contains::Contains;
use geo::algorithm::closest_point::ClosestPoint;
//use rand::Rng;
//use std::error::Error;
//use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::f64::consts::PI;

const Q: f64 = 1.602E-19;
const AMU: f64 = 1.66E-27;
const ANGSTROM: f64 = 1E-10;
const MICRON: f64 = 1E-6;
const EPS0: f64 = 8.85E-12;
const A0: f64 = 0.52918E-10;
const K: f64 = 1.11265E-10;
const ME: f64 = 9.11E-31;
const SQRTPI: f64 = 1.77245385;
const SQRT2PI: f64 = 2.506628274631;
const C: f64 = 300000000.;

pub struct Vector {
    x: f64,
    y: f64,
    z: f64,
}

impl Vector {
    fn magnitude(&self) -> f64 {
        return (self.x*self.x + self.y*self.y + self.z*self.z).sqrt();
    }
    fn assign(&mut self, other: &Vector) {
        self.x = other.x;
        self.y = other.y;
        self.z = other.z;
    }
}

pub struct Vector4 {
    E: f64,
    x: f64,
    y: f64,
    z: f64,
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
}

impl Particle {
    fn add_trajectory(&mut self) {
        self.trajectory.push(Vector4 {E: self.E/Q, x: self.pos.x, y: self.pos.y, z: self.pos.z});
    }
}

pub struct Material {
    n: f64,
    m: f64,
    Z: f64,
    Eb: f64,
    Es: f64,
    surface: Polygon<f64>,
    energy_surface: Polygon<f64>,
    simulation_surface: Polygon<f64>,
}

impl Material {
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

        //Guess for large reduced energy from Mendenhall and Weller, 1991
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

        //Scattering integral quadrature from Mendenhall and Weller, 2005
        let lambda_0 = (0.5 + beta*beta/x0/x0/2. - dphi(x0)/2_f64/reduced_energy).powf(-0.5);
        let alpha = 1_f64/12_f64*(1_f64 + lambda_0 + 5_f64*(0.4206_f64*f(x0/0.9072_f64, beta, reduced_energy) + 0.9072_f64*f(x0/0.4206, beta, reduced_energy)));
        let theta = PI*(1_f64 - beta*alpha/x0);

        //See Eckstein 1991 for details on center of mass and lab frame angles
        let t = x0*a*(theta/2.).sin();
        let psi = (theta.sin()).atan2(Ma/Mb + theta.cos());
        let T = 4.*(Ma*Mb)/(Ma + Mb).powf(2.)*E0*((theta/2.).sin()).powf(2.);

        return (theta, psi, T, t);
    }
}

fn update_coordinates(particle_1: &mut Particle, particle_2: &mut Particle, material: &Material, phi_azimuthal: f64, theta: f64, psi: f64, T: f64, t: f64) {
    let mut mfp = material.mfp(particle_1.pos.x, particle_1.pos.y);

    particle_1.pos_old.assign(&particle_1.pos);

    //TRIDYN-style atomically rough surface - scatters initial collisions by up to 1 mfp
    if particle_1.first_step {
        mfp *= rand::random::<f64>();
        particle_1.first_step = false;
    }

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

        let dir_new_mag: f64 = (ca_new*ca_new + cb_new*cb_new + cg_new*cg_new).sqrt();
        let dir_new = Vector {x: ca_new/dir_new_mag, y: cb_new/dir_new_mag, z: cg_new/dir_new_mag};

        particle_1.dir.assign(&dir_new);
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

        let dir_new_mag: f64 = (ca_new*ca_new + cb_new*cb_new + cg_new*cg_new).sqrt();
        let dir_new = Vector {x: ca_new/dir_new_mag, y: cb_new/dir_new_mag, z: cg_new/dir_new_mag};

        particle_2.dir.assign(&dir_new);
    }

    //Update particle energies
    {
        let Za = particle_1.Z;
        let Ma = particle_1.m;
        let Zb = material.Z;
        let E = particle_1.E;

        let mut stopping_power = 0.;

        //Electronic stopping only inside material
        if material.inside(particle_1.pos.x, particle_1.pos.y) {
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

            let mut B;//Effective empirical shell correction from Biersack and Haggmark
            if Zb < 3. {
                B = 100.*Za/Zb;
            } else {
                B = 5.;
            }

            let prefactor = 8.0735880E-42*Zb*Za*Za/beta/beta; //This is messy, but this prefactor is... not the easiest to compute
            let eb = 2.*ME*v*v/I;
            let S_BB = prefactor*(eb + 1. + B/eb).ln()*n; //Bethe-Bloch stopping, as modified by Biersack and Haggmark (1980) to fit experiment

            //Lindhard-Scharff stopping for low energies
            let S_LS = 1.212*(Za.powf(7./6.)*Zb)/(Za.powf(2./3.) + Zb.powf(2./3.)).powf(3./2.)*(E/Ma*AMU/Q).sqrt()*ANGSTROM*ANGSTROM*Q*n;

            //Biersack-Varelas stopping interpolation scheme
            stopping_power = 1./(1./S_BB + 1./S_LS);
            //println!("S   : {}", stopping_power*6.262E10/8.96);

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

fn pick_collision_partner(particle_1: &Particle, material: &Material) -> (f64, f64, Vector, Vector) {
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

    let dir: Vector = Vector {x: ca, y: cb, z: cg};

    //Collision partner at locus of next collision, displaced by chosen impact parameter and rotated by phi azimuthal
    let x_recoil: f64 = particle_1.pos.x + mfp*ca - impact_parameter*cphi*sa;
    let y_recoil: f64 = particle_1.pos.y + mfp*cb - impact_parameter*(sphi*cg - cphi*cb*ca)/sa;
    let z_recoil: f64 = particle_1.pos.z + mfp*cg + impact_parameter*(sphi*cb - cphi*ca*cg)/sa;

    let pos = Vector {x: x_recoil, y: y_recoil, z: z_recoil};

    return (impact_parameter, phi_azimuthal, dir, pos);
}

fn surface_boundary_condition(particle_1: &mut Particle, material: &Material) -> bool {
    //Planar energy barrier will refract particle and deflect it towards surface
    if !material.inside_energy_barrier(particle_1.pos.x, particle_1.pos.y) & material.inside_energy_barrier(particle_1.pos_old.x, particle_1.pos_old.y) {
        //let p = point!(x: particle_1.pos.x, y: particle_1.pos.y);

        let mut leaving_energy;

        //Find closest point to surface, and use line from particle to point as normal to refract
        if let Closest::SinglePoint(p2) = material.closest_point(particle_1.pos.x, particle_1.pos.y) {
            let dx = p2.x() - particle_1.pos.x;
            let dy = p2.y() - particle_1.pos.y;
            let dz: f64 = 0.;

            let magnitude = (dx*dx + dy*dy + dz*dz).sqrt();

            //Since direction check is already handled in this function, just take abs() here
            leaving_energy = particle_1.E*((dx/magnitude*particle_1.dir.x).abs() + (dy/magnitude*particle_1.dir.y).abs() + (dz/magnitude*particle_1.dir.z).abs());
            //Project energy-scaled velocity vector along direction to nearest point of surface*--*-++++
        } else {
            leaving_energy = particle_1.E;
        }

        if leaving_energy < material.Es {
            //let p = point!(x: particle_1.pos.x, y: particle_1.pos.y);

            //This syntax is a kind of enum matching - closest_point doesn't return the point,
            //Since there are multiple kinds of nearest-point solutions, so it returns an enum
            //That stores the point inside. This code extracts it
            if let Closest::SinglePoint(p2) = material.closest_point(particle_1.pos.x, particle_1.pos.y) {
                let dx = p2.x() - particle_1.pos.x;
                let dy = p2.y() - particle_1.pos.y;

                let axis = Vector {x: dx, y: dy, z: 0.};

                particle_1.dir.x *= -1.;
                particle_1.dir.y *= -1.;

                let new_dir: Vector = rotate_around_axis_by_pi(&particle_1.dir, &axis);

                particle_1.dir.assign(&new_dir);
            }
            return false;
        } else {
            surface_refraction(particle_1, &material);
            return true;
        }
    } else {
        return false;
    }
}

fn surface_refraction(particle_1: &mut Particle, material: &Material) {
    //Deflection code for planar surface refraction
    let Es = material.Es;
    let E0 = particle_1.E;

    let cosx0 = particle_1.dir.x;
    let cosy0 = particle_1.dir.y;
    let cosz0 = particle_1.dir.z;

    //let p = point!(x: particle_1.pos.x, y: particle_1.pos.y);

    if let Closest::SinglePoint(p2) = material.closest_point(particle_1.pos.x, particle_1.pos.y) {
        let dx = particle_1.pos.x - p2.x();
        let dy = particle_1.pos.x - p2.y();
        let dz: f64 = 0.;

        let dot_product = dx*cosx0 + dy*cosy0 + dz*cosz0;
        let sign = dot_product.signum();

        let sign_cosx = cosx0.signum();
        let sign_cosy = cosy0.signum();
        let sign_cosz = cosz0.signum();

        let magnitude = (dx*dx + dy*dy + dz*dz).sqrt();

        let new_cosx = sign_cosx*((cosx0*cosx0*E0 + sign*Es*dx*dx/magnitude/magnitude)/(E0 + sign*Es)).sqrt();
        let new_cosy = sign_cosy*((cosy0*cosy0*E0 + sign*Es*dy*dy/magnitude/magnitude)/(E0 + sign*Es)).sqrt();
        let new_cosz = sign_cosz*(E0*cosz0*cosz0/(E0 + sign*Es)).sqrt();

        particle_1.dir.x = new_cosx;
        particle_1.dir.y = new_cosy;
        particle_1.dir.z = new_cosz;
    }
}

fn bca(N: usize, E0: f64, theta: f64, Ec: f64, Ma: f64, Za: f64, Mb: f64, Zb: f64, Eb: f64, Es: f64, n: f64, track_recoils: bool, track_trajectories: bool, track_recoil_trajectories: bool, write_files: bool, thickness: f64, depth: f64, name: String) {
    let mut particles: Vec<Particle> = vec![];

    //Currently, Material just stores the polygons, which are created on struct initialization
    let dx: f64 = 2_f64*n.powf(-1_f64/3_f64)/SQRT2PI;
    let angles = Array::linspace(0., 2.*PI, 16);
    let r = thickness/2.;
    let re = (thickness + dx)/2.;

    let circ_material = Material {
        n: n,
        m: Mb,
        Z: Zb,
        Eb: Eb,
        Es: Es,
        surface: polygon![
            (x: r*angles[0].cos() + r, y: r*angles[0].sin()),
            (x: r*angles[1].cos() + r, y: r*angles[1].sin()),
            (x: r*angles[2].cos() + r, y: r*angles[2].sin()),
            (x: r*angles[3].cos() + r, y: r*angles[3].sin()),
            (x: r*angles[4].cos() + r, y: r*angles[4].sin()),
            (x: r*angles[5].cos() + r, y: r*angles[5].sin()),
            (x: r*angles[6].cos() + r, y: r*angles[6].sin()),
            (x: r*angles[7].cos() + r, y: r*angles[7].sin()),
            (x: r*angles[8].cos() + r, y: r*angles[8].sin()),
            (x: r*angles[9].cos() + r, y: r*angles[9].sin()),
            (x: r*angles[10].cos() + r, y: r*angles[10].sin()),
            (x: r*angles[11].cos() + r, y: r*angles[11].sin()),
            (x: r*angles[12].cos() + r, y: r*angles[12].sin()),
            (x: r*angles[13].cos() + r, y: r*angles[13].sin()),
            (x: r*angles[14].cos() + r, y: r*angles[14].sin()),
            (x: r*angles[15].cos() + r, y: r*angles[15].sin()),
        ],
        energy_surface: polygon![
            (x: re*angles[0].cos() + re, y: re*angles[0].sin()),
            (x: re*angles[1].cos() + re, y: re*angles[1].sin()),
            (x: re*angles[2].cos() + re, y: re*angles[2].sin()),
            (x: re*angles[3].cos() + re, y: re*angles[3].sin()),
            (x: re*angles[4].cos() + re, y: re*angles[4].sin()),
            (x: re*angles[5].cos() + re, y: re*angles[5].sin()),
            (x: re*angles[6].cos() + re, y: re*angles[6].sin()),
            (x: re*angles[7].cos() + re, y: re*angles[7].sin()),
            (x: re*angles[8].cos() + re, y: re*angles[8].sin()),
            (x: re*angles[9].cos() + re, y: re*angles[9].sin()),
            (x: re*angles[10].cos() + re, y: re*angles[10].sin()),
            (x: re*angles[11].cos() + re, y: re*angles[11].sin()),
            (x: re*angles[12].cos() + re, y: re*angles[12].sin()),
            (x: re*angles[13].cos() + re, y: re*angles[13].sin()),
            (x: re*angles[14].cos() + re, y: re*angles[14].sin()),
            (x: re*angles[15].cos() + re, y: re*angles[15].sin()),
        ],
        simulation_surface: polygon![
            (x: 0.0 - 2.*dx, y: -thickness/2. - 2.*dx),
            (x: depth + 2.*dx, y: -thickness/2. - 2.*dx),
            (x: depth + 2.*dx, y: thickness/2. + 2.*dx),
            (x: 0.0 - 2.*dx, y: thickness/2. + 2.*dx),
        ],
    };

    let material = Material {
        n: n,
        m: Mb,
        Z: Zb,
        Eb: Eb,
        Es: Es,
        surface: polygon![
            (x: 0.0, y: -thickness/2.),
            (x: depth, y: -thickness/2.),
            (x: depth, y: thickness/2.),
            (x: 0.0, y: thickness/2.),
        ],
        energy_surface: polygon![
            (x: 0.0 - dx, y: -thickness/2. - dx),
            (x: depth + dx, y: -thickness/2. - dx),
            (x: depth + dx, y: thickness/2. + dx),
            (x: 0.0 - dx, y: thickness/2. + dx),
        ],
        simulation_surface: polygon![
            (x: 0.0 - 2.*dx, y: -thickness/2. - 2.*dx),
            (x: depth + 2.*dx, y: -thickness/2. - 2.*dx),
            (x: depth + 2.*dx, y: thickness/2. + 2.*dx),
            (x: 0.0 - 2.*dx, y: thickness/2. + 2.*dx),
        ],
    };

    let x0 = -dx;
    let y0 = 0.;
    let z0 = 0.;
    let cosx = (theta*PI/180_f64).cos();
    let sinx = (theta*PI/180_f64).sin();

    //Create N particles, randomly scattered along y
    for i in 0..N {
        let dy = (2.*rand::random::<f64>() - 1.)*thickness/2.;
        particles.push(
            Particle {
                m: Ma,
                Z: Za,
                E: E0,
                pos: Vector {x: x0, y: y0 + dy, z: z0},
                dir: Vector {x: cosx, y: sinx, z: 0.},
                pos_old: Vector  {x: x0, y: y0 + dy, z: z0},
                pos_origin: Vector {x: x0, y: y0 + dy, z: z0},
                t: 0_f64,
                left: false,
                stopped: false,
                incident: true,
                first_step: true,
                track_trajectories: track_trajectories,
                trajectory: vec![],
            }
        )
    }

    let mut particle_index: usize = 0;
    while particle_index < particles.len() {
        if particle_index % 10 == 0 {
            println!("particle {} of {}", particle_index+1, particles.len());
        }
        while !particles[particle_index].stopped & !particles[particle_index].left {

            if !material.inside_simulation_boundary(particles[particle_index].pos.x, particles[particle_index].pos.y) {
                particles[particle_index].left = true;
                //println!("Left! {} {} {}", particles[particle_index].pos.x/ANGSTROM, particles[particle_index].pos.y/ANGSTROM, particles[particle_index].pos.z/ANGSTROM);
                continue;
            }

            if particles[particle_index].E < Ec {
                particles[particle_index].stopped = true;
                continue;
            }

            let (impact_parameter, phi_azimuthal, dir, pos) = pick_collision_partner(&particles[particle_index], &material);

            let mut particle_2 = Particle {
                m: material.m,
                Z: material.Z,
                E: 0.,
                pos: Vector {x: pos.x, y: pos.y, z: pos.z},
                pos_old: Vector {x: pos.x, y: pos.y, z: pos.z},
                pos_origin: Vector {x: pos.x, y: pos.y, z: pos.z},
                dir: Vector {x: dir.x, y: dir.y, z: dir.z},
                t: 0.,
                left: false,
                stopped: false,
                incident: false,
                first_step: false,
                track_trajectories: track_recoil_trajectories,
                trajectory: vec![],
            };

            let (theta, psi, T, t) = binary_collision(&particles[particle_index], &particle_2, &material, impact_parameter, 1E-3, 100);
            update_coordinates(&mut particles[particle_index], &mut particle_2, &material, phi_azimuthal, theta, psi, T, t);
            surface_boundary_condition(&mut particles[particle_index], &material);

            if particles[particle_index].track_trajectories{
                particles[particle_index].add_trajectory();
            }

            if (T > Ec ) & track_recoils {
                particles.push(particle_2);
                let length = particles.len();
                particles[length - 1].add_trajectory();
            }
        }
        particle_index += 1
    }

    //Data output!
    if write_files {
        let mut reflected_file = OpenOptions::new().write(true).create(true).open(format!("{}{}", name, "reflected.dat")).unwrap();
        let mut sputtered_file = OpenOptions::new().write(true).create(true).open(format!("{}{}", name, "sputtered.dat")).unwrap();
        let mut deposited_file = OpenOptions::new().write(true).create(true).open(format!("{}{}", name, "deposited.dat")).unwrap();
        let mut trajectory_file = OpenOptions::new().write(true).create(true).open(format!("{}{}", name, "trajectories.dat")).unwrap();
        let mut trajectory_data = OpenOptions::new().write(true).create(true).open(format!("{}{}", name, "trajectory_data.dat")).unwrap();

        let mut num_sputtered: usize = 0;
        let mut num_reflected: usize = 0;
        let mut num_deposited: usize = 0;

        let mut E_sputtered: f64 = 0.;
        let mut E_reflected: f64 = 0.;
        let mut E_total: f64 = 0.;

        let mut range: f64 = 0.;

        for particle in particles {
            E_total += particle.E;
            if particle.incident & particle.left {
                num_reflected += 1;
                E_reflected += particle.E;
                if write_files {
                    writeln!(reflected_file, "{}, {}, {}, {}, {}, {}, {}, {}", particle.Z, particle.pos.x, particle.pos.y, particle.pos.z, particle.dir.x, particle.dir.y, particle.dir.z, particle.E);
                }
            }

            if particle.incident & particle.stopped {
                num_deposited += 1;
                range += particle.pos.x;
                if write_files {
                    writeln!(deposited_file, "{}, {}, {}, {}", particle.Z, particle.pos.x, particle.pos.y, particle.pos.z);
                }
            }

            if !particle.incident & particle.left {
                num_sputtered += 1;
                E_sputtered += particle.E;
                if write_files {
                    writeln!(sputtered_file, "{}, {}, {}, {}, {}, {}, {}, {}", particle.Z, particle.pos.x, particle.pos.y, particle.pos.z, particle.dir.x, particle.dir.y, particle.dir.z, particle.E);
                }
            }

            if particle.track_trajectories & write_files {
                writeln!(trajectory_data, "{}", particle.trajectory.len());
                for pos in particle.trajectory {
                    writeln!(trajectory_file, "{}, {}, {}, {}, {}", particle.Z, pos.E, pos.x, pos.y, pos.z);
                }
            }
        }
        println!("E: {} Range: {} R: {} E_r: {} Y: {} E_s: {}", E0/Q, range/(num_deposited as f64)/ANGSTROM, (num_reflected as f64)/(N as f64), E_reflected/Q/(num_reflected as f64), (num_sputtered as f64)/(N as f64), E_sputtered/Q/(num_sputtered as f64));
    }
}

fn main() {
    let num_energies = 1;
    let thickness = 1.*MICRON;
    let depth = 1.*MICRON;
    let energies = Array::logspace(10., 6., 4., num_energies);
    for index in 0..num_energies {
        let name: String = index.to_string();
        //track_recoils track_trajectories track_recoil_trajectories write_files
        bca(10000, 1E5*Q, 0.0001, 3.*Q, 4.*AMU, 2., 63.54*AMU, 29., 0.0, 3.52*Q, 8.491E28, true, false, false, true, thickness, depth, name);
    }
}
