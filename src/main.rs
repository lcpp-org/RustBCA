extern crate ndarray;

use ndarray::prelude::*;
use rand::Rng;
use std::error::Error;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;

const Q: f64 = 1.602E-19;
const PI: f64 = 3.14159;
const AMU: f64 = 1.66E-27;
const ANGSTROM: f64 = 1E-10;
const EPS0: f64 = 8.85E-12;
const A0: f64 = 0.52918E-10;
const K: f64 = 1.11265E-10;
const ME: f64 = 9.11E-31;
const SQRTPI: f64 = 1.77245385;
const SQRT2PI: f64 = 2.506628274631;

pub struct Particle {
    m: f64,
    Z: f64,
    E: f64,
    pos: Array1<f64>,
    dir: Array1<f64>,
    pos_old: Array1<f64>,
    pos_origin: Array1<f64>,
    t: f64,
    stopped: bool,
    left: bool,
    incident: bool,
    first_step: bool,
    trajectory: Vec<Array1<f64>>,
    track_trajectories: bool,
}

impl Particle {
    fn add_trajectory(&mut self) {
        self.trajectory.push(array![self.pos[0], self.pos[1], self.pos[2]]);
    }
}

pub struct Material {
    n: f64,
    m: f64,
    Z: f64,
    Eb: f64,
    Es: f64,
    surface_position: f64,
    energy_barrier: f64,
    simulation_boundary: f64,
}

impl Material {
    fn inside(&self, x: f64) -> bool {
        return x > self.surface_position;
    }
    fn mfp(&self, x: f64) -> f64 {
        return self.n.powf(-1_f64/3_f64);
    }
    fn number_density(&self, x: f64) -> f64 {
        return self.n;
    }
    fn inside_energy_barrier(&self, x: f64) -> bool {
        return x > self.energy_barrier;
    }
    fn inside_simulation_boundary(&self, x:f64) -> bool {
        return x > self.simulation_boundary;
    }
}

fn phi(xi: f64) -> f64 {
    return 0.191*(-0.279*xi).exp() + 0.474*(-0.637*xi).exp() + 0.335*(-1.919*xi).exp();
}

fn dphi(xi: f64) -> f64 {
    return -0.279*0.191*(-0.279*xi).exp() - 0.637*0.474*(-0.637*xi).exp() - 0.335*1.919*(-1.919*xi).exp();
}

fn screening_length(Za: f64, Zb: f64) -> f64 {
    return 0.8853*A0/(Za.sqrt() + Zb.sqrt()).powf(2./3.);
}

fn doca_function(x0: f64, beta: f64, reduced_energy: f64) -> f64 {
    return x0 - phi(x0)/reduced_energy - beta*beta/x0;
}

fn diff_doca_function(x0: f64, beta: f64, reduced_energy: f64) -> f64 {
    return beta*beta/x0/x0 - dphi(x0)/reduced_energy + 1.
}

fn f(x: f64, beta: f64, reduced_energy: f64) -> f64 {
    return (1_f64 - phi(x)/x/reduced_energy - beta*beta/x/x).powf(-0.5);
}

fn binary_collision(particle_1: &Particle, particle_2: &Particle,  material: &Material, impact_parameter: f64, tol: f64, max_iter: i32) -> (f64, f64, f64, f64, f64) {
    if !material.inside(particle_2.pos[0]) {
        return (0., 0., 0., 0., 0.);
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

        let mut x0 = 1_f64;
        let mut xn: f64;
        if reduced_energy > 5_f64 {
            let inv_er_2 = 0.5_f64/reduced_energy;
            x0 = inv_er_2 + (inv_er_2.powf(2_f64) + beta.powf(2_f64)).sqrt();
        }

        let mut err: f64;
        for iter in 0..max_iter {
            xn = x0 - doca_function(x0, beta, reduced_energy)/diff_doca_function(x0, beta, reduced_energy);
            err = (xn - x0).abs()/xn;
            x0 = xn;
            if err < tol {
                break;
            }
        }

        let lambda_0 = (0.5 + beta*beta/x0/x0/2. - dphi(x0)/2_f64/reduced_energy).powf(-0.5);
        let alpha = 1_f64/12_f64*(1_f64 + lambda_0 + 5_f64*(0.4206_f64*f(x0/0.9072_f64, beta, reduced_energy) + 0.9072_f64*f(x0/0.4206, beta, reduced_energy)));
        let theta = PI*(1_f64 - beta*alpha/x0);

        let t = x0*a*(theta/2.).sin();
        let psi = (theta.sin()).atan2(Ma/Mb + theta.cos());
        let T = 4.*(Ma*Mb)/(Ma + Mb).powf(2.)*E0*((theta/2.).sin()).powf(2.);

        return (theta, psi, T, t, x0);
    }
}

fn update_coordinates(particle_1: &mut Particle, particle_2: &mut Particle, material: &Material, phi_azimuthal: f64, theta: f64, psi: f64, T: f64, t: f64) {
    let mut mfp = material.mfp(particle_1.pos[0]);

    particle_1.pos_old.assign(&particle_1.pos); //Replace array with other array values

    if particle_1.first_step {
        //let mut rng = rand::thread_rng();
        mfp *= rand::random::<f64>();
        particle_1.first_step = false;
    }

    //println!("{} {} {}", particle_1.pos[0], particle_1.pos[1], particle_1.pos[2]);

    let ffp: f64 = mfp - t + particle_1.t;

    let new_pos: Array1<f64> = array![particle_1.pos[0] + ffp*particle_1.dir[0], particle_1.pos[1] + ffp*particle_1.dir[1], particle_1.pos[2] + ffp*particle_1.dir[2]];
    //let new_pos: Array1<f64> = particle_1.pos + Array1::from_elem(3, ffp)*particle_1.dir;
    particle_1.pos.assign(&new_pos);
    particle_1.t = t;

    let ca: f64 = particle_1.dir[0];
    let cb: f64 = particle_1.dir[1];
    let cg: f64 = particle_1.dir[2];
    let cphi: f64 = phi_azimuthal.cos();
    let sphi: f64 = phi_azimuthal.sin();

    let sa = (1_f64 - ca.powf(2_f64)).sqrt();
    {
        let cpsi: f64 = psi.cos();
        let spsi: f64 = psi.sin();
        let ca_new: f64 = cpsi*ca + spsi*cphi*sa;
        let cb_new: f64 = cpsi*cb - spsi/sa*(cphi*ca*cb - sphi*cg);
        let cg_new: f64 = cpsi*cg - spsi/sa*(cphi*ca*cg + sphi*cb);

        let dir_new_mag: f64 = (ca_new*ca_new + cb_new*cb_new + cg_new*cg_new).sqrt();
        let dir_new: Array1<f64> = array![ca_new/dir_new_mag, cb_new/dir_new_mag, cg_new/dir_new_mag];
        particle_1.dir.assign(&dir_new);
    }

    {
        let psi_b: f64 = (-(theta.sin())).atan2(1. - theta.cos());
        //println!("{}", psi_b);
        let cpsi_b: f64 = psi_b.cos();
        let spsi_b: f64 = psi_b.sin();
        let ca_new: f64 = cpsi_b*ca + spsi_b*cphi*sa;
        let cb_new: f64 = cpsi_b*cb - spsi_b/sa*(cphi*ca*cb - sphi*cg);
        let cg_new: f64 = cpsi_b*cg - spsi_b/sa*(cphi*ca*cg + sphi*cb);

        let dir_new_mag: f64 = (ca_new*ca_new + cb_new*cb_new + cg_new*cg_new).sqrt();
        let dir_new: Array1<f64> = array![ca_new/dir_new_mag, cb_new/dir_new_mag, cg_new/dir_new_mag];
        particle_2.dir.assign(&dir_new);
    }

    {
        let Za = particle_1.Z;
        let Ma = particle_1.m;
        let Mb = material.m;
        let Zb = material.Z;
        let E = particle_1.E;
        let a = screening_length(Za, Zb);

        let mut stopping_factor = 0.;
        if material.inside(particle_1.pos[0]) {
            let Sel = 1.212*(Za.powf(7./6.)*Zb)/((Za.powf(2./3.) + Zb.powf(2./3.)).powf(3./2.))*(E/Ma*AMU/Q).sqrt();
            stopping_factor = material.number_density(particle_1.pos[0])*Sel*ANGSTROM*ANGSTROM*Q;
        }
        let Enl = ffp*stopping_factor;

        particle_1.E = E - T - Enl;
        if particle_1.E < 0_f64 {
            particle_1.E = 0_f64;
        }

        particle_2.E = T - material.Eb;
        if particle_2.E < 0_f64 {
            particle_2.E = 0_f64;
        }
    }
}

fn pick_collision_partner(particle_1: &Particle, material: &Material) -> (f64, f64, Array1<f64>, Array1<f64>) {

    let mfp: f64 = material.mfp(particle_1.pos[0]);
    let pmax : f64 = mfp/SQRTPI;
    let impact_parameter: f64 = pmax*(rand::random::<f64>()).sqrt();
    let phi_azimuthal: f64 = 2.*PI*rand::random::<f64>();

    let sphi: f64 = phi_azimuthal.sin();
    let ca: f64 = particle_1.dir[0];
    let cb: f64 = particle_1.dir[1];
    let cg: f64 = particle_1.dir[2];
    let sa: f64 = (1_f64 - ca*ca).sqrt();
    let cphi: f64 = phi_azimuthal.cos();

    let dir: Array1<f64> = array![ca, cb, cg];

    let x_recoil: f64 = particle_1.pos[0] + mfp*ca - impact_parameter*cphi*sa;
    let y_recoil: f64 = particle_1.pos[1] + mfp*cb - impact_parameter*(sphi*cg - cphi*cb*ca)/sa;
    let z_recoil: f64 = particle_1.pos[2] + mfp*cg + impact_parameter*(sphi*cb - cphi*ca*cg)/sa;

    let pos: Array1<f64> = array![x_recoil, y_recoil, z_recoil];

    return (impact_parameter, phi_azimuthal, dir, pos);
}

fn surface_boundary_condition(particle_1: &mut Particle, material: &Material) -> bool {

    if !material.inside_energy_barrier(particle_1.pos[0]) & material.inside_energy_barrier(particle_1.pos_old[0]) {
        let leaving_energy = particle_1.E*particle_1.dir[0]*particle_1.dir[0];
        //println!("leaving_energy: {}", leaving_energy/Q);

        if leaving_energy < material.Es {
            particle_1.dir[0] = -particle_1.dir[0];
            particle_1.pos[0] = material.energy_barrier;
            return false;
        } else {
            surface_refraction(particle_1, &material, -1_f64);
            return true;
        }
    } else {
        return false;
    }
}

fn surface_refraction(particle_1: &mut Particle, material: &Material, sign: f64) {
    let Es = material.Es;
    let E0 = particle_1.E;
    let cosx0 = particle_1.dir[0];
    let sinx0 = (1_f64 - cosx0*cosx0).sqrt();
    particle_1.dir[0] = ((E0*cosx0*cosx0 + sign*Es)/(E0 +  sign*Es)).sqrt();
    let sinx = (1_f64 - particle_1.dir[0]*particle_1.dir[0]).sqrt();
    particle_1.dir[1] = particle_1.dir[1]*sinx/sinx0;
    particle_1.dir[2] = particle_1.dir[2]*sinx/sinx0;
    particle_1.E = particle_1.E - material.Es;
}

fn bca(N: usize, E0: f64, theta: f64, Ec: f64, Ma: f64, Za: f64, Mb: f64, Zb: f64, Eb: f64, Es: f64, n: f64) {
    let mut particles: Vec<Particle> = vec![];

    println!("Welcome to RustBCA!");
    println!("N particles: {}", N);

    let material = Material {
        n: n,
        m: Mb,
        Z: Zb,
        Eb: Eb,
        Es: Es,
        surface_position: 0_f64,
        energy_barrier: -2_f64*n.powf(-1_f64/3_f64)/SQRT2PI,
        simulation_boundary: -4.*n.powf(-1_f64/3_f64)/SQRT2PI,
    };

    let x0 = material.energy_barrier;
    let y0 = 0.;
    let z0 = 0.;
    let cosx = (theta*PI/180_f64).cos();
    let sinx = (theta*PI/180_f64).sin();

    for i in 0..N {
        particles.push(
            Particle {
                m: Ma,
                Z: Za,
                E: E0,
                pos: array![x0, y0, z0],
                dir: array![cosx, sinx, 0_f64],
                pos_old: array![x0, y0, z0],
                pos_origin: array![x0, y0, z0],
                t: 0_f64,
                left: false,
                stopped: false,
                incident: true,
                first_step: true,
                track_trajectories: true,
                trajectory: vec![],
            }
        )
    }

    let mut particle_index: usize = 0;
    while particle_index < particles.len() {
        if particle_index % 10000 == 0 {
            println!("particle {} of {}", particle_index, particles.len());
        }
        while !particles[particle_index].stopped & !particles[particle_index].left {

            if !material.inside_simulation_boundary(particles[particle_index].pos[0]) {
                particles[particle_index].left = true;
                continue;
            }

            if particles[particle_index].E < Ec {
                particles[particle_index].stopped = true;
                continue;
            }

            let (impact_parameter, phi_azimuthal, dir, pos) = pick_collision_partner(&particles[particle_index], &material);

            let mut particle_2 = Particle {
                m: Mb,
                Z: Zb,
                E: 0.,
                pos: array![pos[0], pos[1], pos[2]],
                pos_old: array![pos[0], pos[1], pos[2]],
                pos_origin: array![pos[0], pos[1], pos[2]],
                dir: array![dir[0], dir[1], dir[2]],
                t: 0.,
                left: false,
                stopped: false,
                incident: false,
                first_step: false,
                track_trajectories: false,
                trajectory: vec![],
            };

            let (theta, psi, T, t, x0) = binary_collision(&particles[particle_index], &particle_2, &material, impact_parameter, 0.00001, 100);
            update_coordinates(&mut particles[particle_index], &mut particle_2, &material, phi_azimuthal, theta, psi, T, t);
            particles[particle_index].left = surface_boundary_condition(&mut particles[particle_index], &material);

            if particles[particle_index].track_trajectories{
                particles[particle_index].add_trajectory();
            }

            if T > Ec {
                particles.push(particle_2);
            }
        }
        particle_index += 1
    }

    let mut file = OpenOptions::new().write(true).create(true).open("output.dat").unwrap();

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
        }

        if particle.incident & particle.stopped {
            num_deposited += 1;
            range += particle.pos[0];
            writeln!(file, "{}, {}, {}", particle.pos[0], particle.pos[1], particle.pos[2]);
        }

        if !particle.incident & particle.left {
            num_sputtered += 1;
            E_sputtered += particle.E;
        }
    }
    println!("Range: {} R: {} E_r: {} Y: {} E_s: {}", range/(num_deposited as f64)/ANGSTROM, (num_reflected as f64)/(N as f64), E_reflected/Q/(num_reflected as f64), (num_sputtered as f64)/(N as f64), E_sputtered/Q/(num_sputtered as f64));
}


fn main() {
    //fn bca(N: usize, E0: f64, theta: f64, Ec: f64, Ma: f64, Za: f64, Mb: f64, Zb: f64, Eb: f64, Es: f64, n: f64) {
    bca(100000, 1000.*Q, 0.0001, 3.*Q, 4.*AMU, 2., 63.54*AMU, 29., 0.0, 3.52*Q, 8.491E28);
}
