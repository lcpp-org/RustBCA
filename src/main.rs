extern crate ndarray;
extern crate geo;

pub use crate::geo::*;
use ndarray::prelude::*;
use geo::algorithm::contains::Contains;
use geo::algorithm::closest_point::ClosestPoint;
//use geo::prelude::*;
//use geo::geo_types::*;
use rand::Rng;
use std::error::Error;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;

const Q: f64 = 1.602E-19;
const PI: f64 = 3.14159;
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
    return beta*beta/x0/x0 - dphi(x0)/reduced_energy + 1.
}

fn f(x: f64, beta: f64, reduced_energy: f64) -> f64 {
    //Function for scattering integral - see Mendenhall and Weller, 1991 & 2005
    return (1_f64 - phi(x)/x/reduced_energy - beta*beta/x/x).powf(-0.5);
}

fn rotate_around_axis_by_pi(vector: &Array1<f64>, axis: &Array1<f64>) -> Array1<f64> {
    //Rotate a vector around an axis by pi (180 degrees) - used for reflection at surface
    let axis_mag = (axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]).sqrt();
    let ux = axis[0]/axis_mag;
    let uy = axis[1]/axis_mag;
    let uz = axis[2]/axis_mag;

    let R = array![
        [2.*ux*ux - 1., 2.*ux*uy, 2.*ux*uz],
        [2.*uy*ux, 2.*uy*uy - 1., 2.*uy*uz],
        [2.*uz*ux, 2.*uz*uy, 2.*uz*uz - 1.]
    ];

    return R.dot(vector);
}

fn binary_collision(particle_1: &Particle, particle_2: &Particle,  material: &Material, impact_parameter: f64, tol: f64, max_iter: i32) -> (f64, f64, f64, f64, f64) {
    if !material.inside(particle_2.pos[0], particle_2.pos[1]) {
        //If the chosen collision partner is not within the target, skip collision
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

        //Guess for large reduced energy from Mendenhall and Weller, 1991
        let mut x0 = 1_f64;
        let mut xn: f64;
        if reduced_energy > 5_f64 {
            let inv_er_2 = 0.5_f64/reduced_energy;
            x0 = inv_er_2 + (inv_er_2.powf(2_f64) + beta.powf(2_f64)).sqrt();
        }

        //Newton-Raphson to determine distance of closest approach
        let mut err: f64;
        for iter in 0..max_iter {
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

        return (theta, psi, T, t, x0);
    }
}

fn update_coordinates(particle_1: &mut Particle, particle_2: &mut Particle, material: &Material, phi_azimuthal: f64, theta: f64, psi: f64, T: f64, t: f64) {
    let mut mfp = material.mfp(particle_1.pos[0], particle_1.pos[1]);

    particle_1.pos_old.assign(&particle_1.pos); //Replace array with other array values

    //TRIDYN-style atomically rough surface - scatters initial collisions by up to 1 mfp
    if particle_1.first_step {
        mfp *= rand::random::<f64>();
        particle_1.first_step = false;
    }

    //path length - must subtract next asymptotic deflection and add previous to obtain correct trajectory
    let ffp: f64 = mfp - t + particle_1.t;

    //There's got to be a better way to do this...
    let new_pos: Array1<f64> = array![particle_1.pos[0] + ffp*particle_1.dir[0], particle_1.pos[1] + ffp*particle_1.dir[1], particle_1.pos[2] + ffp*particle_1.dir[2]];
    particle_1.pos.assign(&new_pos);
    particle_1.t = t;

    let ca: f64 = particle_1.dir[0];
    let cb: f64 = particle_1.dir[1];
    let cg: f64 = particle_1.dir[2];
    let cphi: f64 = phi_azimuthal.cos();
    let sphi: f64 = phi_azimuthal.sin();
    let sa = (1_f64 - ca.powf(2_f64)).sqrt();

    //Update particle 1 direction
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
        let dir_new: Array1<f64> = array![ca_new/dir_new_mag, cb_new/dir_new_mag, cg_new/dir_new_mag];
        particle_2.dir.assign(&dir_new);
    }

    //Update particle energies
    {
        let Za = particle_1.Z;
        let Ma = particle_1.m;
        let Mb = material.m;
        let Zb = material.Z;
        let E = particle_1.E;
        let a = screening_length(Za, Zb);

        let mut stopping_power = 0.;

        //Electronic stopping only inside material
        if material.inside(particle_1.pos[0], particle_1.pos[1]) {
            let n = material.number_density(particle_1.pos[0], particle_1.pos[1]);

            //Bethe-Bloch stopping, as modified by Biersack and Haggmark (1980)
            let v = (2.*E/Ma).sqrt();
            let beta = v/C;

            let mut I0; //Mean excitation potential
            if Zb < 13. {
                I0 = 12. + 7./Zb;
            } else {
                I0 = 9.76 + 58.5*Zb.powf(-1.19);
            }
            let I = Zb*I0*Q;

            let mut B;//Effective empirical shell correction
            if Zb < 3. {
                B = 100.*Za/Zb;
            } else {
                B = 5.;
            }

            let prefactor = 8.0735880E-42*Zb*Za*Za/beta/beta; //This is messy, but this prefactor is... not the easiest to compute
            let eb = 2.*ME*v*v/I;
            let S_BB = prefactor*(eb + 1. + B/eb).ln()*n;

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

fn pick_collision_partner(particle_1: &Particle, material: &Material) -> (f64, f64, Array1<f64>, Array1<f64>) {
    //Liquid model of amorphous material - see Eckstein 1991
    let mfp: f64 = material.mfp(particle_1.pos[0], particle_1.pos[1]);

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
    //Planar energy barrier will refract particle and deflect it towards surface
    if !material.inside_energy_barrier(particle_1.pos[0], particle_1.pos[1]) & material.inside_energy_barrier(particle_1.pos_old[0], particle_1.pos_old[1]) {
        let p = point!(x: particle_1.pos[0], y: particle_1.pos[1]);

        let mut leaving_energy;

        //Find closest point to surface, and use line from particle to point as normal to refract
        if let Closest::SinglePoint(p2) = material.surface.closest_point(&p) {
            let dx = p2.x() - p.x();
            let dy = p2.y() - p.y();
            let dz: f64 = 0.;

            let magnitude = (dx*dx + dy*dy + dz*dz).sqrt();

            //Since direction check is already handled in this function, just take abs() here
            //println!("{}", particle_1.E/Q);
            leaving_energy = (particle_1.E*(dx/magnitude*particle_1.dir[0] + dy/magnitude*particle_1.dir[1] + dz/magnitude*particle_1.dir[2])).abs();
            //println!("{}", leaving_energy/Q)
        } else {
            leaving_energy = particle_1.E;
        }

        if leaving_energy < material.Es {
            let p = point!(x: particle_1.pos[0], y: particle_1.pos[1]);

            //This syntax is a kind of enum matching - closest_point doesn't return the point,
            //Since there are multiple kinds of nearest-point solutions, so it returns an enum
            //That stores the point inside. This code extracts it
            if let Closest::SinglePoint(p2) = material.surface.closest_point(&p) {
                let dx = p2.x() - p.x();
                let dy = p2.y() - p.y();

                let axis = array![dx, dy, 0.];

                particle_1.dir[0] *= -1.;
                particle_1.dir[1] *= -1.;

                let new_dir: Array1<f64> = rotate_around_axis_by_pi(&particle_1.dir, &axis);

                particle_1.dir.assign(&new_dir);
            }
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
    //Deflection code for planar surface refraction
    let Es = material.Es;
    let E0 = particle_1.E;

    let cosx0 = particle_1.dir[0];
    let cosy0 = particle_1.dir[1];
    let cosz0 = particle_1.dir[2];

    let p = point!(x: particle_1.pos[0], y: particle_1.pos[1]);

    if let Closest::SinglePoint(p2) = material.surface.closest_point(&p) {
        let dx = p2.x() - p.x();
        let dy = p2.y() - p.y();
        let dz: f64 = 0.;

        let magnitude = (dx*dx + dy*dy + dz*dz).sqrt();

        let new_cosx = ((cosx0*cosx0 + sign*Es*dx/magnitude)/(E0 + sign*Es)).sqrt();
        let new_cosy = ((cosy0*cosy0 + sign*Es*dy/magnitude)/(E0 + sign*Es)).sqrt();

        let dir_mag = (new_cosx*new_cosx + new_cosy*new_cosy + cosz0*cosz0).sqrt();
        particle_1.dir[0] = new_cosx/dir_mag;
        particle_1.dir[1] = new_cosy/dir_mag;
        particle_1.dir[2] = cosz0/dir_mag;
    }
}

fn bca(N: usize, E0: f64, theta: f64, Ec: f64, Ma: f64, Za: f64, Mb: f64, Zb: f64, Eb: f64, Es: f64, n: f64, track_recoils: bool, track_trajectories: bool, track_recoil_trajectories: bool, write_files: bool, thickness: f64, depth: f64, name: String) {
    let mut particles: Vec<Particle> = vec![];

    //Currently, Material just stores the polygons, which are created on struct initialization
    let dx: f64 = 2_f64*n.powf(-1_f64/3_f64)/SQRT2PI;
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

    for i in 0..N {
        let dy = (2.*rand::random::<f64>() - 1.)*thickness/2.;
        particles.push(
            Particle {
                m: Ma,
                Z: Za,
                E: E0,
                pos: array![x0, y0 + dy , z0],
                dir: array![cosx, sinx, 0_f64],
                pos_old: array![x0, y0 + dy, z0],
                pos_origin: array![x0, y0 + dy, z0],
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
        println!("particle {} of {}", particle_index+1, particles.len());
        while !particles[particle_index].stopped & !particles[particle_index].left {

            if !material.inside_simulation_boundary(particles[particle_index].pos[0], particles[particle_index].pos[1]) {
                particles[particle_index].left = true;
                //println!("Left! {} {} {}", particles[particle_index].pos[0]/ANGSTROM, particles[particle_index].pos[1]/ANGSTROM, particles[particle_index].pos[2]/ANGSTROM);
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
                track_trajectories: track_recoil_trajectories,
                trajectory: vec![],
            };

            let (theta, psi, T, t, x0) = binary_collision(&particles[particle_index], &particle_2, &material, impact_parameter, 1E-3, 100);
            update_coordinates(&mut particles[particle_index], &mut particle_2, &material, phi_azimuthal, theta, psi, T, t);
            surface_boundary_condition(&mut particles[particle_index], &material);

            if particles[particle_index].track_trajectories{
                particles[particle_index].add_trajectory();
            }

            if (T > Ec ) & track_recoils {
                particles.push(particle_2);
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
                    writeln!(reflected_file, "{}, {}, {}, {}, {}, {}, {}, {}", particle.Z, particle.pos[0], particle.pos[1], particle.pos[2], particle.dir[0], particle.dir[1], particle.dir[2], particle.E);
                }
            }

            if particle.incident & particle.stopped {
                num_deposited += 1;
                range += particle.pos[0];
                if write_files {
                    writeln!(deposited_file, "{}, {}, {}, {}", particle.Z, particle.pos[0], particle.pos[1], particle.pos[2]);
                }
            }

            if !particle.incident & particle.left {
                num_sputtered += 1;
                E_sputtered += particle.E;
                if write_files {
                    writeln!(sputtered_file, "{}, {}, {}, {}, {}, {}, {}, {}", particle.Z, particle.pos[0], particle.pos[1], particle.pos[2], particle.dir[0], particle.dir[1], particle.dir[2], particle.E);
                }
            }

            if particle.track_trajectories & write_files {
                for pos in particle.trajectory {
                    writeln!(trajectory_file, "{}, {}, {}, {}", particle.Z, pos[0], pos[1], pos[2]);
                }
            }
        }
        println!("E: {} Range: {} R: {} E_r: {} Y: {} E_s: {}", E0/Q, range/(num_deposited as f64)/ANGSTROM, (num_reflected as f64)/(N as f64), E_reflected/Q/(num_reflected as f64), (num_sputtered as f64)/(N as f64), E_sputtered/Q/(num_sputtered as f64));
    }
}


fn main() {
    let num_energies = 1;
    let thickness = 100.*MICRON;
    let depth = 1000.*MICRON;
    let energies = Array::logspace(10., 6., 4., num_energies);
    for index in 0..num_energies {
        let name: String = index.to_string();
        //track_recoils track_trajectories track_recoil_trajectories write_files
        bca(100000, 1.176E6*Q, 0.0001, 1.*Q, 1.*AMU, 1., 63.54*AMU, 29., 0.0, 3.52*Q, 8.491E28, false, false, false, true, thickness, depth, name);
    }
}
