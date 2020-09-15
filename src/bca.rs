use super::*;
use rcpr::chebyshev::*;

pub struct BinaryCollisionGeometry {
    pub phi_azimuthal: f64,
    pub impact_parameter: f64,
    pub mfp: f64
}

impl BinaryCollisionGeometry {
    pub fn new(phi_azimuthal: f64, impact_parameter: f64, mfp: f64) -> BinaryCollisionGeometry {
        BinaryCollisionGeometry {
            phi_azimuthal: phi_azimuthal,
            impact_parameter: impact_parameter,
            mfp: mfp
        }
    }
}

pub struct BinaryCollisionResult {
    pub theta: f64,
    pub psi: f64,
    pub psi_recoil: f64,
    pub recoil_energy: f64,
    pub asymptotic_deflection: f64,
    pub normalized_distance_of_closest_approach: f64
}

impl BinaryCollisionResult {
    pub fn new(theta: f64, psi: f64, psi_recoil: f64, recoil_energy: f64,
        asymptotic_deflection: f64, normalized_distance_of_closest_approach: f64) -> BinaryCollisionResult {
        BinaryCollisionResult {
            theta,
            psi,
            psi_recoil,
            recoil_energy,
            asymptotic_deflection,
            normalized_distance_of_closest_approach
        }
    }
}

fn scattering_integral_mw(x: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //Function for scattering integral - see Mendenhall and Weller, 1991 & 2005
    return (1. - interactions::phi(x, interaction_potential)/x/reduced_energy - beta*beta/x/x).powf(-0.5);
}

fn scattering_function_gl(u: f64, impact_parameter: f64, r0: f64, relative_energy: f64, interaction_potential: &dyn Fn(f64) -> f64) -> Result<f64, anyhow::Error> {
    let result = 4.*impact_parameter*u/(r0*(1. - interaction_potential(r0/(1. - u*u))/relative_energy - impact_parameter*impact_parameter*(1. - u*u).powf(2.)/r0/r0).sqrt());

    if result.is_nan() {
        Err(anyhow!("Numerical error: Gauss-Mehler scattering integrand complex. Likely incorrect distance of closest approach Er = {}, r0 = {} A, p = {} A - check root-finder.",
            relative_energy/EV, r0/ANGSTROM, impact_parameter/ANGSTROM))
    } else {
        Ok(result)
    }
}

fn scattering_function_gm(u: f64, impact_parameter: f64, r0: f64, relative_energy: f64, interaction_potential: &dyn Fn(f64) -> f64) -> Result<f64, anyhow::Error> {
    let result = impact_parameter/r0/(1. - interaction_potential(r0/u)/relative_energy - (impact_parameter*u/r0).powf(2.)).sqrt();

    if result.is_nan() {
        Err(anyhow!("Numerical error: Gauss-Mehler scattering integrand complex. Likely incorrect distance of closest approach Er = {}, r0 = {} A, p = {} A - check root-finder.",
            relative_energy/EV, r0/ANGSTROM, impact_parameter/ANGSTROM))
    } else {
        Ok(result)
    }
}

fn scattering_integral_gauss_mehler(impact_parameter: f64, relative_energy: f64, r0: f64, interaction_potential: &dyn Fn(f64) -> f64, n_points: usize) -> f64 {
    let x: Vec<f64> = (1..=n_points).map(|i| ((2.*i as f64 - 1.)/4./n_points as f64*PI).cos()).collect();
    let w: Vec<f64> = (1..=n_points).map(|i| PI/n_points as f64*((2.*i as f64 - 1.)/4./n_points as f64*PI).sin()).collect();

    PI - x.iter().zip(w).map(|(&x, w)| w*scattering_function_gm(x, impact_parameter, r0, relative_energy, interaction_potential).unwrap()).sum::<f64>()
}

fn scattering_integral_gauss_legendre(impact_parameter: f64, relative_energy: f64, r0: f64, interaction_potential: &dyn Fn(f64) -> f64) -> f64 {
    let x: Vec<f64> = vec![0., -0.538469, 0.538469, -0.90618, 0.90618].iter().map(|x| x/2. + 1./2.).collect();
    let w: Vec<f64> = vec![0.568889, 0.478629, 0.478629, 0.236927, 0.236927].iter().map(|w| w/2.).collect();

    PI - x.iter().zip(w).map(|(&x, w)| w*scattering_function_gl(x, impact_parameter, r0, relative_energy, interaction_potential).unwrap()).sum::<f64>()
}

pub fn determine_mfp_phi_impact_parameter(particle_1: &mut particle::Particle, material: &material::Material, options: &Options) -> Vec<BinaryCollisionGeometry> {

    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let z = particle_1.pos.z;

    let mut mfp = material.mfp(x, y);

    let mut phis_azimuthal = Vec::with_capacity(options.weak_collision_order + 1);
    let mut binary_collision_geometries = Vec::with_capacity(options.weak_collision_order + 1);

    //Each weak collision gets its own aziumuthal angle in annuli around collision point
    //azimuthal angle randomly selected (0..2pi)
    for k in 0..options.weak_collision_order + 1 {
        phis_azimuthal.push(2.*PI*rand::random::<f64>());
    }

    if options.high_energy_free_flight_paths {

        let Ma: f64 = particle_1.m;
        let Mb: f64  = material.average_mass(x, y);
        let Za: f64  = particle_1.Z;
        let Zb: f64  = material.average_Z(x, y);
        let n: &Vec<f64>  = material.number_densities(x, y);
        let E: f64  = particle_1.E;
        let Ec: f64 = particle_1.Ec;
        let a: f64 = interactions::screening_length(Za, Zb, options.interaction_potential);
        let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E;

        //Minimum energy transfer for generating scattering event set to cutoff energy
        let E_min = Ec*(Ma + Mb)*(Ma + Mb)/4./Ma/Mb;
        let reduced_energy_min: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E_min;

        //Free flight path formulation here described in SRIM textbook chapter 7, and Eckstein 1991 7.5.2
        let ep = (reduced_energy*reduced_energy_min).sqrt();
        let mut pmax = a/(ep + ep.sqrt() + 0.125*ep.powf(0.1));
        let mut ffp = 1./(material.total_number_density(x, y)*pmax*pmax*PI);
        let stopping_powers = material.electronic_stopping_cross_sections(particle_1, INTERPOLATED);

        let delta_energy_electronic = stopping_powers.iter().zip(material.number_densities(x, y)).map(|(&i1, i2)| i1*i2).sum::<f64>()*ffp*material.electronic_stopping_correction_factor;

        //If losing too much energy, scale free-flight-path down
        //5 percent limit set in original TRIM paper, Biersack and Haggmark 1980
        if delta_energy_electronic > 0.05*E {
            ffp = 0.05*E/delta_energy_electronic*ffp;
            pmax = (1./(material.total_number_density(x, y)*PI*ffp)).sqrt()
        }

        //If free-flight-path less than the interatomic spacing, revert to solid model
        //Mentioned in Eckstein 1991, Ziegler, Biersack, and Ziegler 2008 (SRIM textbook 7-8)
        if ffp < mfp {

            ffp = mfp;
            //Cylindrical geometry
            pmax = mfp/SQRT2PI;
            let mut impact_parameter = Vec::with_capacity(1);
            let random_number = rand::random::<f64>();
            let p = pmax*random_number.sqrt();
            impact_parameter.push(p);

            //Atomically rough surface - scatter initial collisions using mfp near interface
            if particle_1.first_step {
                ffp = mfp*rand::random::<f64>();
                particle_1.first_step = false;
            }

            if options.mean_free_path_model == GASEOUS {
                ffp *= -rand::random::<f64>().ln();
            }

            binary_collision_geometries.push(BinaryCollisionGeometry::new(phis_azimuthal[0], impact_parameter[0], ffp));
            return binary_collision_geometries;

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

            if options.mean_free_path_model == GASEOUS {
                ffp *= -rand::random::<f64>().ln();
            }

            binary_collision_geometries.push(BinaryCollisionGeometry::new(phis_azimuthal[0], impact_parameter[0], ffp));
            return binary_collision_geometries;
        }

    } else {

        //If not using free flight paths, use weak collision model
        let pmax = mfp/SQRTPI;

        //Cylindrical geometry
        let mut impact_parameters = Vec::with_capacity(options.weak_collision_order + 1);
        for k in 0..(options.weak_collision_order + 1) {
            let random_number = rand::random::<f64>();
            let p = pmax*(random_number + k as f64).sqrt();
            impact_parameters.push(p)
        }

        //Atomically rough surface - scatter initial collisions
        if particle_1.first_step {
            mfp *= rand::random::<f64>();
            particle_1.first_step = false;
        }

        if options.mean_free_path_model == GASEOUS {
            mfp *= -rand::random::<f64>().ln();
        }

        for k in 0..(options.weak_collision_order + 1) {
            binary_collision_geometries.push(BinaryCollisionGeometry::new(phis_azimuthal[k], impact_parameters[k], mfp))
        }
        return binary_collision_geometries;
    }
}

pub fn choose_collision_partner(particle_1: &particle::Particle, material: &material::Material, binary_collision_geometry: &BinaryCollisionGeometry, options: &Options) -> (usize, particle::Particle) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let z = particle_1.pos.z;

    let impact_parameter = binary_collision_geometry.impact_parameter;
    let mfp = binary_collision_geometry.mfp;
    let phi_azimuthal = binary_collision_geometry.phi_azimuthal;

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
    let (species_index, Z_recoil, M_recoil, Ec_recoil, Es_recoil) = material.choose(x_recoil, y_recoil);

    return (species_index,
        particle::Particle::new(
            M_recoil, Z_recoil, 0., Ec_recoil, Es_recoil,
            x_recoil, y_recoil, z_recoil,
            ca, cb, cg,
            false, options.track_recoil_trajectories
        )
    )
}

fn distance_of_closest_approach(particle_1: &particle::Particle, particle_2: &particle::Particle, binary_collision_geometry: &BinaryCollisionGeometry, options: &Options) -> Result<f64, anyhow::Error> {
    let Za: f64 = particle_1.Z;
    let Zb: f64 = particle_2.Z;
    let Ma: f64 = particle_1.m;
    let Mb: f64 = particle_2.m;
    let E0: f64 = particle_1.E;
    let mu: f64 = Mb/(Ma + Mb);

    //Lindhard screening length and reduced energy
    let a: f64 = interactions::screening_length(Za, Zb, options.interaction_potential);

    let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
    let relative_energy = E0*Mb/(Ma + Mb);
    let beta: f64 = binary_collision_geometry.impact_parameter/a;

    let root_finder = if (relative_energy < interactions::energy_threshold_single_root(options.interaction_potential)) {options.root_finder} else {NEWTON};

    let f = |r: f64| -> f64 {interactions::distance_of_closest_approach_function(r, a, Za, Zb, relative_energy, binary_collision_geometry.impact_parameter, options.interaction_potential)};
    let df = |r: f64| -> f64 {interactions::diff_distance_of_closest_approach_function(r, a, Za, Zb, relative_energy, binary_collision_geometry.impact_parameter, options.interaction_potential)};

    let x0 = match root_finder {
        CPR => {

            let g = |r: f64| -> f64 {interactions::distance_of_closest_approach_function_singularity_free(r, a, Za, Zb, relative_energy, binary_collision_geometry.impact_parameter, options.interaction_potential)*interactions::scaling_function(r, a, options.interaction_potential)};

            let upper_bound = f64::max(10.*a*beta, 10.*a);

            let roots = find_roots_with_newton_polishing(&g, &f, &df, 0., upper_bound,
                options.cpr_n0, options.cpr_epsilon, options.cpr_nmax, options.cpr_complex,
                options.cpr_truncation, options.cpr_interval_limit, options.cpr_far_from_zero);

            let max_root = roots.iter().cloned().fold(f64::NAN, f64::max)/a;

            if roots.is_empty() || max_root.is_nan() {
                return Err(anyhow!("Numerical error: CPR rootfinder failed to find root, {}. E: {}; x, y, z: ({},{},{}); x0: {}, F(a): {}, F(b): {};",
                    options.max_iterations, E0/EV,
                    particle_1.pos.x/ANGSTROM, particle_1.pos.y/ANGSTROM, particle_1.pos.z/ANGSTROM,
                    max_root, g(0.), g(upper_bound)));
            } else {
                return Ok(max_root);
            }
        },
        NEWTON => {
            //Guess for large reduced energy from Mendenhall and Weller 1991
            //For small energies, use pure Newton-Raphson with arbitrary guess of 1
            let mut x0 = beta;
            let mut xn: f64;
            if reduced_energy > 5. {
                let inv_er_2 = 0.5/reduced_energy;
                x0 = inv_er_2 + (inv_er_2*inv_er_2 + beta*beta).sqrt();
            }

            //Newton-Raphson to determine distance of closest approach
            let mut err: f64 = options.tolerance + 1.;
            for k in 0..options.max_iterations {
                xn = x0 - f(x0*a)/df(x0*a);
                err = (xn - x0)*(xn - x0);
                x0 = xn;
                if err < options.tolerance {
                    return Ok(x0);
                }
            }
            return Err(anyhow!("Numerical error: exceeded maximum number of Newton-Raphson iterations, {}. E: {}; x, y, z: ({},{},{}); x0: {}; Error: {}; Tolerance: {}",
                options.max_iterations, E0,
                particle_1.pos.x/ANGSTROM, particle_1.pos.y/ANGSTROM, particle_1.pos.z/ANGSTROM,
                x0, err, options.tolerance));
        },
        _ => panic!("Input error: unimplemented root-finder. Choose 0: Newton 1: Adaptive Chebyshev Proxy Rootfinder with Subdivision")
    };
}

pub fn calculate_binary_collision(particle_1: &particle::Particle, particle_2: &particle::Particle, binary_collision_geometry: &BinaryCollisionGeometry, options: &Options) -> BinaryCollisionResult {
    let Za: f64 = particle_1.Z;
    let Zb: f64 = particle_2.Z;
    let Ma: f64 = particle_1.m;
    let Mb: f64 = particle_2.m;
    let E0: f64 = particle_1.E;
    let mu: f64 = Mb/(Ma + Mb);

    let a: f64 = interactions::screening_length(Za, Zb, options.interaction_potential);
    let x0 = distance_of_closest_approach(particle_1, particle_2, binary_collision_geometry, options).unwrap();
    let r0 = x0*a;

    let theta = match  options.scattering_integral {
        MENDENHALL_WELLER => {
            //Lindhard screening length and reduced energy
            let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
            let beta: f64 = binary_collision_geometry.impact_parameter/a;

            //Scattering integral quadrature from Mendenhall and Weller 2005
            let lambda0 = (0.5 + beta*beta/x0/x0/2. - interactions::dphi(x0,  options.interaction_potential)/2./reduced_energy).powf(-1./2.);
            let alpha = 1./12.*(1. + lambda0 + 5.*(0.4206*scattering_integral_mw(x0/0.9072, beta, reduced_energy,  options.interaction_potential) + 0.9072*scattering_integral_mw(x0/0.4206, beta, reduced_energy,  options.interaction_potential)));
            PI*(1. - beta*alpha/x0)
        },
        GAUSS_MEHLER => {
            let V = |r| {interactions::interaction_potential(r, a, Za, Zb, options.interaction_potential)};
            let relative_energy = E0*Mb/(Ma + Mb);
            let impact_parameter = binary_collision_geometry.impact_parameter;
            scattering_integral_gauss_mehler(impact_parameter, relative_energy, r0, &V, 10)
        },
        GAUSS_LEGENDRE => {
            let V = |r| {interactions::interaction_potential(r, a, Za, Zb, options.interaction_potential)};
            let relative_energy = E0*Mb/(Ma + Mb);
            let impact_parameter = binary_collision_geometry.impact_parameter;
            scattering_integral_gauss_legendre(impact_parameter, relative_energy, r0, &V)
        },
        MAGIC => {
            //MAGIC algorithm
            //Since this is legacy code I don't think I will clean this up
            let C_ = match  options.interaction_potential {
                MOLIERE => vec![ 0.6743, 0.009611, 0.005175, 6.314, 10.0 ],
                KR_C => vec![ 0.7887, 0.01166, 00.006913, 17.16, 10.79 ],
                ZBL => vec![ 0.99229, 0.011615, 0.0071222, 9.3066, 14.813 ],
                TRIDYN => vec![1.0144, 0.235809, 0.126, 69350., 83550.], //Undocumented Tridyn constants
                _ => panic!("Input error: unimplemented interaction potential {} for MAGIC algorithm.",  options.interaction_potential)
            };

            let c = vec![ 0.190945, 0.473674, 0.335381, 0.0 ];
            let d = vec![ -0.278544, -0.637174, -1.919249, 0.0 ];

            let V0 = Za*Zb*Q*Q/4.0/PI/EPS0/a;
            let E_c = E0*Mb/(Ma + Mb);
            let E_r = E0/V0;
            let b = binary_collision_geometry.impact_parameter/a;
            let SQE = E_r.sqrt();
            let R = a*x0;
            //let sum = c[0]*(d[0]*x0).exp() + c[1]*(d[1]*x0).exp() + c[2]*(d[2]*x0).exp();
            let sum = interactions::phi(x0, options.interaction_potential);
            let V = V0*a/R*sum;
            //let sum = d[0]*c[0]*(d[0]*x0).exp() + d[1]*c[1]*(d[1]*x0).exp() + d[2]*c[2]*(d[2]*x0).exp();
            let sum = interactions::dphi(x0, options.interaction_potential);
            let dV = -V/R + V0/R*sum; //1174
            let rho = -2.0*(E_c - V)/dV; //1176
            let D = 2.0*(1.0 + C_[0]/SQE)*E_r*b.powf((C_[1] + SQE)/(C_[2] + SQE)); //1179
            let G = (C_[4] + E_r)/(C_[3] + E_r)*((1. + D*D).sqrt() - D); //F-TRIDYN line 1180
            let delta =  D*G/(1.0 + G)*(x0 - b);
            let ctheta2 = (b + rho/a + delta)/(x0 + rho/a);
            2.*((b + rho/a + delta)/(x0 + rho/a)).acos()
        },
        _ => panic!("Input error: unimplemented scattering integral: {}. Use 0: Mendenhall-Weller Quadrature 1: MAGIC Algorithm",
            options.scattering_integral)
    };

    assert!(!theta.is_nan());

    //See Eckstein 1991 for details on center of mass and lab frame angles
    let asympototic_deflection = x0*a*(theta/2.).sin();
    let psi = (theta.sin().atan2(Ma/Mb + theta.cos())).abs();
    let psi_recoil = (theta.sin().atan2(1. - theta.cos())).abs();
    let recoil_energy = 4.*(Ma*Mb)/(Ma + Mb).powf(2.)*E0*(theta/2.).sin().powf(2.);

    BinaryCollisionResult::new(theta, psi, psi_recoil, recoil_energy, asympototic_deflection, x0)
}

pub fn update_particle_energy(particle_1: &mut particle::Particle, material: &material::Material, distance_traveled: f64,
    recoil_energy: f64, xi: f64, strong_collision_Z: f64, strong_collision_index: usize, options: &Options) {

    //If particle energy  drops below zero before electronic stopping calcualtion, it produces NaNs
    particle_1.E -= recoil_energy;
    if particle_1.E < 0. {
        particle_1.E = 0.;
    }

    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let ck = material.electronic_stopping_correction_factor;

    if material.inside(x, y) {

        let electronic_stopping_powers = material.electronic_stopping_cross_sections(particle_1, options.electronic_stopping_mode);
        let n = material.number_densities(x, y);

        let delta_energy = match options.electronic_stopping_mode {
            INTERPOLATED => electronic_stopping_powers.iter().zip(n).map(|(i1, i2)| i1*i2).collect::<Vec<f64>>().iter().sum::<f64>()*distance_traveled*ck,
            LOW_ENERGY_NONLOCAL => electronic_stopping_powers.iter().zip(n).map(|(i1, i2)| i1*i2).collect::<Vec<f64>>().iter().sum::<f64>()*distance_traveled*ck,
            LOW_ENERGY_LOCAL => {
                let Za: f64  = particle_1.Z;
                let Zb: f64 = strong_collision_Z;

                //Oen-Robinson local electronic stopping power
                let a = interactions::screening_length(Za, Zb, options.interaction_potential);

                //d1 is the first (smallest) interior constant of the screening function
                let d1 = interactions::first_screening_radius(options.interaction_potential);

                d1*d1/2./PI*electronic_stopping_powers[strong_collision_index]*(-d1*xi).exp()/a/a*ck
                },
            LOW_ENERGY_EQUIPARTITION => {
                let Za: f64  = particle_1.Z;
                let Zb: f64 = strong_collision_Z;

                //Oen-Robinson local electronic stopping power
                let a = interactions::screening_length(Za, Zb, options.interaction_potential);

                //d1 is the first (smallest) interior constant of the screening function
                let d1 = interactions::first_screening_radius(options.interaction_potential);
                let delta_energy_local = d1*d1/2./PI*electronic_stopping_powers[strong_collision_index]*(-d1*xi).exp()/a/a;
                let delta_energy_nonlocal = electronic_stopping_powers.iter().zip(n).map(|(i1, i2)| i1*i2).collect::<Vec<f64>>().iter().sum::<f64>()*distance_traveled*ck;

                (0.5*delta_energy_local + 0.5*delta_energy_nonlocal)*ck
            },
            //Panic at unimplemented electronic stopping mode
            _ => panic!("Input error: unimplemented electronic stopping mode. Use 0: Biersack-Varelas 1: Lindhard-Scharff 2: Oen-Robinson 3: Equipartition")
        };

        particle_1.E += -delta_energy;
    }

    //Make sure particle energy doesn't become negative again
    if particle_1.E < 0. {
        particle_1.E = 0.;
    }
}

pub fn single_ion_bca(particle: particle::Particle, material: &material::Material, options: &Options) -> Vec<particle::Particle> {

    let mut particles: Vec<particle::Particle> = Vec::new();
    particles.push(particle);

    let mut particle_output: Vec<particle::Particle> = Vec::new();

    let mut particle_index = particles.len();

    'particle_loop: while particle_index > 0 {
        //Remove particle from top of vector as particle_1
        let mut particle_1 = particles.pop().unwrap();

        //BCA loop
        'trajectory_loop: while !particle_1.stopped & !particle_1.left {

            //Choose impact parameters and azimuthal angles for all collisions, and determine mean free path
            let binary_collision_geometries = bca::determine_mfp_phi_impact_parameter(&mut particle_1, &material, &options);

            let mut total_energy_loss = 0.;
            let mut total_asymptotic_deflection = 0.;
            let mut distance_of_closest_approach = 0.;
            let mut strong_collision_Z = 0.;
            let mut strong_collision_index: usize = 0;

            'collision_loop: for k in 0..options.weak_collision_order + 1 {

                let (species_index, mut particle_2) = bca::choose_collision_partner(&particle_1, &material,
                    &binary_collision_geometries[k], &options);

                //If recoil location is inside, proceed with binary collision loop
                if material.inside(particle_2.pos.x, particle_2.pos.y) & material.inside_energy_barrier(particle_1.pos.x, particle_1.pos.y) {

                    //Determine scattering angle from binary collision
                    let binary_collision_result = bca::calculate_binary_collision(&particle_1,
                        &particle_2, &binary_collision_geometries[k], &options);

                    //Only use 0th order collision for local electronic stopping
                    if k == 0 {
                        distance_of_closest_approach = binary_collision_result.normalized_distance_of_closest_approach;
                        strong_collision_Z = particle_2.Z;
                        strong_collision_index = species_index;
                    }

                    //Energy transfer to recoil
                    particle_2.E = binary_collision_result.recoil_energy - material.average_bulk_binding_energy(particle_2.pos.x, particle_2.pos.y);

                    //Accumulate asymptotic deflections for primary particle
                    total_energy_loss += binary_collision_result.recoil_energy;

                    //total_deflection_angle += psi;
                    total_asymptotic_deflection += binary_collision_result.asymptotic_deflection;

                    //Rotate particle 1, 2 by lab frame scattering angles
                    particle::rotate_particle(&mut particle_1, binary_collision_result.psi,
                        binary_collision_geometries[k].phi_azimuthal);

                    particle::rotate_particle(&mut particle_2, -binary_collision_result.psi_recoil,
                        binary_collision_geometries[k].phi_azimuthal);
                    particle_2.dir_old.x = particle_2.dir.x;
                    particle_2.dir_old.y = particle_2.dir.y;
                    particle_2.dir_old.z = particle_2.dir.z;

                    //Only track number of strong collisions, i.e., k = 0
                    if (binary_collision_result.psi > 0.) & (k == 0) {
                        particle_1.number_collision_events += 1;
                    }

                    //Deep recoil suppression
                    //See Eckstein 1991 7.5.3 for recoil suppression function
                    if options.track_recoils & options.suppress_deep_recoils {
                        let E = particle_1.E;
                        let Za: f64 = particle_1.Z;
                        let Zb: f64 = particle_2.Z;

                        let Ma: f64 = particle_1.m;
                        let Mb: f64 = particle_2.m;

                        let n = material.total_number_density(particle_2.pos.x, particle_2.pos.y);
                        let a: f64 = interactions::screening_length(Za, Zb, options.interaction_potential);
                        let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E;
                        let estimated_range_of_recoils = (reduced_energy.powf(0.3) + 0.1).powf(3.)/n/a/a;

                        if let Closest::SinglePoint(p2) = material.closest_point(particle_2.pos.x, particle_2.pos.y) {
                            let dx = p2.x() - particle_2.pos.x;
                            let dy = p2.y() - particle_2.pos.y;
                            let distance_to_surface = (dx*dx + dy*dy).sqrt();

                            if (distance_to_surface < estimated_range_of_recoils) & (particle_2.E > particle_2.Ec) {
                                particles.push(particle_2);
                            }
                        }
                    //If transferred energy > cutoff energy, add recoil to particle vector
                    } else if options.track_recoils & (particle_2.E > particle_2.Ec) {
                        particles.push(particle_2);
                    }
                }
            }

            //Advance particle in space and track total distance traveled
            let distance_traveled = particle::particle_advance(&mut particle_1,
                binary_collision_geometries[0].mfp, total_asymptotic_deflection);

            //Subtract total energy from all simultaneous collisions and electronic stopping
            bca::update_particle_energy(&mut particle_1, &material, distance_traveled,
                total_energy_loss, distance_of_closest_approach, strong_collision_Z,
                strong_collision_index, &options);

            //Check boundary conditions on leaving and stopping
            material::boundary_condition_2D_planar(&mut particle_1, &material);

            //Set particle index to topmost particle
            particle_index = particles.len();
        }
        particle_output.push(particle_1);
    }
    particle_output
}
