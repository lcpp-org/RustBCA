use super::*;

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
            theta: theta,
            psi: psi,
            psi_recoil: psi_recoil,
            recoil_energy: recoil_energy,
            asymptotic_deflection: asymptotic_deflection,
            normalized_distance_of_closest_approach: normalized_distance_of_closest_approach
        }
    }
}

fn doca_function(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //Transcendental function to _ distance of closest approach
    return x0 - interactions::phi(x0, interaction_potential)/reduced_energy - beta*beta/x0;
}

fn diff_doca_function(x0: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //First differential of distance of closest approach function for N-R solver
    return beta*beta/x0/x0 - interactions::dphi(x0, interaction_potential)/reduced_energy + 1.
}

fn scattering_function(x: f64, beta: f64, reduced_energy: f64, interaction_potential: i32) -> f64 {
    //Function for scattering integral - see Mendenhall and Weller, 1991 & 2005
    return (1. - interactions::phi(x, interaction_potential)/x/reduced_energy - beta*beta/x/x).powf(-0.5);
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
        let n: Vec<f64>  = material.n.clone();
        let E: f64  = particle_1.E;
        let Ec: f64 = particle_1.Ec;
        let a: f64 = interactions::screening_length(Za, Zb, options.interaction_potential);
        let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E;

        //Minimum energy transfer for generating scattering event set to cutoff energy
        let E_min = Ec*(Ma + Mb).powf(2.)/4./Ma/Mb;
        let reduced_energy_min: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E_min;

        //Free flight path formulation here described in SRIM textbook chapter 7, and Eckstein 1991 7.5.2
        let ep = (reduced_energy*reduced_energy_min).sqrt();
        let mut pmax = a/(ep + ep.sqrt() + 0.125*ep.powf(0.1));
        let mut ffp = 1./(material.total_number_density(x, y)*pmax*pmax*PI);
        let stopping_powers = material.electronic_stopping_cross_sections(particle_1, INTERPOLATED);

        let delta_energy_electronic = stopping_powers.iter().zip(material.n.clone()).map(|(&i1, i2)| i1*i2).sum::<f64>()*ffp*material.electronic_stopping_correction_factor;

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
    let (species_index, Z_recoil, M_recoil, Ec_recoil, Es_recoil) = material.choose(x, y);

    return (species_index,
        particle::Particle::new(
            M_recoil, Z_recoil, 0., Ec_recoil, Es_recoil,
            x_recoil, z_recoil, z_recoil,
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
    let beta: f64 = binary_collision_geometry.impact_parameter/a;

    //Guess for large reduced energy from Mendenhall and Weller 1991
    //For small energies, use pure Newton-Raphson with arbitrary guess of 1
    let mut x0 = 1.;
    let mut xn: f64;
    if reduced_energy > 5. {
        let inv_er_2 = 0.5/reduced_energy;
        x0 = inv_er_2 + (inv_er_2*inv_er_2 + beta*beta).sqrt();
    }

    //Newton-Raphson to determine distance of closest approach
    let mut err: f64 = options.tolerance + 1.;
    for k in 0..options.max_iterations {
        let f = doca_function(x0, beta, reduced_energy,  options.interaction_potential);
        let df = diff_doca_function(x0, beta, reduced_energy, options.interaction_potential);
        xn = x0 - f/df;
        err = (xn - x0).powf(2.);
        x0 = xn;
        if err < options.tolerance {
            return Ok(x0);
        }
    }
    return Err(anyhow!("Numerical error: exceeded maximum number of Newton-Raphson iterations, {}. x0: {}; Error: {}; Tolerance: {}",
        options.max_iterations, x0, err, options.tolerance));
}

pub fn calculate_binary_collision(particle_1: &particle::Particle, particle_2: &particle::Particle, binary_collision_geometry: &BinaryCollisionGeometry, options: &Options) -> BinaryCollisionResult {
    let Za: f64 = particle_1.Z;
    let Zb: f64 = particle_2.Z;
    let Ma: f64 = particle_1.m;
    let Mb: f64 = particle_2.m;
    let E0: f64 = particle_1.E;
    let mu: f64 = Mb/(Ma + Mb);

    //Lindhard screening length and reduced energy
    let a: f64 = interactions::screening_length(Za, Zb, options.interaction_potential);
    let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
    let beta: f64 = binary_collision_geometry.impact_parameter/a;

    let x0 = distance_of_closest_approach(particle_1, particle_2, binary_collision_geometry, options).unwrap();

    let theta = match  options.scattering_integral {
        QUADRATURE => {
            //Scattering integral quadrature from Mendenhall and Weller 2005
            let lambda0 = (0.5 + beta*beta/x0/x0/2. - interactions::dphi(x0,  options.interaction_potential)/2./reduced_energy).powf(-1./2.);
            let alpha = 1./12.*(1. + lambda0 + 5.*(0.4206*scattering_function(x0/0.9072, beta, reduced_energy,  options.interaction_potential) + 0.9072*scattering_function(x0/0.4206, beta, reduced_energy,  options.interaction_potential)));
            PI*(1. - beta*alpha/x0)
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
            let sum = c[0]*(d[0]*x0).exp() + c[1]*(d[1]*x0).exp() + c[2]*(d[2]*x0).exp();
            let V = V0*a/R*sum;
            let sum = d[0]*c[0]*(d[0]*x0).exp() + d[1]*c[1]*(d[1]*x0).exp() + d[2]*c[2]*(d[2]*x0).exp();
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
    particle_1.E = particle_1.E - recoil_energy;
    if particle_1.E < 0. {
        particle_1.E = 0.;
    }

    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let ck = material.electronic_stopping_correction_factor;

    //if material.inside_energy_barrier(x, y) {
    if material.inside_energy_barrier(x, y) {

        let electronic_stopping_powers = material.electronic_stopping_cross_sections(particle_1, options.electronic_stopping_mode);
        let n = material.n.clone();

        let delta_energy = match options.electronic_stopping_mode {
            INTERPOLATED => electronic_stopping_powers.iter().zip(n).map(|(i1, i2)| i1*i2).collect::<Vec<f64>>().iter().sum::<f64>()*distance_traveled*ck,

            //v1.iter().zip(v2).map(|(&i1, &i2)| i1 * i2).collect()

            LOW_ENERGY_NONLOCAL => electronic_stopping_powers.iter().zip(n).map(|(i1, i2)| i1*i2).collect::<Vec<f64>>().iter().sum::<f64>()*distance_traveled*ck,
            LOW_ENERGY_LOCAL => {
                let Za: f64  = particle_1.Z;
                let Zb: f64 = strong_collision_Z;

                //Oen-Robinson local electronic stopping power
                let a = interactions::screening_length(Za, Zb, options.interaction_potential);

                //d1 is the first (smallest) interior constant of the screening function
                let d1 = match options.interaction_potential {
                    MOLIERE => 0.3,
                    KR_C => 0.278544,
                    ZBL => 0.20162,
                    LENZ_JENSEN => 0.206,
                    TRIDYN => 0.278544,
                    _ => panic!("Input error: unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
                };

                d1*d1/2./PI*electronic_stopping_powers[strong_collision_index]*(-d1*xi).exp()/a/a*ck
                },
            LOW_ENERGY_EQUIPARTITION => {
                let Za: f64  = particle_1.Z;
                let Zb: f64 = strong_collision_Z;

                //Oen-Robinson local electronic stopping power
                let a = interactions::screening_length(Za, Zb, options.interaction_potential);
                //d1 is the first (smallest) interior constant of the screening function
                let d1 = match options.interaction_potential {
                    MOLIERE => 0.3,
                    KR_C => 0.278544,
                    ZBL => 0.20162,
                    LENZ_JENSEN => 0.206,
                    TRIDYN => 0.278544,
                    _ => panic!("Input error: unimplemented interaction potential. Use 0: MOLIERE 1: KR_C 2: ZBL")
                };
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
