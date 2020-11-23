use super::*;
use geo::{Closest};

#[cfg(any(feature = "cpr_rootfinder_openblas", feature = "cpr_rootfinder_netlib", feature = "cpr_rootfinder_intel_mkl"))]
use rcpr::chebyshev::*;

/// Geometrical quantities of binary collision.
pub struct BinaryCollisionGeometry {
    pub phi_azimuthal: f64,
    pub impact_parameter: f64,
    pub mfp: f64
}

impl BinaryCollisionGeometry {
    /// Constructs a new binary collision geometry object.
    pub fn new(phi_azimuthal: f64, impact_parameter: f64, mfp: f64) -> BinaryCollisionGeometry {
        BinaryCollisionGeometry {
            phi_azimuthal,
            impact_parameter,
            mfp
        }
    }
}

impl fmt::Display for BinaryCollisionGeometry {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Binary collision geometry. \n phi_azimuthal = {} \n p = {} \n mfp = {} \n",
            self.phi_azimuthal, self.impact_parameter, self.mfp)
    }
}

/// Resultant quantities of a binary collision.
pub struct BinaryCollisionResult {
    pub theta: f64,
    pub psi: f64,
    pub psi_recoil: f64,
    pub recoil_energy: f64,
    pub asymptotic_deflection: f64,
    pub normalized_distance_of_closest_approach: f64
}

impl BinaryCollisionResult {

    /// Constructs a new binary collision result object.
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

impl fmt::Display for BinaryCollisionResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Binary collision result. \n theta = {} \n psi = {} \n psi_r = {} \n E_recoil = {} \n tau = {} \n x0 = {} \n",
            self.theta, self.psi, self.psi_recoil, self.recoil_energy, self.asymptotic_deflection, self.normalized_distance_of_closest_approach)
    }
}

/// This function takes a single particle, a material, and an options object and runs a binary-collision-approximation trajectory for that particle in that material, producing a final particle list that consists of the original ion and any material particles displaced thereby.
pub fn single_ion_bca(particle: particle::Particle, material: &material::Material, options: &Options) -> Vec<particle::Particle> {

    let mut particles: Vec<particle::Particle> = Vec::new();
    particles.push(particle);

    let mut particle_output: Vec<particle::Particle> = Vec::new();

    let mut particle_index = particles.len();

    while particle_index > 0 {

        //Remove particle from top of vector as particle_1
        let mut particle_1 = particles.pop().unwrap();

        //BCA loop
        while !particle_1.stopped & !particle_1.left {

            //Choose impact parameters and azimuthal angles for all collisions, and determine mean free path
            let binary_collision_geometries = bca::determine_mfp_phi_impact_parameter(&mut particle_1, &material, &options);

            let mut total_energy_loss = 0.;
            let mut total_asymptotic_deflection = 0.;
            let mut normalized_distance_of_closest_approach = 0.;
            let mut strong_collision_Z = 0.;
            let mut strong_collision_index: usize = 0;

            //Collision loop
            for (k, binary_collision_geometry) in binary_collision_geometries.iter().enumerate().take(options.weak_collision_order + 1) {

                let (species_index, mut particle_2) = bca::choose_collision_partner(&particle_1, &material,
                    &binary_collision_geometry, &options);

                //If recoil location is inside, proceed with binary collision loop
                if material.inside(particle_2.pos.x, particle_2.pos.y) & material.inside_energy_barrier(particle_1.pos.x, particle_1.pos.y) {

                    //Determine scattering angle from binary collision
                    let binary_collision_result = bca::calculate_binary_collision(&particle_1,
                        &particle_2, &binary_collision_geometry, &options)
                        .with_context(|| format!("Numerical error: binary collision calculation failed at x = {} y = {} with {}",
                            particle_1.pos.x, particle_2.pos.x, &binary_collision_geometry))
                        .unwrap();

                    //Only use 0th order collision for local electronic stopping
                    if k == 0 {
                        normalized_distance_of_closest_approach = binary_collision_result.normalized_distance_of_closest_approach;
                        strong_collision_Z = particle_2.Z;
                        strong_collision_index = species_index;
                    }

                    //Energy transfer to recoil
                    particle_2.E = binary_collision_result.recoil_energy - material.average_bulk_binding_energy(particle_2.pos.x, particle_2.pos.y);
                    particle_2.energy_origin = particle_2.E;

                    //Accumulate asymptotic deflections for primary particle
                    total_energy_loss += binary_collision_result.recoil_energy;

                    //total_deflection_angle += psi;
                    total_asymptotic_deflection += binary_collision_result.asymptotic_deflection;

                    //Rotate particle 1, 2 by lab frame scattering angles
                    particle::rotate_particle(&mut particle_1, binary_collision_result.psi,
                        binary_collision_geometry.phi_azimuthal);

                    particle::rotate_particle(&mut particle_2, -binary_collision_result.psi_recoil,
                        binary_collision_geometry.phi_azimuthal);
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
                        //We just need the lindhard screening length here, so the particular potential is not important
                        let a: f64 = interactions::screening_length(Za, Zb, InteractionPotential::MOLIERE);
                        let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E;
                        let estimated_range_of_recoils = (reduced_energy.powf(0.3) + 0.1).powf(3.)/n/a/a;

                        if let Closest::SinglePoint(p2) = material.closest_point(particle_2.pos.x, particle_2.pos.y) {
                            let dx = p2.x() - particle_2.pos.x;
                            let dy = p2.y() - particle_2.pos.y;
                            let distance_to_surface = (dx*dx + dy*dy).sqrt();

                            if (distance_to_surface < estimated_range_of_recoils) & (particle_2.E > particle_2.Ec) {
                                particles.push(particle_2);
                            }
                        } else {
                            panic!("Numerical error: geometry algorithm failed to find distance from particle to surface.")
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
                total_energy_loss, normalized_distance_of_closest_approach, strong_collision_Z,
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

/// For a particle in a material, determine the mean free path and choose the azimuthal angle and impact parameter.
/// The mean free path can be exponentially distributed (gaseous) or constant (amorphous solid/liquid). Azimuthal angles are chosen uniformly. Impact parameters are chosen for collision partners distributed uniformly on a disk of density-dependent radius.
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
        let ck: f64 = material.electronic_stopping_correction_factor(x, y);
        let E: f64  = particle_1.E;
        let Ec: f64 = particle_1.Ec;
        //We just need the Lindhard screening length here, so the particular potential is not important
        let a: f64 = interactions::screening_length(Za, Zb, InteractionPotential::MOLIERE);
        let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E;

        //Minimum energy transfer for generating scattering event set to cutoff energy
        let E_min = Ec*(Ma + Mb)*(Ma + Mb)/4./Ma/Mb;
        let reduced_energy_min: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E_min;

        //Free flight path formulation here described in SRIM textbook chapter 7, and Eckstein 1991 7.5.2
        let ep = (reduced_energy*reduced_energy_min).sqrt();
        let mut pmax = a/(ep + ep.sqrt() + 0.125*ep.powf(0.1));
        let mut ffp = 1./(material.total_number_density(x, y)*pmax*pmax*PI);
        let stopping_powers = material.electronic_stopping_cross_sections(particle_1, ElectronicStoppingMode::INTERPOLATED);

        let delta_energy_electronic = stopping_powers.iter().zip(material.number_densities(x, y)).map(|(&se, number_density)| se*number_density).sum::<f64>()*ffp*ck;

        //If losing too much energy, scale free-flight-path down
        //5 percent limit set in original TRIM paper, Biersack and Haggmark 1980
        if delta_energy_electronic > 0.05*E {
            ffp *= 0.05*E/delta_energy_electronic;
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

            if options.mean_free_path_model == MeanFreePathModel::GASEOUS {
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

            if options.mean_free_path_model == MeanFreePathModel::GASEOUS {
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

        if options.mean_free_path_model == MeanFreePathModel::GASEOUS {
            mfp *= -rand::random::<f64>().ln();
        }

        for k in 0..(options.weak_collision_order + 1) {
            binary_collision_geometries.push(BinaryCollisionGeometry::new(phis_azimuthal[k], impact_parameters[k], mfp))
        }
        return binary_collision_geometries;
    }
}

/// For a particle in a material, and for a particular binary collision geometry, choose a species for the collision partner.
pub fn choose_collision_partner(particle_1: &particle::Particle, material: &material::Material, binary_collision_geometry: &BinaryCollisionGeometry, options: &Options) -> (usize, particle::Particle) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let z = particle_1.pos.z;

    let impact_parameter = binary_collision_geometry.impact_parameter;
    let mfp = binary_collision_geometry.mfp;
    let phi_azimuthal = binary_collision_geometry.phi_azimuthal;

    //Determine cosines and sines
    let sinphi: f64 = phi_azimuthal.sin();
    let cosx: f64 = particle_1.dir.x;
    let cosy: f64 = particle_1.dir.y;
    let cosz: f64 = particle_1.dir.z;
    let sinx: f64 = (1. - cosx*cosx).sqrt();
    let cosphi: f64 = phi_azimuthal.cos();

    //Find recoil location
    let x_recoil: f64 = x + mfp*cosx - impact_parameter*cosphi*sinx;
    let y_recoil: f64 = y + mfp*cosy - impact_parameter*(sinphi*cosz - cosphi*cosy*cosx)/sinx;
    let z_recoil: f64 = z + mfp*cosz + impact_parameter*(sinphi*cosy - cosphi*cosx*cosz)/sinx;

    //Choose recoil Z, M
    let (species_index, Z_recoil, M_recoil, Ec_recoil, Es_recoil, interaction_index) = material.choose(x_recoil, y_recoil);

    return (species_index,
        particle::Particle::new(
            M_recoil, Z_recoil, 0., Ec_recoil, Es_recoil,
            x_recoil, y_recoil, z_recoil,
            cosx, cosy, cosz,
            false, options.track_recoil_trajectories, interaction_index
        )
    )
}

/// Calculate the distance of closest approach of two particles given a particular binary collision geometry.
fn distance_of_closest_approach(particle_1: &particle::Particle, particle_2: &particle::Particle, binary_collision_geometry: &BinaryCollisionGeometry, options: &Options) -> f64 {
    let Za: f64 = particle_1.Z;
    let Zb: f64 = particle_2.Z;
    let Ma: f64 = particle_1.m;
    let Mb: f64 = particle_2.m;
    let E0: f64 = particle_1.E;
    let relative_energy = E0*Mb/(Ma + Mb);
    let p = binary_collision_geometry.impact_parameter;

    let interaction_potential = options.interaction_potential[particle_1.interaction_index][particle_2.interaction_index];

    if let InteractionPotential::COULOMB{Za: Z1, Zb: Z2} = interaction_potential {
        let doca = Z1*Z2*Q*Q/relative_energy/PI/EPS0/8. + (64.*(relative_energy*PI*p*EPS0).powf(2.) + (Z1*Z2*Q*Q).powf(2.)).sqrt()/relative_energy/PI/EPS0/8.;
        return doca/interactions::screening_length(Z1, Z2, interaction_potential);
    }

    let root_finder = if relative_energy < interactions::energy_threshold_single_root(interaction_potential) {
            options.root_finder[particle_1.interaction_index][particle_2.interaction_index]
        } else {Rootfinder::NEWTON{max_iterations: 100, tolerance: 1E-6}};

    #[cfg(any(feature = "cpr_rootfinder_openblas", feature = "cpr_rootfinder_netlib", feature = "cpr_rootfinder_intel_mkl"))]
    match root_finder {
        Rootfinder::POLYNOMIAL{complex_threshold} => polynomial_rootfinder(Za, Zb, Ma, Mb, E0, p, interaction_potential, complex_threshold)
            .with_context(|| format!("Numerical error: CPR rootfinder failed for {} at {} eV with p = {} A.", interaction_potential, E0/EV, p/ANGSTROM))
            .unwrap(),
        Rootfinder::CPR{n0, nmax, epsilon, complex_threshold, truncation_threshold, far_from_zero, interval_limit, derivative_free} =>
            cpr_rootfinder(Za, Zb, Ma, Mb, E0, p, interaction_potential, n0, nmax, epsilon, complex_threshold, truncation_threshold, far_from_zero, interval_limit, derivative_free)
            .with_context(|| format!("Numerical error: CPR rootfinder failed for {} at {} eV with p = {} A.", interaction_potential, E0/EV, p/ANGSTROM))
            .unwrap(),
        Rootfinder::NEWTON{max_iterations, tolerance} => newton_rootfinder(Za, Zb, Ma, Mb, E0, p, interaction_potential, max_iterations, tolerance)
            .with_context(|| format!("Numerical error: CPR rootfinder failed for {} at {} eV with p = {} A.", interaction_potential, E0/EV, p/ANGSTROM))
            .unwrap(),
    }

    #[cfg(not(any(feature = "cpr_rootfinder_openblas", feature = "cpr_rootfinder_netlib", feature = "cpr_rootfinder_intel_mkl")))]
    match root_finder {
        Rootfinder::NEWTON{max_iterations, tolerance} => newton_rootfinder(Za, Zb, Ma, Mb, E0, p, interaction_potential, max_iterations, tolerance)
            .with_context(|| format!("Numerical error: CPR rootfinder failed for {} at {} eV with p = {} A.", interaction_potential, E0/EV, p/ANGSTROM))
            .unwrap(),
        _ => panic!("Input error: unimplemented root-finder. Choose NEWTON or build with cpr_rootfinder to enable CPR and POLYNOMIAL")
    }
}

/// Update a particle's energy following nuclear and electronic interactions of a single BCA step.
pub fn update_particle_energy(particle_1: &mut particle::Particle, material: &material::Material, distance_traveled: f64,
    recoil_energy: f64, x0: f64, strong_collision_Z: f64, strong_collision_index: usize, options: &Options) {

    //If particle energy  drops below zero before electronic stopping calcualtion, it produces NaNs
    particle_1.E -= recoil_energy;
    assert!(!particle_1.E.is_nan(), "Numerical error: particle energy is NaN following collision.");
    if particle_1.E < 0. {
        particle_1.E = 0.;
    }


    let x = particle_1.pos.x;
    let y = particle_1.pos.y;

    if material.inside(x, y) {

        let interaction_potential = options.interaction_potential[particle_1.interaction_index][material.interaction_index[strong_collision_index]];
        let electronic_stopping_powers = material.electronic_stopping_cross_sections(particle_1, options.electronic_stopping_mode);
        let n = material.number_densities(x, y);

        let delta_energy = match options.electronic_stopping_mode {
            ElectronicStoppingMode::INTERPOLATED => electronic_stopping_powers.iter().zip(n).map(|(se, number_density)| se*number_density).collect::<Vec<f64>>().iter().sum::<f64>()*distance_traveled,
            ElectronicStoppingMode::LOW_ENERGY_NONLOCAL => electronic_stopping_powers.iter().zip(n).map(|(se, number_density)| se*number_density).collect::<Vec<f64>>().iter().sum::<f64>()*distance_traveled,
            ElectronicStoppingMode::LOW_ENERGY_LOCAL => oen_robinson_loss(particle_1.Z, strong_collision_Z, electronic_stopping_powers[strong_collision_index], x0, interaction_potential),
            ElectronicStoppingMode::LOW_ENERGY_EQUIPARTITION => {

                let delta_energy_local = oen_robinson_loss(particle_1.Z, strong_collision_Z, electronic_stopping_powers[strong_collision_index], x0, interaction_potential);
                let delta_energy_nonlocal = electronic_stopping_powers.iter().zip(n).map(|(se, number_density)| se*number_density).collect::<Vec<f64>>().iter().sum::<f64>()*distance_traveled;

                (0.5*delta_energy_local + 0.5*delta_energy_nonlocal)
            },
        };

        particle_1.E += -delta_energy;
        //Make sure particle energy doesn't become negative again
        assert!(!particle_1.E.is_nan(), "Numerical error: particle energy is NaN following electronic stopping.");
        if particle_1.E < 0. {
            particle_1.E = 0.;
        }

        particle_1.energy_loss(&options, recoil_energy, delta_energy);
    } else if recoil_energy > 0. {
        particle_1.energy_loss(&options, recoil_energy, 0.);
    }
}

/// Calculate the binary collision result of two particles for a given binary collision geometry. If the calculation fails, return an Error.
pub fn calculate_binary_collision(particle_1: &particle::Particle, particle_2: &particle::Particle, binary_collision_geometry: &BinaryCollisionGeometry, options: &Options) -> Result<BinaryCollisionResult, anyhow::Error> {
    let Za: f64 = particle_1.Z;
    let Zb: f64 = particle_2.Z;
    let Ma: f64 = particle_1.m;
    let Mb: f64 = particle_2.m;
    let E0: f64 = particle_1.E;
    let mu: f64 = Mb/(Ma + Mb);

    let interaction_potential = options.interaction_potential[particle_1.interaction_index][particle_2.interaction_index];
    let scattering_integral = options.scattering_integral[particle_1.interaction_index][particle_2.interaction_index];

    let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
    let x0 = distance_of_closest_approach(particle_1, particle_2, binary_collision_geometry, options);

    let theta: f64 = match  scattering_integral {
        ScatteringIntegral::MENDENHALL_WELLER => mendenhall_weller(Za, Zb, Ma, Mb, E0, binary_collision_geometry.impact_parameter, x0, interaction_potential),
        ScatteringIntegral::GAUSS_MEHLER{n_points} => gauss_mehler(Za, Zb, Ma, Mb, E0, binary_collision_geometry.impact_parameter, x0, interaction_potential, n_points),
        ScatteringIntegral::GAUSS_LEGENDRE => gauss_legendre(Za, Zb, Ma, Mb, E0, binary_collision_geometry.impact_parameter, x0, interaction_potential),
        ScatteringIntegral::MAGIC => magic(Za, Zb, Ma, Mb, E0, binary_collision_geometry.impact_parameter, x0, interaction_potential),
    };

    if theta.is_nan() {
        return Err(anyhow!("Numerical error: CoM deflection angle is NaN for {}. Check input parameters.", binary_collision_geometry));
    }

    //See Eckstein 1991 for details on center of mass and lab frame angles
    let asymptotic_deflection = match interaction_potential {
        InteractionPotential::COULOMB{..} => 0.,
        _ => x0*a*(theta/2.).sin()
    };
    let psi = (theta.sin().atan2(Ma/Mb + theta.cos())).abs();
    let psi_recoil = (theta.sin().atan2(1. - theta.cos())).abs();
    let recoil_energy = 4.*(Ma*Mb)/(Ma + Mb).powf(2.)*E0*(theta/2.).sin().powf(2.);

    Ok(BinaryCollisionResult::new(theta, psi, psi_recoil, recoil_energy, asymptotic_deflection, x0))
}


/// Mendenhall-Weller scattering integrand.
fn scattering_integral_mw(x: f64, beta: f64, reduced_energy: f64, interaction_potential: InteractionPotential) -> f64 {
    //Function for scattering integral - see Mendenhall and Weller, 1991 & 2005
    return (1. - interactions::phi(x, interaction_potential)/x/reduced_energy - beta*beta/x/x).powf(-0.5);
}

/// Gauss-Legendre scattering integrand.
fn scattering_function_gl(u: f64, impact_parameter: f64, r0: f64, relative_energy: f64, interaction_potential: &dyn Fn(f64) -> f64) -> Result<f64, anyhow::Error> {
    let result = 4.*impact_parameter*u/(r0*(1. - interaction_potential(r0/(1. - u*u))/relative_energy - impact_parameter*impact_parameter*(1. - u*u).powf(2.)/r0/r0).sqrt());

    if result.is_nan() {
        Err(anyhow!("Numerical error: Gauss-Legendre scattering integrand complex. Likely incorrect distance of closest approach Er = {}, r0 = {} A, p = {} A - check root-finder.",
            relative_energy/EV, r0/ANGSTROM, impact_parameter/ANGSTROM))
    } else {
        Ok(result)
    }
}

/// Gauss-Mehler scattering integrand.
fn scattering_function_gm(u: f64, impact_parameter: f64, r0: f64, relative_energy: f64, interaction_potential: &dyn Fn(f64) -> f64) -> Result<f64, anyhow::Error> {
    let result = impact_parameter/r0/(1. - interaction_potential(r0/u)/relative_energy - (impact_parameter*u/r0).powf(2.)).sqrt();

    if result.is_nan() {
        Err(anyhow!("Numerical error: Gauss-Mehler scattering integrand complex. Likely incorrect distance of closest approach Er = {} eV r0 = {} A, p = {} A - check root-finder.",
            relative_energy/EV, r0/ANGSTROM, impact_parameter/ANGSTROM))
    } else {
        Ok(result)
    }
}

/// Compute the scattering integral for a given relative energy, distance of closest approach `r0`,  and interaction potential using a Gauss-Mehler, n-point quadrature.
fn scattering_integral_gauss_mehler(impact_parameter: f64, relative_energy: f64, r0: f64, interaction_potential: &dyn Fn(f64) -> f64, n_points: usize) -> f64 {
    let x: Vec<f64> = (1..=n_points).map(|i| ((2.*i as f64 - 1.)/4./n_points as f64*PI).cos()).collect();
    let w: Vec<f64> = (1..=n_points).map(|i| PI/n_points as f64*((2.*i as f64 - 1.)/4./n_points as f64*PI).sin()).collect();

    PI - x.iter().zip(w)
        .map(|(&x, w)| w*scattering_function_gm(x, impact_parameter, r0, relative_energy, interaction_potential)
        .with_context(|| format!("Numerical error: NaN in Gauss-Mehler scattering integral at x = {} with Er = {} eV and p = {} A.", x, relative_energy/EV, impact_parameter/ANGSTROM))
        .unwrap()).sum::<f64>()
}

/// Compute the scattering integral for a given relative energy, distance of closest approach `r0`,  and interaction potential using a Gauss-Legendre, 5-point quadrature.
fn scattering_integral_gauss_legendre(impact_parameter: f64, relative_energy: f64, r0: f64, interaction_potential: &dyn Fn(f64) -> f64) -> f64 {
    let x: Vec<f64> = vec![0., -0.538469, 0.538469, -0.90618, 0.90618].iter().map(|x| x/2. + 1./2.).collect();
    let w: Vec<f64> = vec![0.568889, 0.478629, 0.478629, 0.236927, 0.236927].iter().map(|w| w/2.).collect();

    PI - x.iter().zip(w)
        .map(|(&x, w)| w*scattering_function_gl(x, impact_parameter, r0, relative_energy, interaction_potential)
        .with_context(|| format!("Numerical error: NaN in Gauss-Legendre scattering integral at x = {} with Er = {} eV and p = {} A.", x, relative_energy/EV, impact_parameter/ANGSTROM))
        .unwrap()).sum::<f64>()
}

#[cfg(any(feature = "cpr_rootfinder_openblas", feature = "cpr_rootfinder_netlib", feature = "cpr_rootfinder_intel_mkl"))]
/// Computes the distance of closest approach of two particles with atomic numbers `Za`, `Zb` and masses `Ma`, `Mb` for an inverse-polynomial interaction potential (e.g., Lennard-Jones) for a given impact parameter and incident energy `E0`.
/// Slightly complex roots with an imaginary part smaller than `polynom_complex_threshold` are considered real roots.
pub fn polynomial_rootfinder(Za: f64, Zb: f64, Ma: f64, Mb: f64, E0: f64, impact_parameter: f64,
    interaction_potential: InteractionPotential, polynom_complex_threshold: f64) -> Result<f64, anyhow::Error> {

    let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
    let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
    let relative_energy = E0*Mb/(Ma + Mb);

    let coefficients = interactions::polynomial_coefficients(relative_energy, impact_parameter, interaction_potential);
    let roots = real_polynomial_roots(coefficients.clone(), polynom_complex_threshold).unwrap();
    let max_root = roots.iter().cloned().fold(f64::NAN, f64::max);
    let inverse_transformed_root = interactions::inverse_transform(max_root, interaction_potential);

    if roots.is_empty() || inverse_transformed_root.is_nan() {
        return Err(anyhow!("Numerical error: polynomial rootfinder failed to find root with coefficients {:?}",
            coefficients))
    } else {
        return Ok(inverse_transformed_root/a)
    }
}

#[cfg(any(feature = "cpr_rootfinder_openblas", feature = "cpr_rootfinder_netlib", feature = "cpr_rootfinder_intel_mkl"))]
/// Computes the distance of closest approach of two particles with atomic numbers `Za`, `Zb` and masses `Ma`, `Mb` for an arbitrary interaction potential (e.g., Morse) for a given impact parameter and incident energy `E0` using the Chebyshev-Proxy Root-Finder method.
///
/// # Args:
///
/// `Za`, `Zb`, `Ma`, `Mb`: atomic numbers and masses of the two particles
/// `E0`: incident energy of particle a in the lab frame.
/// `impact_parameter`: the impact parameter between particles a and b.
/// `interaction_potentaial`: the interaction potential between particles a and b.
/// `n0`: initial degree of Chebyshev interpolants.
/// `nmax`: maximum degree of Chebyshev interpolants.
/// `epsilon`: absolute tolerance of Chebyshev interpolant.
/// `complex_threshold`: slightly-complex roots with an imaginary part below this value are considered real.
/// `far_from_zero`: if the distance of closest approach function, evaluated over an interval [a, b] on the Lobatto grid, is always greater than this value, it is assumed that there are no roots in the interval [a, b].
/// `truncation_threshold`: trailing terms of the Chebyshev interpolant with coefficients smaller than this value are ignored.
/// `interval_limit`: if subdivision produces an interval smaller than this value, the root-finder will panic.
/// `derivative_free`: if false, use Newton's method to polish roots from the CPR. If true, use the secant method.
///
/// # Returns
/// Returns the distance of closest approach or an error if the root-finder failed.
pub fn cpr_rootfinder(Za: f64, Zb: f64, Ma: f64, Mb: f64, E0: f64, impact_parameter: f64,
    interaction_potential: InteractionPotential, n0: usize, nmax: usize, epsilon: f64,
    complex_threshold: f64, truncation_threshold: f64, far_from_zero: f64,
    interval_limit: f64, derivative_free: bool) -> Result <f64, anyhow::Error> {

    //Lindhard screening length and reduced energy
    let a = interactions::screening_length(Za, Zb, interaction_potential);
    let reduced_energy = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
    let relative_energy = E0*Mb/(Ma + Mb);
    let p = impact_parameter;

    let f = |r: f64| -> f64 {interactions::distance_of_closest_approach_function(r, a, Za, Zb, relative_energy, impact_parameter, interaction_potential)};
    let g = |r: f64| -> f64 {interactions::distance_of_closest_approach_function_singularity_free(r, a, Za, Zb, relative_energy, impact_parameter, interaction_potential)*
        interactions::scaling_function(r, impact_parameter, interaction_potential)};

    //Using upper bound const, ~10, construct upper bound as a plateau near 0 and linear increase away from that
    //let upper_bound = f64::max(upper_bound_const*p, upper_bound_const*a);
    let upper_bound = impact_parameter + interactions::crossing_point_doca(interaction_potential);

    let roots = match derivative_free {
        true => find_roots_with_secant_polishing(&g, &f, 0., upper_bound,
            n0, epsilon, nmax, complex_threshold,
            truncation_threshold, interval_limit, far_from_zero),

        false => {
            let df = |r: f64| -> f64 {interactions::diff_distance_of_closest_approach_function(r, a, Za, Zb, relative_energy, impact_parameter, interaction_potential)};
            find_roots_with_newton_polishing(&g, &f, &df, 0., upper_bound,
            n0, epsilon, nmax, complex_threshold,
            truncation_threshold, interval_limit, far_from_zero)
        }
    }.with_context(|| format!("Numerical error: CPR Rootfinder failed to converge when calculating distance of closest approach for Er = {} eV p = {} A using {}.",
        relative_energy/EV, impact_parameter/ANGSTROM, interaction_potential))
    .unwrap();

    let max_root = roots.iter().cloned().fold(f64::NAN, f64::max)/a;

    if roots.is_empty() || max_root.is_nan() {
        return Err(anyhow!("Numerical error: CPR rootfinder failed to find root. x0: {}, F(a): {}, F(b): {};", max_root, g(0.), g(upper_bound)));
    } else {
        return Ok(max_root);
    }
}

/// Use Newton's method to find the distance of closest approach between two particles a and b.
pub fn newton_rootfinder(Za: f64, Zb: f64, Ma: f64, Mb: f64, E0: f64, impact_parameter: f64,
    interaction_potential: InteractionPotential, max_iterations: usize, tolerance: f64) -> Result <f64, anyhow::Error> {

    //Lindhard screening length and reduced energy
    let a = interactions::screening_length(Za, Zb, interaction_potential);
    let reduced_energy = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
    let relative_energy = E0*Mb/(Ma + Mb);
    let beta = impact_parameter/a;

    let f = |r: f64| -> f64 {interactions::distance_of_closest_approach_function(r, a, Za, Zb, relative_energy, impact_parameter, interaction_potential)};
    let df = |r: f64| -> f64 {interactions::diff_distance_of_closest_approach_function(r, a, Za, Zb, relative_energy, impact_parameter, interaction_potential)};

    //Guess for large reduced energy from Mendenhall and Weller 1991
    //For small energies, use pure Newton-Raphson with arbitrary guess of 1
    let mut x0 = beta;
    let mut xn: f64;
    if reduced_energy > 5. {
        let inv_er_2 = 0.5/reduced_energy;
        x0 = inv_er_2 + (inv_er_2*inv_er_2 + beta*beta).sqrt();
    }

    //Newton-Raphson to determine distance of closest approach
    let mut err: f64 = tolerance + 1.;
    for k in 0..max_iterations {
        xn = x0 - f(x0*a)/df(x0*a);
        err = (xn - x0)*(xn - x0);
        x0 = xn;
        if err < tolerance {
            return Ok(x0);
        }
    }
    return Err(anyhow!("Numerical error: exceeded maximum number of Newton-Raphson iterations, {}. E: {}; x0: {}; Error: {}; Tolerance: {}",
        max_iterations, E0, x0, err, tolerance));
}

/// Gauss-Mehler quadrature.
pub fn gauss_mehler(Za: f64, Zb: f64, Ma: f64, Mb: f64, E0: f64, impact_parameter: f64, x0: f64, interaction_potential: InteractionPotential, n_points: usize) -> f64 {
    let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
    let r0 = x0*a;
    let V = |r| {interactions::interaction_potential(r, a, Za, Zb, interaction_potential)};
    let relative_energy = E0*Mb/(Ma + Mb);
    scattering_integral_gauss_mehler(impact_parameter, relative_energy, r0, &V, n_points)
}

/// Gauss-Legendre quadrature.
pub fn gauss_legendre(Za: f64, Zb: f64, Ma: f64, Mb: f64, E0: f64, impact_parameter: f64, x0: f64, interaction_potential: InteractionPotential) -> f64 {
    let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
    let r0 = x0*a;
    let V = |r| {interactions::interaction_potential(r, a, Za, Zb, interaction_potential)};
    let relative_energy = E0*Mb/(Ma + Mb);
    scattering_integral_gauss_legendre(impact_parameter, relative_energy, r0, &V)
}

/// Ziegler's MAGIC algorithm for approximating the scattering integral.
pub fn magic(Za: f64, Zb: f64, Ma: f64, Mb: f64, E0: f64, impact_parameter: f64, x0: f64, interaction_potential: InteractionPotential) -> f64 {
    //MAGIC algorithm
    //Since this is legacy code I don't think I will clean this up
    let C_ = match  interaction_potential {
        InteractionPotential::MOLIERE => vec![ 0.6743, 0.009611, 0.005175, 6.314, 10.0 ],
        InteractionPotential::KR_C => vec![ 0.7887, 0.01166, 00.006913, 17.16, 10.79 ],
        InteractionPotential::ZBL => vec![ 0.99229, 0.011615, 0.0071222, 9.3066, 14.813 ],
        InteractionPotential::TRIDYN => vec![1.0144, 0.235809, 0.126, 69350., 83550.], //Undocumented Tridyn constants
        _ => panic!("Input error: unimplemented interaction potential {} for MAGIC algorithm. Use a screened Coulomb potential.",  interaction_potential)
    };
    let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
    let beta: f64 = impact_parameter/a;
    let V0 = Za*Zb*Q*Q/4.0/PI/EPS0/a;
    let relative_energy = E0*Mb/(Ma + Mb);
    let reduced_energy = relative_energy/V0;
    let SQE = reduced_energy.sqrt();
    let r0 = a*x0;
    let V = V0*a/r0*interactions::phi(x0, interaction_potential);
    let dV = -V/r0 + V0/r0*interactions::dphi(x0, interaction_potential);
    let rho = -2.0*(relative_energy - V)/dV;
    let D = 2.0*(1.0 + C_[0]/SQE)*reduced_energy*beta.powf((C_[1] + SQE)/(C_[2] + SQE));
    let G = (C_[4] + reduced_energy)/(C_[3] + reduced_energy)*((1. + D*D).sqrt() - D);
    let delta =  D*G/(1.0 + G)*(x0 - beta);
    let ctheta2 = (beta + rho/a + delta)/(x0 + rho/a);
    2.*((beta + rho/a + delta)/(x0 + rho/a)).acos()
}

/// Mendenhall-Weller quadrature.
pub fn mendenhall_weller(Za: f64, Zb: f64, Ma: f64, Mb: f64, E0: f64, impact_parameter: f64, x0: f64, interaction_potential: InteractionPotential) -> f64 {
    //Lindhard screening length and reduced energy
    let a: f64 = interactions::screening_length(Za, Zb, interaction_potential);
    let reduced_energy: f64 = LINDHARD_REDUCED_ENERGY_PREFACTOR*a*Mb/(Ma+Mb)/Za/Zb*E0;
    let beta: f64 = impact_parameter/a;

    //Scattering integral quadrature from Mendenhall and Weller 2005
    let lambda0 = (0.5 + beta*beta/x0/x0/2. - interactions::dphi(x0,  interaction_potential)/2./reduced_energy).powf(-1./2.);
    let alpha = 1./12.*(1. + lambda0 + 5.*(0.4206*scattering_integral_mw(x0/0.9072, beta, reduced_energy,  interaction_potential) + 0.9072*scattering_integral_mw(x0/0.4206, beta, reduced_energy,  interaction_potential)));
    PI*(1. - beta*alpha/x0)
}

/// Oen-Robinson local electronic energy loss for a collision between particles a and b.
fn oen_robinson_loss(Za: f64, Zb: f64, Se: f64, x0: f64, interaction_potential: InteractionPotential) -> f64 {
    //Oen-Robinson local electronic stopping power
    let a = interactions::screening_length(Za, Zb, interaction_potential);

    //d1 is the first (smallest) interior constant of the screening function
    let d1 = interactions::first_screening_radius(interaction_potential);
    d1*d1/2./PI*Se*(-d1*x0).exp()/a/a
}
