#[cfg(test)]
use super::*;

#[test]
fn test_surface_refraction() {
    let mass = 1.;
    let Z = 1.;
    let E = 10.*EV;
    let Ec = 1.*EV;
    let Es = 5.76*EV;
    let x = 0.;
    let y = 0.;
    let z = 0.;
    let cosx = 1./(2.0_f64).sqrt();
    let cosy = 1./(2.0_f64).sqrt();
    let cosz = 0.;
    let mut particle_1 = particle::Particle::new(mass, Z, E, Ec, Es, x, y, z, cosx, cosy, cosz, false, false);

    //Test particle entering material and gaining energy

    //Eckstein's formulation for particle entering surface
    let cosx_new = ((E*cosx*cosx + Es)/(E + Es)).sqrt();
    let sinx = (1. - cosx*cosx).sqrt();
    let sinx_new = (1. - cosx_new*cosx_new).sqrt();
    let cosy_new = cosy*sinx_new/sinx;
    let cosz_new = cosz*sinx_new/sinx;
    let dir_mag = (cosx_new*cosx_new + cosy_new*cosy_new + cosz_new*cosz_new).sqrt();
    println!("dir_mag: {}", dir_mag);
    let cosx_new = cosx_new/dir_mag;
    let cosy_new = cosy_new/dir_mag;
    let cosz_new = cosz_new/dir_mag;

    let delta_theta = particle::refraction_angle(cosx, E, E + Es);
    particle::rotate_particle(&mut particle_1, delta_theta, 0.);

    println!("{} {} {}", cosx, cosy, cosz);
    println!("{} {} {}", cosx_new, cosy_new, cosz_new);
    println!("{} {} {}", particle_1.dir.x, particle_1.dir.y, particle_1.dir.z);
    println!();

    assert!(approx_eq!(f64, particle_1.dir.x, cosx_new, epsilon=1E-9));
    assert!(approx_eq!(f64, particle_1.dir.y, cosy_new, epsilon=1E-9));
    assert!(approx_eq!(f64, particle_1.dir.z, cosz_new, epsilon=1E-9));

    //Test particle leaving material and losing energy

    let cosx_new = ((particle_1.E*particle_1.dir.x*particle_1.dir.x - Es)/(particle_1.E - Es)).sqrt();
    let sinx = (1. - particle_1.dir.x*particle_1.dir.x).sqrt();
    let sinx_new = (1. - cosx_new*cosx_new).sqrt();
    let cosy_new = particle_1.dir.y*sinx_new/sinx;
    let cosz_new = particle_1.dir.z*sinx_new/sinx;

    println!("{} {} {}", particle_1.dir.x, particle_1.dir.y, particle_1.dir.z);
    let delta_theta = particle::refraction_angle(particle_1.dir.x, particle_1.E, particle_1.E - Es);
    particle::rotate_particle(&mut particle_1, delta_theta, 0.);

    println!("{} {} {}", cosx_new, cosy_new, cosz_new);
    println!("{} {} {}", particle_1.dir.x, particle_1.dir.y, particle_1.dir.z);

    assert!(approx_eq!(f64, particle_1.dir.x, cosx_new, epsilon=1E-9));
    assert!(approx_eq!(f64, particle_1.dir.y, cosy_new, epsilon=1E-9));
    assert!(approx_eq!(f64, particle_1.dir.z, cosz_new, epsilon=1E-9));
}

#[test]
fn test_momentum_conservation() {

    for energy_eV in vec![1., 10., 100., 1000., 10000., 100000., 1000000., 10000000.] {
        //Aluminum
        let m1 = 4.008*AMU;
        let Z1 = 2.;
        let E1 = energy_eV*EV;
        let Ec1 = 1.*EV;
        let Es1 = 1.*EV;
        let x1 = 0.;
        let y1 = 0.;
        let z1 = 0.;

        //Aluminum
        let m2 = 6.941;
        let Z2 = 3.;
        let Ec2 = 1.;
        let Es2 = 1.;

        //Random initial angle
        let theta = 0.974194583091052_f64;
        let cosx = (theta).cos();
        let cosy = (theta).sin();
        let cosz = 0.;

        let material_parameters = material::MaterialParameters{
            energy_unit: "EV".to_string(),
            mass_unit: "AMU".to_string(),
            Eb: vec![0.0],
            Es: vec![Es2],
            Ec: vec![Ec2],
            n: vec![6.026E28],
            Z: vec![Z2],
            m: vec![m2],
            electronic_stopping_correction_factor: 0.0
        };

        //Generate empty geometry
        let geometry = material::Geometry {
            length_unit: "ANGSTROM".to_string(),
            surface: vec![],
            energy_surface: vec![],
            simulation_surface: vec![]
        };

        let thickness: f64 = 1000.;
        let depth: f64 = 1000.;
        let mesh_2d_input = material::Mesh2DInput {
            length_unit: "ANGSTROM".to_string(),
            coordinate_sets: vec![(0., depth, 0., thickness/2., -thickness/2., -thickness/2.), (0., depth, depth, thickness/2., thickness/2., -thickness/2.)],
            data: vec![vec![6.026E28], vec![6.026E28]],
        };

        let material_1 = material::Material::new(material_parameters, geometry, mesh_2d_input);

        for high_energy_free_flight_paths in vec![true, false] {
            for potential in vec![KR_C, MOLIERE, ZBL, TRIDYN] {
                for scattering_integral in vec![QUADRATURE, MAGIC] {

                    println!("Case: {} {} {} {}", energy_eV, high_energy_free_flight_paths, potential, scattering_integral);

                    let mut particle_1 = particle::Particle::new(m1, Z1, E1, Ec1, Es1, x1, y1, z1, cosx, cosy, cosz, false, false);

                    let options = Options {
                        name: "test".to_string(),
                        track_trajectories: true,
                        track_recoils: true,
                        track_recoil_trajectories: true,
                        write_files: true,
                        stream_size: 8000,
                        print: true,
                        print_num: 100,
                        weak_collision_order: 0,
                        suppress_deep_recoils: false,
                        high_energy_free_flight_paths: high_energy_free_flight_paths,
                        electronic_stopping_mode: INTERPOLATED,
                        mean_free_path_model: LIQUID,
                        interaction_potential: potential,
                        scattering_integral: scattering_integral,
                        tolerance: 1E-12,
                        max_iterations: 100,
                    };

                    let binary_collision_geometries = bca::determine_mfp_phi_impact_parameter(&mut particle_1, &material_1, &options);

                    println!("Phi: {} rad p: {} Angstrom mfp: {} Angstrom", binary_collision_geometries[0].phi_azimuthal,
                        binary_collision_geometries[0].impact_parameter/ANGSTROM,
                        binary_collision_geometries[0].mfp/ANGSTROM);

                    let (species_index, mut particle_2) = bca::choose_collision_partner(&mut particle_1, &material_1,
                        &binary_collision_geometries[0], &options);

                    let mom1_0 = particle_1.get_momentum();
                    let mom2_0 = particle_2.get_momentum();

                    let initial_momentum = mom1_0.add(&mom2_0);

                    let binary_collision_result = bca::calculate_binary_collision(&particle_1,
                        &particle_2, &binary_collision_geometries[0], &options);

                    println!("E_recoil: {} eV Psi: {} rad Psi_recoil: {} rad", binary_collision_result.recoil_energy/EV,
                        binary_collision_result.psi,
                        binary_collision_result.psi_recoil);

                    println!("Initial Energies: {} eV {} eV", particle_1.E/EV, particle_2.E/EV);

                    //Energy transfer to recoil
                    particle_2.E = binary_collision_result.recoil_energy - material_1.average_bulk_binding_energy(particle_2.pos.x, particle_2.pos.y);

                    //Rotate particle 1, 2 by lab frame scattering angles
                    particle::rotate_particle(&mut particle_1, binary_collision_result.psi,
                        binary_collision_geometries[0].phi_azimuthal);

                    particle::rotate_particle(&mut particle_2, -binary_collision_result.psi_recoil,
                        binary_collision_geometries[0].phi_azimuthal);

                    //Subtract total energy from all simultaneous collisions and electronic stopping
                    bca::update_particle_energy(&mut particle_1, &material_1, 0.,
                        binary_collision_result.recoil_energy, 0., particle_2.Z, species_index, &options);

                    let mom1_1 = particle_1.get_momentum();
                    let mom2_1 = particle_2.get_momentum();

                    let final_momentum = mom1_1.add(&mom2_1);

                    println!("Final Energies: {} eV {} eV", particle_1.E/EV, particle_2.E/EV);
                    println!("X: {} {} {}% Error", initial_momentum.x/ANGSTROM/AMU, final_momentum.x/ANGSTROM/AMU, 100.*(final_momentum.x - initial_momentum.x)/initial_momentum.magnitude());
                    println!("Y: {} {} {}% Error", initial_momentum.y/ANGSTROM/AMU, final_momentum.y/ANGSTROM/AMU, 100.*(final_momentum.y - initial_momentum.y)/initial_momentum.magnitude());
                    println!("Z: {} {} {}% Error", initial_momentum.z/ANGSTROM/AMU, final_momentum.z/ANGSTROM/AMU, 100.*(final_momentum.z - initial_momentum.z)/initial_momentum.magnitude());
                    println!();

                    assert!(approx_eq!(f64, initial_momentum.x, final_momentum.x, epsilon = 1E-12));
                    assert!(approx_eq!(f64, initial_momentum.y, final_momentum.y, epsilon = 1E-12));
                    assert!(approx_eq!(f64, initial_momentum.z, final_momentum.z, epsilon = 1E-12));

                    assert!(!particle_1.E.is_nan());
                    assert!(!particle_2.E.is_nan());
                    assert!(!initial_momentum.x.is_nan());
                    assert!(!initial_momentum.x.is_nan());
                    assert!(!initial_momentum.x.is_nan());
                }
            }
        }
    }
}

#[test]
fn test_rotate_particle() {
    let mass = 1.;
    let Z = 1.;
    let E = 1.;
    let Ec = 1.;
    let Es = 1.;
    let x = 0.;
    let y = 0.;
    let z = 0.;
    let cosx = (PI/4.).cos();
    let cosy = (PI/4.).sin();
    let cosz = 0.;
    let psi = -PI/4.;
    let phi = 0.;

    let mut particle = particle::Particle::new(mass, Z, E, Ec, Es, x, y, z, cosx, cosy, cosz, false, false);

    //Check that rotation in 2D works
    particle::rotate_particle(&mut particle, psi, phi);
    assert!(approx_eq!(f64, particle.dir.x, 0., epsilon = 1E-12), "particle.dir.x: {} Should be ~0.", particle.dir.x);
    assert!(approx_eq!(f64, particle.dir.y, 1., epsilon = 1E-12), "particle.dir.y: {} Should be ~1.", particle.dir.y);

    //Check that rotating back by negative psi returns to the previous values
    particle::rotate_particle(&mut particle, -psi, phi);
    assert!(approx_eq!(f64, particle.dir.x, cosx, epsilon = 1E-12), "particle.dir.x: {} Should be ~{}", particle.dir.x, cosx);
    assert!(approx_eq!(f64, particle.dir.y, cosy, epsilon = 1E-12), "particle.dir.y: {} Should be ~{}", particle.dir.y, cosy);

    //Check that azimuthal rotation by 180 degrees works correctly
    let phi = PI;
    particle::rotate_particle(&mut particle, psi, phi);
    assert!(approx_eq!(f64, particle.dir.x, 1., epsilon = 1E-12), "particle.dir.x: {} Should be ~1.", particle.dir.x);
    assert!(approx_eq!(f64, particle.dir.y, 0., epsilon = 1E-12), "particle.dir.y: {} Should be ~0.", particle.dir.y);

    //Check that particle direction vector remains normalized following rotations
    assert!(approx_eq!(f64, particle.dir.x.powf(2.) + particle.dir.y.powf(2.) + particle.dir.z.powf(2.), 1.), "Particle direction not normalized.");

}

#[test]
fn test_particle_advance() {
    let mass = 1.;
    let Z = 1.;
    let E = 1.;
    let Ec = 1.;
    let Es = 1.;
    let x = 0.;
    let y = 0.;
    let z = 0.;
    let cosx = (PI/4.).cos();
    let cosy = (PI/4.).sin();
    let cosz = 0.;
    let mfp = 1.;
    let asymptotic_deflection = 0.5;

    let mut particle = particle::Particle::new(mass, Z, E, Ec, Es, x, y, z, cosx, cosy, cosz, false, false);

    let distance_traveled = particle::particle_advance(&mut particle, mfp, asymptotic_deflection);

    assert_eq!(particle.pos.x, (1. - 0.5)*cosx);
    assert_eq!(particle.pos.y, (1. - 0.5)*cosy);
    assert_eq!(particle.pos.z, 0.);
    assert_eq!(distance_traveled, mfp - asymptotic_deflection);
}
