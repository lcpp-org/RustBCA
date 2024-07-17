#[cfg(test)]
use super::*;
#[cfg(test)]
use float_cmp::*;

use geo::algorithm::contains::Contains;
use geo::{Polygon, LineString, Point, point, Closest};
use geo::algorithm::closest_point::ClosestPoint;

use parry2d_f64::shape::Polyline;
use parry2d_f64::query::{PointQuery, PointProjection, Ray, RayCast};
use parry2d_f64::math::{Isometry};
use parry2d_f64::math::Point as Point2d;
use parry2d_f64::bounding_volume::aabb;

/*
        let number_boundary_points = boundary_points_converted.clone().len() as u32;
        let boundary = Polygon::new(LineString::from(boundary_points_converted), vec![]);
        let mut linked_points = (0..number_boundary_points).zip(1..number_boundary_points).map(|(x, y)| [x, y]).collect::<Vec<[u32; 2]>>();
        let boundary2 = Polyline::new(boundary_points_converted2, Some(linked_points));
*/

#[test]
fn test_parry2d() {
    let points: Vec<Point2d<f64>> = vec![Point2d::new(-1.0, -1.0), Point2d::new(-1.0, 1.0), Point2d::new(1.0, 1.0), Point2d::new(1.0, -1.0)];
    //let mut indices: Vec<[u32; 2]> = vec![];
    //let mut indices: Vec<[u32; 2]> = vec![[0, 1], [1, 2], [2, 3], [3, 0]];
    let mut indices: Vec<[u32; 2]> = vec![[0, 3], [3, 2], [2, 1], [1, 0]];
    //indices.reverse();

    //(0.0000006261797114005236, 0.000006009591179670447)
    let query_point = Point2d::new(0.0, 0.0);

    let polyline = Polyline::new(points, Some(indices));

    let (point_projection, (_, _)) = polyline.project_local_point_assuming_solid_interior_ccw(
        query_point
    );

    assert!(point_projection.is_inside);

    assert!(polyline.aabb(&Isometry::identity()).contains_local_point(&query_point));
    //assert!(point_projection.is_inside);

    let test_geometry = vec![(0.25, 0.2), (0.75, 0.2), (0.8, 0.8), (0.5, 0.5), (0.2, 0.8)];
    let num_segments = 5;

    let geometry_input_homogeneous_2D = geometry::HomogeneousMesh2DInput {
        length_unit: "M".to_string(),
        points: test_geometry.clone(),
        densities: vec![0.03, 0.03],
        simulation_boundary_points: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        electronic_stopping_correction_factor: 1.0,
    };

    let material_parameters = material::MaterialParameters{
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![0.0, 0.0],
        Es: vec![2.0, 4.0],
        Ec: vec![1.0, 1.0],
        Ed: vec![0.0, 0.0],
        Z: vec![29., 1.],
        m: vec![63.54, 1.0008],
        interaction_index: vec![0, 0],
        surface_binding_model: SurfaceBindingModel::PLANAR{calculation: SurfaceBindingCalculation::TARGET},
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let material_homogeneous_2D: material::Material<geometry::ParryHomogeneousMesh2D> = material::Material::<geometry::ParryHomogeneousMesh2D>::new(&material_parameters, &geometry_input_homogeneous_2D);

    let query_points_outside = vec![[0.5, 0.1, 0.0], [0.9, 0.5, 0.0], [0.6, 0.7, 0.0], [0.4, 0.7, 0.0], [0.1, 0.5, 0.0]];
    let query_points_inside = vec![[0.5, 0.21, 0.0], [0.75, 0.21, 0.0], [0.51, 0.5, 0.0], [0.49, 0.5, 0.0], [0.25, 0.21, 0.0]];

    for query_point in query_points_outside.clone() {
        assert!(!material_homogeneous_2D.inside(query_point[0], query_point[1], query_point[2]), "Point ({}, {}, {}) failed outsidedness check.", query_point[0], query_point[1], query_point[2]);
    }

    for query_point in query_points_inside.clone() {
        assert!(material_homogeneous_2D.inside(query_point[0], query_point[1], query_point[2]), "Point ({}, {}, {}) failed outsidedness check.", query_point[0], query_point[1], query_point[2]);
    }


    let normals_outside = query_points_outside
        .iter()
        .map( |query_point| {
            (material_homogeneous_2D.geometry.nearest_normal_vector(query_point[0], query_point[1], query_point[2]))
        }
    ).collect::<Vec<(f64, f64, f64)>>();

    let normals_inside = query_points_inside
    .iter()
    .map( |query_point| {
        (material_homogeneous_2D.geometry.nearest_normal_vector(query_point[0], query_point[1], query_point[2]))
    }
    ).collect::<Vec<(f64, f64, f64)>>();

    let test_normals = (0..num_segments)
        .zip(1..num_segments + 1)
        .map(|(i, j)| {
            let dx = test_geometry[i % num_segments].0 - test_geometry[j % num_segments].0;
            let dy = test_geometry[i % num_segments].1 - test_geometry[j % num_segments].1;
            let mag = (dx*dx + dy*dy).sqrt();
            (-dy/mag, dx/mag, 0.0)
        }).collect::<Vec<(f64, f64, f64)>>();

    for (test_normal, normal) in test_normals.iter().zip(normals_outside) {
        assert!(approx_eq!(f64, normal.0, test_normal.0, epsilon=1E-9), "Normal vector check failed. Calculated Normal: ({}, {}, {}) Test: ({}, {}, {})", normal.0, normal.1, normal.2, test_normal.0, test_normal.1, test_normal.2);
        assert!(approx_eq!(f64, normal.1, test_normal.1, epsilon=1E-9));
        assert!(approx_eq!(f64, normal.2, test_normal.2, epsilon=1E-9));
    }

    for (test_normal, normal) in test_normals.iter().zip(normals_inside) {
        assert!(approx_eq!(f64, normal.0, test_normal.0, epsilon=1E-9), "Normal vector check failed. Calculated Normal: ({}, {}, {}) Test: ({}, {}, {})", normal.0, normal.1, normal.2, test_normal.0, test_normal.1, test_normal.2);
        assert!(approx_eq!(f64, normal.1, test_normal.1, epsilon=1E-9));
        assert!(approx_eq!(f64, normal.2, test_normal.2, epsilon=1E-9));
    }

}


#[test]
#[cfg(feature = "cpr_rootfinder")]
fn test_polynom() {
    use rcpr::chebyshev::*;
    let interaction_potential = InteractionPotential::FOUR_EIGHT{alpha: 1., beta: 1.};
    let coefficients = interactions::polynomial_coefficients(1., 1., interaction_potential);
    let roots = real_polynomial_roots(coefficients.clone(), 1e-9).unwrap();

    let max_root = roots.iter().cloned().fold(f64::NAN, f64::max);
    println!("{}", max_root);
    let inverse_transformed_root = interactions::inverse_transform(max_root, interaction_potential);
}

#[test]
#[cfg(feature = "parry3d")]
fn test_parry_cuboid() {
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
    let particle_1 = particle::Particle::new(mass, Z, E, Ec, Es, 0.0, x, y, z, cosx, cosy, cosz, false, false, 0);

    let material_parameters = material::MaterialParameters{
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![0.0, 0.0],
        Es: vec![2.0, 4.0],
        Ec: vec![1.0, 1.0],
        Ed: vec![0.0, 0.0],
        Z: vec![29., 1.],
        m: vec![63.54, 1.0008],
        interaction_index: vec![0, 0],
        surface_binding_model: SurfaceBindingModel::PLANAR{calculation: SurfaceBindingCalculation::TARGET},
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let radius = 1000.*ANGSTROM;

    //cuboid centered at (0, 0, 0) with l=w=h=1 Angstrom
    let geometry_input = parry::ParryTriMeshInput {
        length_unit: "ANGSTROM".to_string(),
        densities: vec![0.03, 0.03],
        electronic_stopping_correction_factor: 1.0,
        vertices: vec![
            [-0.5, -0.5, 0.5],
            [0.5, -0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [0.5, 0.5, 0.5],
            [-0.5, 0.5, -0.5],
            [0.5, 0.5, -0.5],
            [-0.5, -0.5, -0.5],
            [0.5, -0.5, -0.5],
        ],
        indices: vec![
            [0, 1, 3],
            [0, 3, 2],
            [2, 3, 5],
            [2, 5, 4],
            [4, 5, 7],
            [4, 7, 6],
            [6, 7, 1],
            [6, 1, 0],
            [1, 7, 5],
            [1, 5, 3],
            [6, 0, 2],
            [6, 2, 4]
        ]
    };

    let mut material_cuboid: material::Material<parry::ParryTriMesh> = material::Material::<parry::ParryTriMesh>::new(&material_parameters, &geometry_input);
    material_cuboid.geometry.energy_barrier_thickness = 10.*ANGSTROM;

    let surface_binding_energy = material_cuboid.actual_surface_binding_energy(&particle_1, particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);
    assert!(approx_eq!(f64, surface_binding_energy/EV, (2. + 4.)/2., epsilon=1E-24));

    //Test origin is inside
    assert!(material_cuboid.geometry.inside(0., 0., 0.));

    //Test all outside areas for containment
    assert!(!material_cuboid.geometry.inside(0., 1., 0.));
    assert!(!material_cuboid.geometry.inside(0., -1., 0.));
    assert!(!material_cuboid.geometry.inside(1., 0., 0.));
    assert!(!material_cuboid.geometry.inside(-1., 0., 0.));
    assert!(!material_cuboid.geometry.inside(0., 0., -1.));
    assert!(!material_cuboid.geometry.inside(0., 0., -1.));

    //distance to origin
    let (x, y, z) = material_cuboid.geometry.closest_point(0., 0., 0.);
    assert!(approx_eq!(f64, (x.powi(2) + y.powi(2) + z.powi(2)).sqrt(), 0.5*ANGSTROM, epsilon=1E-24));

    //distance to wall should be zero at wall
    let (x, y, z) = material_cuboid.geometry.closest_point(0.5*ANGSTROM, 0., 0.);
    assert!(approx_eq!(f64, ((x - 0.5*ANGSTROM).powi(2) + y.powi(2) + z.powi(2)).sqrt(), 0.0, epsilon=1E-24));
}

#[test]
#[cfg(feature = "parry3d")]
fn test_parry_sphere() {
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
    let mut particle_1 = particle::Particle::new(mass, Z, E, Ec, Es, 0.0, x, y, z, cosx, cosy, cosz, false, false, 0);

    let material_parameters = material::MaterialParameters{
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![0.0, 0.0],
        Es: vec![2.0, 4.0],
        Ec: vec![1.0, 1.0],
        Ed: vec![0.0, 0.0],
        Z: vec![29., 1.],
        m: vec![63.54, 1.0008],
        interaction_index: vec![0, 0],
        surface_binding_model: SurfaceBindingModel::PLANAR{calculation: SurfaceBindingCalculation::TARGET},
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let radius = 1000.*ANGSTROM;

    let geometry_input_sphere = parry::ParryBallInput {
        length_unit: "ANGSTROM".to_string(),
        radius: 1000.0,
        densities: vec![0.03, 0.03],
        electronic_stopping_correction_factor: 1.0
    };

    let mut material_sphere: material::Material<parry::ParryBall> = material::Material::<parry::ParryBall>::new(&material_parameters, &geometry_input_sphere);
    material_sphere.geometry.energy_barrier_thickness = 10.*ANGSTROM;

    particle_1.pos.x = 500.*ANGSTROM;
    particle_1.pos.y = 0.;

    particle_1.pos_old.x = -1500.*ANGSTROM;
    particle_1.pos_old.y = 0.;

    let inside_sphere = material_sphere.inside(particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);
    let inside_old_sphere = material_sphere.inside(particle_1.pos_old.x, particle_1.pos_old.y, particle_1.pos_old.z);

    //println!("{} {}", inside, inside_old);
    assert!(inside_sphere);
    assert!(!inside_old_sphere);

    //Test concentration-dependent surface binding energy
    let surface_binding_energy = material_sphere.actual_surface_binding_energy(&particle_1, particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);
    assert!(approx_eq!(f64, surface_binding_energy/EV, (2. + 4.)/2., epsilon=1E-24));

    assert!(material_sphere.inside_energy_barrier(0., 0., 0.));
    assert!(material_sphere.inside_simulation_boundary(0., 0., 0.));
    assert!(material_sphere.inside(0., 0., 0.));

    let ux = 1214.0*ANGSTROM;
    let uy = 123123.0*ANGSTROM;
    let uz = 1239.0*ANGSTROM;
    let r = (ux.powi(2) + uy.powi(2) + uz.powi(2)).sqrt();
    let R = 1001.0*ANGSTROM;
    let u = (ux/r*R, uy/r*R, uz/r*R);
    assert!(!material_sphere.inside(u.0, u.1, u.2));
    assert!(material_sphere.inside_energy_barrier(u.0, u.1, u.2));
    assert!(material_sphere.inside_simulation_boundary(u.0, u.1, u.2));

    assert!(material_sphere.inside_energy_barrier(0., 0., 0.));
    assert!(material_sphere.inside_simulation_boundary(0., 0., 0.));
    assert!(material_sphere.inside(0., 0., 0.));

    assert!(material_sphere.inside_energy_barrier(-1000.0*ANGSTROM - 1.0*ANGSTROM, 0., 0.));
    assert!(!material_sphere.inside_energy_barrier(-1000.0*ANGSTROM - 11.0*ANGSTROM, 0., 0.));

    assert!(material_sphere.inside_energy_barrier(0., -1000.0*ANGSTROM - 1.0*ANGSTROM, 0.));
    assert!(!material_sphere.inside_energy_barrier(0., -1000.0*ANGSTROM - 11.0*ANGSTROM, 0.));

    assert!(material_sphere.inside_energy_barrier(0., 0., -1000.0*ANGSTROM - 1.0*ANGSTROM));
    assert!(!material_sphere.inside_energy_barrier(0., 0., -1000.0*ANGSTROM - 11.0*ANGSTROM));

    assert!(material_sphere.inside(-1000.0*ANGSTROM + 1.0*ANGSTROM, 0., 0.));
    assert!(!material_sphere.inside(-1000.0*ANGSTROM - 1.0*ANGSTROM, 0., 0.));

    assert!(material_sphere.inside(0., -1000.0*ANGSTROM + 1.0*ANGSTROM, 0.));
    assert!(!material_sphere.inside(0., -1000.0*ANGSTROM - 1.0*ANGSTROM, 0.));

    assert!(material_sphere.inside(0., 0., -1000.0*ANGSTROM + 1.0*ANGSTROM));
    assert!(!material_sphere.inside(0., 0., -1000.0*ANGSTROM - 1.0*ANGSTROM));

    assert!(approx_eq!(f64, material_sphere.closest_point(-2000.*ANGSTROM, 0., 0.).0, (-1000.*ANGSTROM, 0., 0.).0, epsilon=1E-12));
    assert!(approx_eq!(f64, material_sphere.closest_point(0., -2000.*ANGSTROM, 0.).1, (0., -1000.*ANGSTROM, 0.).1, epsilon=1E-12));
    assert!(approx_eq!(f64, material_sphere.closest_point(0., 0., -2000.*ANGSTROM).2, (0., 0., -1000.*ANGSTROM).2, epsilon=1E-12));
}

#[test]
#[cfg(feature = "distributions")]
fn test_distributions() {

    let options = Options {
        name: "test".to_string(),
        track_trajectories: false,
        track_recoils: false,
        track_recoil_trajectories: false,
        write_buffer_size: 8000,
        weak_collision_order: 0,
        suppress_deep_recoils: false,
        high_energy_free_flight_paths: false,
        electronic_stopping_mode: ElectronicStoppingMode::INTERPOLATED,
        mean_free_path_model: MeanFreePathModel::LIQUID,
        interaction_potential: vec![vec![InteractionPotential::KR_C]],
        scattering_integral: vec![vec![ScatteringIntegral::MENDENHALL_WELLER]],
        num_threads: 1,
        num_chunks: 1,
        use_hdf5: false,
        root_finder: vec![vec![Rootfinder::NEWTON{max_iterations: 100, tolerance: 1E-3}]],
        track_displacements: false,
        track_energy_losses: true,
        energy_min: 0.0,
        energy_max: 10.0,
        energy_num: 11,
        angle_min: 0.0,
        angle_max: 90.0,
        angle_num: 11,
        x_min: 0.0,
        y_min: -10.0,
        z_min: -10.0,
        x_max: 10.0,
        y_max: 10.0,
        z_max: 10.0,
        x_num: 11,
        y_num: 11,
        z_num: 11,
    };

    let output_units = OutputUnits {
        length_unit: 1.0,
        energy_unit: 1.0,
        mass_unit: 1.0,
    };

    let mass = 1.0;
    let Z = 1.0;
    let E = 1.5;
    let Ec = 0.0;
    let Es = 0.0;
    let x = 0.0;
    let y = 0.0;
    let z = 0.0;
    let cosx = 0.0;
    let cosy = 0.0;
    let cosz = 0.0;
    let mut particle = particle::Particle::new(mass, Z, E, Ec, Es, 0.0, x, y, z, cosx, cosy, cosz, false, false, 0);

    let mut distributions = output::Distributions::new(&options);
    assert_eq!(distributions.x_range[0], 0.);
    assert_eq!(distributions.x_range[distributions.x_range.len() - 1], 10.);
    assert_eq!(distributions.y_range[0], -10.);
    assert_eq!(distributions.y_range[distributions.y_range.len() - 1], 10.);
    assert_eq!(distributions.z_range[0], -10.);
    assert_eq!(distributions.z_range[distributions.z_range.len() - 1], 10.);
    assert_eq!(distributions.x_range.len(), 11);
    assert_eq!(distributions.y_range.len(), 11);
    assert_eq!(distributions.z_range.len(), 11);

    particle.incident = true;
    particle.stopped = true;
    particle.pos.x = 0.0;
    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.implanted_x[0], 1);

    particle.pos.x = 10.0;
    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.implanted_x[10], 1);

    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.implanted_y[5], 3);
    assert_eq!(distributions.implanted_z[5], 3);

    particle.incident = false;
    particle.stopped = false;
    particle.left = true;
    particle.E = 0.0;
    particle.dir.x = -0.999;
    particle.dir.y = 0.001;
    particle.dir.z = 0.0;
    particle.dir.normalize();

    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.sputtered_ead[[0, 0]], 1);

    particle.E = -0.1;
    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.sputtered_ead[[0, 0]], 1);

    particle.E = 10.0;
    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.sputtered_ead[[10, 0]], 1);

    particle.dir = Vector::new(-0.707, 0.707, 0.0);
    particle.dir.normalize();
    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.sputtered_ead[[10, 5]], 1);

    particle.incident = true;
    particle.E = 0.0;
    particle.dir.x = -0.999;
    particle.dir.y = 0.001;
    particle.dir.z = 0.0;
    particle.dir.normalize();

    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.reflected_ead[[0, 0]], 1);

    particle.E = -0.1;
    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.reflected_ead[[0, 0]], 1);

    particle.E = 10.0;
    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.reflected_ead[[10, 0]], 1);

    particle.dir = Vector::new(-0.707, 0.707, 0.0);
    particle.dir.normalize();
    distributions.update(&particle, &output_units, &options, 1);
    assert_eq!(distributions.reflected_ead[[10, 5]], 1);
}

#[test]
fn test_spherical_geometry() {
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
    let mut particle_1 = particle::Particle::new(mass, Z, E, Ec, Es, 0.0, x, y, z, cosx, cosy, cosz, false, false, 0);

    let material_parameters = material::MaterialParameters{
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![0.0, 0.0],
        Es: vec![2.0, 4.0],
        Ec: vec![1.0, 1.0],
        Ed: vec![0.0, 0.0],
        Z: vec![29., 1.],
        m: vec![63.54, 1.0008],
        interaction_index: vec![0, 0],
        surface_binding_model: SurfaceBindingModel::PLANAR{calculation: SurfaceBindingCalculation::TARGET},
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let radius = 1000.*ANGSTROM;

    let geometry_input_sphere = sphere::SphereInput {
        length_unit: "ANGSTROM".to_string(),
        radius: 1000.0,
        densities: vec![0.03, 0.03],
        electronic_stopping_correction_factor: 1.0
    };

    let mut material_sphere: material::Material<sphere::Sphere> = material::Material::<sphere::Sphere>::new(&material_parameters, &geometry_input_sphere);
    material_sphere.geometry.energy_barrier_thickness = 10.*ANGSTROM;

    particle_1.pos.x = 500.*ANGSTROM;
    particle_1.pos.y = 0.;

    particle_1.pos_old.x = -1500.*ANGSTROM;
    particle_1.pos_old.y = 0.;

    let inside_sphere = material_sphere.inside(particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);
    let inside_old_sphere = material_sphere.inside(particle_1.pos_old.x, particle_1.pos_old.y, particle_1.pos_old.z);

    //println!("{} {}", inside, inside_old);
    assert!(inside_sphere);
    assert!(!inside_old_sphere);

    //Test concentration-dependent surface binding energy
    let surface_binding_energy = material_sphere.actual_surface_binding_energy(&particle_1, particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);
    assert!(approx_eq!(f64, surface_binding_energy/EV, (2. + 4.)/2., epsilon=1E-24));

    assert!(material_sphere.inside_energy_barrier(0., 0., 0.));
    assert!(material_sphere.inside_simulation_boundary(0., 0., 0.));
    assert!(material_sphere.inside(0., 0., 0.));

    let ux = 1214.0*ANGSTROM;
    let uy = 123123.0*ANGSTROM;
    let uz = 1239.0*ANGSTROM;
    let r = (ux.powi(2) + uy.powi(2) + uz.powi(2)).sqrt();
    let R = 1001.0*ANGSTROM;
    let u = (ux/r*R, uy/r*R, uz/r*R);
    assert!(!material_sphere.inside(u.0, u.1, u.2));
    assert!(material_sphere.inside_energy_barrier(u.0, u.1, u.2));
    assert!(material_sphere.inside_simulation_boundary(u.0, u.1, u.2));

    assert!(material_sphere.inside_energy_barrier(0., 0., 0.));
    assert!(material_sphere.inside_simulation_boundary(0., 0., 0.));
    assert!(material_sphere.inside(0., 0., 0.));

    assert!(material_sphere.inside_energy_barrier(-1000.0*ANGSTROM - 1.0*ANGSTROM, 0., 0.));
    assert!(!material_sphere.inside_energy_barrier(-1000.0*ANGSTROM - 11.0*ANGSTROM, 0., 0.));

    assert!(material_sphere.inside_energy_barrier(0., -1000.0*ANGSTROM - 1.0*ANGSTROM, 0.));
    assert!(!material_sphere.inside_energy_barrier(0., -1000.0*ANGSTROM - 11.0*ANGSTROM, 0.));

    assert!(material_sphere.inside_energy_barrier(0., 0., -1000.0*ANGSTROM - 1.0*ANGSTROM));
    assert!(!material_sphere.inside_energy_barrier(0., 0., -1000.0*ANGSTROM - 11.0*ANGSTROM));

    assert!(material_sphere.inside(-1000.0*ANGSTROM + 1.0*ANGSTROM, 0., 0.));
    assert!(!material_sphere.inside(-1000.0*ANGSTROM - 1.0*ANGSTROM, 0., 0.));

    assert!(material_sphere.inside(0., -1000.0*ANGSTROM + 1.0*ANGSTROM, 0.));
    assert!(!material_sphere.inside(0., -1000.0*ANGSTROM - 1.0*ANGSTROM, 0.));

    assert!(material_sphere.inside(0., 0., -1000.0*ANGSTROM + 1.0*ANGSTROM));
    assert!(!material_sphere.inside(0., 0., -1000.0*ANGSTROM - 1.0*ANGSTROM));

    assert!(approx_eq!(f64, material_sphere.closest_point(-2000.*ANGSTROM, 0., 0.).0, (-1000.*ANGSTROM, 0., 0.).0, epsilon=1E-12));
    assert!(approx_eq!(f64, material_sphere.closest_point(0., -2000.*ANGSTROM, 0.).1, (0., -1000.*ANGSTROM, 0.).1, epsilon=1E-12));
    assert!(approx_eq!(f64, material_sphere.closest_point(0., 0., -2000.*ANGSTROM).2, (0., 0., -1000.*ANGSTROM).2, epsilon=1E-12));

}

#[test]
fn extended_test_2D_geometry() {
    let mass = 1.;
    let Z = 1.;
    let E = 10.*EV;
    let Ec = 1.*EV;
    let Es = 5.76*EV;
    let x = 1.5e-6;
    let y = 10.0e-6;
    let z = 0.;
    let cosx = 0.0;
    let cosy = -1.0;
    let cosz = 0.;

    let mut particle_1 = particle::Particle::new(mass, Z, E, Ec, Es, 0.0, x, y, z, cosx, cosy, cosz, false, false, 0);

    let material_parameters = material::MaterialParameters{
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![0.0, 0.0],
        Es: vec![2.0, 4.0],
        Ec: vec![1.0, 1.0],
        Ed: vec![0.0, 0.0],
        Z: vec![29., 1.],
        m: vec![63.54, 1.0008],
        interaction_index: vec![0, 0],
        surface_binding_model: SurfaceBindingModel::PLANAR{calculation: SurfaceBindingCalculation::TARGET},
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let thickness: f64 = 1000.;
    let depth: f64 = 1000.;

    let geometry_input_2D = geometry::Mesh2DInput {
        length_unit: "MICRON".to_string(),
        points: vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0), (4.0, 1.0), (3.0, 1.0), (3.0, 10.0), (2.0, 10.0), (2.0, 1.0), (1.0, 1.0), (1.0, 10.0), (0.0, 10.0), (0.0, 1.0)],
        triangles: vec![(13, 12, 11), (13, 11, 10), (0, 13, 10), (0, 10, 1), (1, 10, 9), (1, 9, 2), (9, 8, 7), (9, 7, 6), (2, 9, 6), (2, 6, 3), (3, 6, 5), (3, 5, 4)],
        densities: vec![vec![3e10, 3e10]; 12],
        boundary: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0],
        simulation_boundary_points: vec![(-0.1, -0.1), (4.1, -0.1), (4.1, 10.1), (-0.1, 10.1), (-0.1, -0.1)],
        energy_barrier_thickness: 4e-4,
        electronic_stopping_correction_factors: vec![1.0; 12],
    };

    let material_2D: material::Material<geometry::Mesh2D> = material::Material::<geometry::Mesh2D>::new(&material_parameters, &geometry_input_2D);

    if let Closest::Intersection(p) = material_2D.geometry.boundary.closest_point(&point!(x: 1.5e-6, y: 1.0001e-6)) {
        println!("({} {})", p.x(), p.y());
    } else if let Closest::SinglePoint(p) = material_2D.geometry.boundary.closest_point(&point!(x: 1.5e-6, y: 1.0001e-6)) {
        println!("({} {})", p.x(), p.y());
    } else {
        panic!();
    }

    assert!(material_2D.inside(1.5e-6, 0.999e-6, 0.0));
    assert!(!material_2D.inside(1.5e-6, 1.0001e-6, 0.0));

    material::surface_binding_energy(&mut particle_1, &material_2D);

    material::boundary_condition_planar(&mut particle_1, &material_2D);
}

#[test]
fn test_geometry() {
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
    let mut particle_1 = particle::Particle::new(mass, Z, E, Ec, Es, 0.0, x, y, z, cosx, cosy, cosz, false, false, 0);

    let material_parameters = material::MaterialParameters{
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![0.0, 0.0],
        Es: vec![2.0, 4.0],
        Ec: vec![1.0, 1.0],
        Ed: vec![0.0, 0.0],
        Z: vec![29., 1.],
        m: vec![63.54, 1.0008],
        interaction_index: vec![0, 0],
        surface_binding_model: SurfaceBindingModel::PLANAR{calculation: SurfaceBindingCalculation::TARGET},
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let thickness: f64 = 1000.;
    let depth: f64 = 1000.;

    let geometry_input_2D = geometry::Mesh2DInput {
        length_unit: "ANGSTROM".to_string(),
        points: vec![(0., -thickness/2.), (depth, -thickness/2.), (depth, thickness/2.), (0., thickness/2.)],
        triangles: vec![(0, 1, 2), (0, 2, 3)],
        densities: vec![vec![0.03, 0.03], vec![0.03, 0.03]],
        boundary: vec![0, 1, 2, 3],
        simulation_boundary_points: vec![(0., 1.1*thickness/2.), (depth, 1.1*thickness/2.), (depth, -1.1*thickness/2.), (0., -1.1*thickness/2.), (0., 1.1*thickness/2.)],
        energy_barrier_thickness: 10.,
        electronic_stopping_correction_factors: vec![1.0, 1.0],
    };

    let geometry_input_homogeneous_2D = geometry::HomogeneousMesh2DInput {
        length_unit: "ANGSTROM".to_string(),
        points: vec![(0., -thickness/2.), (depth, -thickness/2.), (depth, thickness/2.), (0., thickness/2.)],
        densities: vec![0.03, 0.03],
        simulation_boundary_points: vec![(0., 1.1*thickness/2.), (depth, 1.1*thickness/2.), (depth, -1.1*thickness/2.), (0., -1.1*thickness/2.), (0., 1.1*thickness/2.)],
        electronic_stopping_correction_factor: 1.0,
    };

    let geometry_input_1D = geometry::Mesh1DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: vec![vec![0.03, 0.03], vec![0.03, 0.03]],
        layer_thicknesses: vec![thickness/2., thickness/2.],
        electronic_stopping_correction_factors: vec![1.0, 1.0, 1.0]

    };

    let geometry_input_0D = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: vec![0.03, 0.03],
        electronic_stopping_correction_factor: 1.0
    };

    let material_2D: material::Material<geometry::Mesh2D> = material::Material::<geometry::Mesh2D>::new(&material_parameters, &geometry_input_2D);
    let material_homogeneous_2D: material::Material<geometry::HomogeneousMesh2D> = material::Material::<geometry::HomogeneousMesh2D>::new(&material_parameters, &geometry_input_homogeneous_2D);
    let material_1D: material::Material<geometry::Mesh1D> = material::Material::<geometry::Mesh1D>::new(&material_parameters, &geometry_input_1D);
    let mut material_0D: material::Material<geometry::Mesh0D> = material::Material::<geometry::Mesh0D>::new(&material_parameters, &geometry_input_0D);
    material_0D.geometry.energy_barrier_thickness = 10.*ANGSTROM;

    particle_1.pos.x = 500.*ANGSTROM;
    particle_1.pos.y = 0.;

    particle_1.pos_old.x = -500.*ANGSTROM;
    particle_1.pos_old.y = 0.;

    assert!(material_2D.inside(particle_1.pos.x, particle_1.pos.y, particle_1.pos.z));
    assert!(!material_2D.inside(particle_1.pos_old.x, particle_1.pos_old.y, particle_1.pos_old.z));

    assert!(material_homogeneous_2D.inside(particle_1.pos.x, particle_1.pos.y, particle_1.pos.z));
    assert!(!material_homogeneous_2D.inside(particle_1.pos_old.x, particle_1.pos_old.y, particle_1.pos_old.z));

    assert!(material_0D.inside(particle_1.pos.x, particle_1.pos.y, particle_1.pos.z));
    assert!(!material_0D.inside(particle_1.pos_old.x, particle_1.pos_old.y, particle_1.pos_old.z));

    assert_eq!(material_2D.geometry.get_ck(0., 0., 0.), material_0D.geometry.get_ck(0., 0., 0.));
    assert_eq!(material_2D.geometry.get_densities(0., 0., 0.), material_0D.geometry.get_densities(0., 0., 0.));
    assert_eq!(material_2D.geometry.get_total_density(0., 0., 0.), material_0D.geometry.get_total_density(0., 0., 0.));
    assert_eq!(material_2D.geometry.get_concentrations(0., 0., 0.), material_0D.geometry.get_concentrations(0., 0., 0.));
    assert_eq!(material_2D.geometry.closest_point(-10., 0., 5.), material_0D.geometry.closest_point(-10., 0., 5.));
    assert_eq!(material_2D.geometry.get_densities(-10., 0., 5.), material_0D.geometry.get_densities(-10., 0., 5.));
    assert_eq!(material_2D.geometry.get_ck(-10., 0., 5.), material_0D.geometry.get_ck(-10., 0., 5.));

    assert_eq!(material_2D.geometry.get_ck(0., 0., 0.), material_homogeneous_2D.geometry.get_ck(0., 0., 0.));
    assert_eq!(material_2D.geometry.get_densities(0., 0., 0.), material_homogeneous_2D.geometry.get_densities(0., 0., 0.));
    assert_eq!(material_2D.geometry.get_total_density(0., 0., 0.), material_homogeneous_2D.geometry.get_total_density(0., 0., 0.));
    assert_eq!(material_2D.geometry.get_concentrations(0., 0., 0.), material_homogeneous_2D.geometry.get_concentrations(0., 0., 0.));
    assert_eq!(material_2D.geometry.closest_point(-10., 0., 5.), material_homogeneous_2D.geometry.closest_point(-10., 0., 5.));
    assert_eq!(material_2D.geometry.get_densities(-10., 0., 5.), material_homogeneous_2D.geometry.get_densities(-10., 0., 5.));
    assert_eq!(material_2D.geometry.get_ck(-10., 0., 5.), material_homogeneous_2D.geometry.get_ck(-10., 0., 5.));
}

#[test]
fn test_surface_binding_energy_barrier() {
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
    let mut particle_1 = particle::Particle::new(mass, Z, E, Ec, Es, 0.0, x, y, z, cosx, cosy, cosz, false, false, 0);

    let material_parameters = material::MaterialParameters{
        energy_unit: "EV".to_string(),
        mass_unit: "AMU".to_string(),
        Eb: vec![0.0, 0.0],
        Es: vec![2.0, 4.0],
        Ec: vec![1.0, 1.0],
        Ed: vec![0.0, 0.0],
        Z: vec![29., 1.],
        m: vec![63.54, 1.0008],
        interaction_index: vec![0, 0],
        surface_binding_model: SurfaceBindingModel::PLANAR{calculation: SurfaceBindingCalculation::TARGET},
        bulk_binding_model: BulkBindingModel::INDIVIDUAL,
    };

    let thickness: f64 = 1000.;
    let depth: f64 = 1000.;

    let geometry_input_2D = geometry::Mesh2DInput {
        length_unit: "ANGSTROM".to_string(),
        points: vec![(0., -thickness/2.), (depth, -thickness/2.), (depth, thickness/2.), (0., thickness/2.)],
        triangles: vec![(0, 1, 2), (0, 2, 3)],
        densities: vec![vec![0.03, 0.03], vec![0.03, 0.03]],
        boundary: vec![0, 1, 2, 3],
        simulation_boundary_points: vec![(0., 1.1*thickness/2.), (depth, 1.1*thickness/2.), (depth, -1.1*thickness/2.), (0., -1.1*thickness/2.), (0., 1.1*thickness/2.)],
        energy_barrier_thickness: 10.,
        electronic_stopping_correction_factors: vec![1.0, 1.0],
    };

    let geometry_input_0D = geometry::Mesh0DInput {
        length_unit: "ANGSTROM".to_string(),
        densities: vec![0.03, 0.03],
        electronic_stopping_correction_factor: 1.0
    };

    let material_2D: material::Material<geometry::Mesh2D> = material::Material::<geometry::Mesh2D>::new(&material_parameters, &geometry_input_2D);
    let mut material_0D: material::Material<geometry::Mesh0D> = material::Material::<geometry::Mesh0D>::new(&material_parameters, &geometry_input_0D);
    material_0D.geometry.energy_barrier_thickness = 10.*ANGSTROM;

    particle_1.pos.x = 500.*ANGSTROM;
    particle_1.pos.y = 0.;

    particle_1.pos_old.x = -500.*ANGSTROM;
    particle_1.pos_old.y = 0.;

    let inside_0D = material_0D.inside(particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);
    let inside_old_0D = material_0D.inside(particle_1.pos_old.x, particle_1.pos_old.y, particle_1.pos_old.z);

    let inside_2D = material_2D.inside(particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);
    let inside_old_2D = material_2D.inside(particle_1.pos_old.x, particle_1.pos_old.y, particle_1.pos_old.z);

    //println!("{} {}", inside, inside_old);
    assert!(inside_0D);
    assert!(!inside_old_0D);

    assert!(inside_2D);
    assert!(!inside_old_2D);

    //Test concentration-dependent surface binding energy
    let surface_binding_energy_2D = material_2D.actual_surface_binding_energy(&particle_1, particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);
    let surface_binding_energy_0D = material_0D.actual_surface_binding_energy(&particle_1, particle_1.pos.x, particle_1.pos.y, particle_1.pos.z);

    assert!(approx_eq!(f64, surface_binding_energy_0D/EV, (2. + 4.)/2., epsilon=1E-24));
    assert!(approx_eq!(f64, surface_binding_energy_2D/EV, (2. + 4.)/2., epsilon=1E-24));

    println!("sbv 0D: {}", surface_binding_energy_0D/EV);
    println!("sbv 2D: {}", surface_binding_energy_2D/EV);

    //Test leftmost boundary
    assert!(material_2D.inside_energy_barrier(500.*ANGSTROM, 0., 0.));
    assert!(material_2D.inside_energy_barrier(-5.*ANGSTROM, 0., 0.));
    assert!(!material_2D.inside_energy_barrier(-15.*ANGSTROM, 0., 0.));

    assert!(material_0D.inside_energy_barrier(500.*ANGSTROM, 0., 0.));
    assert!(material_0D.inside_energy_barrier(-5.*ANGSTROM, 0., 0.));
    assert!(!material_0D.inside_energy_barrier(-15.*ANGSTROM, 0., 0.));

    assert!(material_0D.inside_energy_barrier(500.*ANGSTROM, 0., 0.));
    assert!(!material_0D.inside_energy_barrier(-500.*ANGSTROM, 0., 0.));
    assert!(material_0D.inside_energy_barrier(-9.*ANGSTROM, 0., 0.));
    assert!(!material_0D.inside_energy_barrier(-11.*ANGSTROM, 0., 0.));

    //Test top boundary
    assert!(material_2D.inside_energy_barrier(500.*ANGSTROM, 0., 0.));
    assert!(material_2D.inside_energy_barrier(500.*ANGSTROM, 505.*ANGSTROM, 0.));
    assert!(!material_2D.inside_energy_barrier(500.*ANGSTROM, 515.*ANGSTROM, 0.));

    assert!(material_0D.inside_energy_barrier(500.*ANGSTROM, 0., 0.));
    assert!(material_0D.inside_energy_barrier(500.*ANGSTROM, 505.*ANGSTROM, 0.));
    assert!(material_0D.inside_energy_barrier(500.*ANGSTROM, 515.*ANGSTROM, 0.));

    //Test bottom boundary
    assert!(material_2D.inside_energy_barrier(500.*ANGSTROM, -505.*ANGSTROM, 0.));
    assert!(!material_2D.inside_energy_barrier(500.*ANGSTROM, -515.*ANGSTROM, 0.));

    assert!(material_0D.inside_energy_barrier(500.*ANGSTROM, -505.*ANGSTROM, 0.));
    assert!(material_0D.inside_energy_barrier(500.*ANGSTROM, -515.*ANGSTROM, 0.));

    //Test rightmost boundary
    assert!(material_2D.inside_energy_barrier(1005.*ANGSTROM, 0., 0.));
    assert!(!material_2D.inside_energy_barrier(1015.*ANGSTROM, 0., 0.));

    assert!(material_0D.inside_energy_barrier(1005.*ANGSTROM, 0., 0.));
    assert!(material_0D.inside_energy_barrier(1015.*ANGSTROM, 0., 0.));
}

#[test]
fn test_triangle_contains() {
    let triangle_1 = geometry::Triangle2D::new((0., 2., 0., 2., 0., 0.));
    assert!(triangle_1.contains(0.5, 0.5));
    assert!(!triangle_1.contains(2., 2.));

    let triangle_2 = geometry::Triangle2D::new((-2., 0., 0., 0., 0., -2.));
    assert!(triangle_2.contains(-0.5, -0.5));
    assert!(!triangle_2.contains(0.5, 0.5));
    assert!(!triangle_2.contains(-2., -2.));
}

#[test]
fn test_triangle_distance_to() {
    let triangle_1 = geometry::Triangle2D::new((0., 2., 0., 2., 0., 0.));
    assert!(approx_eq!(f64, triangle_1.distance_to(-2., 0.), 2., epsilon=1E-12), "{}", triangle_1.distance_to(-2., 0.));

    assert!(approx_eq!(f64, triangle_1.distance_to(2., 2.), (2.0_f64).sqrt(), epsilon=1E-12), "{}", triangle_1.distance_to(2., 2.));

    assert!(approx_eq!(f64, triangle_1.distance_to(0., 0.), 0., epsilon=1E-12), "{}", triangle_1.distance_to(0., 0.));
    assert!(approx_eq!(f64, triangle_1.distance_to(2., 0.), 0., epsilon=1E-12), "{}", triangle_1.distance_to(2., 0.));
    assert!(approx_eq!(f64, triangle_1.distance_to(0., 2.), 0., epsilon=1E-12), "{}", triangle_1.distance_to(0., 2.));
}

#[test]
fn test_surface_refraction() {
    let print_output = false;

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
    let mut particle_1 = particle::Particle::new(mass, Z, E, Ec, Es, 0.0, x, y, z, cosx, cosy, cosz, false, false, 0);

    //Test particle entering material and gaining energy

    //Eckstein's formulation for particle entering surface
    let cosx_new = ((E*cosx*cosx + Es)/(E + Es)).sqrt();
    let sinx = (1. - cosx*cosx).sqrt();
    let sinx_new = (1. - cosx_new*cosx_new).sqrt();
    let cosy_new = cosy*sinx_new/sinx;
    let cosz_new = cosz*sinx_new/sinx;
    let dir_mag = (cosx_new*cosx_new + cosy_new*cosy_new + cosz_new*cosz_new).sqrt();
    let cosx_new = cosx_new/dir_mag;
    let cosy_new = cosy_new/dir_mag;
    let cosz_new = cosz_new/dir_mag;

    let normal = Vector::new(1.0, 0.0, 0.0);
    particle::surface_refraction(&mut particle_1, normal, Es);

    if print_output {
        println!("dir_mag: {}", dir_mag);
        println!("{} {} {}", cosx, cosy, cosz);
        println!("{} {} {}", cosx_new, cosy_new, cosz_new);
        println!("{} {} {}", particle_1.dir.x, particle_1.dir.y, particle_1.dir.z);
        println!();
    }

    assert!(approx_eq!(f64, particle_1.dir.x, cosx_new, epsilon=1E-12));
    assert!(approx_eq!(f64, particle_1.dir.y, cosy_new, epsilon=1E-12));
    assert!(approx_eq!(f64, particle_1.dir.z, cosz_new, epsilon=1E-12));

    //Test particle leaving material and losing energy

    let cosx_new = ((particle_1.E*particle_1.dir.x*particle_1.dir.x - Es)/(particle_1.E - Es)).sqrt();
    let sinx = (1. - particle_1.dir.x*particle_1.dir.x).sqrt();
    let sinx_new = (1. - cosx_new*cosx_new).sqrt();
    let cosy_new = particle_1.dir.y*sinx_new/sinx;
    let cosz_new = particle_1.dir.z*sinx_new/sinx;

    if print_output {
        println!("{} {} {}", particle_1.dir.x, particle_1.dir.y, particle_1.dir.z);
    }

    let normal = Vector::new(1.0, 0.0, 0.0);
    particle::surface_refraction(&mut particle_1, normal, -Es);


    if print_output {
        println!("{} {} {}", cosx_new, cosy_new, cosz_new);
        println!("{} {} {}", particle_1.dir.x, particle_1.dir.y, particle_1.dir.z);
    }

    assert!(approx_eq!(f64, particle_1.dir.x, cosx_new, epsilon=1E-9));
    assert!(approx_eq!(f64, particle_1.dir.y, cosy_new, epsilon=1E-9));
    assert!(approx_eq!(f64, particle_1.dir.z, cosz_new, epsilon=1E-9));

}

#[test]
fn test_momentum_conservation() {

    for energy_eV in vec![1., 10., 100., 1000., 1000.] {
        //Aluminum
        let m1 = 183.8*AMU;
        let Z1 = 74.;
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

        //Arbitrary initial angle
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
            Ed: vec![0.0],
            Z: vec![Z2],
            m: vec![m2],
            interaction_index: vec![0],
            surface_binding_model: SurfaceBindingModel::PLANAR{calculation: SurfaceBindingCalculation::TARGET},
            bulk_binding_model: BulkBindingModel::INDIVIDUAL,
        };

        let thickness: f64 = 1000.;
        let depth: f64 = 1000.;
        let geometry_input = geometry::Mesh2DInput {
            length_unit: "ANGSTROM".to_string(),
            triangles: vec![(0, 1, 2), (0, 2, 3)],
            points: vec![(0., -thickness/2.), (depth, -thickness/2.), (depth, thickness/2.), (0., thickness/2.)],
            densities: vec![vec![0.06306], vec![0.06306]],
            boundary: vec![0, 1, 2, 3],
            simulation_boundary_points: vec![(0., 1.1*thickness/2.), (depth, 1.1*thickness/2.), (depth, -1.1*thickness/2.), (0., -1.1*thickness/2.), (0., 1.1*thickness/2.)],
            electronic_stopping_correction_factors: vec![0.0, 0.0],
            energy_barrier_thickness: 0.,
        };

        let material_1: material::Material<geometry::Mesh2D> = material::Material::<geometry::Mesh2D>::new(&material_parameters, &geometry_input);

        for high_energy_free_flight_paths in vec![true, false] {
            for potential in vec![InteractionPotential::KR_C, InteractionPotential::MOLIERE, InteractionPotential::ZBL, InteractionPotential::LENZ_JENSEN] {
                for scattering_integral in vec![ScatteringIntegral::MENDENHALL_WELLER, ScatteringIntegral::GAUSS_MEHLER{n_points: 10}, ScatteringIntegral::GAUSS_LEGENDRE] {
                    for root_finder in vec![Rootfinder::NEWTON{max_iterations: 100, tolerance: 1E-3}] {

                        //println!("Case: {} {} {} {}", energy_eV, high_energy_free_flight_paths, potential, scattering_integral);

                        let mut particle_1 = particle::Particle::new(m1, Z1, E1, Ec1, Es1, 0.0, x1, y1, z1, cosx, cosy, cosz, false, false, 0);

                        #[cfg(not(feature = "distributions"))]
                        let options = Options {
                            name: "test".to_string(),
                            track_trajectories: false,
                            track_recoils: true,
                            track_recoil_trajectories: false,
                            write_buffer_size: 8000,
                            weak_collision_order: 0,
                            suppress_deep_recoils: false,
                            high_energy_free_flight_paths: high_energy_free_flight_paths,
                            electronic_stopping_mode: ElectronicStoppingMode::INTERPOLATED,
                            mean_free_path_model: MeanFreePathModel::LIQUID,
                            interaction_potential: vec![vec![potential]],
                            scattering_integral: vec![vec![scattering_integral]],
                            num_threads: 1,
                            num_chunks: 1,
                            use_hdf5: false,
                            root_finder: vec![vec![root_finder]],
                            track_displacements: false,
                            track_energy_losses: false,
                        };

                        #[cfg(feature = "distributions")]
                        let options = Options {
                            name: "test".to_string(),
                            track_trajectories: false,
                            track_recoils: true,
                            track_recoil_trajectories: false,
                            write_buffer_size: 8000,
                            weak_collision_order: 0,
                            suppress_deep_recoils: false,
                            high_energy_free_flight_paths: high_energy_free_flight_paths,
                            electronic_stopping_mode: ElectronicStoppingMode::INTERPOLATED,
                            mean_free_path_model: MeanFreePathModel::LIQUID,
                            interaction_potential: vec![vec![potential]],
                            scattering_integral: vec![vec![scattering_integral]],
                            num_threads: 1,
                            num_chunks: 1,
                            use_hdf5: false,
                            root_finder: vec![vec![root_finder]],
                            track_displacements: false,
                            track_energy_losses: false,
                            energy_min: 0.0,
                            energy_max: 10.0,
                            energy_num: 11,
                            angle_min: 0.0,
                            angle_max: 90.0,
                            angle_num: 11,
                            x_min: 0.0,
                            y_min: -10.0,
                            z_min: -10.0,
                            x_max: 10.0,
                            y_max: 10.0,
                            z_max: 10.0,
                            x_num: 11,
                            y_num: 11,
                            z_num: 11,
                        };

                        let binary_collision_geometries = bca::determine_mfp_phi_impact_parameter(&mut particle_1, &material_1, &options);

                        /*
                        println!("Phi: {} rad p: {} Angstrom mfp: {} Angstrom", binary_collision_geometries[0].phi_azimuthal,
                            binary_collision_geometries[0].impact_parameter/ANGSTROM,
                            binary_collision_geometries[0].mfp/ANGSTROM);
                        */

                        let (species_index, mut particle_2) = bca::choose_collision_partner(&mut particle_1, &material_1,
                            &binary_collision_geometries[0], &options);

                        let mom1_0 = particle_1.get_momentum();
                        let mom2_0 = particle_2.get_momentum();

                        let initial_momentum = mom1_0.add(&mom2_0);

                        let binary_collision_result = bca::calculate_binary_collision(&particle_1,
                            &particle_2, &binary_collision_geometries[0], &options).unwrap();

                        /*
                        println!("E_recoil: {} eV Psi: {} rad Psi_recoil: {} rad", binary_collision_result.recoil_energy/EV,
                            binary_collision_result.psi,
                            binary_collision_result.psi_recoil);

                        println!("Initial Energies: {} eV {} eV", particle_1.E/EV, particle_2.E/EV);
                        */
                        //Energy transfer to recoil
                        particle_2.E = binary_collision_result.recoil_energy - material_1.average_bulk_binding_energy(particle_2.pos.x, particle_2.pos.y, particle_2.pos.z);

                        //Rotate particle 1, 2 by lab frame scattering angles
                        particle_1.rotate(binary_collision_result.psi,
                            binary_collision_geometries[0].phi_azimuthal);

                        particle_2.rotate(-binary_collision_result.psi_recoil,
                            binary_collision_geometries[0].phi_azimuthal);

                        //Subtract total energy from all simultaneous collisions and electronic stopping
                        particle_1.E += -binary_collision_result.recoil_energy;
                        bca::subtract_electronic_stopping_energy(&mut particle_1, &material_1, 0.,
                            0., particle_2.Z, species_index, &options);

                        let mom1_1 = particle_1.get_momentum();
                        let mom2_1 = particle_2.get_momentum();

                        let final_momentum = mom1_1.add(&mom2_1);

                        /*
                        println!("Final Energies: {} eV {} eV", particle_1.E/EV, particle_2.E/EV);
                        println!("X: {} {} {}% Error", initial_momentum.x/ANGSTROM/AMU, final_momentum.x/ANGSTROM/AMU, 100.*(final_momentum.x - initial_momentum.x)/initial_momentum.magnitude());
                        println!("Y: {} {} {}% Error", initial_momentum.y/ANGSTROM/AMU, final_momentum.y/ANGSTROM/AMU, 100.*(final_momentum.y - initial_momentum.y)/initial_momentum.magnitude());
                        println!("Z: {} {} {}% Error", initial_momentum.z/ANGSTROM/AMU, final_momentum.z/ANGSTROM/AMU, 100.*(final_momentum.z - initial_momentum.z)/initial_momentum.magnitude());
                        println!();
                        */

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
}

#[test]
fn test_rotate() {
    let mass = 1.;
    let Z = 1.;
    let E = 1.;
    let Ec = 1.;
    let Es = 1.;
    let Ed = 0.0;
    let x = 0.;
    let y = 0.;
    let z = 0.;
    let cosx = (PI/4.).cos();
    let cosy = (PI/4.).sin();
    let cosz = 0.;
    let psi = -PI/4.;
    let phi = 0.;

    let mut particle = particle::Particle::new(mass, Z, E, Ec, Es, Ed, x, y, z, cosx, cosy, cosz, false, false, 0);

    //Check that rotation in 2D works
    particle.rotate(psi, phi);
    assert!(approx_eq!(f64, particle.dir.x, 0., epsilon = 1E-12), "particle.dir.x: {} Should be ~0.", particle.dir.x);
    assert!(approx_eq!(f64, particle.dir.y, 1., epsilon = 1E-12), "particle.dir.y: {} Should be ~1.", particle.dir.y);

    //Check that rotating back by negative psi returns to the previous values
    particle.rotate(-psi, phi);
    assert!(approx_eq!(f64, particle.dir.x, cosx, epsilon = 1E-12), "particle.dir.x: {} Should be ~{}", particle.dir.x, cosx);
    assert!(approx_eq!(f64, particle.dir.y, cosy, epsilon = 1E-12), "particle.dir.y: {} Should be ~{}", particle.dir.y, cosy);

    //Check that azimuthal rotation by 180 degrees works correctly
    let phi = PI;
    particle.rotate(psi, phi);
    assert!(approx_eq!(f64, particle.dir.x, 1., epsilon = 1E-12), "particle.dir.x: {} Should be ~1.", particle.dir.x);
    assert!(approx_eq!(f64, particle.dir.y, 0., epsilon = 1E-12), "particle.dir.y: {} Should be ~0.", particle.dir.y);

    //Check that particle direction vector remains normalized following rotations
    assert!(approx_eq!(f64, particle.dir.x.powi(2) + particle.dir.y.powi(2) + particle.dir.z.powi(2), 1.), "Particle direction not normalized.");

}

#[test]
fn test_particle_advance() {
    let mass = 1.;
    let Z = 1.;
    let E = 1.;
    let Ec = 1.;
    let Es = 1.;
    let Ed = 0.0;
    let x = 0.;
    let y = 0.;
    let z = 0.;
    let cosx = (PI/4.).cos();
    let cosy = (PI/4.).sin();
    let cosz = 0.;
    let mfp = 1.;
    let asymptotic_deflection = 0.5;

    let mut particle = particle::Particle::new(mass, Z, E, Ec, Es, Ed, x, y, z, cosx, cosy, cosz, false, false, 0);

    let distance_traveled = particle.advance(mfp, asymptotic_deflection);

    assert_eq!(particle.pos.x, (1. - 0.5)*cosx);
    assert_eq!(particle.pos.y, (1. - 0.5)*cosy);
    assert_eq!(particle.pos.z, 0.);
    assert_eq!(distance_traveled, mfp - asymptotic_deflection);
}

#[test]
fn test_quadrature() {
    let Za = 1.;
    let Zb = 13.;
    let Ma = 1.008;
    let Mb = 26.9815385;
    let E0 = 10.*EV;
    let p = 1.*ANGSTROM;
    let a = interactions::screening_length(Za, Zb, InteractionPotential::KR_C);

    #[cfg(not(feature = "distributions"))]
    let options = Options {
        name: "test".to_string(),
        track_trajectories: false,
        track_recoils: true,
        track_recoil_trajectories: false,
        write_buffer_size: 8000,
        weak_collision_order: 0,
        suppress_deep_recoils: false,
        high_energy_free_flight_paths: false,
        electronic_stopping_mode: ElectronicStoppingMode::INTERPOLATED,
        mean_free_path_model: MeanFreePathModel::LIQUID,
        interaction_potential:  vec![vec![InteractionPotential::KR_C]],
        scattering_integral: vec![vec![ScatteringIntegral::MENDENHALL_WELLER]],
        num_threads: 1,
        num_chunks: 1,
        use_hdf5: false,
        root_finder: vec![vec![Rootfinder::NEWTON{max_iterations: 100, tolerance: 1E-14}]],
        track_displacements: false,
        track_energy_losses: false,
    };

    #[cfg(feature = "distributions")]
    let options = Options {
        name: "test".to_string(),
        track_trajectories: false,
        track_recoils: true,
        track_recoil_trajectories: false,
        write_buffer_size: 8000,
        weak_collision_order: 0,
        suppress_deep_recoils: false,
        high_energy_free_flight_paths: false,
        electronic_stopping_mode: ElectronicStoppingMode::INTERPOLATED,
        mean_free_path_model: MeanFreePathModel::LIQUID,
        interaction_potential:  vec![vec![InteractionPotential::KR_C]],
        scattering_integral: vec![vec![ScatteringIntegral::MENDENHALL_WELLER]],
        num_threads: 1,
        num_chunks: 1,
        use_hdf5: false,
        root_finder: vec![vec![Rootfinder::NEWTON{max_iterations: 100, tolerance: 1E-14}]],
        track_displacements: false,
        track_energy_losses: false,
        energy_min: 0.0,
        energy_max: 10.0,
        energy_num: 11,
        angle_min: 0.0,
        angle_max: 90.0,
        angle_num: 11,
        x_min: 0.0,
        y_min: -10.0,
        z_min: -10.0,
        x_max: 10.0,
        y_max: 10.0,
        z_max: 10.0,
        x_num: 11,
        y_num: 11,
        z_num: 11,
    };

    let x0_newton = bca::newton_rootfinder(Za, Zb, Ma, Mb, E0, p, InteractionPotential::KR_C, 100, 1E-12).unwrap();

    //If cpr_rootfinder is enabled, compare Newton to CPR - they should be nearly identical
    #[cfg(feature = "cpr_rootfinder")]
    if let Ok(x0_cpr) = bca::cpr_rootfinder(Za, Zb, Ma, Mb, E0, p, InteractionPotential::KR_C, 2, 10000, 1E-6, 1E-6, 1E-9, 1E9, 1E-13, true) {
        println!("CPR: {} Newton: {}", x0_cpr, x0_newton);
        assert!(approx_eq!(f64, x0_newton, x0_cpr, epsilon=1E-3));
    };

    //Compute center of mass deflection angle with each algorithm
    let theta_gm = bca::gauss_mehler(Za, Zb, Ma, Mb, E0, p, x0_newton, InteractionPotential::KR_C, 10);
    let theta_gl = bca::gauss_legendre(Za, Zb, Ma, Mb, E0, p, x0_newton, InteractionPotential::KR_C);
    let theta_mw = bca::mendenhall_weller(Za, Zb, Ma, Mb, E0, p, x0_newton, InteractionPotential::KR_C);
    let theta_magic = bca::magic(Za, Zb, Ma, Mb, E0, p, x0_newton, InteractionPotential::KR_C);

    //Gauss-Mehler and Gauss-Legendre should be very close to each other
    assert!(approx_eq!(f64, theta_gm, theta_gl, epsilon=0.001));
    assert!(approx_eq!(f64, theta_gm, theta_mw, epsilon=0.001));
    assert!(approx_eq!(f64, theta_gm, theta_magic, epsilon=0.15));

    println!("Gauss-Mehler: {} Gauss-Legendre: {} Mendenhall-Weller: {} MAGIC: {}",
        theta_gm, theta_gl, theta_mw, theta_magic);
}