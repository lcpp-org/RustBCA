use super::*;

use itertools::Itertools;

use parry2d_f64::shape::{Polyline, TriMesh, FeatureId::*};
use parry2d_f64::query::PointQuery;
use parry2d_f64::math::Point as Point2d;
use parry2d_f64::math::Isometry;

///Trait for a Geometry object - all forms of geometry must implement these traits to be used
pub trait Geometry {

    type InputFileFormat: InputFile + Clone;

    fn new(input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> Self;
    fn get_densities(&self,  x: f64, y: f64, z: f64) -> &Vec<f64>;
    fn get_ck(&self,  x: f64, y: f64, z: f64) -> f64;
    fn get_total_density(&self,  x: f64, y: f64, z: f64) -> f64;
    fn get_concentrations(&self, x: f64, y: f64, z: f64) -> &Vec<f64>;
    fn inside(&self, x: f64, y: f64, z: f64) -> bool;
    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool;
    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool;
    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64);
    fn nearest_normal_vector(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64);
}

pub trait GeometryElement: Clone {
    fn contains(&self, x: f64, y: f64, z: f64) -> bool;
    fn distance_to(&self, x: f64, y: f64, z: f64) -> f64;
    fn get_densities(&self) -> &Vec<f64>;
    fn get_concentrations(&self) -> &Vec<f64>;
    fn get_electronic_stopping_correction_factor(&self) -> f64;
}

#[derive(Deserialize, Clone)]
pub struct Mesh0DInput {
    pub length_unit: String,
    pub densities: Vec<f64>,
    pub electronic_stopping_correction_factor: f64,
}

#[derive(Clone)]
pub struct Mesh0D {
    pub densities: Vec<f64>,
    pub concentrations: Vec<f64>,
    pub electronic_stopping_correction_factor: f64,
    pub energy_barrier_thickness: f64,
}

impl Geometry for Mesh0D {

    type InputFileFormat = Input0D;

    fn new(input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> Mesh0D {

        let length_unit: f64 = match input.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "MM" => MM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => input.length_unit.parse()
                .expect(format!(
                        "Input errror: could nor parse length unit {}. Use a valid float or one of ANGSTROM, NM, MICRON, CM, MM, M",
                        &input.length_unit.as_str()
                    ).as_str()),
        };

        let electronic_stopping_correction_factor = input.electronic_stopping_correction_factor;

        let densities: Vec<f64> = input.densities.iter().map(|element| element/(length_unit).powi(3)).collect();
        assert!(densities.len() > 0, "Input Error: density list empty.");

        let total_density: f64 = densities.iter().sum();

        let energy_barrier_thickness = total_density.powf(-1./3.)/SQRTPI*2.;

        let concentrations: Vec<f64> = densities.iter().map(|&density| density/total_density).collect::<Vec<f64>>();

        Mesh0D {
            densities,
            concentrations,
            electronic_stopping_correction_factor,
            energy_barrier_thickness
        }
    }

    fn get_densities(&self,  x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.densities
    }
    fn get_ck(&self,  x: f64, y: f64, z: f64) -> f64 {
        self.electronic_stopping_correction_factor
    }
    fn get_total_density(&self,  x: f64, y: f64, z: f64) -> f64{
        self.densities.iter().sum()
    }
    fn get_concentrations(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.concentrations
    }
    fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        x > 0.0
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        x > -10.*self.energy_barrier_thickness
    }

    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        x > -self.energy_barrier_thickness
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        (0., y, z)
    }

    fn nearest_normal_vector(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        (-1., 0., 0.)
    }
}

#[derive(Deserialize, Clone)]
pub struct Mesh1DInput {
    pub length_unit: String,
    pub layer_thicknesses: Vec<f64>,
    pub densities: Vec<Vec<f64>>,
    pub electronic_stopping_correction_factors: Vec<f64>,
}

#[derive(Clone)]
pub struct Mesh1D {
    layers: Vec<Layer1D>,
    top: f64,
    bottom: f64,
    pub top_energy_barrier_thickness: f64,
    pub bottom_energy_barrier_thickness: f64,
}

impl Geometry for Mesh1D {

    type InputFileFormat = input::Input1D;

    fn new(geometry_input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> Self {

        let layer_thicknesses = geometry_input.layer_thicknesses.clone();
        let electronic_stopping_correction_factors = geometry_input.electronic_stopping_correction_factors.clone();
        assert!(electronic_stopping_correction_factors.len() > 0, "Input Error: Electronic stopping correction factor list empty.");
        let n = layer_thicknesses.len();

        let mut layers: Vec<Layer1D> =  Vec::with_capacity(n);

        //Multiply all coordinates by value of geometry unit.
        let length_unit: f64 = match geometry_input.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "MM" => MM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => geometry_input.length_unit.parse()
                .expect(format!(
                        "Input errror: could nor parse length unit {}. Use a valid float or one of ANGSTROM, NM, MICRON, CM, MM, M", &geometry_input.length_unit.as_str()
                    ).as_str()),
        };

        let densities: Vec<Vec<f64>> = geometry_input.densities
            .iter()
            .map( |row| row.iter().map(|element| element/(length_unit).powi(3)).collect() ).collect();

        //Assert all layer density lists are equal length
        assert!(
            densities.iter()
            .all(|density_list| density_list.len() == densities[0].len()),
            "Input error: all layer density lists must be the same size."
        );
        assert_eq!(layer_thicknesses.len(), densities.len(), "Input error: coordinates and data of unequal length.");

        let mut layer_top = 0.;
        let mut layer_bottom = 0.;
        for ((layer_thickness, densities), ck) in layer_thicknesses.iter().zip(densities).zip(electronic_stopping_correction_factors) {

            layer_bottom += layer_thickness*length_unit;

            let total_density: f64 = densities.iter().sum();
            let concentrations: Vec<f64> = densities.iter().map(|&density| density/total_density).collect::<Vec<f64>>();

            layers.push(Layer1D::new(layer_top, layer_bottom, densities, concentrations, ck));
            layer_top = layer_bottom;
        }

        let top_energy_barrier_thickness = layers[0].densities.iter().sum::<f64>().powf(-1./3.)/SQRTPI*2.;
        let bottom_energy_barrier_thickness = layers[layers.len() - 1].densities.iter().sum::<f64>().powf(-1./3.)/SQRTPI*2.;

        Mesh1D {
            layers,
            top: 0.0,
            bottom: layer_bottom,
            top_energy_barrier_thickness,
            bottom_energy_barrier_thickness,
        }
    }
    fn get_densities(&self,  x: f64, y: f64, z: f64) -> &Vec<f64> {
        if self.inside(x, y, z) {
            for layer in &self.layers {
                if layer.contains(x, y, z) {
                    return &layer.densities
                }
            }
            panic!("Geometry error: method inside() is returning true for points outside all layers. Check geometry.")
        } else if x < self.top {
            &self.layers[0].densities
        } else {
            &self.layers[self.layers.len() - 1].densities
        }
    }
    fn get_ck(&self,  x: f64, y: f64, z: f64) -> f64 {
        if self.inside(x, y, z) {
            for layer in &self.layers {
                if layer.contains(x, y, z) {
                    return layer.electronic_stopping_correction_factor
                }
            }
            panic!("Geometry error: method inside() is returning true for points outside all layers. Check geometry.")
        } else if x < self.top {
            self.layers[0].electronic_stopping_correction_factor
        } else {
            self.layers[self.layers.len() - 1].electronic_stopping_correction_factor
        }
    }
    fn get_total_density(&self,  x: f64, y: f64, z: f64) -> f64 {
        if self.inside(x, y, z) {
            for layer in &self.layers {
                if layer.contains(x, y, z) {
                    return layer.densities.iter().sum()
                }
            }
            panic!("Geometry error: method inside() is returning true for points outside all layers. Check geometry.")
        } else if x < self.top {
            self.layers[0].densities.iter().sum()
        } else {
            self.layers[self.layers.len() - 1].densities.iter().sum()
        }
    }
    fn get_concentrations(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        if self.inside(x, y, z) {
            for layer in &self.layers {
                if layer.contains(x, y, z) {
                    return &layer.concentrations
                }
            }
            panic!("Geometry error: method inside() is returning true for points outside all layers. Check geometry.")
        } else if x < self.top {
            &self.layers[0].concentrations
        } else {
            &self.layers[self.layers.len() - 1].concentrations
        }
    }
    fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        (x > self.top) & (x < self.bottom)
    }
    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        (x > self.top - 10.*self.top_energy_barrier_thickness) & (x < self.bottom + 10.*self.bottom_energy_barrier_thickness)
    }
    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        (x > self.top - self.top_energy_barrier_thickness) & (x < self.bottom + self.bottom_energy_barrier_thickness)
    }
    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let distance_from_bottom = (x - self.bottom).abs();
        let distance_from_top = (x - self.top).abs();
        if distance_from_bottom < distance_from_top {
            (self.bottom, y, z)
        } else {
            (self.top, y, z)
        }
    }

    fn nearest_normal_vector(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        (-1., 0., 0.)
    }
}

#[derive(Deserialize, Clone)]
pub struct ParryHomogeneousMesh2DInput {
    pub length_unit: String,
    pub points: Vec<(f64, f64)>,
    pub simulation_boundary_points: Vec<(f64, f64)>,
    pub densities: Vec<f64>,
    pub electronic_stopping_correction_factor: f64
}

#[derive(Clone)]
pub struct ParryHomogeneousMesh2D {
    pub boundary: Polyline,
    pub simulation_boundary: Polyline,
    pub energy_barrier_thickness: f64,
    pub densities: Vec<f64>,
    pub electronic_stopping_correction_factor: f64,
    pub concentrations: Vec<f64>
}


impl Geometry for ParryHomogeneousMesh2D {

    type InputFileFormat = InputHomogeneous2D;

    fn new(input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> Self {
        let length_unit: f64 = match input.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "MM" => MM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => input.length_unit.parse()
                .expect(format!(
                        "Input errror: could nor parse length unit {}. Use a valid float or one of ANGSTROM, NM, MICRON, CM, MM, M",
                        &input.length_unit.as_str()
                    ).as_str()),
        };

        let mut boundary_points_converted: Vec<Point2d<f64>> = input.points.iter().map(|(x, y)| Point2d::new(x*length_unit, y*length_unit)).collect();
        let number_boundary_points = boundary_points_converted.len() as u32;

        for p in boundary_points_converted.iter().combinations(2) {
            assert!(p[0] != p[1], "Input error: duplicate vertices in boundary geometry input. Boundary must be defined ccw with no duplicate points (terminal segment will span from final to initial point).")
        }

        let mut simulation_boundary_points_converted: Vec<Point2d<f64>> = input.simulation_boundary_points.iter().map(|(x, y)| Point2d::new(x*length_unit, y*length_unit)).collect();
        let number_simulation_boundary_points = simulation_boundary_points_converted.len() as u32;

        let test_simulation_boundary_ccw = (0..number_simulation_boundary_points as usize)
        .map(|i| (simulation_boundary_points_converted[(i + 1) % number_simulation_boundary_points as usize].x - simulation_boundary_points_converted[i].x)*(simulation_boundary_points_converted[i].y + simulation_boundary_points_converted[(i + 1) % number_simulation_boundary_points as usize].y))
        .sum::<f64>();

        if test_simulation_boundary_ccw > 0.0 {
            simulation_boundary_points_converted = simulation_boundary_points_converted.clone().into_iter().rev().collect::<Vec<Point2d<f64>>>();
        }

        for p in simulation_boundary_points_converted.iter().combinations(2) {
            assert!(p[0] != p[1], "Input error: duplicate vertices in simulation boundary geometry input. Boundary must be defined ccw with no duplicate points (terminal segment will span from final to initial point).")
        }

        let test_ccw = (0..number_boundary_points as usize)
                .map(|i| (boundary_points_converted[(i + 1) % number_boundary_points as usize].x - boundary_points_converted[i].x)*(boundary_points_converted[i].y + boundary_points_converted[(i + 1) % number_boundary_points as usize].y))
                .sum::<f64>();

        if test_ccw > 0.0 {
            boundary_points_converted = boundary_points_converted.clone().into_iter().rev().collect::<Vec<Point2d<f64>>>();
        }

        let electronic_stopping_correction_factor = input.electronic_stopping_correction_factor;

        let densities: Vec<f64> = input.densities.iter().map(|element| element/(length_unit).powi(3)).collect();

        let total_density: f64 = densities.iter().sum();

        let energy_barrier_thickness = total_density.powf(-1./3.)/SQRTPI*2.;

        let concentrations: Vec<f64> = densities.iter().map(|&density| density/total_density).collect::<Vec<f64>>();

        let mut linked_boundary_points = (0..number_boundary_points).zip(1..number_boundary_points).map(|(x, y)| [x, y]).collect::<Vec<[u32; 2]>>();
        linked_boundary_points.push([number_boundary_points - 1, 0]);
        let boundary2 = Polyline::new(boundary_points_converted, Some(linked_boundary_points));

        let mut linked_simulation_boundary_points = (0..number_simulation_boundary_points).zip(1..number_simulation_boundary_points).map(|(x, y)| [x, y]).collect::<Vec<[u32; 2]>>();
        linked_simulation_boundary_points.push([number_simulation_boundary_points - 1, 0]);
        let simulation_boundary2 = Polyline::new(simulation_boundary_points_converted, Some(linked_simulation_boundary_points));

        ParryHomogeneousMesh2D {
            densities,
            simulation_boundary: simulation_boundary2,
            boundary: boundary2,
            electronic_stopping_correction_factor,
            energy_barrier_thickness,
            concentrations
        }
    }

    fn get_densities(&self,  x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.densities
    }

    fn get_ck(&self,  x: f64, y: f64, z: f64) -> f64 {
        self.electronic_stopping_correction_factor
    }

    fn get_total_density(&self,  x: f64, y: f64, z: f64) -> f64 {
        self.densities.iter().sum::<f64>()
    }

    fn get_concentrations(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.concentrations
    }

    fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point2d::new(x, y);
        if self.boundary.aabb(&Isometry::identity()).contains_local_point(&p) {
            let (point_projection, (_, _)) = self.boundary.project_local_point_assuming_solid_interior_ccw(p);
            point_projection.is_inside
        } else {
            false
        }
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point2d::new(x, y);

        let (point_projection, (_, _)) = self.simulation_boundary.project_local_point_assuming_solid_interior_ccw(p);
        point_projection.is_inside

    }
    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        if self.inside(x, y, z) {
            true
        } else {
            let p = Point2d::new(x, y);
            (self.boundary.distance_to_local_point(&p, true) as f64) < self.energy_barrier_thickness
        }
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let p = Point2d::new(x, y);
        let (point_projection, (_, _)) = self.boundary.project_local_point_assuming_solid_interior_ccw(p);
        let (x_, y_) = (point_projection.point.x, point_projection.point.y);
        (x_ as f64, y_ as f64, z)
    }

    fn nearest_normal_vector(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let p = Point2d::new(x, y);

        let (point_projection, (_, _)) = self.boundary.project_local_point_assuming_solid_interior_ccw(p);
        let (x_intersect, y_intersect, z_intersect) = self.closest_point(x, y, z);

        let dx = x_intersect - x;
        let dy = y_intersect - y;
        let dz = z_intersect - z;
        let mag = (dx*dx + dy*dy + dz*dz).sqrt();

        if point_projection.is_inside {
            (dx/mag, dy/mag, dz/mag)
        } else {
            (-dx/mag, -dy/mag, -dz/mag)
        }
    }
}

#[derive(Deserialize, Clone)]
pub struct HomogeneousMesh2DInput {
    pub length_unit: String,
    pub points: Vec<(f64, f64)>,
    pub simulation_boundary_points: Vec<(f64, f64)>,
    pub densities: Vec<f64>,
    pub electronic_stopping_correction_factor: f64
}

/*
#[derive(Clone)]
pub struct HomogeneousMesh2D {
    pub boundary: Polygon<f64>,
    pub simulation_boundary: Polygon<f64>,
    pub energy_barrier_thickness: f64,
    pub densities: Vec<f64>,
    pub electronic_stopping_correction_factor: f64,
    pub concentrations: Vec<f64>
}

impl Geometry for HomogeneousMesh2D {

    type InputFileFormat = InputHomogeneous2D;

    fn new(input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> Self {
        let length_unit: f64 = match input.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "MM" => MM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => input.length_unit.parse()
                .expect(format!(
                        "Input errror: could nor parse length unit {}. Use a valid float or one of ANGSTROM, NM, MICRON, CM, MM, M",
                        &input.length_unit.as_str()
                    ).as_str()),
        };

        let boundary_points_converted: Vec<(f64, f64)> = input.points.iter().map(|(x, y)| (x*length_unit, y*length_unit)).collect();

        let simulation_boundary_points_converted: Vec<(f64, f64)> = input.simulation_boundary_points.iter().map(|(x, y)| (x*length_unit, y*length_unit)).collect();

        let electronic_stopping_correction_factor = input.electronic_stopping_correction_factor;

        let densities: Vec<f64> = input.densities.iter().map(|element| element/(length_unit).powi(3)).collect();

        let total_density: f64 = densities.iter().sum();

        let energy_barrier_thickness = total_density.powf(-1./3.)/SQRTPI*2.;

        let concentrations: Vec<f64> = densities.iter().map(|&density| density/total_density).collect::<Vec<f64>>();

        let boundary = Polygon::new(LineString::from(boundary_points_converted), vec![]);

        let simulation_boundary =  Polygon::new(LineString::from(simulation_boundary_points_converted), vec![]);

        HomogeneousMesh2D {
            densities,
            simulation_boundary,
            boundary,
            electronic_stopping_correction_factor,
            energy_barrier_thickness,
            concentrations
        }
    }

    fn get_densities(&self,  x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.densities
    }

    fn get_ck(&self,  x: f64, y: f64, z: f64) -> f64 {
        self.electronic_stopping_correction_factor
    }

    fn get_total_density(&self,  x: f64, y: f64, z: f64) -> f64 {
        self.densities.iter().sum::<f64>()
    }

    fn get_concentrations(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.concentrations
    }

    fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        self.boundary.contains(&point!(x: x, y: y))
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        self.simulation_boundary.contains(&point!(x: x, y: y))
    }
    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        if self.inside(x, y, z) {
            true
        } else {
            if let Closest::SinglePoint(p) = self.boundary.closest_point(&point!(x: x, y: y)) {
                let distance = ((x - p.x()).powf(2.) +  (y - p.y()).powf(2.)).sqrt();
                distance < self.energy_barrier_thickness
            } else if let Closest::Intersection(p) = self.boundary.closest_point(&point!(x: x, y: y)) {
                true
            } else {
                panic!("Geometry error: closest point routine failed to find single closest point to ({}, {}, {}).", x, y, z);
            }
        }
    }
    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        if let Closest::SinglePoint(p) = self.boundary.closest_point(&point!(x: x, y: y)) {
            (p.x(), p.y(), z)
        } else if let Closest::Intersection(p) = self.boundary.closest_point(&point!(x: x, y: y)) {
            return (p.x(), p.y(), z)
        } else {
            panic!("Geometry error: closest point routine failed to find single closest point to ({}, {}, {}).", x, y, z);
        }
    }

    fn nearest_normal_vector(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        panic!("Not implemented.")
    }
}
*/

/// Object that contains raw mesh input data.
#[derive(Deserialize, Clone)]
pub struct Mesh2DInput {
    pub length_unit: String,
    pub triangles: Vec<(usize, usize, usize)>,
    pub points: Vec<(f64, f64)>,
    pub boundary: Vec<usize>,
    pub simulation_boundary_points: Vec<(f64, f64)>,
    pub densities: Vec<Vec<f64>>,
    pub electronic_stopping_correction_factors: Vec<f64>,
    pub energy_barrier_thickness: f64,
}

#[derive(Clone)]
pub struct ParryMesh2D {
    trimesh: TriMesh,
    densities: Vec<Vec<f64>>,
    electronic_stopping_correction_factors: Vec<f64>,
    concentrations: Vec<Vec<f64>>,
    pub boundary: Polyline,
    pub simulation_boundary: Polyline,
    pub energy_barrier_thickness: f64
}

impl ParryMesh2D {
    fn get_id(&self, x: f64, y: f64, z: f64) -> usize {
        let p = Point2d::new(x, y);
        let (point_projection, feature_id) = self.trimesh.project_local_point_and_get_feature(&p);
        match feature_id {
            Vertex(vertex_id) => {
                for (triangle_id, triangle) in self.trimesh.triangles().enumerate() {
                    if triangle.contains_local_point(&self.trimesh.vertices()[vertex_id as usize]) {
                        return triangle_id as usize
                    }
                }
                panic!("Geometry error: point ({}, {}) located on a vertex that is not associated with any triangle. Check mesh input.", x, y);
            },
            Face(triangle_id) => {
                triangle_id as usize
            }
            Unknown => {
                panic!("Geometry error: unknown feature detected near ({}, {}). Check mesh input.", x, y);
            }
        }
    }
}

impl Geometry for ParryMesh2D {
    type InputFileFormat = Input2D;

    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        if self.inside(x, y, z) {
            true
        } else {
            let p = Point2d::new(x, y);
            if self.trimesh.distance_to_local_point(&p, true) < self.energy_barrier_thickness {
                true
            } else {
                false
            }
        }
    }

    fn inside_simulation_boundary (&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point2d::new(x, y);

        let (point_projection, (_, _)) = self.simulation_boundary.project_local_point_assuming_solid_interior_ccw(p);

        point_projection.is_inside
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let p = Point2d::new(x, y);
        let (point_projection, (_, _)) = self.boundary.project_local_point_assuming_solid_interior_ccw(p);
        let (x_, y_) = (point_projection.point.x, point_projection.point.y);
        (x_ as f64, y_ as f64, z)
    }

    fn get_densities(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.densities[self.get_id(x, y, z)]
    }

    fn get_ck(&self, x: f64, y: f64, z: f64) -> f64 {
        self.electronic_stopping_correction_factors[self.get_id(x, y, z)]
    }

    /// Determine the total number density of the triangle that contains or is nearest to (x, y).
    fn get_total_density(&self, x: f64, y: f64, z: f64) -> f64 {
        self.get_densities(x, y, z).iter().sum::<f64>()
    }

    /// Find the concentrations of the triangle that contains or is nearest to (x, y).
    fn get_concentrations(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.concentrations[self.get_id(x, y, z)]
    }

    /// Determines whether the point (x, y) is inside the mesh.
    fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point2d::new(x, y);
        if self.boundary.aabb(&Isometry::identity()).contains_local_point(&p) {
            let (point_projection, (_, _)) = self.boundary.project_local_point_assuming_solid_interior_ccw(p);
            point_projection.is_inside
        } else {
            false
        }
    }

    fn nearest_normal_vector(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let p = Point2d::new(x, y);

        let (point_projection, (_, _)) = self.boundary.project_local_point_assuming_solid_interior_ccw(p);
        let (x_intersect, y_intersect, z_intersect) = self.closest_point(x, y, z);

        let dx = x_intersect - x;
        let dy = y_intersect - y;
        let dz = z_intersect - z;
        let mag = (dx*dx + dy*dy + dz*dz).sqrt();

        if point_projection.is_inside {
            (dx/mag, dy/mag, dz/mag)
        } else {
            (-dx/mag, -dy/mag, -dz/mag)
        }
    }

    fn new(geometry_input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> ParryMesh2D {
        let length_unit: f64 = match geometry_input.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "MM" => MM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => geometry_input.length_unit.parse()
                .expect(format!(
                        "Input errror: could nor parse length unit {}. Use a valid float or one of ANGSTROM, NM, MICRON, CM, MM, M",
                        &geometry_input.length_unit.as_str()
                    ).as_str()),
        };

        let mut boundary_points_converted: Vec<Point2d<f64>> = geometry_input.points.iter().map(|(x, y)| Point2d::new(x*length_unit, y*length_unit)).collect();
        let number_boundary_points = boundary_points_converted.len() as u32;

        let mut simulation_boundary_points_converted: Vec<Point2d<f64>> = geometry_input.simulation_boundary_points.iter().map(|(x, y)| Point2d::new(x*length_unit, y*length_unit)).collect();
        let number_simulation_boundary_points = simulation_boundary_points_converted.len() as u32;

        let test_simulation_boundary_ccw = (0..number_simulation_boundary_points as usize)
        .map(|i| (simulation_boundary_points_converted[(i + 1) % number_simulation_boundary_points as usize].x - simulation_boundary_points_converted[i].x)*(simulation_boundary_points_converted[i].y + simulation_boundary_points_converted[(i + 1) % number_simulation_boundary_points as usize].y))
        .sum::<f64>();

        dbg!(test_simulation_boundary_ccw > 0.0);

        dbg!(&simulation_boundary_points_converted);

        if test_simulation_boundary_ccw > 0.0 {
            simulation_boundary_points_converted = simulation_boundary_points_converted.clone().into_iter().rev().collect::<Vec<Point2d<f64>>>();
        }

        dbg!(&simulation_boundary_points_converted);

        for p in simulation_boundary_points_converted.iter().combinations(2) {
            assert!(p[0] != p[1], "Input error: duplicate vertices in simulation boundary geometry input. Boundary must be defined ccw with no duplicate points (terminal segment will span from final to initial point).")
        }

        let test_ccw = (0..number_boundary_points as usize)
                .map(|i| (boundary_points_converted[(i + 1) % number_boundary_points as usize].x - boundary_points_converted[i].x)*(boundary_points_converted[i].y + boundary_points_converted[(i + 1) % number_boundary_points as usize].y))
                .sum::<f64>();

        if test_ccw > 0.0 {
            boundary_points_converted = boundary_points_converted.clone().into_iter().rev().collect::<Vec<Point2d<f64>>>();
        }

        for p in boundary_points_converted.iter().combinations(2) {
            assert!(p[0] != p[1], "Input error: duplicate vertices in boundary geometry input. Boundary must be defined ccw with no duplicate points (terminal segment will span from final to initial point).")
        }

        let triangles = geometry_input.triangles.iter().map(|(i, j, k)| [*i as u32, *j as u32, *k as u32]).collect::<Vec<[u32; 3]>>();
        let points_converted = geometry_input.points.iter().map(|(x, y)| Point2d::new(x*length_unit, y*length_unit)).collect::<Vec<Point2d<f64>>>();

        let trimesh = TriMesh::new(points_converted, triangles);

        let electronic_stopping_correction_factors = geometry_input.electronic_stopping_correction_factors.clone();

        let densities: Vec<Vec<f64>> = geometry_input.densities.iter().map(|density_list| density_list.iter().map(|density| density / (length_unit).powi(3)).collect::<Vec<f64>>()).collect::<Vec<Vec<f64>>>();


        let energy_barrier_thickness = geometry_input.energy_barrier_thickness*length_unit;

        let concentrations: Vec<Vec<f64>> = densities.iter().map(|density_list| density_list.iter().map(|density| density / density_list.iter().sum::<f64>() ).collect::<Vec<f64>>()).collect::<Vec<Vec<f64>>>();

        for concentration in &concentrations {
            assert!((concentration.iter().sum::<f64>() - 1.0).abs() < 1e-6);
        }

        let mut linked_boundary_points = (0..number_boundary_points).zip(1..number_boundary_points).map(|(x, y)| [x, y]).collect::<Vec<[u32; 2]>>();
        linked_boundary_points.push([number_boundary_points - 1, 0]);
        let boundary2 = Polyline::new(boundary_points_converted, Some(linked_boundary_points));

        let number_simulation_boundary_points = simulation_boundary_points_converted.clone().len() as u32;
        let mut linked_simulation_boundary_points = (0..number_simulation_boundary_points).zip(1..number_simulation_boundary_points).map(|(x, y)| [x, y]).collect::<Vec<[u32; 2]>>();
        linked_simulation_boundary_points.push([number_simulation_boundary_points - 1, 0]);
        let simulation_boundary2 = Polyline::new(simulation_boundary_points_converted, Some(linked_simulation_boundary_points));

        ParryMesh2D {
            trimesh.expect(""),
            densities,
            simulation_boundary: simulation_boundary2,
            boundary: boundary2,
            electronic_stopping_correction_factors,
            energy_barrier_thickness,
            concentrations
        }
    }
}

/*
/// Triangular mesh for rustbca.
#[derive(Clone)]
pub struct Mesh2D {
    mesh: Vec<Cell2D>,
    pub boundary: Polygon<f64>,
    pub simulation_boundary: Polygon<f64>,
    pub energy_barrier_thickness: f64
}

impl Mesh2D {
    /// Finds the cell that is nearest to (x, y).
    fn nearest_to(&self, x: f64, y: f64, z: f64) -> &Cell2D {

        let mut min_distance: f64 = std::f64::MAX;
        let mut index: usize = 0;

        for (cell_index, cell) in self.mesh.iter().enumerate() {
            let distance_to = cell.distance_to(x, y, z);
            if distance_to < min_distance {
                min_distance = distance_to;
                index = cell_index;
            }
        }

        return &self.mesh[index];
    }
}

impl Geometry for Mesh2D {

    type InputFileFormat = Input2D;

    /// Constructor for Mesh2D object from geometry_input.
    fn new(geometry_input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> Mesh2D {

        let triangles = geometry_input.triangles.clone();
        let points = geometry_input.points.clone();
        let material_boundary_point_indices = geometry_input.boundary.clone();

        let simulation_boundary_points = geometry_input.simulation_boundary_points.clone();
        let electronic_stopping_correction_factors = geometry_input.electronic_stopping_correction_factors.clone();
        assert!(electronic_stopping_correction_factors.len() > 0, "Input Error: Electronic stopping correction factor list empty.");
        let n = triangles.len();

        let mut cells: Vec<Cell2D> =  Vec::with_capacity(n);

        //Multiply all coordinates by value of geometry unit.
        let length_unit: f64 = match geometry_input.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "MM" => MM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => geometry_input.length_unit.parse()
                .expect(format!(
                        "Input errror: could nor parse length unit {}. Use a valid float or one of ANGSTROM, NM, MICRON, CM, MM, M", &geometry_input.length_unit.as_str()
                    ).as_str()),
        };

        let densities: Vec<Vec<f64>> = geometry_input.densities
            .iter()
            .map( |row| row.iter().map(|element| element/(length_unit).powi(3)).collect() ).collect();

        //Assert all triangle density lists are equal length
        assert!(
            densities.iter()
            .all(|density_list| density_list.len() == densities[0].len()),
            "Input error: all triangle density lists must be the same size."
        );

        assert_eq!(triangles.len(), densities.len(), "Input error: coordinates and data of unequal length.");


        for ((triangle, densities), ck) in triangles.iter().zip(densities).zip(electronic_stopping_correction_factors) {

            let x1 = points[triangle.0].0;
            let x2 = points[triangle.1].0;
            let x3 = points[triangle.2].0;
            let y1 = points[triangle.0].1;
            let y2 = points[triangle.1].1;
            let y3 = points[triangle.2].1;

            let coordinate_set_converted = (
                x1*length_unit,
                x2*length_unit,
                x3*length_unit,
                y1*length_unit,
                y2*length_unit,
                y3*length_unit,
            );
            let total_density: f64 = densities.iter().sum();
            let concentrations: Vec<f64> = densities.iter().map(|&density| density/total_density).collect::<Vec<f64>>();

            cells.push(Cell2D::new(coordinate_set_converted, densities, concentrations, ck));
        }

        let mut boundary_points_converted = Vec::with_capacity(material_boundary_point_indices.len());
        for index in material_boundary_point_indices.iter() {
            boundary_points_converted.push((points[*index].0*length_unit, points[*index].1*length_unit));
        }
        let boundary: Polygon<f64> = Polygon::new(LineString::from(boundary_points_converted), vec![]);

        let mut simulation_boundary_points_converted = Vec::with_capacity(simulation_boundary_points.len());
        for (x, y) in simulation_boundary_points.iter() {
            simulation_boundary_points_converted.push((x*length_unit, y*length_unit));
        }
        let simulation_boundary: Polygon<f64> = Polygon::new(LineString::from(simulation_boundary_points_converted), vec![]);

        let energy_barrier_thickness = geometry_input.energy_barrier_thickness*length_unit;

        Mesh2D {
            mesh: cells,
            boundary: boundary,
            simulation_boundary: simulation_boundary,
            energy_barrier_thickness: energy_barrier_thickness,
        }
    }

    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        if self.inside(x, y, z) {
            true
        } else {
            let nearest_cell = self.nearest_to(x, y, z);
            let distance = nearest_cell.distance_to(x, y, z);
            distance < self.energy_barrier_thickness
        }
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        self.simulation_boundary.contains(&point!(x: x, y: y))
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        if let Closest::SinglePoint(p) = self.boundary.closest_point(&point!(x: x, y: y)) {
            if p.x() == x && p.y() == y {
                dbg!(&self.boundary);
                println!("single point");
                panic!("({} {} {}) ({} {} {})", p.x(), p.y(), z, x, y, z);
            }
            (p.x(), p.y(), z)
        } else if let Closest::Intersection(p) = self.boundary.closest_point(&point!(x: x, y: y)) {
            if p.x() == x && p.y() == y {
                dbg!(&self.boundary);
                println!("intersection");
                panic!("({} {} {}) ({} {} {})", p.x(), p.y(), z, x, y, z);
            }
            (p.x(), p.y(), z)
        } else {
            panic!("Geometry error: closest point routine failed to find single closest point to ({}, {}, {}).", x, y, z);
        }
    }

    /// Find the number densities of the triangle that contains or is nearest to (x, y).
    fn get_densities(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        if self.inside(x, y, z) {
            for cell in &self.mesh {
                if cell.contains(x, y, z) {
                    return &cell.densities;
                }
            }
            panic!("Geometry error: point ({}, {}) not found in any cell of the mesh.", x, y);
        } else {
            self.nearest_to(x, y, z).get_densities()
        }
    }

    /// Find the number densities of the triangle that contains or is nearest to (x, y).
    fn get_ck(&self, x: f64, y: f64, z: f64) -> f64 {
        if self.inside(x, y, z) {
            for cell in &self.mesh {
                if cell.contains(x, y, z) {
                    return cell.electronic_stopping_correction_factor;
                }
            }
            panic!("Geometry error: point ({}, {}) not found in any cell of the mesh.", x, y);
        } else {
            return self.nearest_to(x, y, z).electronic_stopping_correction_factor;
        }
    }

    /// Determine the total number density of the triangle that contains or is nearest to (x, y).
    fn get_total_density(&self, x: f64, y: f64, z: f64) -> f64 {
        if self.inside(x, y, z) {
            for cell in &self.mesh {
                if cell.contains(x, y, z) {
                    return cell.densities.iter().sum::<f64>();
                }
            }
            panic!("Geometry error: point ({}, {}) not found in any cell of the mesh.", x, y);
        } else {
            return self.nearest_to(x, y, z).densities.iter().sum::<f64>();
        }
    }

    /// Find the concentrations of the triangle that contains or is nearest to (x, y).
    fn get_concentrations(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        if self.inside(x, y, z) {
            for cell in &self.mesh {
                if cell.contains(x, y, z) {
                    return &cell.concentrations;
                }
            }
            panic!("Geometry error: method inside() is returning true for points outside all cells. Check boundary points.")
        } else {
            return &self.nearest_to(x, y, z).concentrations;
        }
    }

    /// Determines whether the point (x, y) is inside the mesh.
    fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        self.boundary.contains(&Point::new(x, y))
    }

    fn nearest_normal_vector(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        panic!("Not implemented.");
    }
}

/// A mesh cell that contains a triangle and the local number densities and concentrations.
#[derive(Clone)]
pub struct Cell2D {
    triangle: Triangle2D,
    pub densities: Vec<f64>,
    pub concentrations: Vec<f64>,
    pub electronic_stopping_correction_factor: f64
}
impl Cell2D {
    /// Constructor for Cell2D from a set of coordinates and a list of densities and concentrations.
    pub fn new(coordinate_set: (f64, f64, f64, f64, f64, f64), densities: Vec<f64>, concentrations: Vec<f64>, electronic_stopping_correction_factor: f64) -> Cell2D {
        Cell2D {
            triangle: Triangle2D::new(coordinate_set),
            densities,
            concentrations,
            electronic_stopping_correction_factor
        }
    }
}

impl GeometryElement for Cell2D {
    /// Determines whether this cell contains the point (x, y).
    fn contains(&self, x: f64, y: f64, z: f64) -> bool {
        self.triangle.contains(x, y)
    }

    /// Computes the shortest distance between this cell and the point (x, y).
    fn distance_to(&self, x: f64, y: f64, z: f64) -> f64 {
        self.triangle.distance_to(x, y)
    }

    fn get_densities(&self) -> &Vec<f64> {
        &self.densities
    }

    fn get_concentrations(&self) -> &Vec<f64> {
        &self.concentrations
    }

    fn get_electronic_stopping_correction_factor(&self) -> f64 {
        self.electronic_stopping_correction_factor
    }
}
*/

#[derive(Clone)]
pub struct Layer1D {
    top: f64,
    bottom: f64,
    densities: Vec<f64>,
    concentrations: Vec<f64>,
    electronic_stopping_correction_factor: f64
}

impl Layer1D {
    pub fn new(top: f64, bottom: f64, densities: Vec<f64>, concentrations: Vec<f64>, electronic_stopping_correction_factor: f64) -> Layer1D {

        assert!(top < bottom, "Layer cannot have negative thickness.");

        Layer1D{
            top,
            bottom,
            densities,
            concentrations,
            electronic_stopping_correction_factor
        }
    }
}

impl GeometryElement for Layer1D {
    fn contains(&self, x: f64, y: f64, z: f64) -> bool {
        (x > self.top) & (x <= self.bottom)
    }
    fn distance_to(&self, x: f64, y: f64, z: f64) -> f64 {
        let distance_to_top = (x - self.top).abs();
        let distance_to_bottom = (x - self.bottom).abs();
        if distance_to_top < distance_to_bottom {
            distance_to_top
        } else {
            distance_to_bottom
        }
    }
    fn get_densities(&self) -> &Vec<f64> {
        &self.densities
    }
    fn get_concentrations(&self) -> &Vec<f64> {
        &self.concentrations
    }
    fn get_electronic_stopping_correction_factor(&self) -> f64 {
        self.electronic_stopping_correction_factor
    }
}

/*
/// A triangle in 2D, with points (x1, y1), (x2, y2), (x3, y3), and the three line segments bewtween them.
#[derive(Clone)]
pub struct Triangle2D {
    x1: f64,
    x2: f64,
    x3: f64,
    y1: f64,
    y2: f64,
    y3: f64,
    segments: Vec<(f64, f64, f64, f64)>
}
impl Triangle2D {
    /// Constructor for a triangle from a list of coordinates (x1, x2, x3, y1, y2, y3).
    pub fn new(coordinate_set: (f64, f64, f64, f64, f64, f64)) -> Triangle2D {

        let x1 = coordinate_set.0;
        let x2 = coordinate_set.1;
        let x3 = coordinate_set.2;

        let y1 = coordinate_set.3;
        let y2 = coordinate_set.4;
        let y3 = coordinate_set.5;

        Triangle2D {
            x1,
            x2,
            x3,
            y1,
            y2,
            y3,
            segments: vec![(x1, y1, x2, y2), (x2, y2, x3, y3), (x3, y3, x1, y1)]
        }
    }

    /// Determines whether the point (x, y) is inside this triangle.
    pub fn contains(&self, x: f64, y: f64) -> bool {
        let x1 = self.x1;
        let x2 = self.x2;
        let x3 = self.x3;
        let y1 = self.y1;
        let y2 = self.y2;
        let y3 = self.y3;

        let a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
        let b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
        let c = 1. - a - b;

         (0. <= a) & (a <= 1.) & (0. <= b) & (b <= 1.) & (0. <= c) & (c <= 1.)
    }

    /// Returns a point (x, y) that is the centroid of the triangle.
    pub fn centroid(&self) -> (f64, f64) {
        ((self.x1 + self.x2 + self.x3)/3., (self.y1 + self.y2 + self.y3)/3.)
    }

    /// Calculates the shortest distance from this triangle to the point (x, y).
    pub fn distance_to(&self, x: f64, y: f64) -> f64 {
        let mut distance_to = std::f64::MAX;

        for segment in &self.segments {
            let length_2 = (segment.2 - segment.0)*(segment.2 - segment.0) + (segment.3 - segment.1)*(segment.3 - segment.1);

            assert!(length_2 != 0., "Geometry error: mesh contains triangle with zero-length side. (x1, y1), (x2, y2) = ({}, {}) ({}, {})", segment.0, segment.2, segment.1, segment.3);

            let u = ((x - segment.0)*(segment.2 - segment.0) + (y - segment.1)*(segment.3 - segment.1))/length_2;

            let distance = if u < 0. {
                ((x - segment.0)*(x - segment.0) + (y - segment.1)*(y - segment.1)).sqrt()
            } else if u > 1. {
                ((x - segment.2)*(x - segment.2) + (y - segment.3)*(y - segment.3)).sqrt()
            } else {
                let x_intersect = segment.0 + u*(segment.2 - segment.0);
                let y_intersect = segment.1 + u*(segment.3 - segment.1);

                ((x - x_intersect)*(x - x_intersect) + (y - y_intersect)*(y - y_intersect)).sqrt()
            };

            if distance < distance_to {
                distance_to = distance;
            }
        }
        distance_to
    }

    /// Calculates the distance between the center of this triangle and the point (x, y).
    pub fn distance_to_center(&self, x: f64, y: f64) -> f64 {
        let centroid = self.centroid();
        ((x - centroid.0)*(x - centroid.0) + (y - centroid.1)*(y - centroid.1)).sqrt()
    }
}
*/