use super::*;
use parry3d_f64::shape::{Ball, TriMesh};
use parry3d_f64::query::{PointQuery, Ray, RayCast};
use parry3d_f64::math::{Isometry, Point, Vector};
use parry3d_f64::bounding_volume::AABB;

#[derive(Deserialize, Clone)]
pub struct InputParryBall {
    pub options: Options,
    pub material_parameters: material::MaterialParameters,
    pub particle_parameters: particle::ParticleParameters,
    pub geometry_input: ParryBallInput,
}

impl InputFile for InputParryBall {

    fn new(string: &str) -> InputParryBall {
        toml::from_str(string).context("Could not parse TOML file. Be sure you are using the correct input file mode (e.g., ./RustBCA SPHERE sphere.toml or RustBCA.exe 0D mesh_0d.toml).").unwrap()
    }

    fn get_options(&self) -> &Options{
        &self.options
    }
    fn get_material_parameters(&self) -> &material::MaterialParameters{
        &self.material_parameters
    }
    fn get_particle_parameters(&self) ->& particle::ParticleParameters{
        &self.particle_parameters
    }
    fn get_geometry_input(&self) -> &Self::GeometryInput{
        &self.geometry_input
    }
}

#[derive(Deserialize, Clone)]
pub struct ParryBallInput {
    pub length_unit: String,
    pub radius: f64,
    pub densities: Vec<f64>,
    pub electronic_stopping_correction_factor: f64,
}

#[derive(Clone)]
pub struct ParryBall {
    pub densities: Vec<f64>,
    pub concentrations: Vec<f64>,
    pub radius: f64,
    pub electronic_stopping_correction_factor: f64,
    pub energy_barrier_thickness: f64,
    pub ball: Ball,
}

impl GeometryInput for InputParryBall {
    type GeometryInput = ParryBallInput;
}

impl Geometry for ParryBall {

    type InputFileFormat = InputParryBall;

    fn new(input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> ParryBall {

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
        let total_density: f64 = densities.iter().sum();
        let energy_barrier_thickness = total_density.powf(-1./3.)/SQRTPI*2.;
        let concentrations: Vec<f64> = densities.iter().map(|&density| density/total_density).collect::<Vec<f64>>();
        let radius = input.radius*length_unit;

        ParryBall {
            densities,
            concentrations,
            radius,
            electronic_stopping_correction_factor,
            energy_barrier_thickness,
            ball: Ball::new(radius)
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
        let p = Point::new(x , y , z );
        self.ball.contains_local_point(&p)
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point::new(x , y , z );
        (self.ball.distance_to_local_point(&p, true) as f64) < 10.*self.energy_barrier_thickness
    }

    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point::new(x , y , z );
        (self.ball.distance_to_local_point(&p, true) as f64) < self.energy_barrier_thickness
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let p = Point::new(x , y , z );
        let point_projection = self.ball.project_point(&Isometry::identity(), &p, false);
        let (x_, y_, z_) = (point_projection.point.x, point_projection.point.y, point_projection.point.z);
        (x_ as f64, y_ as f64, z_ as f64)
    }
}


#[derive(Deserialize, Clone)]
pub struct InputParryTriMesh {
    pub options: Options,
    pub material_parameters: material::MaterialParameters,
    pub particle_parameters: particle::ParticleParameters,
    pub geometry_input: ParryTriMeshInput,
}

impl InputFile for InputParryTriMesh {

    fn new(string: &str) -> InputParryTriMesh {
        toml::from_str(string).context("Could not parse TOML file. Be sure you are using the correct input file mode (e.g., ./RustBCA SPHERE sphere.toml or RustBCA.exe 0D mesh_0d.toml).").unwrap()
    }

    fn get_options(&self) -> &Options{
        &self.options
    }
    fn get_material_parameters(&self) -> &material::MaterialParameters{
        &self.material_parameters
    }
    fn get_particle_parameters(&self) ->& particle::ParticleParameters{
        &self.particle_parameters
    }
    fn get_geometry_input(&self) -> &Self::GeometryInput{
        &self.geometry_input
    }
}

#[derive(Deserialize, Clone)]
pub struct ParryTriMeshInput {
    pub length_unit: String,
    pub densities: Vec<f64>,
    pub electronic_stopping_correction_factor: f64,
    pub vertices: Vec<[f64; 3]>,
    pub indices: Vec<[u32; 3]>,
}

#[derive(Clone)]
pub struct ParryTriMesh {
    pub densities: Vec<f64>,
    pub concentrations: Vec<f64>,
    pub electronic_stopping_correction_factor: f64,
    pub energy_barrier_thickness: f64,
    pub trimesh: TriMesh,
    pub boundary: AABB,
}

impl GeometryInput for InputParryTriMesh {
    type GeometryInput = ParryTriMeshInput;
}

impl Geometry for ParryTriMesh {

    type InputFileFormat = InputParryTriMesh;

    fn new(input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> ParryTriMesh {

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
        let total_density: f64 = densities.iter().sum();
        let energy_barrier_thickness = total_density.powf(-1./3.)/SQRTPI*2.;
        let concentrations: Vec<f64> = densities.iter().map(|&density| density/total_density).collect::<Vec<f64>>();
        let points = input.vertices.iter().map(|p| Point::new(p[0]*length_unit , p[1]*length_unit , p[2]*length_unit)).collect();
        let trimesh = TriMesh::new(points, input.indices.clone());
        let boundary = trimesh.aabb(&Isometry::identity());

        ParryTriMesh {
            densities,
            concentrations,
            electronic_stopping_correction_factor,
            energy_barrier_thickness,
            trimesh,
            boundary,
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
        inside_trimesh(&self.trimesh, x, y, z)
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        if self.inside(x, y, z) {
            true
        } else {
            let p = Point::new(x , y , z );
            let distance = self.boundary.distance_to_local_point(&p, true);
            //dbg!(distance);
            distance < 10.*self.energy_barrier_thickness
        }
    }

    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        if self.inside(x, y, z) {
            true
        } else {
            let p = Point::new(x , y , z );
            let distance = self.trimesh.distance_to_local_point(&p, true);
            //dbg!(distance/ANGSTROM);
            //dbg!(self.energy_barrier_thickness/ANGSTROM);
            //dbg!(distance < self.energy_barrier_thickness);
            distance < self.energy_barrier_thickness
        }
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let p = Point::new(x , y , z );
        let point_projection = self.trimesh.project_local_point(&p, false);
        let (x_, y_, z_) = (point_projection.point.x, point_projection.point.y, point_projection.point.z);
        (x_ as f64, y_ as f64, z_ as f64)
    }
}

fn inside_trimesh(trimesh: &TriMesh, x: f64, y: f64, z: f64) -> bool {
    let p = Point::new(x , y , z );

    if !trimesh.aabb(&Isometry::identity())
        .contains_local_point(&p) {
        return false;
    }

    let dir = Vector::new(1.0, 0.0, 0.0);

    let result = trimesh.cast_ray_and_get_normal(
        &Isometry::identity(),
        &Ray::new(p, dir),
        1.0,
        false
    );

    //dbg!(result);

    match result {
        Some(intersection) => {
            //dbg!(intersection.feature);
            trimesh.is_backface(intersection.feature)
        },
        None => false,
    }
}
