use super::*;
use parry3d::shape::{Ball};
use parry3d::query::PointQuery;
use parry3d::math::{Isometry, Point};

#[derive(Deserialize, Clone)]
pub struct InputParryBall {
    pub options: Options,
    pub material_parameters: material::MaterialParameters,
    pub particle_parameters: particle::ParticleParameters,
    pub geometry_input: ParryBallInput,
}

impl GeometryInput for InputParryBall {
    type GeometryInput = ParryBallInput;
}

impl InputFile for InputParryBall {

    fn new(string: &str) -> InputParryBall {
        toml::from_str(string).expect("Could not parse TOML file.")
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

impl GeometryInput for ParryBall {
    type GeometryInput = ParryBallInput;
}

impl Geometry for ParryBall {

    type InputFileFormat = InputParryBall;

    fn new(input: &<Self as GeometryInput>::GeometryInput) -> ParryBall {

        let length_unit: f64 = match input.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
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
            ball: Ball::new(radius as f32)
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
    fn get_densities_nearest_to(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        &self.densities
    }
    fn get_ck_nearest_to(&self, x: f64, y: f64, z: f64) -> f64 {
        self.electronic_stopping_correction_factor
    }
    fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point::new(x as f32, y as f32, z as f32);
        self.ball.contains_local_point(&p)
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point::new(x as f32, y as f32, z as f32);
        (self.ball.distance_to_local_point(&p, true) as f64) < 10.*self.energy_barrier_thickness
    }

    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        let p = Point::new(x as f32, y as f32, z as f32);
        (self.ball.distance_to_local_point(&p, true) as f64) < self.energy_barrier_thickness
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let p = Point::new(x as f32, y as f32, z as f32);
        let point_projection = self.ball.project_point(&Isometry::identity(), &p, false);
        let (x_, y_, z_) = (point_projection.point.x, point_projection.point.y, point_projection.point.z);
        (x_ as f64, y_ as f64, z_ as f64)
    }
}
