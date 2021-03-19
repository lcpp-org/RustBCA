use super::*;

#[derive(Deserialize, Clone)]
pub struct InputSphere {
    pub options: Options,
    pub material_parameters: material::MaterialParameters,
    pub particle_parameters: particle::ParticleParameters,
    pub geometry_input: sphere::SphereInput,
}

impl InputFile for InputSphere {

    fn new(string: &str) -> InputSphere {
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
pub struct SphereInput {
    pub length_unit: String,
    pub radius: f64,
    pub densities: Vec<f64>,
    pub electronic_stopping_correction_factor: f64,
}

#[derive(Clone)]
pub struct Sphere {
    pub densities: Vec<f64>,
    pub concentrations: Vec<f64>,
    pub radius: f64,
    pub electronic_stopping_correction_factor: f64,
    pub energy_barrier_thickness: f64,
}

impl GeometryInput for InputSphere {
    type GeometryInput = SphereInput;
}

impl Geometry for Sphere {

    type InputFileFormat = InputSphere;

    fn new(input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> Sphere {

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

        Sphere {
            densities,
            concentrations,
            radius,
            electronic_stopping_correction_factor,
            energy_barrier_thickness,
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
        let r = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        r < self.radius
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        let r = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        r < (10.*self.energy_barrier_thickness + self.radius)
    }

    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        let r = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        r < (self.energy_barrier_thickness + self.radius)
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let r = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        let R = self.radius;
        let ux = x/r;
        let uy = y/r;
        let uz = z/r;
        (ux*R, uy*R, uz*R)
    }
}
