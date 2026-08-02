use super::*;

#[derive(Deserialize, Clone)]
pub struct InputSphereInCuboid {
    pub options: Options,
    pub material_parameters: material::MaterialParameters,
    pub particle_parameters: particle::ParticleParameters,
    pub geometry_input: sphereincuboid::SphereInCuboidInput,
}

impl InputFile for InputSphereInCuboid {

    fn new(string: &str) -> InputSphereInCuboid {
        toml::from_str(string).context("Could not parse TOML file. Be sure you are using the correct input file mode (e.g., ./RustBCA SPHERE sphere.toml or RustBCA.exe 0D mesh_0d.toml).").unwrap()
    }

    fn get_options(&self) -> &Options{
        &self.options
    }
    fn get_material_parameters(&self) -> &material::MaterialParameters{
        &self.material_parameters
    }
    fn get_particle_parameters(&self) -> &particle::ParticleParameters{
        &self.particle_parameters
    }
    fn get_geometry_input(&self) -> &Self::GeometryInput{
        &self.geometry_input
    }
}

#[derive(Deserialize, Clone)]
pub struct SphereInCuboidInput {
    pub length_unit: String,

    pub sphere_radius: f64,
    pub sphere_densities: Vec<f64>,
    pub sphere_electronic_stopping_correction_factor: f64,

    pub cuboid_corners: Vec<(f64, f64, f64)>,
    pub cuboid_densities: Vec<f64>,
    pub cuboid_electronic_stopping_correction_factor: f64,
}

#[derive(Clone)]
pub struct SphereInCuboid {
    pub sphere_radius: f64,
    pub sphere_densities: Vec<f64>,
    pub sphere_concentrations: Vec<f64>,
    pub sphere_electronic_stopping_correction_factor: f64,

    pub cuboid_corners: Vec<(f64, f64, f64)>,
    pub cuboid_densities: Vec<f64>,
    pub cuboid_concentrations: Vec<f64>,
    pub cuboid_electronic_stopping_correction_factor: f64,

    pub energy_barrier_thickness: f64,
}

impl GeometryInput for InputSphereInCuboid {
    type GeometryInput = SphereInCuboidInput;
}

impl Geometry for SphereInCuboid {

    type InputFileFormat = InputSphereInCuboid;

    fn new(input: &<<Self as Geometry>::InputFileFormat as GeometryInput>::GeometryInput) -> SphereInCuboid {

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

        let sphere_radius: f64 = input.sphere_radius*length_unit;
        let sphere_densities: Vec<f64> = input.sphere_densities
            .iter().map(|element| element/(length_unit).powi(3)).collect();
        let sphere_total_density: f64 = sphere_densities.iter().sum();
        let sphere_concentrations: Vec<f64> = sphere_densities
            .iter().map(|&density| density/sphere_total_density).collect::<Vec<f64>>();
        let sphere_electronic_stopping_correction_factor: f64 = input.sphere_electronic_stopping_correction_factor;

        let cuboid_corners: Vec<(f64, f64, f64)> = input.cuboid_corners
            .iter().map(
                |&(x, y, z)| (x * length_unit, y*length_unit, z*length_unit)
            ).collect();
        let cuboid_densities: Vec<f64> = input.cuboid_densities
            .iter().map(|element| element/(length_unit).powi(3)).collect();
        let cuboid_total_density: f64 = cuboid_densities.iter().sum();
        let cuboid_concentrations: Vec<f64> = cuboid_densities
            .iter().map(|&density| density/cuboid_total_density).collect::<Vec<f64>>();
        let cuboid_electronic_stopping_correction_factor: f64 = input.cuboid_electronic_stopping_correction_factor;

        let energy_barrier_thickness: f64 = cuboid_total_density.powf(-1./3.)/SQRTPI*2.;

        SphereInCuboid {
            sphere_radius,
            sphere_densities,
            sphere_concentrations,
            sphere_electronic_stopping_correction_factor,

            cuboid_corners,
            cuboid_densities,
            cuboid_concentrations,
            cuboid_electronic_stopping_correction_factor,

            energy_barrier_thickness,
        }
    }

    fn get_densities(&self,  x: f64, y: f64, z: f64) -> &Vec<f64> {
        let r: f64 = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        if r > self.sphere_radius {
            &self.cuboid_densities
        } else {
            &self.sphere_densities
        }
    }

    fn get_ck(&self,  x: f64, y: f64, z: f64) -> f64 {
        let r: f64 = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        if r > self.sphere_radius {
            self.cuboid_electronic_stopping_correction_factor
        } else {
            self.sphere_electronic_stopping_correction_factor
        }
    }

    fn get_total_density(&self,  x: f64, y: f64, z: f64) -> f64{
        let r: f64 = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        if r > self.sphere_radius {
            self.cuboid_densities.iter().sum()
        } else {
            self.sphere_densities.iter().sum()
        }
    }

    fn get_concentrations(&self, x: f64, y: f64, z: f64) -> &Vec<f64> {
        let r: f64 = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        if r > self.sphere_radius {
            &self.cuboid_concentrations
        } else {
            &self.sphere_concentrations
        }
    }

    fn inside(&self, x: f64, y: f64, z: f64) -> bool {
        let (xlo, ylo, zlo): (f64, f64, f64) = self.cuboid_corners[0];
        let (xhi, yhi, zhi): (f64, f64, f64) = self.cuboid_corners[1];

        xlo <= x && x <= xhi &&
        ylo <= y && y <= yhi &&
        zlo <= z && z <= zhi 
    }

    fn inside_simulation_boundary(&self, x: f64, y: f64, z: f64) -> bool {
        let (xlo, ylo, zlo): (f64, f64, f64) = self.cuboid_corners[0];
        let (xhi, yhi, zhi): (f64, f64, f64) = self.cuboid_corners[1];
        let ebt: f64 = self.energy_barrier_thickness;

        (xlo - 10.0 * ebt) <= x && x <= (xhi + 10.0 * ebt) &&
        (ylo - 10.0 * ebt) <= y && y <= (yhi + 10.0 * ebt) &&
        (zlo - 10.0 * ebt) <= z && z <= (zhi + 10.0 * ebt)
    }

    fn inside_energy_barrier(&self, x: f64, y: f64, z: f64) -> bool {
        let (xlo, ylo, zlo): (f64, f64, f64) = self.cuboid_corners[0];
        let (xhi, yhi, zhi): (f64, f64, f64) = self.cuboid_corners[1];
        let ebt: f64 = self.energy_barrier_thickness;

        (xlo - ebt) <= x && x <= (xhi + ebt) &&
        (ylo - ebt) <= y && y <= (yhi + ebt) &&
        (zlo - ebt) <= z && z <= (zhi + ebt)
    }

    fn closest_point(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let (xlo, ylo, zlo): (f64, f64, f64) = self.cuboid_corners[0];
        let (xhi, yhi, zhi): (f64, f64, f64) = self.cuboid_corners[1];

        if self.inside(x, y, z) {
            let faces: [((f64, f64, f64), f64); 6] = [
                ((xlo, y, z), x - xlo),
                ((xhi, y, z), xhi - x),
                ((x, ylo, z), y - ylo),
                ((x, yhi, z), yhi - y),
                ((x, y, zlo), z - zlo),
                ((x, y, zhi), zhi - z),
            ];

            let closest: ((f64, f64, f64), f64) = *faces
                .iter()
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .unwrap();
            closest.0
        } else {
            (x.max(xlo).min(xhi), y.max(ylo).min(yhi), z.max(zlo).min(zhi))
        }
    }
}
