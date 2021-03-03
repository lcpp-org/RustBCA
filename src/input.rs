use super::*;

pub trait InputFile: GeometryInput {
    fn new(string: &str) -> Self;
    fn get_options(&self) -> &Options;
    fn get_material_parameters(&self) -> &material::MaterialParameters;
    fn get_particle_parameters(&self) -> &particle::ParticleParameters;
    fn get_geometry_input(&self) -> &<Self as GeometryInput>::GeometryInput;
}

pub trait GeometryInput {
    type GeometryInput;
}

/// Rustbca's internal representation of an input file.
#[derive(Deserialize, Clone)]
pub struct Input2D {
    pub options: Options,
    pub material_parameters: material::MaterialParameters,
    pub particle_parameters: particle::ParticleParameters,
    pub geometry_input: geometry::Mesh2DInput,
}

impl GeometryInput for Input2D {
    type GeometryInput = geometry::Mesh2DInput;
}

impl InputFile for Input2D {

    fn new(string: &str) -> Input2D {
        toml::from_str(string).unwrap()
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
pub struct Input0D {
    pub options: Options,
    pub material_parameters: material::MaterialParameters,
    pub particle_parameters: particle::ParticleParameters,
    pub geometry_input: geometry::Mesh0DInput,
}

impl GeometryInput for Input0D {
    type GeometryInput = geometry::Mesh0DInput;
}

impl InputFile for Input0D {

    fn new(string: &str) -> Input0D {
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
pub struct Input1D {
    pub options: Options,
    pub material_parameters: material::MaterialParameters,
    pub particle_parameters: particle::ParticleParameters,
    pub geometry_input: geometry::Mesh1DInput,
}

impl GeometryInput for Input1D {
    type GeometryInput = geometry::Mesh1DInput;
}

impl InputFile for Input1D {

    fn new(string: &str) -> Input1D {
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

/// Rustbca's internal representation of the simulation-level options.
#[cfg(not(feature = "distributions"))]
#[derive(Deserialize, Clone)]
pub struct Options {
    pub name: String,
    pub track_trajectories: bool,
    pub track_recoils: bool,
    pub track_recoil_trajectories: bool,
    pub write_buffer_size: usize,
    pub weak_collision_order: usize,
    pub suppress_deep_recoils: bool,
    pub high_energy_free_flight_paths: bool,
    pub electronic_stopping_mode: ElectronicStoppingMode,
    pub mean_free_path_model: MeanFreePathModel,
    pub interaction_potential: Vec<Vec<InteractionPotential>>,
    pub scattering_integral: Vec<Vec<ScatteringIntegral>>,
    pub root_finder: Vec<Vec<Rootfinder>>,
    pub num_threads: usize,
    pub num_chunks: u64,
    pub use_hdf5: bool,
    pub track_displacements: bool,
    pub track_energy_losses: bool,
}

#[cfg(feature = "distributions")]
#[derive(Deserialize, Clone)]
pub struct Options {
    pub name: String,
    pub track_trajectories: bool,
    pub track_recoils: bool,
    pub track_recoil_trajectories: bool,
    pub write_buffer_size: usize,
    pub weak_collision_order: usize,
    pub suppress_deep_recoils: bool,
    pub high_energy_free_flight_paths: bool,
    pub electronic_stopping_mode: ElectronicStoppingMode,
    pub mean_free_path_model: MeanFreePathModel,
    pub interaction_potential: Vec<Vec<InteractionPotential>>,
    pub scattering_integral: Vec<Vec<ScatteringIntegral>>,
    pub root_finder: Vec<Vec<Rootfinder>>,
    pub num_threads: usize,
    pub num_chunks: u64,
    pub use_hdf5: bool,
    pub track_displacements: bool,
    pub track_energy_losses: bool,
    pub energy_min: f64,
    pub energy_max: f64,
    pub energy_num: usize,
    pub angle_min: f64,
    pub angle_max: f64,
    pub angle_num: usize,
    pub x_min: f64,
    pub y_min: f64,
    pub z_min: f64,
    pub x_max: f64,
    pub y_max: f64,
    pub z_max: f64,
    pub x_num: usize,
    pub y_num: usize,
    pub z_num: usize,
}

pub fn input<T: Geometry>(input_file: String) -> (Vec<particle::ParticleInput>, material::Material<T>, Options, OutputUnits)
where <T as Geometry>::InputFileFormat: Deserialize<'static> + 'static, <T as GeometryInput>::GeometryInput: Clone {

    //Read input file, convert to string, and open with toml
    let mut input_toml = String::new();
    let mut file = OpenOptions::new()
        .read(true)
        .write(false)
        .create(false)
        .open(&input_file)
        .expect(format!("Input errror: could not open input file {}.", &input_file).as_str());
    file.read_to_string(&mut input_toml).context("Could not convert TOML file to string.").unwrap();

    let input: <T as Geometry>::InputFileFormat = InputFile::new(&input_toml);

    //Unpack toml information into structs
    let options = (*input.get_options()).clone();
    let particle_parameters = (*input.get_particle_parameters()).clone();
    let material_parameters = (*input.get_material_parameters()).clone();
    let material: material::Material<T> = material::Material::<T>::new(&material_parameters, input.get_geometry_input());

    //Ensure nonsensical threads/chunks options crash on input
    assert!(options.num_threads > 0, "Input error: num_threads must be greater than zero.");
    assert!(options.num_chunks > 0, "Input error: num_chunks must be greater than zero.");

    //Check that all material arrays are of equal length.
    assert!(material.m.len() == material.Z.len(), "Input error: material input arrays of unequal length.");
    assert!(material.m.len() == material.Eb.len(), "Input error: material input arrays of unequal length.");
    assert!(material.m.len() == material.Es.len(), "Input error: material input arrays of unequal length.");
    assert!(material.m.len() == material.interaction_index.len(), "Input error: material input arrays of unequal length.");

    //Check that incompatible options are not on simultaneously
    if options.high_energy_free_flight_paths {
        assert!(options.electronic_stopping_mode == ElectronicStoppingMode::INTERPOLATED,
            "Input error: High energy free flight paths used with low energy stoppping power. Change to INTERPOLATED.");
        assert!(options.weak_collision_order == 0,
            "Input error: Cannot use weak collisions with free flight paths. Set weak_collision_order = 0.");
    }

    if options.mean_free_path_model == MeanFreePathModel::GASEOUS {
        assert!(options.weak_collision_order == 0,
            "Input error: Cannot use weak collisions with gaseous mean free path model. Set weak_collision_order = 0.");
    }

    for i in 0..options.interaction_potential.len() {
        assert!(&options.interaction_potential.len() == &options.interaction_potential[i].len(),
            "Input error: interaction matrix not square.");

        assert!(&options.scattering_integral.len() == &options.scattering_integral[i].len(),
            "Input error: scattering intergral matrix not square.");

        assert!(&options.root_finder.len() == &options.root_finder[i].len(),
            "Input error: rootfinder matrix not square.");
    }

    for ((interaction_potentials, scattering_integrals), root_finders) in options.interaction_potential.clone().iter().zip(options.scattering_integral.clone()).zip(options.root_finder.clone()) {
        for ((interaction_potential, scattering_integral), root_finder) in interaction_potentials.iter().zip(scattering_integrals).zip(root_finders) {

            if cfg!(not(any(feature="cpr_rootfinder_openblas", feature="cpr_rootfinder_netlib", feature="cpr_rootfinder_intel_mkl",))) {
                assert!( match root_finder {
                    Rootfinder::POLYNOMIAL{..} => false,
                    Rootfinder::CPR{..} => false,
                    _ => true,
                },
                "Input error: CPR rootfinder not enabled. Build with --features cpr_rootfinder");
            }

            assert!(
                match (interaction_potential, root_finder) {
                    (InteractionPotential::LENNARD_JONES_12_6{..}, Rootfinder::CPR{..}) => true,
                    (InteractionPotential::LENNARD_JONES_12_6{..}, Rootfinder::POLYNOMIAL{..}) => true,
                    (InteractionPotential::LENNARD_JONES_12_6{..}, _) => false,
                    (InteractionPotential::LENNARD_JONES_65_6{..}, Rootfinder::CPR{..}) => true,
                    (InteractionPotential::LENNARD_JONES_65_6{..}, Rootfinder::POLYNOMIAL{..}) => true,
                    (InteractionPotential::LENNARD_JONES_65_6{..}, _) => false,
                    (InteractionPotential::MORSE{..}, Rootfinder::CPR{..}) => true,
                    (InteractionPotential::MORSE{..}, _) => false,
                    (_, Rootfinder::POLYNOMIAL{..}) => false,
                    (_, _) => true,
                },
            "Input error: cannot use {} with {}. Try switching to a different rootfinder.", interaction_potential, root_finder);
        }
    }

    //Check that particle arrays are equal length
    assert_eq!(particle_parameters.Z.len(), particle_parameters.m.len(),
        "Input error: particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.E.len(),
        "Input error: particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.pos.len(),
        "Input error: particle input arrays of unequal length.");
    assert_eq!(particle_parameters.Z.len(), particle_parameters.dir.len(),
        "Input error: particle input arrays of unequal length.");

    //Check that interaction indices are all within interaction matrices
    assert!(material.interaction_index.iter().max().unwrap() < &options.interaction_potential.len(),
        "Input error: interaction matrix too small for material interaction indices.");
    assert!(particle_parameters.interaction_index.iter().max().unwrap() < &options.interaction_potential.len(),
        "Input error: interaction matrix too small for particle interaction indices.");

    //N is the number of distinct particles.
    let N = particle_parameters.Z.len();

    //Determine the length, energy, and mass units for particle input
    let length_unit: f64 = match particle_parameters.length_unit.as_str() {
        "MICRON" => MICRON,
        "CM" => CM,
        "ANGSTROM" => ANGSTROM,
        "NM" => NM,
        "M" => 1.,
        _ => particle_parameters.length_unit.parse()
            .expect(format!(
                    "Input errror: could nor parse length unit {}. Use a valid float or one of ANGSTROM, NM, MICRON, CM, MM, M", &particle_parameters.length_unit.as_str()
                ).as_str()),
    };

    let energy_unit: f64 = match particle_parameters.energy_unit.as_str() {
        "EV" => EV,
        "J"  => 1.,
        "KEV" => EV*1E3,
        "MEV" => EV*1E6,
        _ => particle_parameters.energy_unit.parse()
            .expect(format!(
                    "Input errror: could nor parse energy unit {}. Use a valid float or one of EV, J, KEV, MEV", &particle_parameters.energy_unit.as_str()
                ).as_str()),
    };

    let mass_unit: f64 = match particle_parameters.mass_unit.as_str() {
        "AMU" => AMU,
        "KG" => 1.0,
        _ => particle_parameters.mass_unit.parse()
            .expect(format!(
                    "Input errror: could nor parse mass unit {}. Use a valid float or one of AMU, KG", &particle_parameters.mass_unit.as_str()
                ).as_str()),
    };

    //HDF5
    #[cfg(feature = "hdf5_input")]
    let particle_input_array: Vec<particle::ParticleInput> = {
        if options.use_hdf5 {
            let particle_input_filename = particle_parameters.particle_input_filename.as_str();
            let _e = hdf5::silence_errors();
            let particle_input_file = hdf5::File::open(particle_input_filename)
                .context("Input error: cannot open HDF5 file.")
                .unwrap();
            let particle_input = particle_input_file.dataset("particles")
                .context("Input error: cannot read from HDF5 file.")
                .unwrap();
            particle_input.read_raw::<particle::ParticleInput>().unwrap()

        } else {
            let mut particle_input: Vec<particle::ParticleInput> = Vec::new();

            for particle_index in 0..N {
                let N_ = particle_parameters.N[particle_index];
                let m = particle_parameters.m[particle_index];
                let Z = particle_parameters.Z[particle_index];
                let E = particle_parameters.E[particle_index];
                let Ec = particle_parameters.Ec[particle_index];
                let Es = particle_parameters.Es[particle_index];
                let interaction_index = particle_parameters.interaction_index[particle_index];
                let (x, y, z) = particle_parameters.pos[particle_index];
                let (cosx, cosy, cosz) = particle_parameters.dir[particle_index];
                assert!(cosx < 1.,
                    "Input error: particle x-direction cannot be exactly equal to 1 to avoid numerical gimbal lock.");
                for sub_particle_index in 0..N_ {
                    //Add new particle to particle vector
                    particle_input.push(
                        particle::ParticleInput{
                            m: m*mass_unit,
                            Z: Z,
                            E: E*energy_unit,
                            Ec: Ec*energy_unit,
                            Es: Es*energy_unit,
                            x: x*length_unit,
                            y: y*length_unit,
                            z: z*length_unit,
                            ux: cosx,
                            uy: cosy,
                            uz: cosz,
                            interaction_index: interaction_index
                        }
                    );
                }
            }
            particle_input
        }
    };

    #[cfg(not(feature = "hdf5_input"))]
    let particle_input_array: Vec<particle::ParticleInput> = {
        if options.use_hdf5 {
            panic!("HDF5 particle input not enabled. Enable with: cargo build --features hdf5_input")
        } else {
            let mut particle_input: Vec<particle::ParticleInput> = Vec::new();

            for particle_index in 0..N {
                let N_ = particle_parameters.N[particle_index];
                let m = particle_parameters.m[particle_index];
                let Z = particle_parameters.Z[particle_index];
                let E = particle_parameters.E[particle_index];
                let Ec = particle_parameters.Ec[particle_index];
                let Es = particle_parameters.Es[particle_index];
                let interaction_index = particle_parameters.interaction_index[particle_index];
                let (x, y, z) = particle_parameters.pos[particle_index];
                let (cosx, cosy, cosz) = particle_parameters.dir[particle_index];
                assert!(cosx < 1.,
                    "Input error: particle x-direction cannot be exactly equal to 1 to avoid numerical gimbal lock.");
                for sub_particle_index in 0..N_ {

                    //Add new particle to particle vector
                    particle_input.push(
                        particle::ParticleInput{
                            m: m*mass_unit,
                            Z: Z,
                            E: E*energy_unit,
                            Ec: Ec*energy_unit,
                            Es: Es*energy_unit,
                            x: x*length_unit,
                            y: y*length_unit,
                            z: z*length_unit,
                            ux: cosx,
                            uy: cosy,
                            uz: cosz,
                            interaction_index
                        }
                    );
                }
            }
            particle_input
        }
    };
    (particle_input_array, material, options, OutputUnits {length_unit, energy_unit, mass_unit})
}
