use super::*;
use std::fs::File;

pub struct OutputUnits {
    pub length_unit: f64,
    pub energy_unit: f64,
    pub mass_unit: f64
}

/// Converts from 6D particle coordinates to energy, angle coordinates w.r.t. negative x-axis
pub fn energy_angle_from_particle(particle: &particle::Particle, units: &OutputUnits) -> (f64, f64) {
    let energy = particle.E/units.energy_unit;
    let ux = particle.dir.x;
    let uy = particle.dir.y;
    let uz = particle.dir.z;

    let vyz = ((uy).powi(2) + (uz).powi(2)).sqrt();
    let angle = vyz.atan2(-ux) * 180.0 / PI;

    (energy, angle)
}

#[cfg(feature = "distributions")]
extern crate ndarray;

#[cfg(feature = "distributions")]
use ndarray::prelude::*;

/// Distribution tracker for tracking EADs and implantation distributions
#[derive(Serialize)]
#[cfg(feature = "distributions")]
pub struct Distributions {
    pub energies: Array1<f64>,
    pub angles: Array1<f64>,
    pub x_range: Array1<f64>,
    pub y_range: Array1<f64>,
    pub z_range: Array1<f64>,
    pub reflected_ead: Array2<usize>,
    pub sputtered_ead: Array2<usize>,
    pub implanted_x: Array1<usize>,
    pub implanted_y: Array1<usize>,
    pub implanted_z: Array1<usize>,
    pub electronic_energy_loss_x: Array1<f64>,
    pub nuclear_energy_loss_x: Array1<f64>,
}

#[cfg(feature = "distributions")]
impl Distributions {
    pub fn new(options: &Options) -> Distributions {
        Distributions {
            energies: Array::linspace(options.energy_min, options.energy_max, options.energy_num),
            angles: Array::linspace(options.angle_min, options.angle_max, options.angle_num),
            x_range: Array::linspace(options.x_min, options.x_max, options.x_num),
            y_range: Array::linspace(options.y_min, options.y_max, options.y_num),
            z_range: Array::linspace(options.z_min, options.z_max, options.z_num),
            reflected_ead: Array::zeros((options.energy_num, options.angle_num)),
            sputtered_ead: Array::zeros((options.energy_num, options.angle_num)),
            implanted_x: Array::zeros(options.x_num),
            implanted_y: Array::zeros(options.y_num),
            implanted_z: Array::zeros(options.z_num),
            electronic_energy_loss_x: Array::zeros(options.x_num),
            nuclear_energy_loss_x: Array::zeros(options.x_num),
        }
    }

    /// Write distributions to toml
    pub fn print(&self, options: &Options) {

        let distribution_output_file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(format!("{}{}", options.name, "distributions.toml"))
            .context("Could not open distributions output file.")
            .unwrap();
        let mut distribution_file_stream = BufWriter::with_capacity(8000, distribution_output_file);
        let toml = toml::to_string(&self).unwrap();
        writeln!(distribution_file_stream, "{}", toml).unwrap();

    }

    /// Updates distributions with a single particle
    pub fn update(&mut self, particle: &particle::Particle, units: &OutputUnits, options: &Options, num_incident: usize) {
        let (energy, angle) = energy_angle_from_particle(particle, units);

        let delta_energy = self.energies[1] - self.energies[0];
        let energy_index_left: i32 = ((energy - self.energies[0])/delta_energy).floor() as i32;

        let delta_angle = self.angles[1] - self.angles[0];
        let angle_index_left: i32 = ((angle - self.angles[0])/delta_angle).floor() as i32;

        let inside_energy = (energy_index_left >= 0) & (energy_index_left < self.energies.len() as i32);
        let inside_angle = (angle_index_left >= 0) & (angle_index_left < self.angles.len() as i32);

        if particle.incident & particle.left {
            if inside_energy & inside_angle {
                self.reflected_ead[[energy_index_left as usize, angle_index_left as usize]] += 1;
            }
        }

        if !particle.incident & particle.left {
            if inside_energy & inside_angle {
                self.sputtered_ead[[energy_index_left as usize, angle_index_left as usize]] += 1;
            }
        }

        if particle.incident & particle.stopped {
            let x = particle.pos.x/units.length_unit;
            let y = particle.pos.y/units.length_unit;
            let z = particle.pos.z/units.length_unit;

            let delta_x = self.x_range[1] - self.x_range[0];
            let delta_y = self.y_range[1] - self.y_range[0];
            let delta_z = self.z_range[1] - self.z_range[0];

            let x_index_left: i32 = ((x - self.x_range[0])/delta_x) as i32;
            let y_index_left: i32 = ((y - self.y_range[0])/delta_y) as i32;
            let z_index_left: i32 = ((z - self.z_range[0])/delta_z) as i32;

            let inside_x = (x_index_left >= 0) & (x_index_left < self.x_range.len() as i32);
            let inside_y = (y_index_left >= 0) & (y_index_left < self.y_range.len() as i32);
            let inside_z = (z_index_left >= 0) & (z_index_left < self.z_range.len() as i32);

            if inside_x {
                self.implanted_x[x_index_left as usize] += 1;
            }
            if inside_y {
                self.implanted_y[y_index_left as usize] += 1;
            }
            if inside_z {
                self.implanted_z[z_index_left as usize] += 1;
            }
        }

        if options.track_energy_losses & particle.incident {

            for energy_loss in particle.energies.iter() {
                let x = energy_loss.x/units.length_unit;
                let delta_x = self.x_range[1] - self.x_range[0];

                let electronic_energy_loss = energy_loss.Ee/units.energy_unit/(num_incident as f64)/delta_x;
                let nuclear_energy_loss = energy_loss.En/units.energy_unit/(num_incident as f64)/delta_x;

                let x_index_left: i32 = ((x - self.x_range[0])/delta_x) as i32;
                let inside_x = (x_index_left >= 0) & (x_index_left < self.x_range.len() as i32);

                if inside_x {
                    self.electronic_energy_loss_x[x_index_left as usize] += electronic_energy_loss;
                    self.nuclear_energy_loss_x[x_index_left as usize] += nuclear_energy_loss;
                }
            }
        }
    }
}

/// File streams for list output
pub struct OutputListStreams {
    reflected_file_stream: BufWriter<File>,
    sputtered_file_stream: BufWriter<File>,
    deposited_file_stream: BufWriter<File>,
    trajectory_file_stream: BufWriter<File>,
    trajectory_data_stream: BufWriter<File>,
    displacements_file_stream: BufWriter<File>,
    energy_loss_file_stream: BufWriter<File>,
}

/// Simulation-wide summary output tracker.
pub struct Summary {
    pub num_incident: u64,
    pub num_sputtered: u64,
    pub num_reflected: u64,
}

/// Summary tracker of sputtering and reflection
impl Summary {
    pub fn new(num_incident: u64) -> Summary {
        Summary {
            num_incident,
            num_sputtered: 0,
            num_reflected: 0,
        }
    }

    pub fn add(&mut self, particle: &particle::Particle) {
        if particle.incident & particle.left {
            self.num_reflected += 1;
        }
        if !particle.incident & particle.left {
            self.num_sputtered += 1;
        }
    }
}

pub struct SummaryPerSpecies {
    pub m: Vec<f64>,
    pub sputtered: Vec<usize>,
    pub reflected: Vec<usize>,
    pub deposited: Vec<usize>,
    pub summary_stream_file: BufWriter<File>,
}

impl SummaryPerSpecies {
    pub fn new(options: &Options) -> SummaryPerSpecies {

        let summary_output_file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(format!("{}{}", options.name, "summary.output"))
            .context("Could not open output file.")
            .unwrap();
        let writer = BufWriter::with_capacity(8000, summary_output_file);

        SummaryPerSpecies {
            m: vec![],
            sputtered: vec![],
            reflected: vec![],
            deposited: vec![],
            summary_stream_file: writer,
        }
    }

    pub fn print(&mut self, options: &Options, output_units: &OutputUnits) {
        //Write to summary file
        writeln!(self.summary_stream_file, "mass, reflected, sputtered, deposited")
            .expect(format!("Output error: could not write to {}summary.output.", options.name).as_str());

        for (mass, reflected, sputtered, deposited) in izip!(&self.m, &self.reflected, &self.sputtered, &self.deposited) {
            writeln!(self.summary_stream_file, "{}, {}, {}, {},", mass/output_units.mass_unit, reflected, sputtered, deposited)
                .expect(format!("Output error: could not write to {}summary.output.", options.name).as_str());
        }
        self.summary_stream_file.flush().unwrap();
    }

    pub fn update(&mut self, particle: &particle::Particle) {

        if self.m.contains(&(particle.m)) {

            let index = self.m.iter().position(|m| *m == particle.m).unwrap();

            match (particle.incident, particle.left) {
                (true, true) => self.reflected[index] += 1,
                (true, false) => self.deposited[index] += 1,
                (false, true) => self.sputtered[index] += 1,
                _ => (),
            }
        } else {
            self.m.push(particle.m);
            match (particle.incident, particle.left) {
                (true, true) => {self.reflected.push(1); self.deposited.push(0); self.sputtered.push(0);},
                (true, false) => {self.reflected.push(0); self.deposited.push(1); self.sputtered.push(0);},
                (false, true) => {self.reflected.push(0); self.deposited.push(0); self.sputtered.push(1);},
                _ => {self.reflected.push(0); self.deposited.push(0); self.sputtered.push(0);},
            }
        }
    }
}

/// Open list output files for streaming write
pub fn open_output_lists(options: &Options) -> OutputListStreams {
    //Open output files for streaming output
    let reflected_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "reflected.output"))
        .context("Could not open output file.")
        .unwrap();
    let reflected_file_stream = BufWriter::with_capacity(options.write_buffer_size, reflected_file);

    let sputtered_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "sputtered.output"))
        .context("Could not open output file.")
        .unwrap();
    let sputtered_file_stream = BufWriter::with_capacity(options.write_buffer_size, sputtered_file);

    let deposited_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "deposited.output"))
        .context("Could not open output file.")
        .unwrap();
    let deposited_file_stream = BufWriter::with_capacity(options.write_buffer_size, deposited_file);

    let trajectory_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "trajectories.output"))
        .context("Could not open output file.")
        .unwrap();
    let trajectory_file_stream = BufWriter::with_capacity(options.write_buffer_size, trajectory_file);

    let trajectory_data = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "trajectory_data.output"))
        .context("Could not open output file.")
        .unwrap();
    let trajectory_data_stream = BufWriter::with_capacity(options.write_buffer_size, trajectory_data);

    let displacements_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "displacements.output"))
        .context("Could not open output file.")
        .unwrap();
    let displacements_file_stream = BufWriter::with_capacity(options.write_buffer_size, displacements_file);

    let energy_loss_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "energy_loss.output"))
        .context("Could not open output file.")
        .unwrap();
    let energy_loss_file_stream = BufWriter::with_capacity(options.write_buffer_size, energy_loss_file);

    OutputListStreams {
        reflected_file_stream,
        sputtered_file_stream,
        deposited_file_stream,
        trajectory_file_stream,
        trajectory_data_stream,
        displacements_file_stream,
        energy_loss_file_stream,
    }
}

/// Write output lists
pub fn output_lists(output_list_streams: &mut OutputListStreams, particle: particle::Particle, options: &Options, output_units: &OutputUnits) {

    let length_unit = output_units.length_unit;
    let energy_unit = output_units.energy_unit;
    let mass_unit = output_units.mass_unit;

    if !particle.incident & options.track_displacements {
        writeln!(
            output_list_streams.displacements_file_stream, "{},{},{},{},{},{},{},{},{}",
            particle.m/mass_unit, particle.Z, particle.energy_origin/energy_unit,
            particle.pos_origin.x/length_unit, particle.pos_origin.y/length_unit, particle.pos_origin.z/length_unit,
            particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit
        ).expect(format!("Output error: could not write to {}displacements.output.", options.name).as_str());
    }

    //Incident particle, left simulation: reflected
    if particle.incident & particle.left {
        writeln!(
            output_list_streams.reflected_file_stream, "{},{},{},{},{},{},{},{},{},{}",
            particle.m/mass_unit, particle.Z, particle.E/energy_unit,
            particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit,
            particle.dir.x, particle.dir.y, particle.dir.z,
            particle.number_collision_events
        ).expect(format!("Output error: could not write to {}reflected.output.", options.name).as_str());
    }

    //Incident particle, stopped in material: deposited
    if particle.incident & particle.stopped {
        writeln!(
            output_list_streams.deposited_file_stream, "{},{},{},{},{},{}",
            particle.m/mass_unit, particle.Z,
            particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit,
            particle.number_collision_events
        ).expect(format!("Output error: could not write to {}deposited.output.", options.name).as_str());
    }

    //Not an incident particle, left material: sputtered
    if !particle.incident & particle.left {
        writeln!(
            output_list_streams.sputtered_file_stream, "{},{},{},{},{},{},{},{},{},{},{},{},{}",
            particle.m/mass_unit, particle.Z, particle.E/energy_unit,
            particle.pos.x/length_unit, particle.pos.y/length_unit, particle.pos.z/length_unit,
            particle.dir.x, particle.dir.y, particle.dir.z,
            particle.number_collision_events,
            particle.pos_origin.x/length_unit, particle.pos_origin.y/length_unit, particle.pos_origin.z/length_unit
        ).expect(format!("Output error: could not write to {}sputtered.output.", options.name).as_str());
    }

    //Trajectory output
    if particle.track_trajectories {
        writeln!(output_list_streams.trajectory_data_stream, "{}", particle.trajectory.len())
            .expect(format!("Output error: could not write to {}trajectory_data.output.", options.name).as_str());

        for pos in particle.trajectory {
            writeln!(
                output_list_streams.trajectory_file_stream, "{},{},{},{},{},{}",
                particle.m/mass_unit, particle.Z, pos.E/energy_unit,
                pos.x/length_unit, pos.y/length_unit, pos.z/length_unit,
            ).expect(format!("Output error: could not write to {}trajectories.output.", options.name).as_str());
        }
    }

    if particle.incident & options.track_energy_losses {
        for energy_loss in particle.energies {
            writeln!(
                output_list_streams.energy_loss_file_stream, "{},{},{},{},{},{},{}",
                particle.m/mass_unit, particle.Z,
                energy_loss.En/energy_unit, energy_loss.Ee/energy_unit,
                energy_loss.x/length_unit, energy_loss.y/length_unit, energy_loss.z/length_unit,
            ).expect(format!("Output error: could not write to {}energy_loss.output.", options.name).as_str());
        }
    }
}

/// Flush output list streams
pub fn output_list_flush(output_list_streams: &mut OutputListStreams) {
    output_list_streams.displacements_file_stream.flush().unwrap();
    output_list_streams.reflected_file_stream.flush().unwrap();
    output_list_streams.deposited_file_stream.flush().unwrap();
    output_list_streams.sputtered_file_stream.flush().unwrap();
    output_list_streams.trajectory_data_stream.flush().unwrap();
    output_list_streams.trajectory_file_stream.flush().unwrap();
}
