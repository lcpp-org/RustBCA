use super::*;
use std::fs::File;

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

pub fn open_output_summary(options: &Options) -> BufWriter<File> {
    let summary_output_file = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(format!("{}{}", options.name, "summary.output"))
        .context("Could not open output file.")
        .unwrap();
    BufWriter::with_capacity(options.write_buffer_size, summary_output_file)
}

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

pub fn output_list_flush(output_list_streams: &mut OutputListStreams) {
    output_list_streams.reflected_file_stream.flush().unwrap();
    output_list_streams.deposited_file_stream.flush().unwrap();
    output_list_streams.sputtered_file_stream.flush().unwrap();
    output_list_streams.trajectory_data_stream.flush().unwrap();
    output_list_streams.trajectory_file_stream.flush().unwrap();
}
