use super::*;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::iter::{IndexedParallelIterator, ParallelExtend, IntoParallelIterator, ParallelIterator};

pub fn silent_physics_loop<T: Geometry + Sync>(particle_input_array: Vec<particle::ParticleInput>, material: material::Material<T>, options: Options, output_units: OutputUnits) -> Vec<particle::Particle> {

    let mut finished_particles: Vec<particle::Particle> = Vec::new();

    finished_particles.par_extend(
        particle_input_array.into_par_iter()
        .enumerate()
        .map_init(
            || if options.seed < 0 { 
                ChaCha8Rng::seed_from_u64(rand::random())
            } else {
                ChaCha8Rng::seed_from_u64(u64::try_from(options.seed).expect("Value Error: seed not u64."))
            },
            | rng, (particle_index, particle_input)| {
                rng.set_stream((particle_index) as u64);
                bca::single_ion_bca(particle::Particle::from_input(particle_input, &options), &material, &options, rng)
        }).flatten()
    );
    finished_particles
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FinishedParticlesContainer {
    sputtered: Vec<bool>,
    implanted: Vec<bool>,
    atomic_number: Vec<usize>,
    mass: Vec<f64>,
    energy: Vec<f64>,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    ux: Vec<f64>,
    uy: Vec<f64>,
    uz: Vec<f64>,
}
impl FinishedParticlesContainer {
    pub fn new() -> FinishedParticlesContainer {
        FinishedParticlesContainer {
            sputtered: Vec::new(),
            implanted: Vec::new(),
            atomic_number: Vec::new(),
            mass: Vec::new(),
            energy: Vec::new(),
            x: Vec::new(),
            y: Vec::new(),
            z: Vec::new(),
            ux: Vec::new(),
            uy: Vec::new(),
            uz: Vec::new(),
        }
    }
}

pub fn process_finished_particles_to_arrays(finished_particles: Vec<particle::Particle>, output_units: OutputUnits) -> FinishedParticlesContainer {
    let mut finished_particles_container = FinishedParticlesContainer::new();
    for particle in finished_particles {

        let implanted = particle.incident & !particle.left;
        let sputtered = !particle.incident & particle.left;
        let reflected = particle.incident & particle.left;

        if implanted | sputtered | reflected {
            finished_particles_container.sputtered.push(sputtered);
            finished_particles_container.implanted.push(implanted);

            finished_particles_container.atomic_number.push(particle.Z as usize);
            finished_particles_container.mass.push(particle.m/output_units.mass_unit);
            finished_particles_container.energy.push(particle.E/output_units.energy_unit);
            finished_particles_container.x.push(particle.pos.x/output_units.length_unit);
            finished_particles_container.y.push(particle.pos.y/output_units.length_unit);
            finished_particles_container.z.push(particle.pos.z/output_units.length_unit);
            finished_particles_container.ux.push(particle.dir.x);
            finished_particles_container.uy.push(particle.dir.y);
            finished_particles_container.uz.push(particle.dir.z);
        }
    }
    finished_particles_container
}

pub fn physics_loop<T: Geometry + Sync>(particle_input_array: Vec<particle::ParticleInput>, material: material::Material<T>, options: Options, output_units: OutputUnits) {

        let total_count: u64 = particle_input_array.len() as u64;
        assert!(total_count/options.num_chunks > 0, "Input error: chunk size == 0 - reduce num_chunks or increase particle count.");

        #[cfg(not(feature = "no_list_output"))]
        let mut output_list_streams = output::open_output_lists(&options);

        let mut summary = output::SummaryPerSpecies::new(&options);

        #[cfg(feature = "distributions")]
        let mut distributions = output::Distributions::new(&options);

        //Create and configure progress bar
        let bar: ProgressBar = ProgressBar::new(total_count);
        bar.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}][{bar:40.cyan/blue}][{eta_precise}] {percent}%").expect("")
            .progress_chars("#>-"));

        //Main loop
        let chunk_size = (total_count/options.num_chunks) as usize;
        for (chunk_index, particle_input_chunk) in particle_input_array.chunks(chunk_size).enumerate() {

            let mut finished_particles: Vec<particle::Particle> = Vec::new();

            // BCA loop is implemented as parallelized extension of a per-chunk initially empty
            // finished particle array via map from particle -> finished particles via BCA
            finished_particles.par_extend(
                particle_input_chunk.into_par_iter()
                .enumerate()
                .map_init(
                    || if options.seed < 0 { 
                        ChaCha8Rng::seed_from_u64(rand::random())
                    } else {
                        ChaCha8Rng::seed_from_u64(u64::try_from(options.seed).expect("Value Error: seed not u64."))
                    },
                    | rng, (particle_index, particle_input)| {
                        rng.set_stream((chunk_index * chunk_size + particle_index) as u64);
                        bar.tick();
                        bar.inc(1);
                        bca::single_ion_bca(particle::Particle::from_input(*particle_input, &options), &material, &options, rng)
                }).flatten()
            );

            // Process this chunk of finished particles for output
            for particle in finished_particles {

                summary.update(&particle);

                #[cfg(feature = "distributions")]
                distributions.update(&particle, &output_units, &options, total_count as usize);

                #[cfg(not(feature = "no_list_output"))]
                output::output_lists(&mut output_list_streams, particle, &options, &output_units);

            }
            //Flush all file streams before dropping to ensure all data is written
            #[cfg(not(feature = "no_list_output"))]
            output::output_list_flush(&mut output_list_streams);
        }

        summary.print(&options, &output_units);

        //Write distributions to file
        #[cfg(feature = "distributions")]
        distributions.print(&options);

        bar.finish();
        println!("Finished!");
}
