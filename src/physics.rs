use super::*;

pub fn physics_loop<T: Geometry + Sync>(particle_input_array: Vec<particle::ParticleInput>, material: material::Material<T>, options: Options, output_units: OutputUnits) {

        println!("Processing {} ions...", particle_input_array.len());

        static SEED: u64 = 0;

        let total_count: u64 = particle_input_array.len() as u64;
        assert!(total_count/options.num_chunks > 0, "Input error: chunk size == 0 - reduce num_chunks or increase particle count.");

        #[cfg(not(feature = "no_list_output"))]
        let mut output_list_streams = output::open_output_lists(&options);

        let mut summary = output::SummaryPerSpecies::new(&options);

        #[cfg(feature = "distributions")]
        let mut distributions = output::Distributions::new(&options);

        //Initialize threads with rayon
        println!("Initializing with {} threads...", options.num_threads);
        let pool = rayon::ThreadPoolBuilder::new().num_threads(options.num_threads).build_global().unwrap();

        //Create and configure progress bar
        let bar: ProgressBar = ProgressBar::new(total_count);
        bar.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}][{bar:40.cyan/blue}][{eta_precise}] {percent}%")
            .progress_chars("#>-"));

        //Main loop
        for (chunk_index, particle_input_chunk) in particle_input_array.chunks((total_count/options.num_chunks) as usize).enumerate() {

            let mut finished_particles: Vec<particle::Particle> = Vec::new();

            // BCA loop is implemented as parallelized extension of a per-chunk initially empty
            // finished particle array via map from particle -> finished particles via BCA
            finished_particles.par_extend(
                particle_input_chunk.into_par_iter()
                .enumerate()
                .map_init(
                    || ChaCha8Rng::seed_from_u64(SEED),
                    | rng, (particle_index, particle_input)| {
                        rng.set_stream(particle_index as u64);
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
