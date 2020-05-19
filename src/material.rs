use super::*;

#[derive(Deserialize)]
pub struct MaterialParameters {
    pub energy_unit: String,
    pub mass_unit: String,
    pub Eb: Vec<f64>,
    pub Es: Vec<f64>,
    pub Ec: Vec<f64>,
    pub n: Vec<f64>,
    pub Z: Vec<f64>,
    pub m: Vec<f64>,
    pub electronic_stopping_correction_factor: f64
}

#[derive(Deserialize)]
pub struct Geometry {
    pub length_unit: String,
    pub surface: Vec<(f64, f64)>,
    pub energy_surface: Vec<(f64, f64)>,
    pub simulation_surface: Vec<(f64, f64)>
}

pub struct Material {
    pub n: Vec<f64>,
    pub m: Vec<f64>,
    pub Z: Vec<f64>,
    pub Eb: Vec<f64>,
    pub Es: Vec<f64>,
    pub Ec: Vec<f64>,
    pub surface: Polygon<f64>,
    pub energy_surface: Polygon<f64>,
    pub simulation_surface: Polygon<f64>,
    pub electronic_stopping_correction_factor: f64
}
impl Material {
    pub fn new(material_parameters: MaterialParameters, geometry: Geometry) -> Material {

        //Multiply all coordinates by value of geometry unit.
        let length_unit: f64 = match geometry.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => panic!("Incorrect unit {} in input file. Choose one of: MICRON, CM, ANGSTROM, NM, M",
                geometry.length_unit.as_str())
        };
        let energy_unit: f64 = match material_parameters.energy_unit.as_str() {
            "EV" => EV,
            "J"  => 1.,
            "KEV" => EV*1E3,
            "MEV" => EV*1E6,
            _ => panic!("Incorrect unit {} in input file. Choose one of: EV, J, KEV, MEV",
                material_parameters.energy_unit.as_str())
        };
        let mass_unit: f64 = match material_parameters.mass_unit.as_str() {
            "AMU" => AMU,
            "KG" => 1.,
            _ => panic!("Incorrect unit {} in input file. Choose one of: AMU, KG",
                material_parameters.mass_unit.as_str())
        };

        let mut unit_coords = geometry.surface.clone();
        for pair in &mut unit_coords {
            pair.0 *= length_unit;
            pair.1 *= length_unit;
        }
        let mut unit_e_coords = geometry.energy_surface.clone();
        for pair in &mut unit_e_coords {
            pair.0 *= length_unit;
            pair.1 *= length_unit;
        }
        let mut unit_b_coords = geometry.simulation_surface.clone();
        for pair in &mut unit_b_coords {
            pair.0 *= length_unit;
            pair.1 *= length_unit;
        }

        Material {
            n: material_parameters.n,
            m: material_parameters.m.iter().map(|&i| i*mass_unit).collect(),
            Z: material_parameters.Z,
            Eb: material_parameters.Eb.iter().map(|&i| i*energy_unit).collect(),
            Es: material_parameters.Es.iter().map(|&i| i*energy_unit).collect(),
            Ec: material_parameters.Ec.iter().map(|&i| i*energy_unit).collect(),
            electronic_stopping_correction_factor: material_parameters.electronic_stopping_correction_factor,
            surface: Polygon::new(LineString::from(unit_coords), vec![]),
            energy_surface: Polygon::new(LineString::from(unit_e_coords), vec![]),
            simulation_surface: Polygon::new(LineString::from(unit_b_coords), vec![])
        }
    }

    pub fn get_concentrations(&self) -> Vec<f64> {
        let total_number_density: f64 = self.n.iter().sum();
        return self.n.iter().map(|&i| i / total_number_density).collect();
    }

    pub fn get_cumulative_concentrations(&self) -> Vec<f64> {
        let mut sum: f64 = 0.;
        let concentrations = self.get_concentrations();
        let mut cumulative_concentrations = Vec::with_capacity(concentrations.len());
        for concentration in concentrations {
            sum += concentration;
            cumulative_concentrations.push(sum);
        }
        return cumulative_concentrations;
    }

    pub fn inside(&self, x: f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.surface.contains(&p);
    }

    pub fn mfp(&self, x: f64, y: f64) -> f64 {
        return self.total_number_density(x, y).powf(-1./3.);
    }

    pub fn total_number_density(&self, x: f64, y: f64) -> f64 {
        return self.n.iter().sum::<f64>();
    }

    pub fn inside_energy_barrier(&self, x: f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.energy_surface.contains(&p);
    }

    pub fn inside_simulation_boundary(&self, x:f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.simulation_surface.contains(&p);
    }

    pub fn closest_point(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.surface.closest_point(&p)
    }

    pub fn closest_point_on_energy_barrier(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.energy_surface.closest_point(&p)
    }

    pub fn closest_point_on_simulation_surface(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.simulation_surface.closest_point(&p)
    }

    pub fn average_Z(&self, x: f64, y: f64) -> f64 {
        return self.Z.iter().sum::<f64>()/self.Z.len() as f64;
    }

    pub fn average_mass(&self, x: f64, y: f64) -> f64 {
        return self.m.iter().sum::<f64>()/self.m.len() as f64;
    }

    pub fn average_bulk_binding_energy(&self, x: f64, y: f64) -> f64 {
        //returns average bulk binding energy
        return self.Eb.iter().sum::<f64>()/self.Eb.len() as f64;
    }

    pub fn minimum_cutoff_energy(&self) -> f64 {
        let mut min_Ec = self.Ec.iter().sum::<f64>();
        for Ec in self.Ec.iter() {
            if min_Ec > *Ec {
                min_Ec = *Ec;
            }
        }
        return min_Ec;
    }

    pub fn choose(&self, x: f64, y: f64) -> (usize, f64, f64, f64, f64) {
        let random_number: f64 = rand::random::<f64>();
        let cumulative_concentrations = self.get_cumulative_concentrations();

        for (component_index, cumulative_concentration) in cumulative_concentrations.iter().enumerate() {
            if random_number < *cumulative_concentration {
                return (component_index, self.Z[component_index], self.m[component_index], self.Ec[component_index], self.Es[component_index]);
            }
        }
        panic!("Species choice operation failed. Check concentrations.");
    }

    pub fn electronic_stopping_cross_sections(&self, particle_1: &super::particle::Particle, electronic_stopping_mode: i32) -> Vec<f64> {

        let E = particle_1.E;
        let Ma = particle_1.m;
        let Za = particle_1.Z;

        let mut stopping_powers = Vec::with_capacity(self.n.len());

        //Bragg's rule: total stopping power is sum of stopping powers of individual atoms
        //Therefore, compute each stopping power separately, and add them up
        for (n, Zb) in self.n.iter().zip(self.Z.clone()) {
            //let n = self.number_density(particle_1.pos.x, particle_1.pos.y);
            //let Zb = self.Z_eff(particle_1.pos.x, particle_1.pos.y);

            let beta = (1. - (1. + E/Ma/C.powf(2.)).powf(-2.)).sqrt();
            let v = beta*C;

            let I0 = match Zb < 13. {
                true => 12. + 7./Zb,
                false => 9.76 + 58.5*Zb.powf(-1.19),
            };
            let I = Zb*I0*Q;

            //See Biersack and Haggmark - this looks like an empirical shell correction
            let B = match Zb < 3. {
                true => 100.*Za/Zb,
                false => 5.
            };

        //Bethe stopping modified by Biersack and Varelas
            let prefactor = BETHE_BLOCH_PREFACTOR*Zb*Za*Za/beta/beta;
            let eb = 2.*ME*v*v/I;
            let S_high = prefactor*(eb + 1. + B/eb).ln();

            //Lindhard-Scharff electronic stopping
            let S_low = LINDHARD_SCHARFF_PREFACTOR*(Za.powf(7./6.)*Zb)/(Za.powf(2./3.) + Zb.powf(2./3.)).powf(3./2.)*(E/Q/Ma*AMU).sqrt();

            let stopping_power = match electronic_stopping_mode {
                //Biersack-Varelas Interpolation
                INTERPOLATED => 1./(1./S_high + 1./S_low),
                //Oen-Robinson
                LOW_ENERGY_LOCAL => S_low,
                //Lindhard-Scharff
                LOW_ENERGY_NONLOCAL => S_low,
                //Lindhard-Scharff and Oen-Robinson, using Lindhard Equipartition
                LOW_ENERGY_EQUIPARTITION => S_low,
                //Panic at unimplemented electronic stopping mode
                _ => panic!("Unimplemented electronic stopping mode. Use 0: Biersack-Varelas 1: Lindhard-Scharff 2: Oen-Robinson 3: Equipartition")
            };

            stopping_powers.push(stopping_power);
        }
        return stopping_powers;
    }
}

pub fn boundary_condition_2D_planar(particle_1: &mut particle::Particle, material: &material::Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let cosx = particle_1.dir.x;
    let cosy = particle_1.dir.y;
    let cosz = particle_1.dir.z;
    let E = particle_1.E;
    let Es = particle_1.Es;
    let Ec = particle_1.Ec;

    if !material.inside_energy_barrier(x, y) {
        if let Closest::SinglePoint(p2) = material.closest_point(x, y) {
            let dx = p2.x() - x;
            let dy = p2.y() - y;
            let mag = (dx*dx + dy*dy).sqrt();

            let costheta = dx*cosx/mag + dy*cosy/mag;
            let leaving_energy = E*costheta*costheta;

            if costheta < 0. {
                if leaving_energy > Es {

                    particle_1.left = true;
                    particle_1.E += -Es;

                    //Surface refraction via Snell's law
                    let sintheta0 = (1. - costheta*costheta).sqrt();
                    let sintheta1 = sintheta0*(E/(E - Es)).sqrt();
                    let delta_theta = sintheta1.asin() - sintheta0.asin();

                    super::particle::rotate_particle(particle_1, delta_theta, 0.);
                    particle_1.add_trajectory();

                } else {
                    //Specular reflection at local surface normal
                    particle_1.dir.x = -2.*(costheta)*dx/mag + cosx;
                    particle_1.dir.y = -2.*(costheta)*dy/mag + cosy;

                    particle_1.add_trajectory();
                }
            }
        } else {
            panic!("Surface boundary algorithm encountered an error. Check geometry.");
        }
    }

    if (E < Ec) & !particle_1.left {
        particle_1.stopped = true;
        particle_1.add_trajectory();
    }
}
