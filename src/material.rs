use super::*;
use geo::algorithm::contains::Contains;
use geo::algorithm::closest_point::ClosestPoint;
use geo::{point, Closest};


/// Holds material input parameters from [material_params].
#[derive(Deserialize)]
pub struct MaterialParameters {
    pub energy_unit: String,
    pub mass_unit: String,
    pub Eb: Vec<f64>,
    pub Es: Vec<f64>,
    pub Ec: Vec<f64>,
    pub Z: Vec<f64>,
    pub m: Vec<f64>,
    pub interaction_index: Vec<usize>,
    pub surface_binding_model: SurfaceBindingModel
}

/// Material in rustbca. Includes the material properties and the mesh that defines the material geometry.
pub struct Material {
    pub m: Vec<f64>,
    pub Z: Vec<f64>,
    pub Eb: Vec<f64>,
    pub Es: Vec<f64>,
    pub Ec: Vec<f64>,
    pub interaction_index: Vec<usize>,
    pub mesh_2d: mesh::Mesh2D,
    pub surface_binding_model: SurfaceBindingModel

}
impl Material {
    /// Constructs a new material object from a material parameters object and a mesh_2d_input object.
    pub fn new(material_parameters: MaterialParameters, mesh_2d_input: mesh::Mesh2DInput) -> Material {

        //Multiply all coordinates by value of geometry unit.
        let length_unit: f64 = match mesh_2d_input.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => panic!("Input error: incorrect unit {} in input file. Choose one of: MICRON, CM, ANGSTROM, NM, M",
                mesh_2d_input.length_unit.as_str())
        };

        let energy_unit: f64 = match material_parameters.energy_unit.as_str() {
            "EV" => EV,
            "J"  => 1.,
            "KEV" => EV*1E3,
            "MEV" => EV*1E6,
            _ => panic!("Input error: incorrect unit {} in input file. Choose one of: EV, J, KEV, MEV",
                material_parameters.energy_unit.as_str())
        };

        let mass_unit: f64 = match material_parameters.mass_unit.as_str() {
            "AMU" => AMU,
            "KG" => 1.,
            _ => panic!("Input error: incorrect unit {} in input file. Choose one of: AMU, KG",
                material_parameters.mass_unit.as_str())
        };

        Material {
            m: material_parameters.m.iter().map(|&i| i*mass_unit).collect(),
            Z: material_parameters.Z,
            Eb: material_parameters.Eb.iter().map(|&i| i*energy_unit).collect(),
            Es: material_parameters.Es.iter().map(|&i| i*energy_unit).collect(),
            Ec: material_parameters.Ec.iter().map(|&i| i*energy_unit).collect(),
            interaction_index: material_parameters.interaction_index,
            mesh_2d: mesh::Mesh2D::new(mesh_2d_input),
            surface_binding_model: material_parameters.surface_binding_model,
        }
    }

    /// Gets concentrations of triangle that contains or is nearest to (x, y)
    pub fn get_concentrations(&self, x: f64, y: f64) -> Vec<f64> {

        let total_number_density: f64 = if self.inside(x, y) {
            self.mesh_2d.get_total_density(x,y)
        } else {
            self.mesh_2d.nearest_cell_to(x, y).densities.iter().sum::<f64>()
        };

        return self.mesh_2d.get_densities(x, y).iter().map(|&i| i / total_number_density).collect();
    }

    /// Gets cumulative concentrations of triangle that contains or is nearest to (x, y).
    /// Used for weighting the random choice of one of the constituent species.
    pub fn get_cumulative_concentrations(&self, x: f64, y: f64) -> Vec<f64> {
        let mut sum: f64 = 0.;
        let concentrations = self.mesh_2d.get_concentrations(x, y);
        let mut cumulative_concentrations = Vec::with_capacity(concentrations.len());
        for concentration in concentrations {
            sum += concentration;
            cumulative_concentrations.push(sum);
        }
        return cumulative_concentrations;
    }

    /// Determines whether (x, y) is inside the material.
    pub fn inside(&self, x: f64, y: f64) -> bool {
        return self.mesh_2d.inside(x, y);
    }

    /// Gets electronic stopping correction factor for LS and OR
    pub fn electronic_stopping_correction_factor(&self, x: f64, y: f64) -> f64 {
        if self.inside(x, y) {
            return self.mesh_2d.get_ck(x, y);
        } else {
            return self.mesh_2d.nearest_cell_to(x, y).electronic_stopping_correction_factor;
        }
    }

    /// Determines the local mean free path from the formula sum(n(x, y))^(-1/3)
    pub fn mfp(&self, x: f64, y: f64) -> f64 {
        return self.total_number_density(x, y).powf(-1./3.);
    }

    /// Total number density of triangle that contains or is nearest to (x, y).
    pub fn total_number_density(&self, x: f64, y: f64) -> f64 {
        if self.inside(x, y) {
            return self.mesh_2d.get_densities(x, y).iter().sum::<f64>();
        } else {
            return self.mesh_2d.nearest_cell_to(x, y).densities.iter().sum::<f64>();
        }
        //return self.n.iter().sum::<f64>();
    }

    /// Lists number density of each species of triangle that contains or is nearest to (x, y).
    pub fn number_densities(&self, x: f64, y: f64) -> &Vec<f64> {
        if self.inside(x, y) {
            return &self.mesh_2d.get_densities(x, y);
        } else {
            return &self.mesh_2d.nearest_cell_to(x, y).densities;
        }
        //return self.n.iter().sum::<f64>();
    }

    /// Determines whether a point (x, y) is inside the energy barrier of the material.
    pub fn inside_energy_barrier(&self, x: f64, y: f64) -> bool {

        if self.mesh_2d.inside(x, y) {
            true
        } else {
            let nearest_cell = self.mesh_2d.nearest_cell_to(x, y);
            let distance = nearest_cell.distance_to(x, y);
            distance < self.mesh_2d.energy_barrier_thickness
        }
        //let p = point!(x: x, y: y);
        //return self.energy_surface.contains(&p);
    }

    /// Determines whether a point (x, y) is inside the simulation boundary.
    pub fn inside_simulation_boundary(&self, x:f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.mesh_2d.simulation_boundary.contains(&p);
    }

    /// Finds the closest point on the material boundary to the point (x, y).
    pub fn closest_point(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.mesh_2d.boundary.closest_point(&p);
        //return self.surface.closest_point(&p)
    }

    /// Finds the closest point on the simulation boundary to the point (x, y).
    pub fn closest_point_on_simulation_surface(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.mesh_2d.simulation_boundary.closest_point(&p)
    }

    /// Finds the average, concentration-weighted atomic number, Z_effective, of the triangle that contains or is nearest to (x, y).
    pub fn average_Z(&self, x: f64, y: f64) -> f64 {
        let concentrations = self.mesh_2d.get_concentrations(x, y);
        return self.Z.iter().zip(concentrations).map(|(charge, concentration)| charge*concentration).collect::<Vec<f64>>().iter().sum();
        //return self.Z.iter().sum::<f64>()/self.Z.len() as f64;
    }

    /// Finds the average, concentration-weighted atomic mass, m_effective, of the triangle that contains or is nearest to (x, y).
    pub fn average_mass(&self, x: f64, y: f64) -> f64 {
        let concentrations = self.mesh_2d.get_concentrations(x, y);
        return self.m.iter().zip(concentrations).map(|(mass, concentration)| mass*concentration).collect::<Vec<f64>>().iter().sum();
        //return self.m.iter().sum::<f64>()/self.m.len() as f64;
    }

    /// Finds the average, concentration-weighted bulk binding energy of the triangle that contains or is nearest to (x, y).
    pub fn average_bulk_binding_energy(&self, x: f64, y: f64) -> f64 {
        //returns average bulk binding energy
        let concentrations = self.mesh_2d.get_concentrations(x, y);
        return self.Eb.iter().zip(concentrations).map(|(bulk_binding_energy, concentration)| bulk_binding_energy*concentration).collect::<Vec<f64>>().iter().sum();
        //return self.Eb.iter().sum::<f64>()/self.Eb.len() as f64;
    }

    /// Finds the concentration-dependent surface binding energy of the triangle that contains or is nearest to (x, y).
    /// The surface binding energy is calculated using one of three methods:
    /// 1. INDIVIDUAL: the surface binding energies are set individually for each species, as Es.
    /// 2. TARGET: the surface binding energy is calculated as the local concentration-weighted average of the target surface binding energies, unless the particle has Es = 0, in which case it is 0.
    /// 3. AVERAGE: the surface binding energy is the average of the particle surface binding energy and the TARGET surface binding energy, unless either is 0 in which case it is 0.
    pub fn actual_surface_binding_energy(&self, particle: &particle::Particle, x: f64, y: f64) -> f64 {
        let concentrations = self.mesh_2d.get_concentrations(x, y);

        match self.surface_binding_model {
            SurfaceBindingModel::INDIVIDUAL => particle.Es,

            SurfaceBindingModel::TARGET => {
                if particle.Es == 0. {
                    0.
                } else {
                    self.Es.iter().zip(concentrations).map(|(surface_binding_energy, concentration)| surface_binding_energy*concentration).collect::<Vec<f64>>().iter().sum()
                }
            },

            SurfaceBindingModel::AVERAGE => {
                if (particle.Es == 0.) | (self.Es.iter().sum::<f64>() == 0.) {
                    0.
                } else {
                    0.5*(particle.Es + self.Es.iter().zip(concentrations).map(|(surface_binding_energy, concentration)| surface_binding_energy*concentration).collect::<Vec<f64>>().iter().sum::<f64>())
                }
            },

        }
    }

    ///The minimum cutoff energy of all species that make up the material.
    pub fn minimum_cutoff_energy(&self) -> f64 {
        let mut min_Ec = self.Ec.iter().sum::<f64>();
        for Ec in self.Ec.iter() {
            if min_Ec > *Ec {
                min_Ec = *Ec;
            }
        }
        return min_Ec;
    }

    ///Choose the parameters of a target atom as a concentration-weighted random draw from the species in the triangle that contains or is nearest to (x, y).
    pub fn choose(&self, x: f64, y: f64) -> (usize, f64, f64, f64, f64, usize) {
        let random_number: f64 = rand::random::<f64>();
        let cumulative_concentrations = self.get_cumulative_concentrations(x, y);

        for (component_index, cumulative_concentration) in cumulative_concentrations.iter().enumerate() {
            if random_number < *cumulative_concentration {
                return (component_index, self.Z[component_index], self.m[component_index], self.Ec[component_index], self.Es[component_index], self.interaction_index[component_index]);
            }
        }
        panic!("Input error: method choose() operation failed to choose a valid species. Check densities.");
    }

    /// Calculate the electronic stopping cross-sections using the mode set in [options].
    pub fn electronic_stopping_cross_sections(&self, particle_1: &super::particle::Particle, electronic_stopping_mode: ElectronicStoppingMode) -> Vec<f64> {

        let E = particle_1.E;
        let Ma = particle_1.m;
        let Za = particle_1.Z;

        let mut stopping_powers = Vec::with_capacity(self.Z.len());

        //Bragg's rule: total stopping power is sum of stopping powers of individual atoms
        //Therefore, compute each stopping power separately, and add them up

        let x = particle_1.pos.x;
        let y = particle_1.pos.y;
        let ck = self.electronic_stopping_correction_factor(x, y);

        for (n, Zb) in self.number_densities(x, y).iter().zip(&self.Z) {
            //let n = self.number_density(particle_1.pos.x, particle_1.pos.y);
            //let Zb = self.Z_eff(particle_1.pos.x, particle_1.pos.y);

            let beta = (1. - (1. + E/Ma/C.powf(2.)).powf(-2.)).sqrt();
            let v = beta*C;

            let I0 = match *Zb < 13. {
                true => 12. + 7./Zb,
                false => 9.76 + 58.5*Zb.powf(-1.19),
            };
            let I = Zb*I0*Q;

            //See Biersack and Haggmark - this looks like an empirical shell correction
            let B = match *Zb < 3. {
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
                ElectronicStoppingMode::INTERPOLATED => 1./(1./S_high + 1./S_low*ck),
                //Oen-Robinson
                ElectronicStoppingMode::LOW_ENERGY_LOCAL => S_low*ck,
                //Lindhard-Scharff
                ElectronicStoppingMode::LOW_ENERGY_NONLOCAL => S_low*ck,
                //Lindhard-Scharff and Oen-Robinson, using Lindhard Equipartition
                ElectronicStoppingMode::LOW_ENERGY_EQUIPARTITION => S_low*ck,
            };

            stopping_powers.push(stopping_power);
        }
        return stopping_powers;
    }
}

/// Calculate the effects of the planar surface binding potential of a material on a particle.
/// These effects include surface reflection and refraction of particles with non-zero surface binding energies.
pub fn surface_binding_energy(particle_1: &mut particle::Particle, material: &material::Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let x_old = particle_1.pos_old.x;
    let y_old = particle_1.pos_old.y;
    let cosx = particle_1.dir.x;
    let cosy = particle_1.dir.y;
    let cosz = particle_1.dir.z;
    let E = particle_1.E;

    //Actual surface binding energies
    let Es = material.actual_surface_binding_energy(particle_1, x_old, y_old);
    let Ec = particle_1.Ec;

    let inside_now = material.inside_energy_barrier(x, y);
    let inside_old = material.inside_energy_barrier(x_old, y_old);

    let leaving = !inside_now & inside_old;
    let entering = inside_now & !inside_old;

    if entering {
        if particle_1.backreflected {
            particle_1.backreflected = false;
        } else if let Closest::SinglePoint(p2) = material.closest_point(x, y) {
            let dx = p2.x() - x;
            let dy = p2.y() - y;
            let mag = (dx*dx + dy*dy).sqrt();

            let costheta = dx*cosx/mag + dy*cosy/mag;
            particle_1.E += Es;

            //Surface refraction via Snell's law
            let delta_theta = particle::refraction_angle(costheta, E, E + Es);
            particle::rotate_particle(particle_1, delta_theta, 0.);

            //particle_1.add_trajectory();
        } else {
            panic!("Numerical error: surface boundary algorithm encountered an error. Check geometry.");
        }
    }

    if leaving {
        if let Closest::SinglePoint(p2) = material.closest_point(x, y) {

            let dx = p2.x() - x;
            let dy = p2.y() - y;
            let mag = (dx*dx + dy*dy).sqrt();

            let costheta = dx*cosx/mag + dy*cosy/mag;
            let leaving_energy = E*costheta*costheta;

            if costheta < 0. {
                if leaving_energy > Es {

                    //particle_1.left = true;
                    particle_1.E += -Es;

                    //Surface refraction via Snell's law
                    let delta_theta = particle::refraction_angle(costheta, E, E - Es);
                    particle::rotate_particle(particle_1, delta_theta, 0.);

                    //particle_1.add_trajectory();

                } else {

                    //Specular reflection at local surface normal
                    particle_1.dir.x = -2.*(costheta)*dx/mag + cosx;
                    particle_1.dir.y = -2.*(costheta)*dy/mag + cosy;

                    particle_1.backreflected = true;
                    //particle_1.add_trajectory();
                }
            }
        } else {
            panic!("Numerical error: surface boundary algorithm encountered an error. Check geometry.");
        }
    }
}

/// Apply the boundary conditions of a material on a particle, including stopping, leaving, and reflection/refraction by/through the surface binding potential.
pub fn boundary_condition_2D_planar(particle_1: &mut particle::Particle, material: &material::Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let E = particle_1.E;
    let Ec = particle_1.Ec;

    if !material.inside_simulation_boundary(x, y) {
        particle_1.left = true;
        particle_1.add_trajectory();
    }

    if (E < Ec) & !particle_1.left & material.inside_energy_barrier(x, y) {
        particle_1.stopped = true;
        particle_1.add_trajectory();
    }

    if !particle_1.stopped & !particle_1.left {
        surface_binding_energy(particle_1, material);
    }
}
