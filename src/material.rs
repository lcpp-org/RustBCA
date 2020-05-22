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

#[derive(Deserialize)]
pub struct Mesh2DInput {
    length_unit: String,
    coordinate_sets: Vec<(f64, f64, f64, f64, f64, f64)>,
    data: Vec<Vec<f64>>,
}

pub struct Mesh2D {
    mesh: Vec<Cell2D>,
    boundary: Polygon<f64>
}
impl Mesh2D {
    pub fn new(mesh_2d_input: Mesh2DInput) -> Mesh2D {

        let coordinate_sets = mesh_2d_input.coordinate_sets;
        let data = mesh_2d_input.data;

        assert_eq!(coordinate_sets.len(), data.len(), "Coordinates and data of unequal length.");
        let n = coordinate_sets.len();

        let mut cells: Vec<Cell2D> =  Vec::with_capacity(n);
        let mut points: Vec<(f64, f64)> = Vec::with_capacity(n);

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

        for (coordinate_set, densities) in coordinate_sets.iter().zip(data) {
            cells.push(Cell2D::new(coordinate_set.clone(), densities));

            points.push((coordinate_set.0*length_unit, coordinate_set.3*length_unit));
            points.push((coordinate_set.1*length_unit, coordinate_set.4*length_unit));
            points.push((coordinate_set.2*length_unit, coordinate_set.5*length_unit));
        }
        let multi_point: MultiPoint<f64> = points.into();
        let boundary: Polygon<f64> = multi_point.convex_hull();

        Mesh2D {
            mesh: cells,
            boundary: boundary
        }
    }

    pub fn get_densities(&self, x: f64, y: f64) -> Vec<f64> {
        for cell in &self.mesh {
            if cell.contains(x, y) {
                return cell.data.clone();
            }
        }
        panic!("Point ({}, {}) not found in any cell of the mesh.", x, y);
    }

    pub fn get_total_density(&self, x: f64, y: f64) -> f64 {
        for cell in &self.mesh {
            if cell.contains(x, y) {
                return cell.data.clone().iter().sum::<f64>();
            }
        }
        panic!("Point ({}, {}) not found in any cell of the mesh.", x, y);
    }

    pub fn get_concentrations(&self, x: f64, y: f64) -> Vec<f64> {
        if self.inside(x, y) {
            for cell in &self.mesh {
                if cell.contains(x, y) {
                    let densities = cell.data.clone();
                    let total_density: f64 = cell.data.clone().iter().sum();
                    return densities.iter().map(|&i| i/total_density).collect::<Vec<f64>>()
                }
            }
            panic!("Geometry error.")
        } else {
            let densities = self.nearest_cell_to(x, y).data.clone();
            let total_density: f64 = densities.clone().iter().sum();
            return densities.iter().map(|&i| i/total_density).collect::<Vec<f64>>()
        }
    }

    pub fn inside(&self, x: f64, y: f64) -> bool {
        self.boundary.contains(&Point::new(x, y))
    }

    pub fn nearest_cell_to(&self, x: f64, y: f64) -> &Cell2D {

        let mut min_distance: f64 = std::f64::MAX;
        let mut index: usize = 0;

        for (cell_index, cell) in self.mesh.iter().enumerate() {
            let distance_to = cell.distance_to(x, y);
            if distance_to < min_distance {
                min_distance = distance_to;
                index = cell_index;
            }
        }

        return &self.mesh[index];
    }
}

pub struct Cell2D {
    triangle: Triangle2D,
    pub data: Vec<f64>
}
impl Cell2D {
    pub fn new(coordinate_set: (f64, f64, f64, f64, f64, f64), data: Vec<f64>) -> Cell2D {
        Cell2D {
            triangle: Triangle2D::new(coordinate_set),
            data: data
        }
    }

    pub fn contains(&self, x: f64, y: f64) -> bool {
        self.triangle.contains(x, y)
    }

    pub fn distance_to(&self, x: f64, y: f64) -> f64 {
        self.triangle.distance_to(x, y)
    }
}

pub struct Triangle2D {
    x1: f64,
    x2: f64,
    x3: f64,
    y1: f64,
    y2: f64,
    y3: f64,
}

impl Triangle2D {
    pub fn new(coordinate_set: (f64, f64, f64, f64, f64, f64)) -> Triangle2D {

        let x1 = coordinate_set.0;
        let x2 = coordinate_set.1;
        let x3 = coordinate_set.2;

        let y1 = coordinate_set.3;
        let y2 = coordinate_set.4;
        let y3 = coordinate_set.5;

        Triangle2D {
            x1: x1,
            x2: x2,
            x3: x3,
            y1: y1,
            y2: y2,
            y3: y3
        }
    }

    pub fn contains(&self, x: f64, y: f64) -> bool {
        let x1 = self.x1;
        let x2 = self.x2;
        let x3 = self.x3;
        let y1 = self.y1;
        let y2 = self.y2;
        let y3 = self.y3;

        let a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
        let b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
        let c = 1. - a - b;

         (0. <= a) & (a <= 1.) & (0. <= b) & (b <= 1.) & (0. <= c) & (c <= 1.)
    }

    pub fn centroid(&self) -> (f64, f64) {
        ((self.x1 + self.x2 + self.x3)/3., (self.y1 + self.y2 + self.y3)/3.)
    }

    pub fn distance_to(&self, x: f64, y: f64) -> f64 {
        let centroid = self.centroid();
        ((x - centroid.0).powf(2.) + (y - centroid.1).powf(2.)).sqrt()
    }
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
    pub electronic_stopping_correction_factor: f64,
    pub mesh_2d: Mesh2D

}
impl Material {
    pub fn new(material_parameters: MaterialParameters, geometry: Geometry, mesh_2d_input: Mesh2DInput) -> Material {

        //Multiply all coordinates by value of geometry unit.
        let length_unit: f64 = match geometry.length_unit.as_str() {
            "MICRON" => MICRON,
            "CM" => CM,
            "ANGSTROM" => ANGSTROM,
            "NM" => NM,
            "M" => 1.,
            _ => panic!("Input error: incorrect unit {} in input file. Choose one of: MICRON, CM, ANGSTROM, NM, M",
                geometry.length_unit.as_str())
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
            simulation_surface: Polygon::new(LineString::from(unit_b_coords), vec![]),
            mesh_2d: Mesh2D::new(mesh_2d_input),
        }
    }

    pub fn get_concentrations(&self, x: f64, y: f64) -> Vec<f64> {

        let total_number_density: f64 = if self.inside(x, y) {
            self.mesh_2d.get_total_density(x,y)
        } else {
            self.mesh_2d.nearest_cell_to(x, y).data.iter().sum::<f64>()
        };

        return self.mesh_2d.get_densities(x, y).iter().map(|&i| i / total_number_density).collect();
    }

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

    pub fn inside(&self, x: f64, y: f64) -> bool {
        return self.mesh_2d.inside(x, y);
        //let p = point!(x: x, y: y);
        //return self.surface.contains(&p);
    }

    pub fn mfp(&self, x: f64, y: f64) -> f64 {
        return self.total_number_density(x, y).powf(-1./3.);
    }

    pub fn total_number_density(&self, x: f64, y: f64) -> f64 {
        if self.inside(x, y) {
            return self.mesh_2d.get_densities(x, y).iter().sum::<f64>();
        } else {
            return self.mesh_2d.nearest_cell_to(x, y).data.iter().sum::<f64>();
        }
        //return self.n.iter().sum::<f64>();
    }

    pub fn inside_energy_barrier(&self, x: f64, y: f64) -> bool {
        return self.mesh_2d.inside(x, y);
        //let p = point!(x: x, y: y);
        //return self.energy_surface.contains(&p);
    }

    pub fn inside_simulation_boundary(&self, x:f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.simulation_surface.contains(&p);
    }

    pub fn closest_point(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.mesh_2d.boundary.closest_point(&p);
        //return self.surface.closest_point(&p)
    }

    pub fn closest_point_on_energy_barrier(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.mesh_2d.boundary.closest_point(&p);
        //return self.energy_surface.closest_point(&p)
    }

    pub fn closest_point_on_simulation_surface(&self, x: f64, y: f64) -> Closest<f64> {
        let p = point!(x: x, y: y);
        return self.simulation_surface.closest_point(&p)
    }

    pub fn average_Z(&self, x: f64, y: f64) -> f64 {
        let concentrations = self.mesh_2d.get_concentrations(x, y);
        return self.Z.iter().zip(concentrations).map(|(i1, i2)| i1*i2).collect::<Vec<f64>>().iter().sum();
        //return self.Z.iter().sum::<f64>()/self.Z.len() as f64;
    }

    pub fn average_mass(&self, x: f64, y: f64) -> f64 {
        let concentrations = self.mesh_2d.get_concentrations(x, y);
        return self.m.iter().zip(concentrations).map(|(i1, i2)| i1*i2).collect::<Vec<f64>>().iter().sum();
        //return self.m.iter().sum::<f64>()/self.m.len() as f64;
    }

    pub fn average_bulk_binding_energy(&self, x: f64, y: f64) -> f64 {
        //returns average bulk binding energy
        let concentrations = self.mesh_2d.get_concentrations(x, y);
        return self.Eb.iter().zip(concentrations).map(|(i1, i2)| i1*i2).collect::<Vec<f64>>().iter().sum();
        //return self.Eb.iter().sum::<f64>()/self.Eb.len() as f64;
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
        let cumulative_concentrations = self.get_cumulative_concentrations(x, y);

        for (component_index, cumulative_concentration) in cumulative_concentrations.iter().enumerate() {
            if random_number < *cumulative_concentration {
                return (component_index, self.Z[component_index], self.m[component_index], self.Ec[component_index], self.Es[component_index]);
            }
        }
        panic!("Input error: species choice operation failed. Check densities.");
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
                _ => panic!("Input error: unimplemented electronic stopping mode. Use 0: Biersack-Varelas 1: Lindhard-Scharff 2: Oen-Robinson 3: Equipartition")
            };

            stopping_powers.push(stopping_power);
        }
        return stopping_powers;
    }
}

pub fn surface_binding_energy(particle_1: &mut particle::Particle, material: &material::Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let x_old = particle_1.pos_old.x;
    let y_old = particle_1.pos_old.y;
    let cosx = particle_1.dir.x;
    let cosy = particle_1.dir.y;
    let cosz = particle_1.dir.z;
    let E = particle_1.E;
    let Es = particle_1.Es;
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

            particle_1.add_trajectory();
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

                    particle_1.add_trajectory();

                } else {

                    //Specular reflection at local surface normal
                    particle_1.dir.x = -2.*(costheta)*dx/mag + cosx;
                    particle_1.dir.y = -2.*(costheta)*dy/mag + cosy;

                    particle_1.backreflected = true;
                    particle_1.add_trajectory();
                }
            }
        } else {
            panic!("Numerical error: surface boundary algorithm encountered an error. Check geometry.");
        }
    }
}

pub fn boundary_condition_2D_planar(particle_1: &mut particle::Particle, material: &material::Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let x_old = particle_1.pos_old.x;
    let y_old = particle_1.pos_old.y;
    let cosx = particle_1.dir.x;
    let cosy = particle_1.dir.y;
    let cosz = particle_1.dir.z;
    let E = particle_1.E;
    let Es = particle_1.Es;
    let Ec = particle_1.Ec;

    surface_binding_energy(particle_1, material);

    if !material.inside_simulation_boundary(x, y) {
        particle_1.left = true;
        particle_1.add_trajectory();
    }

    if (E < Ec) & !particle_1.left {
        if material.inside_energy_barrier(x, y) {
            particle_1.stopped = true;
            particle_1.add_trajectory();
        }
    }
}
