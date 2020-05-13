use super::*;

#[derive(Deserialize)]
pub struct MaterialParameters {
    pub energy_unit: String,
    pub mass_unit: String,
    pub Eb: f64,
    pub Es: f64,
    pub Ec: f64,
    pub n: f64,
    pub Z: f64,
    pub m: f64,
    pub electronic_stopping_correction_factor: f64
}

#[derive(Deserialize)]
pub struct Geometry {
    length_unit: String,
    surface: Vec<(f64, f64)>,
    energy_surface: Vec<(f64, f64)>,
    simulation_surface: Vec<(f64, f64)>
}

pub struct Material {
    pub n: f64,
    pub m: f64,
    pub Z: f64,
    pub Eb: f64,
    pub Es: f64,
    pub Ec: f64,
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
            m: material_parameters.m*mass_unit,
            Z: material_parameters.Z,
            Eb: material_parameters.Eb*energy_unit,
            Es: material_parameters.Es*energy_unit,
            Ec: material_parameters.Ec*energy_unit,
            electronic_stopping_correction_factor: material_parameters.electronic_stopping_correction_factor,
            surface: Polygon::new(LineString::from(unit_coords), vec![]),
            energy_surface: Polygon::new(LineString::from(unit_e_coords), vec![]),
            simulation_surface: Polygon::new(LineString::from(unit_b_coords), vec![])
        }
    }

    pub fn inside(&self, x: f64, y: f64) -> bool {
        let p = point!(x: x, y: y);
        return self.surface.contains(&p);
    }

    pub fn inside_1D(&self, x: f64) -> bool {
        return x > 0.;
    }

    pub fn inside_energy_barrier_1D(&self, x: f64) -> bool {
        let dx = -2.*self.n.powf(-1./3.)/SQRT2PI;
        return x > dx;
    }

    pub fn inside_simulation_boundary_1D(&self, x: f64) -> bool {
        let dx = -2.*self.n.powf(-1./3.)/SQRT2PI;
        return x > 2.*dx;
    }

    pub fn mfp(&self, x: f64, y: f64) -> f64 {
        return self.n.powf(-1./3.);
    }

    pub fn number_density(&self, x: f64, y: f64) -> f64 {
        return self.n;
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

    pub fn Z_eff(&self, x: f64, y: f64) -> f64 {
        return self.Z;
    }

    pub fn m_eff(&self, x: f64, y: f64) -> f64 {
        return self.m;
    }

    pub fn Eb(&self, x: f64, y: f64) -> f64 {
        return self.Eb;
    }

    pub fn choose(&self, x: f64, y: f64) -> (f64, f64, f64, f64) {
        return (self.Z, self.m, self.Ec, self.Es);
    }

    pub fn electronic_stopping_power(&self, particle_1: &super::particle::Particle, electronic_stopping_mode: i32) -> f64 {

        let n = self.number_density(particle_1.pos.x, particle_1.pos.y);
        let E = particle_1.E;
        let Ma = particle_1.m;
        let Za = particle_1.Z;
        let Zb = self.Z_eff(particle_1.pos.x, particle_1.pos.y);

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

        return stopping_power;
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

            if (costheta < 0.) {
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

pub fn boundary_condition_1D_planar(particle_1: &mut particle::Particle, material: &material::Material) {
    let x = particle_1.pos.x;
    let y = particle_1.pos.y;
    let cosx = particle_1.dir.x;
    let cosy = particle_1.dir.y;
    let cosz = particle_1.dir.z;
    let E = particle_1.E;
    let xc = -2.*material.number_density(x, y).powf(-1./3.)/SQRT2PI;
    let Es = particle_1.Es;
    let Ec = particle_1.Ec;

    if (x < xc) & (cosx < 0.) {
        let leaving_energy = E*cosx*cosx;

        if leaving_energy > Es {
            particle_1.left = true;
            particle_1.E += -Es;

            let cosx_new = ((E*cosx*cosx - Es)/(E - Es)).sqrt();
            let sinx = (1. - cosx*cosx).sqrt();
            let sinx_new = (1. - cosx_new*cosx_new).sqrt();
            particle_1.dir.x = -cosx_new;
            particle_1.dir.y *= sinx_new/sinx;
            particle_1.dir.z *= sinx_new/sinx;
            particle_1.dir.normalize();
            particle_1.add_trajectory();

        } else {
            particle_1.dir.x = cosx.abs();
        }
    }

    if (E < Ec) & !particle_1.left {
        particle_1.stopped = true;
        particle_1.add_trajectory();
    }
}
