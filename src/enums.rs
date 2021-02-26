use super::*;

/// Mode of electronic stopping to use.
#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum ElectronicStoppingMode {
    /// Biersack-Varelas interpolated electronic stopping. Valid for ~eV/nucleon to ~GeV/nucleon.
    INTERPOLATED,
    /// Oen-Robinson Firsov-type local electronic stopping. Valid up to ~25 keV/nucleon.
    LOW_ENERGY_LOCAL,
    /// Lindhard-Scharff nonlocal electronic stopping. Valid up to ~25 keV/nucleon.
    LOW_ENERGY_NONLOCAL,
    /// Equipartition between Oen-Robinson and Lindhard-Scharff electronic stopping formulas.
    LOW_ENERGY_EQUIPARTITION,
}

impl fmt::Display for ElectronicStoppingMode {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            ElectronicStoppingMode::INTERPOLATED => write!(f, "Biersack-Varelas electronic stopping"),
            ElectronicStoppingMode::LOW_ENERGY_NONLOCAL => write!(f, "Lindhard-Scharff electronic stopping"),
            ElectronicStoppingMode::LOW_ENERGY_LOCAL => write!(f, "Oen-Robinson electronic stopping"),
            ElectronicStoppingMode::LOW_ENERGY_EQUIPARTITION => write!(f, "Equipartition with Lindhard-Scharff and Oen-Robinson"),
        }
    }
}

/// Mode of surface binding energy calculation.
#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum SurfaceBindingModel {
    /// Surface binding energies will be determined individually depending only on the particle's `Es`.
    INDIVIDUAL,
    /// Surface binding energies will be a concentration-weighted average of material surface-binding energies, unless `particle.Es == 0` in which case it will be zero.
    TARGET,
    /// Surface binding energies will be the average of the particle and TARGET, unless either is zero in which case it will be zero.
    AVERAGE,
}

impl fmt::Display for SurfaceBindingModel {
    fn fmt (&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            SurfaceBindingModel::INDIVIDUAL => write!(f,
                "Individual surface binding energies."),
            SurfaceBindingModel::TARGET => write!(f,
                "Concentration-dependent linear combinaion of target binding energies."),
            SurfaceBindingModel::AVERAGE => write!(f,
                "Average between particle and concentration-dependent linear combination of target binding energies."),
        }
    }
}

/// Mean-free-path model.
#[derive(Deserialize, PartialEq, Clone, Copy)]
pub enum MeanFreePathModel {
    /// Constant mean-free-path for liquids and amorphous solids.
    LIQUID,
    /// Exponentially-distributed mean-free-paths for gases.
    GASEOUS,
}

impl fmt::Display for MeanFreePathModel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            MeanFreePathModel::LIQUID => write!(f, "Amorphous Solid/Liquid Model"),
            MeanFreePathModel::GASEOUS => write!(f, "Gaseous Model"),
        }
    }
}

/// Interatomic potentials between particles in rustbca.
#[derive(Deserialize, Clone, Copy)]
pub enum InteractionPotential {
    /// TRIDYN-style Kr-C. Equivalent to KR_C, except for the MAGIC constants.
    TRIDYN,
    /// Moliere's approximation to the Thomas-Fermi interatomic potential.
    MOLIERE,
    /// Krypton-Carbon "universal" interatomic potential.
    KR_C,
    /// Ziegler-Biersack-Littmark "unviversal" semi-empirical interatomic potential.
    ZBL,
    /// Lenz-Jensen screened Coulomb potential.
    LENZ_JENSEN,
    /// Lennard-Jones 12-6 potential, with user-defined sigma and epsilon.
    LENNARD_JONES_12_6 {sigma: f64, epsilon: f64},
    /// Lennard-Jones 6.5-6 potential, with user-defined sigma and epsilon.
    LENNARD_JONES_65_6 {sigma: f64, epsilon: f64},
    /// Morse potential, with user-defined D, alpha, and r0.
    MORSE{D: f64, alpha: f64, r0: f64},
    /// Tungsten-tungsten cubic spline potential (Following Bjorkas et al.)
    WW,
    /// Unscreened Coulombic interatomic potential between ions with charges Za and Zb.
    COULOMB{Za: f64, Zb: f64}
}

impl fmt::Display for InteractionPotential {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            InteractionPotential::TRIDYN => write!(f, "TRIDYN-style Kr-C (Different MAGIC constants)"),
            InteractionPotential::MOLIERE => write!(f, "Moliere Potential"),
            InteractionPotential::KR_C => write!(f, "Kr-C Potential"),
            InteractionPotential::ZBL => write!(f, "ZBL Potential"),
            InteractionPotential::LENZ_JENSEN => write!(f, "Lenz-Jensen Potential"),
            InteractionPotential::LENNARD_JONES_12_6{sigma, epsilon} => write!(f, "Lennard-Jones 12-6 Potential with sigma = {} A, epsilon = {} eV", sigma/ANGSTROM, epsilon/EV),
            InteractionPotential::LENNARD_JONES_65_6{sigma, epsilon} => write!(f, "Lennard-Jones 6.5-6 Potential with sigma = {} A, epsilon = {} eV", sigma/ANGSTROM, epsilon/EV),
            InteractionPotential::MORSE{D, alpha, r0} => write!(f, "Morse potential with D = {} eV, alpha = {} 1/A, and r0 = {} A", D/EV, alpha*ANGSTROM, r0/ANGSTROM),
            InteractionPotential::WW => write!(f, "W-W cubic spline interaction potential."),
            InteractionPotential::COULOMB{Za, Zb} => write!(f, "Coulombic interaction with Za = {} and Zb = {}", Za, Zb)
        }
    }
}

impl PartialEq for InteractionPotential {
    fn eq(&self, other: &Self) -> bool {
        discriminant(self) == discriminant(other)
    }
}

/// Method for solving the scattering integral.
#[derive(Deserialize, Clone, Copy)]
pub enum ScatteringIntegral {
    /// Mendenhall-Weller Gauss-Lobatto 4-point quadrature.
    MENDENHALL_WELLER,
    /// Ziegler's MAGIC algorithm.
    MAGIC,
    /// Gauss-Mehler n-point quadrature.
    GAUSS_MEHLER{n_points: usize},
    /// Gauss-Legendre 5-point quadrature.
    GAUSS_LEGENDRE
}

impl fmt::Display for ScatteringIntegral {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            ScatteringIntegral::MENDENHALL_WELLER => write!(f, "Mendenhall-Weller 4-Point Lobatto Quadrature"),
            ScatteringIntegral::MAGIC => write!(f, "MAGIC Algorithm"),
            ScatteringIntegral::GAUSS_MEHLER{n_points} => write!(f, "Gauss-Mehler {}-point Quadrature", n_points),
            ScatteringIntegral::GAUSS_LEGENDRE => write!(f, "Gauss-Legendre 5-point Quadrature"),
        }
    }
}

impl PartialEq for ScatteringIntegral {
    fn eq(&self, other: &Self) -> bool {
        discriminant(self) == discriminant(other)
    }
}

/// Root-finding algorithm.
#[derive(Deserialize, Clone, Copy)]
pub enum Rootfinder {
    /// Newton root-finder with user-defined `max_iterations` and `tolerance`.
    NEWTON{max_iterations: usize, tolerance: f64},
    CPR{n0: usize, nmax: usize, epsilon: f64, complex_threshold: f64, truncation_threshold: f64,
        far_from_zero: f64, interval_limit: f64, derivative_free: bool},
    POLYNOMIAL{complex_threshold: f64},
}

impl fmt::Display for Rootfinder {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Rootfinder::NEWTON{max_iterations, tolerance} => write!(f, "Newton-Raphson Rootfinder with maximum {} iterations and toleance = {}", max_iterations, tolerance),
            Rootfinder::CPR{n0, nmax, epsilon, complex_threshold, truncation_threshold, far_from_zero, interval_limit, derivative_free} =>
                write!(f, "Chebyshev-Proxy Rootfinder with {}-polishing", match derivative_free { true => "Secant", false => "Newton"}),
            Rootfinder::POLYNOMIAL{complex_threshold} => write!(f, "Frobenius Companion Matrix Polynomial Real Rootfinder with a complex tolerance of {}", complex_threshold),
        }
    }
}
impl PartialEq for Rootfinder {
    fn eq(&self, other: &Self) -> bool {
        discriminant(self) == discriminant(other)
    }
}
