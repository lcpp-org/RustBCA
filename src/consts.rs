use super::*;

//Physical constants
///Fundamental charge in Coulombs.
pub const Q: f64 = 1.602176634E-19;
/// One electron-volt in Joules.
pub const EV: f64 = Q;
/// One atomic mass unit in kilograms.
pub const AMU: f64 = 1.66053906660E-27;
/// One Angstrom in meters.
pub const ANGSTROM: f64 = 1E-10;
/// One micron in meters.
pub const MICRON: f64 = 1E-6;
/// One nanometer in meters.
pub const NM: f64 = 1E-9;
/// One centimeter in meters.
pub const CM: f64 = 1E-2;
/// Vacuum permitivity in Farads/meter.
pub const EPS0: f64 = 8.8541878128E-12;
/// Bohr radius in meters.
pub const A0: f64 = 5.29177210903E-11;
/// Electron mass in kilograms.
pub const ME: f64 = 9.1093837015E-31;
/// sqrt(pi).
pub const SQRTPI: f64 = 2. / FRAC_2_SQRT_PI;
/// sqrt(2 * pi).
pub const SQRT2PI: f64 = 2. * SQRT_2 / FRAC_2_SQRT_PI;
/// Speed of light in meters/second.
pub const C: f64 = 299792458.;
/// Bethe-Bloch electronic stopping prefactor, in SI units.
pub const BETHE_BLOCH_PREFACTOR: f64 = 4.*PI*(Q*Q/(4.*PI*EPS0))*(Q*Q/(4.*PI*EPS0))/ME/C/C;
/// Lindhard-Scharff electronic stopping prefactor, in SI units.
pub const LINDHARD_SCHARFF_PREFACTOR: f64 = 1.212*ANGSTROM*ANGSTROM*Q;
/// Lindhard reduced energy prefactor, in SI units.
pub const LINDHARD_REDUCED_ENERGY_PREFACTOR: f64 = 4.*PI*EPS0/Q/Q;
