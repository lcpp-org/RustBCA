#include <cstdarg>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <new>
#include <cmath>

static const double PI = M_PI;

static const double FRAC_2_SQRT_PI = M_2_PI;

static const double SQRT_2 = M_SQRT2;

///Fundamental charge in Coulombs.
static const double Q = 1.602176634e-19;

/// One electron-volt in Joules.
static const double EV = Q;

/// One atomic mass unit in kilograms.
static const double AMU = 1.66053906660e-27;

/// One Angstrom in meters.
static const double ANGSTROM = 1e-10;

/// One micron in meters.
static const double MICRON = 1e-6;

/// One nanometer in meters.
static const double NM = 1e-9;

/// One centimeter in meters.
static const double CM = 1e-2;

/// Vacuum permitivity in Farads/meter.
static const double EPS0 = 8.8541878128e-12;

/// Bohr radius in meters.
static const double A0 = 5.29177210903e-11;

/// Electron mass in kilograms.
static const double ME = 9.1093837015e-31;

/// sqrt(pi).
static const double SQRTPI = (2. / FRAC_2_SQRT_PI);

/// sqrt(2 * pi).
static const double SQRT2PI = ((2. * SQRT_2) / FRAC_2_SQRT_PI);

/// Speed of light in meters/second.
static const double C = 299792458.;

/// Bethe-Bloch electronic stopping prefactor, in SI units.
static const double BETHE_BLOCH_PREFACTOR = ((((((4. * PI) * ((Q * Q) / ((4. * PI) * EPS0))) * ((Q * Q) / ((4. * PI) * EPS0))) / ME) / C) / C);

/// Lindhard-Scharff electronic stopping prefactor, in SI units.
static const double LINDHARD_SCHARFF_PREFACTOR = (((1.212 * ANGSTROM) * ANGSTROM) * Q);

/// Lindhard reduced energy prefactor, in SI units.
static const double LINDHARD_REDUCED_ENERGY_PREFACTOR = ((((4. * PI) * EPS0) / Q) / Q);

struct OutputBCA {
  uintptr_t len;
  double (*particles)[9];
};

struct InputSimpleBCA {
  uintptr_t len;
  /// vx, vy, vz
  double (*velocities)[3];
  double Z1;
  double m1;
  double Ec1;
  double Es1;
  double Z2;
  double m2;
  double n2;
  double Ec2;
  double Es2;
  double Eb2;
};

struct InputCompoundBCA {
  uintptr_t len;
  /// vx, vy, vz
  double (*velocities)[3];
  double Z1;
  double m1;
  double Ec1;
  double Es1;
  uintptr_t num_species_target;
  double *Z2;
  double *m2;
  double *n2;
  double *Ec2;
  double *Es2;
  double *Eb2;
};

struct OutputTaggedBCA {
  uintptr_t len;
  double (*particles)[9];
  double *weights;
  int32_t *tags;
  bool *incident;
};

struct InputTaggedBCA {
  uintptr_t len;
  /// x y z
  double (*positions)[3];
  /// vx, vy, vz
  double (*velocities)[3];
  double Z1;
  double m1;
  double Ec1;
  double Es1;
  uintptr_t num_species_target;
  double *Z2;
  double *m2;
  double *n2;
  double *Ec2;
  double *Es2;
  double *Eb2;
  int32_t *tags;
  double *weights;
};

extern "C" {

OutputTaggedBCA compound_tagged_bca_list_c(InputTaggedBCA input);

OutputBCA simple_bca_list_c(InputSimpleBCA input);

OutputBCA compound_bca_list_c(InputCompoundBCA input);

OutputBCA simple_bca_c(double x,
                       double y,
                       double z,
                       double ux,
                       double uy,
                       double uz,
                       double E1,
                       double Z1,
                       double m1,
                       double Ec1,
                       double Es1,
                       double Z2,
                       double m2,
                       double Ec2,
                       double Es2,
                       double n2,
                       double Eb2);

} // extern "C"
