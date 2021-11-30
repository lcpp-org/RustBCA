#include "RustBCA.h"
#include <iostream>
#include <vector>

int main(int argc, char * argv[]) {
  OutputTaggedBCA output;
  double velocities[2][3] = {{500000.0, 0.1, 0.0}, {500000.0, 0.1, 0.0}};
  double positions[2][3] = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
  int tags[2] = {0, 1};
  double weights[2] = {1.0, 1.0};
  double Z[3] = {74.0, 74.0};
  double m[2] = {184.0, 184.0};
  double n[2] = {0.06306, 0.06306};
  double Ec[2] = {1.0, 1.0};
  double Es[2] = {8.79, 8.79};
  double Eb[2] = {0.0, 0.0};

  InputTaggedBCA input = {
    2,
    positions,
    velocities,
    1.0,
    1.0,
    1.0,
    1.0,
    2,
    Z,
    m,
    n,
    Ec,
    Es,
    Eb,
    tags,
    weights
  };

  //output = simple_bca_c(0., 0., 0., 0.5, 0.5, 0.00, 2000.0, 2.0, 4.0, 1.0, 0.0, 74.0, 184.0, 1.0, 8.79, 0.06306, 0.0);
  //output = compound_bca_list_c(input);
  output = compound_tagged_bca_list_c(input);

  std::cout << "Particle 1 Z: ";
  std::cout << output.particles[0][0];
  std::cout << std::endl;
  std::cout << "Particle 1 E [eV]: ";
  std::cout << output.particles[0][2];
  std::cout << std::endl;
  std::cout << "Particle 2 Z: ";
  std::cout << output.particles[1][0];
  std::cout << std::endl;
  std::cout << "Particle 2 E [eV]: ";
  std::cout << output.particles[1][2];
  std::cout << std::endl;
  return 0;
}
