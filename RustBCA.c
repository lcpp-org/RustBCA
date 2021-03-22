#include "RustBCA.h"
#include <iostream>
#include <vector>

int main() {
  OutputBCA output;

  output = simple_bca_c(0., 0., 0., 0.5, 0.5, 0.00, 2000.0, 2.0, 4.0, 1.0, 0.0, 74.0, 184.0, 1.0, 8.79, 0.06306, 0.0);

  std::cout << output.len;
  std::cout << output.particles[1][0];

  return 0.;
}
