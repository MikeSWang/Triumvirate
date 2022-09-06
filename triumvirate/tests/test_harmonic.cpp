#include <cstdlib>
#include <iostream>

#include "harmonic.hpp"

/**
 * Test the evaluation of normalised spherical harmonics.
 */
int main(int argc, char* argv[]) {
  /// Read the degree, order and vector components from command-line input.
  int ell = atoi(argv[1]);
  int m = atoi(argv[2]);
  double x = atof(argv[3]);
  double y = atof(argv[4]);
  double z = atof(argv[5]);

  double pos[3] = {x, y, z};

  /// Initialise the tool class.
  trv::maths::SphericalHarmonicCalculator ylm;

  /// Display the numerical result.
  std::cout
    << "y_" << ell << m << "(" << x << ", " << y << ", " << z << ") = "
    << ylm.calc_reduced_spherical_harmonic(ell, m, pos)
    << std::endl;

  return 0;
}
