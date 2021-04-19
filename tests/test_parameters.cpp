#include <cstdlib>
#include <fstream>

#include <sys/stat.h>

#include "parameters.hpp"

/**
 * Test initial parameter input.
 */
int main(int argc, char* argv[]) {
  ParameterSet parameters_ini;
  parameters_ini.read_parameters(argv);
  return 0;
}
