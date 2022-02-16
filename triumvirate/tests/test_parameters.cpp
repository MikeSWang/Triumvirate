#include <cstdlib>
#include <fstream>

#include <sys/stat.h>

#include "parameters.hpp"

/**
 * Test initial parameter input.
 */
int main(int argc, char* argv[]) {
  Parameters params_ini;
  params_ini.read(argv);
  return 0;
}
