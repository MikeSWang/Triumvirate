#include "parameters.hpp"

/// TODO: Write test.
/**
 * Test initial parameter input.
 */
int main(int argc, char* argv[]) {
  trv::scheme::ParameterSet params_test;
  params_test.read_from_file(argv);
  params_test.printout();
  return 0;
}
