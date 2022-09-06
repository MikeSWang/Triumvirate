#include "parameters.hpp"

/// TODO: Write test.
/**
 * Test initial parameter input.
 */
int main(int argc, char* argv[]) {
  trv::ParameterSet params_test;
  params_test.read_from_file(argv[1]);
  params_test.validate();
  params_test.printout();
  return 0;
}
