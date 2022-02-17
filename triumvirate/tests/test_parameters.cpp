#include "parameters.hpp"

/// TODO: Write test.
/**
 * Test initial parameter input.
 */
int main(int argc, char* argv[]) {
  Parameters params_test;
  params_test.read_from_file(argv);
  params_test.print_to_file();
  return 0;
}
