#include <fstream>

#include "fftlog.hpp"

int main() {
  const int ell = 0;
  const int m = 2;

  /// Load test data.
  std::ifstream fin;
  std::string fline;

  char test_fname[] = "triumvirate/tests/test_input/clustats/test_pk0.dat";

  fin.open(test_fname, std::ios::in);
  int lineno = 0;
  while (std::getline(fin, fline)) {lineno++;}
  fin.close();

  int Nk = lineno;
  double k[Nk], pk[Nk];

  fin.open(test_fname, std::ios::in);
  lineno = 0;
  while (std::getline(fin, fline)) {
    sscanf(fline.data(), "%lf %lf", &k[lineno], &pk[lineno]);
    lineno++;
  }
  fin.close();

  /// Perform transform.
  double r[Nk], xi[Nk];

  trv::maths::sj_transform(ell, m, Nk, k, pk, r, xi);
  // trv::transform_powspec_to_corrfunc_multipole(ell, Nk, k, pk, r, xi);

  /// Save test results.
  char test_fname_out[]= "triumvirate/tests/test_output/xi0.dat";
  FILE *test_file_out = std::fopen(test_fname_out, "w");

  for (int i = 0; i < Nk; i++) {
    std::fprintf(test_file_out, "%.9e %.9e\n", r[i], xi[i]);
  }

  std::fclose(test_file_out);

  return 0;
}
