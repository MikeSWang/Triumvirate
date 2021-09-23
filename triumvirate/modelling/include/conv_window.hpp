#ifndef CONVOLUTION_H_INCLUDED_
#define CONVOLUTION_FFTLOG_H_INCLUDED_

#include <complex>
#include <vector>

#include <gsl/gsl_sf_coupling.h>

/// Provide standalone functions.
#define wigner_3j(j1, j2, j3, m1, m2, m3) ( \
  gsl_sf_coupling_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3) \
)

#define wigner_9j(j1, j2, j3, j4, j5, j6, j7, j8, j9) ( \
  gsl_sf_coupling_9j(2*j1, 2*j2, 2*j3, 2*j4, 2*j5, 2*j6, 2*j7, 2*j8, 2*j9) \
)

#define N_factor(ell1, ell2, ELL) ( \
  (2*ell1 + 1) * (2*ell2 + 1) * (2*ELL + 1) \
)

#define H_factor(ell1, ell2, ELL) ( \
  wigner_3j(ell1, ell2, ELL, 0, 0, 0) \
)

const double PI = 3.14159265358979323846;
const std::complex<double> I(0., 1.);

/**
 * Convolve three-point correlation function and window function by
 * wide-angle correction orders before Hankel transform.
 *
 * @param[in] ell1 First degree of the multipole.
 * @param[in] ell2 Second degree of the multipole.
 * @param[in] ELL Third, main degree of the multipole.
 * @param[in] ell_max Max degree of the multipoles to be summed over to.
 * @param[in] win_multipoles Wide-angle window correction multipoles.
 * @param[in] corr_multipoles Three-point correlation function wide-angle correction multipoles.
 * @param[out] array_out Output array prepared for double Hankel transform.
 */
int pre_hankel_win_conv_by_wide_angle_order(
  int ell1, int ell2, int ELL, int ell_max,
  std::vector<std::vector< std::complex<double> >> win_multipoles,
  std::vector<std::vector< std::complex<double> >> corr_multipoles,
  std::vector<std::vector< std::complex<double> >> array_out
) {
  /*
  std::complex<double> prefactor = (16 * PI * PI) * pow(I, -ell1 - ell2)
    * double(N_factor(ell1, ell2, ELL) * H_factor(ell1, ell2, ELL));
  */

  /// Get array dimensions and initialise output array.
  int ndim = sizeof(win_multipoles) / sizeof(win_multipoles[0]);

  for (int idx_row = 0; idx_row < ndim; idx_row++) {
    for (int idx_col = 0; idx_col < ndim; idx_col++) {
      array_out[idx_row][idx_col] = 0.;
    }
  }

  /// Perform convolution.
  for (int ell3 = 0; ell3 <= ell_max; ell3++) {
    for (int ell4 = 0; ell4 <= ell_max; ell4++) {
      for (int ell5 = abs(ell3 - ell4); ell5 <= ell3 + ell4; ell5++) {
        for (int ell6 = abs(ell1 - ell3); ell5 <= ell1 + ell3; ell6++) {
          for (int ell7 = abs(ell2 - ell4); ell7 <= ell2 + ell4; ell7++) {
            for (
              int ell8 = std::max(abs(ELL - ell5), abs(ell6 - ell7));
              ell8 <= std::min(ELL + ell5, ell6 + ell7);
              ell8++
            ) {
              /// Check for trivial terms.
              double H_136 = H_factor(ell1, ell3, ell6);
              if (H_136 < 1.e-10) continue;
              double H_247 = H_factor(ell2, ell4, ell7);
              if (H_247 < 1.e-10) continue;
              double H_L58 = H_factor(ELL, ell5, ell8);
              if (H_L58 < 1.e-10) continue;
              double H_345 = H_factor(ell3, ell4, ell5);
              if (H_345 < 1.e-10) continue;
              double H_678 = H_factor(ell6, ell7, ell8);
              if (H_678 < 1.e-10) continue;

              double nine_j = wigner_9j(
                ell1, ell2, ELL, ell3, ell4, ell5, ell6, ell7, ell8
              );
              if (nine_j < 1.e-10) continue;

              for (int idx_row = 0; idx_row < ndim; idx_row++) {
                for (int idx_col = 0; idx_col < ndim; idx_col++) {
                  array_out[idx_row][idx_col] += nine_j
                    * H_136 * H_247 * H_L58 / H_345 / H_678
                    * win_multipoles[idx_row][idx_col]
                    * corr_multipoles[idx_row][idx_col];
                }
              }
            }
          }
        }
      }
    }
  }

  return 0;
}

#endif
