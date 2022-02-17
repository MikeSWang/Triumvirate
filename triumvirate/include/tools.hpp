#ifndef TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_

#include <cmath>

#include <gsl/gsl_sf_coupling.h>

#include "parameters.hpp"

/**
 * Calculate Wigner 3-j symbol.
 */
#define wigner_3j(j1, j2, j3, m1, m2, m3) ( \
  gsl_sf_coupling_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3) \
)

/**
 * Line-of-sight vector.
 *
 */
struct LineOfSight{
  double pos[3];   ///< 3-d position vector
};

/**
 * Binning scheme.
 *
 * This sets up isotropic wavenumber or separation bins
 * in configuration or Fourier space.
 *
 */
class BinScheme {
 public:
  /**
   * Set wavenumber bins.
   *
   * @param[in] params Parameter set.
   * @param[out] kbin_out Wavenumber bins.
   * @returns Exit status.
   */
  static int set_kbin(Parameters& params, double* kbin_out) {
    double dk = (params.kmax - params.kmin) / double(params.num_kbin - 1);

    for (int i = 0; i < params.num_kbin; i++) {
      kbin_out[i] = params.kmin + dk * i;
    }

    return 0;
  }

  /**
   * Set separation bins.
   *
   * @param[in] params Parameter set.
   * @param[out] rbin_out Separation bins.
   * @returns Exit status.
   */
  static int set_rbin(Parameters& params, double* rbin_out) {
    if (0) {
    } else if (params.binning == "custom") {
      /// RFE: Insert customised binning code here.
    } else if (params.binning == "logpad") {
      /// CAVEAT: Discretionary choices.
      int nbins_custom = 5;

      for (int idx_bin = 0; idx_bin < nbins_custom; idx_bin++) {
        rbin_out[idx_bin] = 10. * idx_bin;
      }

      double rmin = 10. * nbins_custom;

      double dlnr = (log(params.rmax) - log(rmin))
        / double((params.num_rbin - nbins_custom) - 1);

      for (int idx_bin = nbins_custom; idx_bin < params.num_rbin; idx_bin++) {
        rbin_out[idx_bin] = rmin * exp(dlnr * (idx_bin - nbins_custom));
      }
    } else if (params.binning == "linpad") {
      /// CAVEAT: Discretionary choices.
      int nbins_custom = 5;

      for (int idx_bin = 0; idx_bin < nbins_custom; idx_bin++) {
        rbin_out[idx_bin] = 10. * idx_bin;
      }

      double rmin = 10. * nbins_custom;

      double dr = (params.rmax - rmin)
        / double((params.num_rbin - nbins_custom) - 1);

      for (int idx_bin = nbins_custom; idx_bin < params.num_rbin; idx_bin++) {
        rbin_out[idx_bin] = rmin + dr * (idx_bin - nbins_custom);
      }
    } else if (params.binning == "log")  {
      double rmin;
      if (params.rmin == 0.) {
        /// CAVEAT: Discretionary choice.
        rmin = 1.;
      } else {
        rmin = params.rmin;
      }

      double dlnr = (log(params.rmax) - log(rmin))
        / double(params.num_rbin - 1);

      for (int idx_bin = 0; idx_bin < params.num_rbin; idx_bin++) {
        rbin_out[idx_bin] = rmin * exp(dlnr * idx_bin);
      }
    } else {  // by default ``params.binning == "lin"``
      double dr = (params.rmax - params.rmin) / double(params.num_rbin - 1);

      for (int idx_bin = 0; idx_bin < params.num_rbin; idx_bin++) {
        rbin_out[idx_bin] = params.rmin + dr * idx_bin;
      }
    }

    return 0;
  }
};

#endif
