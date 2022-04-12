/**
 * @file tools.hpp
 * @brief Miscellaneous tools.
 *
 */

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

const std::complex<double> M_I(0., 1.);  ///< imaginary unit

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
   */
  static void set_kbin(ParameterSet& params, double* kbin_out) {
    double dk = (params.kmax - params.kmin) / double(params.num_kbin - 1);

    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      kbin_out[ibin] = params.kmin + dk * ibin;
    }
  }

  /**
   * Set separation bins.
   *
   * @param[in] params Parameter set.
   * @param[out] rbin_out Separation bins.
   */
  static void set_rbin(ParameterSet& params, double* rbin_out) {
    if (
      params.binning == "custom"
    ) {
      /// RFE: Insert customised binning code here.
    } else if (
      params.binning == "logpad"
    ) {
      /// CAVEAT: Discretionary choices.
      int nbin_custom = 5;

      for (int ibin = 0; ibin < nbin_custom; ibin++) {
        rbin_out[ibin] = 10. * ibin;
      }

      double rmin = 10. * nbin_custom;

      double dlnr = (log(params.rmax) - log(rmin))
        / double((params.num_rbin - nbin_custom) - 1);

      for (int ibin = nbin_custom; ibin < params.num_rbin; ibin++) {
        rbin_out[ibin] = rmin * exp(dlnr * (ibin - nbin_custom));
      }
    } else if (
      params.binning == "linpad"
    ) {
      /// CAVEAT: Discretionary choices.
      int nbin_custom = 5;

      for (int ibin = 0; ibin < nbin_custom; ibin++) {
        rbin_out[ibin] = 10. * ibin;
      }

      double rmin = 10. * nbin_custom;

      double dr = (params.rmax - rmin)
        / double((params.num_rbin - nbin_custom) - 1);

      for (int ibin = nbin_custom; ibin < params.num_rbin; ibin++) {
        rbin_out[ibin] = rmin + dr * (ibin - nbin_custom);
      }
    } else if (
      params.binning == "log"
    )  {
      double rmin;
      if (params.rmin == 0.) {
        /// CAVEAT: Discretionary choice.
        rmin = 1.;
      } else {
        rmin = params.rmin;
      }

      double dlnr = (log(params.rmax) - log(rmin))
        / double(params.num_rbin - 1);

      for (int ibin = 0; ibin < params.num_rbin; ibin++) {
        rbin_out[ibin] = rmin * exp(dlnr * ibin);
      }
    } else {  // by default ``params.binning == "lin"``
      double dr = (params.rmax - params.rmin) / double(params.num_rbin - 1);

      for (int ibin = 0; ibin < params.num_rbin; ibin++) {
        rbin_out[ibin] = params.rmin + dr * ibin;
      }
    }
  }
};

#endif  // TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_
