#ifndef TRIUM_TOOLS_H_INCLUDED_
#define TRIUM_TOOLS_H_INCLUDED_

#ifndef TRIUM_PARAMETERS_H_INCLUDED_
#include "parameters.hpp"
#endif

/**
 * Collection of tools.
 *
 * This includes the set up of isotropic wavenumber or separation bins
 * in configuration or Fourier space.
 *
 */
class ToolCollection {
 public:
  /**
   * Set wavenumber bins.
   *
   * @param[in] params (Reference to) the input parameter set.
   * @param[out] kbin_out (Pointer to) set wavenumber bins.
   * @returns Exit status.
   */
  static int set_kbin(ParameterSet& params, double* kbin_out) {
    double dk = (params.kmax - params.kmin) / double(params.num_kbin - 1);
    for (int i = 0; i < params.num_kbin; i++) {
      kbin_out[i] = params.kmin + dk * i;
    }
    return 0;
  }

  /**
   * Set separation bins.
   *
   * @param[in] params (Reference to) the input parameter set.
   * @param[out] rbin_out (Pointer to) set separation bins.
   * @returns Exit status.
   */
  static int set_rbin(ParameterSet& params, double* rbin_out) {
    if (0) {
    } else if (params.binning == "log")  {
      double rmin;
      if (params.rmin == 0.) {
        rmin = 1.;
      } else {
        rmin = params.rmin;
      }
      double dlnr = (log(params.rmax) - log(rmin))
        / double(params.num_rbin - 1);

      for (int i = 0; i < params.num_rbin; i++) {
        rbin_out[i] = rmin * exp(dlnr * i);
      }
    } else if (params.binning == "linpad") {
      rbin_out[0] = 0.;
      rbin_out[1] = 10.;
      rbin_out[2] = 20.;
      rbin_out[3] = 30.;
      rbin_out[4] = 40.;
      double rmin = 50.;
      double dr = (params.rmax - rmin) / double((params.num_rbin - 5) - 1);

      for (int i = 5; i < params.num_rbin; i++) {
        rbin_out[i] = rmin + dr * (i - 5);
      }
    } else if (params.binning == "logpad") {
      rbin_out[0] = 0.;
      rbin_out[1] = 10.;
      rbin_out[2] = 20.;
      rbin_out[3] = 30.;
      rbin_out[4] = 40.;
      double rmin = 50.;
      double dlnr = (log(params.rmax) - log(rmin))
        / double((params.num_rbin - 5) - 1);

      for (int i = 5; i < params.num_rbin; i++) {
        rbin_out[i] = rmin * exp(dlnr * (i - 5));
      }
    } else {  // generically ``params.binning == "lin"``
      double dr = (params.rmax - params.rmin) / double(params.num_rbin - 1);
      for (int i = 0; i < params.num_rbin; i++) {
        rbin_out[i] = params.rmin + dr * i;
      }
    }

    return 0;
  }
};

#endif
