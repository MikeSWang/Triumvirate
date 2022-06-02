/**
 * @file tools.hpp
 * @brief Miscellaneous scientific tools.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_

#include <cmath>
#include <complex>
#include <cstdarg>
#include <stdexcept>
#include <vector>

#include "parameters.hpp"

const std::complex<double> M_I(0., 1.);  ///< imaginary unit

namespace trv {

/// //////////////////////////////////////////////////////////////////////
/// Mathematics
/// //////////////////////////////////////////////////////////////////////

namespace maths {

/**
 * Evaluate a complex number f@$ r e^{i \theta} f@$ in the polar form.
 *
 * @param r Modulus of the complex number.
 * @param theta Argument of the complex number.
 * @returns Value of the complex number.
 */
std::complex<double> eval_polar(double r, double theta) {
  return r * (std::cos(theta) + M_I * std::sin(theta));
}

}  // trv::maths::

namespace runtime {

/// //////////////////////////////////////////////////////////////////////
/// Arrays
/// //////////////////////////////////////////////////////////////////////

/**
 * Exception raised when an extrapolation error occurs.
 *
 */
class ExtrapError: public std::runtime_error {
 public:
  std::string err_mesg;

  ExtrapError(const char* fmt_string, ...): std::runtime_error(
    "Extrapolation error."  // default error message; essential
  ) {
    std::va_list args;
    char err_mesg_buf[4096];

    va_start(args, fmt_string);
    std::vsprintf(err_mesg_buf, fmt_string, args);
    va_end(args);

    this->err_mesg = std::string(err_mesg_buf);
  }

  virtual const char* what() const noexcept {
    return err_mesg.c_str();
  }
};

}  // trv::runtime::

namespace ops {

/**
 * Extrapolate sample series exponentially (i.e. log-linearly).
 *
 * @param[in] a Samples series.
 * @param[in] N Sample number.
 * @param[in] N_ext Extrapolation number on either end.
 * @param[out] a_ext Extrapolated sample series.
 */
void extrap_loglin(double* a, int N, int N_ext, double* a_ext) {
  double dlna_left = std::log(a[1] / a[0]);
  if (std::isnan(dlna_left)) {
    throw trv::runtime::ExtrapError(
      "[ERRO] Sign change or zero detected in log-linear extrapolation "
      "at the lower end."
    );
  }

  double dlna_right = std::log(a[N - 1] / a[N - 2]);
  if (std::isnan(dlna_left)) {
    throw trv::runtime::ExtrapError(
      "[ERRO] Sign change or zero detected in log-linear extrapolation "
      "at the upper end."
    );
  }

  /// Extrapolate the lower end.
  if (a[0] > 0) {
    for (int i = 0; i < N_ext; i++) {
      a_ext[i] = a[0] * std::exp((i - N_ext) * dlna_left);
    }
  } else {
    for (int i = 0; i < N_ext; i++) {
      a_ext[i] = 0.;
    }
  }

  /// Fill in the middle part.
  for (int i = N_ext; i < N_ext + N; i++) {
    a_ext[i] = a[i - N_ext];
  }

  /// Extrapolate the upper end.
  if (a[N - 1] > 0) {
    for (int i = N_ext + N; i < N + 2*N_ext; i++) {
      a_ext[i] = a[N - 1] * std::exp((i - N_ext - (N - 1)) * dlna_right);
    }
  } else {
    for (int i = N_ext + N; i < N + 2*N_ext; i++) {
      a_ext[i] = 0.;
    }
  }
}

/**
 * Extrapolate sample bi-series exponentially (i.e. log-bilinearly).
 *
 * @param[in] a Samples bi-series.
 * @param[in] N Sample number (in both dimensions).
 * @param[in] N_ext Extrapolation number on either end
 *                  (in both dimensions).
 * @param[out] a_ext Extrapolated sample bi-series.
 */
void extrap2d_logbilin(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
) {
  double dlna_left, dlna_right;
  for (int i = N_ext; i < N_ext + N; i++) {
    dlna_left = std::log(a[i - N_ext][1] / a[i - N_ext][0]);
    if (std::isnan(dlna_left)) {
      throw trv::runtime::ExtrapError(
        "[ERRO] Sign change or zero detected in log-linear extrapolation "
        "at the left end."
      );
    }

    dlna_right = std::log(a[i - N_ext][N - 1] / a[i - N_ext][N - 2]);
    if (std::isnan(dlna_right)) {
      throw trv::runtime::ExtrapError(
        "[ERRO] Sign change or zero detected in log-linear extrapolation "
        "at the right end."
      );
    }

    /// Extrapolate the middle left edge.
    for (int j = 0; j < N_ext; j++) {
      if (a[i - N_ext][0] > 0) {
        a_ext[i][j] = a[i - N_ext][0] * std::exp((j - N_ext) * dlna_left);
      } else {
        a_ext[i][j] = 0.;
      }
    }

    /// Fill in the central part.
    for (int j = N_ext; j < N_ext + N; j++) {
      a_ext[i][j] = a[i - N_ext][j - N_ext];
    }

    /// Extrapolate the middle right edge.
    for (int j = N_ext + N; j < N + 2*N_ext; j++) {
      if (a[i - N_ext][N - 1] > 0) {
        a_ext[i][j] = a[i - N_ext][N - 1] * std::exp(
          (j - N_ext - (N - 1)) * dlna_right
        );
      } else {
        a_ext[i][j] = 0.;
      }
    }
  }

  double dlna_up, dlna_down;
  for (int j = 0; j < N + 2*N_ext; j++) {
    dlna_up = std::log(a_ext[N_ext + 1][j] / a_ext[N_ext][j]);
    if (std::isnan(dlna_up)) {
      throw trv::runtime::ExtrapError(
        "[ERRO] Sign change or zero detected in log-linear extrapolation "
        "at the top end."
      );
    }

    dlna_down = std::log(a_ext[N_ext + N - 1][j] / a_ext[N_ext + N - 2][j]);
    if (std::isnan(dlna_down)) {
      throw trv::runtime::ExtrapError(
        "[ERRO] Sign change or zero detected in log-linear extrapolation "
        "at the bottom end."
      );
    }

    /// Extrapolate the upper edge.
    for (int i = 0; i < N_ext; i++) {
      if (a_ext[N_ext][j] > 0) {
        a_ext[i][j] = a_ext[N_ext][j] * std::exp((i - N_ext) * dlna_up);
      } else {
        a_ext[i][j] = 0.;
      }
    }

    /// Extrapolate the lower edge.
    for (int i = N_ext + N; i < N + 2*N_ext; i++) {
      if (a_ext[N_ext + N - 1][j] > 0) {
        a_ext[i][j] = a_ext[N_ext + N - 1][j] * std::exp(
          (i - N_ext - (N - 1)) * dlna_down
        );
      } else {
        a_ext[i][j] = 0.;
      }
    }
  }
}

/**
 * Extrapolate sample bi-series linearly (i.e. bilinearly).
 *
 * @param[in] a Samples bi-series.
 * @param[in] N Sample number (in both dimensions).
 * @param[in] N_ext Extrapolation number on either end
 *                  (in both dimensions).
 * @param[out] a_ext Extrapolated sample bi-series.
 */
void extrap2d_bilin(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
) {
  double da_left, da_right;
  for (int i = N_ext; i < N + N_ext; i++) {
    da_left = a[i - N_ext][1] -  a[i - N_ext][0];
    da_right = a[i - N_ext][N - 1] - a[i - N_ext][N - 2];

    /// Extrapolate the middle left edge.
    for (int j = 0; j < N_ext; j++) {
      a_ext[i][j] = a[i - N_ext][0] + (j - N_ext) * da_left;
    }

    /// Fill in the central part.
    for (int j = N_ext; j < N_ext + N; j++) {
      a_ext[i][j] = a[i - N_ext][j - N_ext];
    }

    /// Extrapolate the middle right edge.
    for (int j = N_ext + N; j < N + 2*N_ext; j++) {
      a_ext[i][j] = a[i - N_ext][N-1] + (j - N_ext - (N - 1)) * da_right;
    }
  }

  double da_up, da_down;
  for (int j = 0; j < N + 2*N_ext; j++) {
    da_up = a_ext[N_ext + 1][j] - a_ext[N_ext][j];
    da_down = a_ext[N_ext + N - 1][j] - a_ext[N_ext + N - 2][j];

    /// Extrapolate the upper edge.
    for (int i =  0; i < N_ext; i++) {
      a_ext[i][j] = a_ext[N_ext][j] + (i - N_ext) * da_up;
    }

    /// Extrapolate the lower edge.
    for (int i =  N_ext + N; i < 2*N_ext + N; i++) {
      a_ext[i][j] = a_ext[N_ext + N - 1][j] + (i - N_ext - (N - 1)) * da_down;
    }
  }
}

/**
 * Extrapolate sample bi-series by zero padding.
 *
 * @param[in] a Samples bi-series.
 * @param[in] N Sample number (in both dimensions).
 * @param[in] N_ext Extrapolation number on either end
 *                  (in both dimensions).
 * @param[out] a_ext Extrapolated sample bi-series.
 */
void extrap2d_bizeros(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
) {
  for (int i = N_ext; i < N + N_ext; i++) {
    for (int j = 0; j < N_ext; j++) {
      a_ext[i][j] = 0.;
    }
    for (int j = N_ext; j < N_ext + N; j++) {
      a_ext[i][j] = a[i - N_ext][j - N_ext];
    }
    for (int j = N_ext + N; j < N + 2*N_ext; j++) {
      a_ext[i][j] = 0.;
    }
  }

  for (int j = 0; j < N + 2*N_ext; j++) {
    for (int i = 0; i < N_ext; i++) {
      a_ext[i][j] = 0.;
    }
    for (int i = N_ext + N; i < N + 2*N_ext; i++) {
      a_ext[i][j] = 0.;
    }
  }
}

}  // trv::ops::

namespace scheme {

/// //////////////////////////////////////////////////////////////////////
/// Misc
/// //////////////////////////////////////////////////////////////////////

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
  static void set_kbin(trv::scheme::ParameterSet& params, double* kbin_out) {
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
  static void set_rbin(trv::scheme::ParameterSet& params, double* rbin_out) {
    if (params.binning == "custom") {
      /// RFE: Insert customised binning code here.
    } else
    if (params.binning == "logpad") {
      /// CAVEAT: Discretionary choices.
      int nbin_custom = 5;

      for (int ibin = 0; ibin < nbin_custom; ibin++) {
        rbin_out[ibin] = 10. * ibin;
      }

      double rmin = 10. * nbin_custom;

      double dlnr = (std::log(params.rmax) - std::log(rmin))
        / double((params.num_rbin - nbin_custom) - 1);

      for (int ibin = nbin_custom; ibin < params.num_rbin; ibin++) {
        rbin_out[ibin] = rmin * std::exp(dlnr * (ibin - nbin_custom));
      }
    } else
    if (params.binning == "linpad") {
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
    } else
    if (params.binning == "log")  {
      double rmin;
      if (params.rmin == 0.) {
        /// CAVEAT: Discretionary choice.
        rmin = 1.;
      } else {
        rmin = params.rmin;
      }

      double dlnr = (std::log(params.rmax) - std::log(rmin))
        / double(params.num_rbin - 1);

      for (int ibin = 0; ibin < params.num_rbin; ibin++) {
        rbin_out[ibin] = rmin * std::exp(dlnr * ibin);
      }
    } else {  // by default ``params.binning == "lin"``
      double dr = (params.rmax - params.rmin) / double(params.num_rbin - 1);

      for (int ibin = 0; ibin < params.num_rbin; ibin++) {
        rbin_out[ibin] = params.rmin + dr * ibin;
      }
    }
  }
};

}  // trv::scheme::

}  // trv::

#endif  // TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_
