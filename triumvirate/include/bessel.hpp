#ifndef TRIUMVIRATE_INCLUDE_BESSEL_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_BESSEL_HPP_INCLUDED_

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>

/**
 * Interpolated spherical Bessel function @f$ j_\ell(x) @f$
 * of the first kind.
 *
 */
class SphericalBesselCalculator {
 public:
  /**
   * Initialise the interpolated function.
   *
   * @param ell Order @f$ \ell @f$ of the spherical Bessel function.
   */
  SphericalBesselCalculator(const int ell) {
    /// Set up sampling range and number.
    /// CAVEAT: Discretionary choices.
    double xmin = 0.;
    double xmax = 10000.;
    int nsample = 1000000;

    /// Initialise and evaluate sample points.
    double dx = (xmax - xmin) / (nsample - 1);

    double* x = new double[nsample];
    double* j_ell = new double[nsample];
    for (int i = 0; i < nsample; i++) {
      x[i] = xmin + dx * double(i);  // FIXME: remove double conversion
      j_ell[i] = gsl_sf_bessel_jl(ell, x[i]);
    }

    /// Initialise the interpolator.
    this->spline = gsl_spline_alloc(gsl_interp_cspline, nsample);
      // use cubic spline
    this->accel = gsl_interp_accel_alloc();
      // store state variables for faster lookup

    gsl_spline_init(this->spline, x, j_ell, nsample);

    /// Delete sample points.
    delete[] x; x = NULL;
    delete[] j_ell; j_ell = NULL;
  }

  /**
   * Destruct the interpolated function.
   */
  ~SphericalBesselCalculator() {
    if (this->accel != NULL) {
      gsl_interp_accel_free(this->accel); this->accel = NULL;
    }

    if (this->spline != NULL) {
      gsl_spline_free(this->spline); this->spline = NULL;
    }
  }

  /**
   * Evaluate the interpolated function.
   *
   * @param x Argument of the function.
   * @returns Value of the function.
   */
  double eval(double x) {
    return gsl_spline_eval(this->spline, x, this->accel);
  }

 private:
  gsl_spline* spline;
  gsl_interp_accel* accel;
};

#endif
