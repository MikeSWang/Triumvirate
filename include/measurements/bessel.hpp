#ifndef TRIUM_BESSEL_H_INCLUDED_
#define TRIUM_BESSEL_H_INCLUDED_

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
  SphericalBesselCalculator(int ell) {
    /// Set up sampling range and number.
    double xmin = 0.;
    double xmax = 10000.;  // NOTE: discretionary choice
    int nsample = 1000000;  // NOTE: discretionary choice

    /// Initialise and evaluate sample points.
    double* x = new double[nsample];
    double* j_ell = new double[nsample];

    double dx = (xmax - xmin) / (nsample - 1);
    for (int i = 0; i < nsample; i++) {
      x[i] = xmin + dx * double(i);
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
