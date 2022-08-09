/**
 * @file twopt.hpp
 * @brief Two-point correlator computations.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_

#include "monitor.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "particles.hpp"
#include "field.hpp"

using namespace trv::obj;

const double EPS_COUPLING_2PT = 1.e-10;  /**< zero-tolerance for two-point
                                              coupling coefficients */
                                         // CAVEAT: discretionary choice

namespace trv {
namespace algo {

/**
 * Print out the header to the measurement output file.
 *
 * @param save_fileptr Saved file pointer.
 * @param params Parameter set.
 * @param catalogue_data (Data-source) catalogue.
 * @param catalogue_rand (Random-source) catalogue.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (only printed out
 *                 if non-zero).
 * @param space Either 'config'(-uration) space or 'fourier' space.
 *
 * @overload
 */
void print_2pt_meas_file_header(
  std::FILE* save_fileptr, trv::scheme::ParameterSet& params,
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  float norm, float norm_alt, std::string space
) {
  std::fprintf(
    save_fileptr,
    "# Data catalogue: %d particles of total weight %.3f\n",
    catalogue_data.ntotal, catalogue_data.wtotal
  );
  std::fprintf(
    save_fileptr,
    "# Random catalogue: %d particles of total weight %.3f\n",
    catalogue_rand.ntotal, catalogue_rand.wtotal
  );
  std::fprintf(
    save_fileptr,
    "# Box size: (%.3f, %.3f, %.3f)\n",
    params.boxsize[0], params.boxsize[1], params.boxsize[2]
  );
  std::fprintf(
    save_fileptr,
    "# Mesh number: (%d, %d, %d)\n",
    params.ngrid[0], params.ngrid[1], params.ngrid[2]
  );
  std::fprintf(
    save_fileptr,
    "# Particle extents: ([%.3f, %.3f], [%.3f, %.3f], [%.3f, %.3f])\n",
    catalogue_data.pos_min[0], catalogue_data.pos_max[0],
    catalogue_data.pos_min[1], catalogue_data.pos_max[1],
    catalogue_data.pos_min[2], catalogue_data.pos_max[2]
  );
  std::fprintf(
    save_fileptr,
    "# Alignment and assignment: %s, %s\n",
    params.alignment.c_str(), params.assignment.c_str()
  );
  std::fprintf(
    save_fileptr,
    "# Interlacing: %s\n", params.interlace.c_str()
  );
  std::fprintf(
    save_fileptr,
    "# Normalisation: %.9e, %s-based\n", norm, params.norm_convention.c_str()
  );
  if (norm_alt != 0.) {
    std::fprintf(save_fileptr, "# Normalisation alternative: %.9e\n", norm_alt);
  }
  if (space == "config") {
    std::fprintf(
      save_fileptr,
      "# [0] r_central, [1] r_eff, [2] npair, "
      "[3] xi%d.real, [4] xi%d.imag",
      params.ELL, params.ELL
    );
  } else
  if (space == "fourier") {
    std::fprintf(
      save_fileptr,
      "# [0] k_central, [1] k_eff, [2] nmode, "
      "[3] pk%d_raw.real, [4] pk%d_raw.imag, "
      "[5] pk_shot.real, [6] pk_shot.imag\n",
      params.ELL, params.ELL
    );
  } else {
    throw trv::runtime::InvalidParameter(
      "[%s ERRO] `space` must be either 'config' or 'fourier': %s\n",
      trv::runtime::show_timestamp().c_str(), space.c_str()
    );
  }
}

/**
 * Print out the header to the measurement output file.
 *
 * @param save_fileptr Saved file pointer.
 * @param params Parameter set.
 * @param catalogue Catalogue.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (only printed out
 *                 if non-zero).
 * @param space Either 'config'(-uration) space or 'fourier' space.
 */
void print_2pt_meas_file_header(
  std::FILE* save_fileptr,trv::scheme::ParameterSet& params,
  ParticleCatalogue& catalogue, float norm, float norm_alt, std::string space
) {
  std::fprintf(
    save_fileptr,
    "# Catalogue: %d particles of total weight %.3f\n",
    catalogue.ntotal, catalogue.wtotal
  );
  std::fprintf(
    save_fileptr,
    "# Box size: (%.3f, %.3f, %.3f)\n",
    params.boxsize[0], params.boxsize[1], params.boxsize[2]
  );
  std::fprintf(
    save_fileptr,
    "# Mesh number: (%d, %d, %d)\n",
    params.ngrid[0], params.ngrid[1], params.ngrid[2]
  );
  std::fprintf(
    save_fileptr,
    "# Particle extents: ([%.3f, %.3f], [%.3f, %.3f], [%.3f, %.3f])\n",
    catalogue.pos_min[0], catalogue.pos_max[0],
    catalogue.pos_min[1], catalogue.pos_max[1],
    catalogue.pos_min[2], catalogue.pos_max[2]
  );
  std::fprintf(
    save_fileptr,
    "# Alignment and assignment: %s, %s\n",
    params.alignment.c_str(), params.assignment.c_str()
  );
  std::fprintf(
    save_fileptr,
    "# Interlacing: %s\n", params.interlace.c_str()
  );
  std::fprintf(
    save_fileptr,
    "# Normalisation: %.9e, %s-based\n", norm, params.norm_convention.c_str()
  );
  if (norm_alt != 0.) {
    std::fprintf(save_fileptr, "# Normalisation alternative: %.9e\n", norm_alt);
  }
  if (space == "config") {
    std::fprintf(
      save_fileptr,
      "# [0] r_central, [1] r_eff, [2] npair, "
      "[3] xi%d.real, [4] xi%d.imag",
      params.ELL, params.ELL
    );
  } else
  if (space == "fourier") {
    std::fprintf(
      save_fileptr,
      "# [0] k_central, [1] k_eff, [2] nmode, "
      "[3] pk%d_raw.real, [4] pk%d_raw.imag, "
      "[5] pk_shot.real, [6] pk_shot.imag\n",
      params.ELL, params.ELL
    );
  } else {
    throw trv::runtime::InvalidParameter(
      "[%s ERRO] `space` must be either 'config' or 'fourier': %s\n",
      trv::runtime::show_timestamp().c_str(), space.c_str()
    );
  }
}

/**
 * Power spectrum measurements.
 *
 */
struct PowspecMeasurements {
  std::vector<double> kbin;  ///< central wavenumber in bins
  std::vector<double> keff;  ///< effective wavenumber in bins
  std::vector<int> nmode;  ///< number of contributing wavevector modes
  std::vector< std::complex<double> > pk_raw;  /**< power spectrum
                                                    raw measurements */
  std::vector< std::complex<double> > pk_shot;  ///< power spectrum shot noise
};

/**
 * Two-point correlation function measurements.
 *
 */
struct CorrfuncMeasurements {
  std::vector<double> rbin;  ///< central separation in bins
  std::vector<double> reff;  ///< effective separation in bins
  std::vector<int> npair;  ///< number of contributing separation pairs
  std::vector< std::complex<double> > xi;  /**< two-point correlation function
                                                measurements */
};

/**
 * Power spectrum window measurements.
 *
 */
struct PowspecWindowMeasurements {
  std::vector<double> kbin;  ///< central wavenumber in bins
  std::vector<double> keff;  ///< effective wavenumber in bins
  std::vector<int> nmode;  ///< number of contributing wavevector modes
  std::vector< std::complex<double> > pk;  /**< power spectrum
                                                window measurements */
};

/**
 * Two-point correlation function window measurements.
 *
 */
struct CorrfuncWindowMeasurements {
  std::vector<double> rbin;  ///< central separation in bins
  std::vector<double> reff;  ///< effective separation in bins
  std::vector<int> npair;  ///< number of contributing separation pairs
  std::vector< std::complex<double> > xi;  /**< two-point correlation function
                                                window measurements */
};

/**
 * Calculate mesh-based power spectrum normalisation.
 *
 * @param catalogue Particle catalogue.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @returns norm_factor Power spectrum normalisation factor.
 */
double calc_powspec_normalisation_from_mesh(
  ParticleCatalogue& catalogue,
  trv::scheme::ParameterSet& params,
  double alpha=1.
) {
  PseudoDensityField<ParticleCatalogue> catalogue_mesh(params);

  double norm_factor = catalogue_mesh._calc_wgt_sq_volume_norm(catalogue)
    / std::pow(alpha, 2);

  catalogue_mesh.finalise_density_field();

  return norm_factor;
}

/**
 * Calculate particle-based power spectrum normalisation.
 *
 * @param catalogue Particle catalogue.
 * @param alpha Alpha ratio.
 * @returns norm_factor Power spectrum normalisation factor.
 */
double calc_powspec_normalisation_from_particles(
  ParticleCatalogue& catalogue,
  double alpha=1.
) {
  double norm_factor = catalogue._calc_powspec_normalisation() / alpha;

  return norm_factor;
}

/// NOBUG: Standard naming convention is not always followed for
/// intermediary quantities in the functions below.

/**
 * Compute power spectrum from paired catalogues and
 * optionally save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param particles_rand (Random-source) particle container.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param kbin Wavenumber bins.
 * @param alpha Alpha ratio.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (default is 0.) which
 *                 may be printed out in the header of the
 *                 measurement output file.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns powspec_out Output power spectrum measurements.
 */
PowspecMeasurements compute_powspec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::scheme::ParameterSet& params,
  double* kbin,
  double alpha,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: power spectrum from "
      "data and random catalogues.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up input.
  int ell1 = params.ELL;

  /// Set up output.
  int* nmode_save = new int[params.num_kbin];
  double* k_save = new double[params.num_kbin];
  std::complex<double>* pk_save = new std::complex<double>[params.num_kbin];
  std::complex<double>* sn_save = new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    nmode_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute power spectrum.
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    PseudoDensityField<ParticleCatalogue> dn_LM(params);
    dn_LM.compute_ylm_wgtd_fluctuation(
      particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    Pseudo2ptStats<ParticleCatalogue> stats2pt(params);
    std::complex<double> sn_amp = stats2pt.calc_ylm_wgtd_shotnoise_for_powspec(
      particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, M_
    );

    /// Compute quantity equivalent to (-1)^m_1 δ_{m_1, -M} which, after
    /// being summed over m_1, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = (2*params.ELL + 1) * (2*ell1 + 1)
        * trv::maths::wigner_3j(ell1, 0, params.ELL, 0, 0, 0)
        * trv::maths::wigner_3j(ell1, 0, params.ELL, m1, 0, M_);

      if (std::fabs(coupling) < EPS_COUPLING_2PT) {continue;}

      stats2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
        dn_LM, dn_00, sn_amp, kbin, ell1, m1
      );

      for (int ibin = 0; ibin < params.num_kbin; ibin++) {
        pk_save[ibin] += coupling * stats2pt.pk[ibin];
        sn_save[ibin] += coupling * stats2pt.sn[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_kbin; ibin++) {
          nmode_save[ibin] = stats2pt.nmode[ibin];
          k_save[ibin] = stats2pt.k[ibin];
        }
      }
    }

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Power spectrum term at order M = %d computed.\n",
        trv::runtime::show_timestamp().c_str(),
        M_
      );
    }
  }

  /// Compute shot noise if the approach is particle-based.
  std::complex<double> sn_particle =
    particles_data._calc_powspec_shotnoise()
    + alpha * alpha * particles_rand._calc_powspec_shotnoise();

  /* * Output ************************************************************** */

  /// Fill in output struct.
  PowspecMeasurements powspec_out;
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    powspec_out.kbin.push_back(kbin[ibin]);
    powspec_out.keff.push_back(k_save[ibin]);
    powspec_out.nmode.push_back(nmode_save[ibin]);
    powspec_out.pk_raw.push_back(norm * pk_save[ibin]);
    if (params.shotnoise_convention == "mesh") {
      powspec_out.pk_shot.push_back(norm * sn_save[ibin]);
    } else
    if (params.shotnoise_convention == "particle") {
      powspec_out.pk_shot.push_back(norm * sn_particle);
    }
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    std::sprintf(
      save_filepath, "%s/pk%d%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    print_2pt_meas_file_header(
      save_fileptr, params, particles_data, particles_rand, norm, norm_alt,
      "fourier"
    );
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        powspec_out.kbin[ibin], powspec_out.keff[ibin],
        powspec_out.nmode[ibin],
        powspec_out.pk_raw[ibin].real(), powspec_out.pk_raw[ibin].imag(),
        powspec_out.pk_shot[ibin].real(), powspec_out.pk_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::runtime::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] nmode_save; delete[] k_save; delete[] pk_save; delete[] sn_save;

  return powspec_out;
}

/**
 * Compute two-point correlation function from paired catalogues and
 * optionally save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param particles_rand (Random-source) particle container.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbin Separation bins.
 * @param alpha Alpha ratio.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (default is 0.) which
 *                 may be printed out in the header of the
 *                 measurement output file.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns corrfunc_out Output two-point correlation function
 *                       measurements.
 */
CorrfuncMeasurements compute_corrfunc(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::scheme::ParameterSet& params,
  double* rbin,
  double alpha,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: two-point correlation function "
      "from data and random catalogues.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up input.
  int ell1 = params.ELL;

  /// Set up output.
  int* npair_save = new int[params.num_rbin];
  double* r_save = new double[params.num_rbin];
  std::complex<double>* xi_save = new std::complex<double>[params.num_rbin];
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    npair_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute 2PCF.
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    PseudoDensityField<ParticleCatalogue> dn_LM(params);
    dn_LM.compute_ylm_wgtd_fluctuation(
      particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    Pseudo2ptStats<ParticleCatalogue> stats2pt(params);
    std::complex<double> sn_amp = stats2pt.calc_ylm_wgtd_shotnoise_for_powspec(
      particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, M_
    );

    /// Compute quantity equivalent to (-1)^m_1 δ_{m_1, -M} which, after
    /// being summed over m_1, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = (2*params.ELL + 1) * (2*ell1 + 1)
        * trv::maths::wigner_3j(ell1, 0, params.ELL, 0, 0, 0)
        * trv::maths::wigner_3j(ell1, 0, params.ELL, m1, 0, M_);

      if (std::fabs(coupling) < EPS_COUPLING_2PT) {continue;}

      stats2pt.compute_ylm_wgtd_2pt_stats_in_config(
        dn_LM, dn_00, sn_amp, rbin, ell1, m1
      );

      for (int ibin = 0; ibin < params.num_rbin; ibin++) {
        xi_save[ibin] += coupling * stats2pt.xi[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_kbin; ibin++) {
          npair_save[ibin] = stats2pt.npair[ibin];
          r_save[ibin] = stats2pt.r[ibin];
        }
      }
    }

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Two-point correlation function term at "
        "order M = %d computed.\n",
        trv::runtime::show_timestamp().c_str(),
        M_
      );
    }
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  CorrfuncMeasurements corrfunc_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    corrfunc_out.rbin.push_back(rbin[ibin]);
    corrfunc_out.reff.push_back(r_save[ibin]);
    corrfunc_out.npair.push_back(npair_save[ibin]);
    corrfunc_out.xi.push_back(norm * xi_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    std::sprintf(
      save_filepath, "%s/xi%d%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    print_2pt_meas_file_header(
      save_fileptr, params, particles_data, particles_rand, norm, norm_alt,
      "config"
    );
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      std::fprintf(
        save_fileptr, "%.9e \t %.9e \t %d \t %.9e \t %.9e\n",
        corrfunc_out.rbin[ibin], corrfunc_out.reff[ibin],
        corrfunc_out.npair[ibin],
        corrfunc_out.xi[ibin].real(), corrfunc_out.xi[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::runtime::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] npair_save; delete[] r_save; delete[] xi_save;

  return corrfunc_out;
}

/**
 * Compute power spectrum in a periodic box and optionally
 * save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param params Parameter set.
 * @param kbin Wavenumber bins.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (default is 0.) which
 *                 may be printed out in the header of the
 *                 measurement output file.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns powspec_out Output power spectrum measurements.
 */
PowspecMeasurements compute_powspec_in_box(
  ParticleCatalogue& particles_data,
  trv::scheme::ParameterSet& params,
  double* kbin,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: power spectrum in a periodic box.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  int* nmode_save = new int[params.num_kbin];
  double* k_save = new double[params.num_kbin];
  std::complex<double>* pk_save = new std::complex<double>[params.num_kbin];
  std::complex<double>* sn_save = new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    nmode_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute power spectrum.
  PseudoDensityField<ParticleCatalogue> dn(params);
  dn.compute_unweighted_fluctuation_insitu(particles_data, params.volume);
  dn.fourier_transform();

  std::complex<double> sn_amp = double(particles_data.ntotal);

  Pseudo2ptStats<ParticleCatalogue> stats2pt(params);
  stats2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
    dn, dn, sn_amp, kbin, params.ELL, 0
  );

  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    nmode_save[ibin] = stats2pt.nmode[ibin];
    k_save[ibin] = stats2pt.k[ibin];
    pk_save[ibin] += double(2*params.ELL + 1) * stats2pt.pk[ibin];
    sn_save[ibin] += double(2*params.ELL + 1) * stats2pt.sn[ibin];
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Power spectrum terms computed.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute shot noise if the approach is particle-based.
  std::complex<double> sn_particle = particles_data._calc_powspec_shotnoise();

  /// Normalisation factor should redue to
  /// ``params.volume / particles_data.ntotal / particles_data.ntotal``
  /// in the case of uniform particle weights.

  /* * Output ************************************************************** */

  /// Fill in output struct.
  PowspecMeasurements powspec_out;
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    powspec_out.kbin.push_back(kbin[ibin]);
    powspec_out.keff.push_back(k_save[ibin]);
    powspec_out.nmode.push_back(nmode_save[ibin]);
    powspec_out.pk_raw.push_back(norm * pk_save[ibin]);
    if (params.shotnoise_convention == "mesh") {
      powspec_out.pk_shot.push_back(norm * sn_save[ibin]);
    } else
    if (params.shotnoise_convention == "particle") {
      powspec_out.pk_shot.push_back(norm * sn_particle);
    }
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    std::sprintf(
      save_filepath, "%s/pk%d%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    print_2pt_meas_file_header(
      save_fileptr, params, particles_data, norm, norm_alt, "fourier"
    );
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        powspec_out.kbin[ibin], powspec_out.keff[ibin],
        powspec_out.nmode[ibin],
        powspec_out.pk_raw[ibin].real(), powspec_out.pk_raw[ibin].imag(),
        powspec_out.pk_shot[ibin].real(), powspec_out.pk_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::runtime::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] nmode_save; delete[] k_save; delete[] pk_save; delete[] sn_save;

  return powspec_out;
}

/**
 * Compute two-point correlation function in a periodic box
 * and optionally save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param params Parameter set.
 * @param rbin Separation bins.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (default is 0.) which
 *                 may be printed out in the header of the
 *                 measurement output file.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns corrfunc_out Output two-point correlation function
 *                       measurements.
 */
CorrfuncMeasurements compute_corrfunc_in_box(
  ParticleCatalogue& particles_data,
  trv::scheme::ParameterSet& params,
  double* rbin,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: two-point correlation function "
      "in a periodic box.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  int* npair_save = new int[params.num_rbin];
  double* r_save = new double[params.num_kbin];
  std::complex<double>* xi_save = new std::complex<double>[params.num_rbin];
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    npair_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute 2PCF.
  PseudoDensityField<ParticleCatalogue> dn(params);
  dn.compute_unweighted_fluctuation_insitu(particles_data, params.volume);
  dn.fourier_transform();

  std::complex<double> sn_amp = double(particles_data.ntotal);

  Pseudo2ptStats<ParticleCatalogue> stats2pt(params);
  stats2pt.compute_ylm_wgtd_2pt_stats_in_config(
    dn, dn, sn_amp, rbin, params.ELL, 0
  );

  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    npair_save[ibin] = stats2pt.npair[ibin];
    r_save[ibin] = stats2pt.r[ibin];
    xi_save[ibin] += double(2*params.ELL + 1) * stats2pt.xi[ibin];
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Two-point correlation function terms computed.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Normalisation factor should redue to
  /// ``params.volume / particles_data.ntotal / particles_data.ntotal``
  /// in the case of uniform particle weights.

  /* * Output ************************************************************** */

  /// Fill in output struct.
  CorrfuncMeasurements corrfunc_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    corrfunc_out.rbin.push_back(rbin[ibin]);
    corrfunc_out.reff.push_back(r_save[ibin]);
    corrfunc_out.npair.push_back(npair_save[ibin]);
    corrfunc_out.xi.push_back(norm * xi_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    std::sprintf(
      save_filepath, "%s/xi%d%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    print_2pt_meas_file_header(
      save_fileptr, params, particles_data, norm, norm_alt, "config"
    );
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      std::fprintf(
        save_fileptr, "%.9e \t %.9e \t %d \t %.9e \t %.9e\n",
        corrfunc_out.rbin[ibin], corrfunc_out.reff[ibin],
        corrfunc_out.npair[ibin],
        corrfunc_out.xi[ibin].real(), corrfunc_out.xi[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::runtime::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] npair_save; delete[] r_save; delete[] xi_save;

  return corrfunc_out;
}

#ifdef TRV_USE_DISABLED_CODE
/**
 * Compute power spectrum window from a random catalogue and
 * optionally save the results.
 *
 * @param particles_rand (Random-source) particle container.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param kbin Wavenumber bins.
 * @param alpha Alpha ratio.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (default is 0.) which
 *                 may be printed out in the header of the
 *                 measurement output file.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns powwin_out Output power spectrum window measurements.
 */
PowspecWindowMeasurements compute_powspec_window(
  ParticleCatalogue& particles_rand,
  LineOfSight* los_rand,
  trv::scheme::ParameterSet& params,
  double* kbin,
  double alpha,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: power spectrum window "
      "from random catalogue.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up output.
  int* nmode_save = new int[params.num_kbin];
  double* k_save = new double[params.num_kbin];
  std::complex<double>* pk_save = new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    nmode_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute power spectrum window.
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_density(particles_rand, los_rand, alpha, 0, 0);
  dn_00.fourier_transform();

  Pseudo2ptStats<ParticleCatalogue> stats2pt(params);
  std::complex<double> sn_amp = stats2pt.calc_ylm_wgtd_shotnoise_for_powspec(
    particles_rand, los_rand, alpha, params.ELL, 0
  );

  stats2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
    dn_00, dn_00, sn_amp, kbin, params.ELL, 0
  );

  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    nmode_save[ibin] = stats2pt.nmode[ibin];
    k_save[ibin] = stats2pt.k[ibin];
    pk_save[ibin] += stats2pt.pk[ibin];
  }

  /// Renormalise the normalisation factor to be dimensionless.
  /// CAVEAT: Discretionary choice to eliminate the physical volume
  /// dimension by physical volume.
  norm /= params.volume;

  /* * Output ************************************************************** */

  /// Fill in output struct.
  PowspecWindowMeasurements powwin_out;
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    powwin_out.kbin.push_back(kbin[ibin]);
    powwin_out.keff.push_back(k_save[ibin]);
    powwin_out.nmode.push_back(nmode_save[ibin]);
    powwin_out.pk.push_back(norm * pk_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    std::sprintf(
      save_filepath, "%s/pk%d_win%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    print_2pt_meas_file_header(
      save_fileptr, params, particles_rand, norm, norm_alt, "fourier"
    );
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      std::fprintf(
        save_fileptr, "%.9e \t %.9e \t %d \t %.9e \t %.9e\n",
        powwin_out.kbin[ibin], powwin_out.keff[ibin],
        powwin_out.nmode[ibin],
        powwin_out.pk[ibin].real(), powwin_out.pk[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::runtime::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] nmode_save; delete[] k_save; delete[] pk_save;

  return powwin_out;
}
#endif  // TRV_USE_DISABLED_CODE

/**
 * Compute two-point correlation function window from a random catalogue
 * and optionally save the results.
 *
 * @param particles_rand (Random-source) particle container.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param kbin Wavenumber bins.
 * @param alpha Alpha ratio.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (default is 0.) which
 *                 may be printed out in the header of the
 *                 measurement output file.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns corrfwin_out Output two-point correlation function
 *                       window measurements.
 */
CorrfuncWindowMeasurements compute_corrfunc_window(
  ParticleCatalogue& particles_rand,
  LineOfSight* los_rand,
  trv::scheme::ParameterSet& params,
  double* rbin,
  double alpha,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: two-point correlation function window "
      "from random catalogue.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up input.
  int ell1 = params.ELL;

  /// Set up output.
  int* npair_save = new int[params.num_rbin];
  double* r_save = new double[params.num_rbin];
  std::complex<double>* xi_save = new std::complex<double>[params.num_rbin];
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    npair_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute 2PCF window.
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_density(particles_rand, los_rand, alpha, 0, 0);
  dn_00.fourier_transform();

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    PseudoDensityField<ParticleCatalogue> dn_LM(params);
    dn_LM.compute_ylm_wgtd_density(
      particles_rand, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    Pseudo2ptStats<ParticleCatalogue> stats2pt(params);
    std::complex<double> sn_amp = stats2pt.calc_ylm_wgtd_shotnoise_for_powspec(
      particles_rand, los_rand, alpha, params.ELL, M_
    );

    /// Calculate equivalent to (-1)^m_1 δ_{m_1, -M} which, after being
    /// summed over m_1, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = (2*params.ELL + 1) * (2*ell1 + 1)
          * trv::maths::wigner_3j(ell1, 0, params.ELL, 0, 0, 0)
          * trv::maths::wigner_3j(ell1, 0, params.ELL, m1, 0, M_);

      if (std::fabs(coupling) < EPS_COUPLING_2PT) {continue;}

      stats2pt.compute_ylm_wgtd_2pt_stats_in_config(
        dn_LM, dn_00, sn_amp, rbin, ell1, m1
      );

      for (int ibin = 0; ibin < params.num_rbin; ibin++) {
        xi_save[ibin] += coupling * stats2pt.xi[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_kbin; ibin++) {
          npair_save[ibin] = stats2pt.npair[ibin];
          r_save[ibin] = stats2pt.r[ibin];
        }
      }
    }

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Two-point correlation function window term "
        "at order M = %d computed.\n",
        trv::runtime::show_timestamp().c_str(),
        M_
      );
    }
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  CorrfuncWindowMeasurements corrfwin_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    corrfwin_out.rbin.push_back(rbin[ibin]);
    corrfwin_out.reff.push_back(r_save[ibin]);
    corrfwin_out.npair.push_back(npair_save[ibin]);
    corrfwin_out.xi.push_back(norm * xi_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    std::sprintf(
      save_filepath, "%s/xi%d_win%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    print_2pt_meas_file_header(
      save_fileptr, params, particles_rand, norm, norm_alt, "config"
    );
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      std::fprintf(
        save_fileptr, "%.9e \t %.9e \t %d \t %.9e \t %.9e\n",
        corrfwin_out.rbin[ibin], corrfwin_out.reff[ibin],
        corrfwin_out.npair[ibin],
        corrfwin_out.xi[ibin].real(), corrfwin_out.xi[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::runtime::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] npair_save; delete[] r_save; delete[] xi_save;

  return corrfwin_out;
}

}  // trv::algo::
}  // trv::

#endif  // TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_
