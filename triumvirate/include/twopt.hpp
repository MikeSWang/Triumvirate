/**
 * @file twopt.hpp
 * @brief Two-point correlator computations.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_

#include "monitor.hpp"
#include "parameters.hpp"
#include "maths.hpp"
#include "particles.hpp"
#include "io.hpp"
#include "dataobjs.hpp"
#include "field.hpp"

using namespace trv::maths;

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
 * @param pk_sn_particle Shot noise calculated based on particles (only
 *                       printed out if non-zero).
 * @param space Either 'config'(-uration) space or 'fourier' space.
 *
 * @overload
 */
void print_2pt_meas_file_header(
  std::FILE* save_fileptr, trv::ParameterSet& params,
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  float norm, float norm_alt, std::complex<double> pk_sn_particle
) {
  trv::print_premeasurement_info(
    save_fileptr, params, catalogue_data, catalogue_rand
  );
  std::fprintf(
    save_fileptr,
    "# Normalisation (%s-based): %.9e\n",
    params.norm_convention.c_str(), norm
  );
  if (norm_alt != 0.) {
    std::fprintf(
      save_fileptr,
      "# Normalisation (alternative): %.9e\n", norm_alt
    );
  }
  std::fprintf(
    save_fileptr,
    "# Shot noise (%s-based): see data below\n",
    params.shotnoise_convention.c_str()
  );
  if (pk_sn_particle != 0.) {
    std::fprintf(
      save_fileptr,
      "# Shot noise (particle-based): %.9e\n",
      pk_sn_particle.real()
    );
  }
  if (params.space == "config") {
    std::fprintf(
      save_fileptr,
      "# [0] r_central, [1] r_eff, [2] npairs, "
      "[3] xi%d.real, [4] xi%d.imag",
      params.ELL, params.ELL
    );
  } else
  if (params.space == "fourier") {
    std::fprintf(
      save_fileptr,
      "# [0] k_central, [1] k_eff, [2] nmodes, "
      "[3] pk%d_raw.real, [4] pk%d_raw.imag, "
      "[5] pk_shot.real, [6] pk_shot.imag\n",
      params.ELL, params.ELL
    );
  } else {
    throw trv::sys::InvalidParameter(
      "[%s ERRO] `space` must be either 'config' or 'fourier': %s\n",
      trv::sys::show_timestamp().c_str(), params.space.c_str()
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
 * @param pk_sn_particle Shot noise calculated based on particles (only
 *                       printed out if non-zero).
 * @param space Either 'config'(-uration) space or 'fourier' space.
 */
void print_2pt_meas_file_header(
  std::FILE* save_fileptr,trv::ParameterSet& params,
  ParticleCatalogue& catalogue, float norm, float norm_alt,
  std::complex<double> pk_sn_particle
) {
  trv::print_premeasurement_info(save_fileptr, params, catalogue);
  std::fprintf(
    save_fileptr,
    "# Normalisation (%s-based): %.9e\n",
    params.norm_convention.c_str(), norm
  );
  if (norm_alt != 0.) {
    std::fprintf(
      save_fileptr,
      "# Normalisation (alternative): %.9e\n", norm_alt
    );
  }
  std::fprintf(
    save_fileptr,
    "# Shot noise (%s-based): see data below\n",
    params.shotnoise_convention.c_str()
  );
  if (pk_sn_particle != 0.) {
    std::fprintf(
      save_fileptr,
      "# Shot noise (particle-based): %.9e\n",
      pk_sn_particle.real()
    );
  }
  if (params.space == "config") {
    std::fprintf(
      save_fileptr,
      "# [0] r_central, [1] r_eff, [2] npairs, "
      "[3] xi%d.real, [4] xi%d.imag",
      params.ELL, params.ELL
    );
  } else
  if (space == "fourier") {
    std::fprintf(
      save_fileptr,
      "# [0] k_central, [1] k_eff, [2] nmodes, "
      "[3] pk%d_raw.real, [4] pk%d_raw.imag, "
      "[5] pk_shot.real, [6] pk_shot.imag\n",
      params.ELL, params.ELL
    );
  } else {
    throw trv::sys::InvalidParameter(
      "[%s ERRO] `space` must be either 'config' or 'fourier': %s\n",
      trv::sys::show_timestamp().c_str(), space.c_str()
    );
  }
}

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
  trv::ParameterSet& params,
  double alpha=1.
) {
  MeshField catalogue_mesh(params);

  double norm_factor = catalogue_mesh.calc_grid_based_powlaw_norm(catalogue, 2)
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
  if (catalogue.pdata == nullptr) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle data are uninitialised.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  double vol_eff_inv = 0.;  // I_2
  for (int pid = 0; pid < catalogue.ntotal; pid++) {
    vol_eff_inv += catalogue[pid].nz
      * catalogue[pid].ws * std::pow(catalogue[pid].wc, 2);
  }

  if (vol_eff_inv == 0.) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle 'nz' values appear to be all zeros. "
        "Check the input catalogue contains valid 'nz' field.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  double norm_factor = 1. / vol_eff_inv / alpha;  // I_2^(-1)

  return norm_factor;
}

/**
 * Calculate particle-based power spectrum shot noise.
 *
 * @param catalogue Particle catalogue.
 * @param alpha Alpha ratio.
 * @returns shotnoise Power spectrum normalisation factor.
 */
double calc_powspec_shotnoise_from_particles(
  ParticleCatalogue& catalogue,
  double alpha=1.
) {
  if (catalogue.pdata == nullptr) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle data are uninitialised.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  double shotnoise = 0.;
  for (int pid = 0; pid < catalogue.ntotal; pid++) {
    shotnoise += std::pow(catalogue[pid].ws, 2)
      * std::pow(catalogue[pid].wc, 2);
  }

  return shotnoise;
}

/// STYLE: Standard naming convention is not always followed for
/// intermediary quantities in the functions below.

/**
 * Compute power spectrum from paired catalogues and
 * optionally save the results.
 *
 * @param particles_data (Data-source) particle catalogue.
 * @param particles_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param kbins Wavenumber bins.
 * @param alpha Alpha ratio.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (default is 0.) which
 *                 may be printed out in the header of the
 *                 measurement output file.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns powspec_out Output power spectrum measurements.
 */
trv::PowspecMeasurements compute_powspec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params,
  trv::Binning& kbins,
  double alpha,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: power spectrum from "
      "data and random catalogues.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up input.
  int ell1 = params.ELL;

  /// Set up output.
  int* nmodes_save = new int[params.num_bins];
  double* k_save = new double[params.num_bins];
  std::complex<double>* pk_save = new std::complex<double>[params.num_bins];
  std::complex<double>* sn_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute power spectrum.
  MeshField dn_00(params);
  dn_00.compute_ylm_wgtd_field(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params);
    dn_LM.compute_ylm_wgtd_field(
      particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    FieldStats stats2pt(params);
    std::complex<double> sn_amp = stats2pt.calc_ylm_wgtd_shotnoise_amp_for_powspec(
      particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, M_
    );

    /// Compute quantity equivalent to (-1)^m_1 δ_{m_1, -M} which, after
    /// being summed over m_1, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = (2*params.ELL + 1) * (2*ell1 + 1)
        * wigner_3j(ell1, 0, params.ELL, 0, 0, 0)
        * wigner_3j(ell1, 0, params.ELL, m1, 0, M_);

      if (std::fabs(coupling) < EPS_COUPLING_2PT) {continue;}

      stats2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
        dn_LM, dn_00, sn_amp, ell1, m1, kbins
      );

      for (int ibin = 0; ibin < params.num_bins; ibin++) {
        pk_save[ibin] += coupling * stats2pt.pk[ibin];
        sn_save[ibin] += coupling * stats2pt.sn[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_bins; ibin++) {
          nmodes_save[ibin] = stats2pt.nmodes[ibin];
          k_save[ibin] = stats2pt.k[ibin];
        }
      }
    }

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Power spectrum term at order M = %d computed.\n",
        trv::sys::show_timestamp().c_str(),
        M_
      );
    }
  }

  /// Compute shot noise if the approach is particle-based.
  std::complex<double> sn_particle =
    calc_powspec_shotnoise_from_particles(particles_data)
    + alpha * alpha * calc_powspec_shotnoise_from_particles(particles_rand);

  /* * Output ************************************************************** */

  /// Fill in output struct.
  trv::PowspecMeasurements powspec_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    powspec_out.kbin.push_back(kbins.bin_centres[ibin]);
    powspec_out.keff.push_back(k_save[ibin]);
    powspec_out.nmodes.push_back(nmodes_save[ibin]);
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
      save_fileptr, params, particles_data, particles_rand,
      norm, norm_alt, norm * sn_particle,
      "fourier"
    );
    for (int ibin = 0; ibin < params.num_bins; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        powspec_out.kbin[ibin], powspec_out.keff[ibin],
        powspec_out.nmodes[ibin],
        powspec_out.pk_raw[ibin].real(), powspec_out.pk_raw[ibin].imag(),
        powspec_out.pk_shot[ibin].real(), powspec_out.pk_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::sys::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] nmodes_save; delete[] k_save; delete[] pk_save; delete[] sn_save;

  return powspec_out;
}

/**
 * Compute two-point correlation function from paired catalogues and
 * optionally save the results.
 *
 * @param particles_data (Data-source) particle catalogue.
 * @param particles_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbins Separation bins.
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
trv::TwoPCFMeasurements compute_corrfunc(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params,
  trv::Binning& rbins,
  double alpha,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: two-point correlation function "
      "from data and random catalogues.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up input.
  int ell1 = params.ELL;

  /// Set up output.
  int* npairs_save = new int[params.num_bins];
  double* r_save = new double[params.num_bins];
  std::complex<double>* xi_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute 2PCF.
  MeshField dn_00(params);
  dn_00.compute_ylm_wgtd_field(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params);
    dn_LM.compute_ylm_wgtd_field(
      particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    FieldStats stats2pt(params);
    std::complex<double> sn_amp = stats2pt.calc_ylm_wgtd_shotnoise_amp_for_powspec(
      particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, M_
    );

    /// Compute quantity equivalent to (-1)^m_1 δ_{m_1, -M} which, after
    /// being summed over m_1, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = (2*params.ELL + 1) * (2*ell1 + 1)
        * wigner_3j(ell1, 0, params.ELL, 0, 0, 0)
        * wigner_3j(ell1, 0, params.ELL, m1, 0, M_);

      if (std::fabs(coupling) < EPS_COUPLING_2PT) {continue;}

      stats2pt.compute_ylm_wgtd_2pt_stats_in_config(
        dn_LM, dn_00, sn_amp, ell1, m1, rbins
      );

      for (int ibin = 0; ibin < params.num_bins; ibin++) {
        xi_save[ibin] += coupling * stats2pt.xi[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_bins; ibin++) {
          npairs_save[ibin] = stats2pt.npairs[ibin];
          r_save[ibin] = stats2pt.r[ibin];
        }
      }
    }

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Two-point correlation function term at "
        "order M = %d computed.\n",
        trv::sys::show_timestamp().c_str(),
        M_
      );
    }
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  trv::TwoPCFMeasurements corrfunc_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    corrfunc_out.rbin.push_back(rbins.bin_centres[ibin]);
    corrfunc_out.reff.push_back(r_save[ibin]);
    corrfunc_out.npairs.push_back(npairs_save[ibin]);
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
      save_fileptr, params, particles_data, particles_rand,
      norm, norm_alt, 0.,
      "config"
    );
    for (int ibin = 0; ibin < params.num_bins; ibin++) {
      std::fprintf(
        save_fileptr, "%.9e \t %.9e \t %d \t %.9e \t %.9e\n",
        corrfunc_out.rbin[ibin], corrfunc_out.reff[ibin],
        corrfunc_out.npairs[ibin],
        corrfunc_out.xi[ibin].real(), corrfunc_out.xi[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::sys::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  return corrfunc_out;
}

/**
 * Compute power spectrum in a periodic box and optionally
 * save the results.
 *
 * @param particles_data (Data-source) particle catalogue.
 * @param params Parameter set.
 * @param kbins Wavenumber bins.
 * @param norm Normalisation factor.
 * @param norm_alt Alternative normalisation factor (default is 0.) which
 *                 may be printed out in the header of the
 *                 measurement output file.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns powspec_out Output power spectrum measurements.
 */
trv::PowspecMeasurements compute_powspec_in_box(
  ParticleCatalogue& particles_data,
  trv::ParameterSet& params,
  trv::Binning kbins,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: power spectrum in a periodic box.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  int* nmodes_save = new int[params.num_bins];
  double* k_save = new double[params.num_bins];
  std::complex<double>* pk_save = new std::complex<double>[params.num_bins];
  std::complex<double>* sn_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute power spectrum.
  MeshField dn(params);
  dn.compute_unweighted_field_fluctuations_insitu(particles_data);
  dn.fourier_transform();

  std::complex<double> sn_amp = double(particles_data.ntotal);

  FieldStats stats2pt(params);
  stats2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
    dn, dn, sn_amp, params.ELL, 0, kbins
  );

  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    nmodes_save[ibin] = stats2pt.nmodes[ibin];
    k_save[ibin] = stats2pt.k[ibin];
    pk_save[ibin] += double(2*params.ELL + 1) * stats2pt.pk[ibin];
    sn_save[ibin] += double(2*params.ELL + 1) * stats2pt.sn[ibin];
  }

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Power spectrum terms computed.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// Compute shot noise if the approach is particle-based.
  std::complex<double> sn_particle = calc_powspec_shotnoise_from_particles(
    particles_data
  );

  /// Normalisation factor should redue to
  /// ``params.volume / particles_data.ntotal / particles_data.ntotal``
  /// in the case of uniform particle weights.

  /* * Output ************************************************************** */

  /// Fill in output struct.
  trv::PowspecMeasurements powspec_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    powspec_out.kbin.push_back(kbins.bin_centres[ibin]);
    powspec_out.keff.push_back(k_save[ibin]);
    powspec_out.nmodes.push_back(nmodes_save[ibin]);
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
      save_fileptr, params, particles_data,
      norm, norm_alt, norm * sn_particle,
      "fourier"
    );
    for (int ibin = 0; ibin < params.num_bins; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        powspec_out.kbin[ibin], powspec_out.keff[ibin],
        powspec_out.nmodes[ibin],
        powspec_out.pk_raw[ibin].real(), powspec_out.pk_raw[ibin].imag(),
        powspec_out.pk_shot[ibin].real(), powspec_out.pk_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::sys::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] nmodes_save; delete[] k_save; delete[] pk_save; delete[] sn_save;

  return powspec_out;
}

/**
 * Compute two-point correlation function in a periodic box
 * and optionally save the results.
 *
 * @param particles_data (Data-source) particle catalogue.
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
trv::TwoPCFMeasurements compute_corrfunc_in_box(
  ParticleCatalogue& particles_data,
  trv::ParameterSet& params,
  trv::Binning& rbins,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: two-point correlation function "
      "in a periodic box.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  int* npairs_save = new int[params.num_bins];
  double* r_save = new double[params.num_bins];
  std::complex<double>* xi_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute 2PCF.
  MeshField dn(params);
  dn.compute_unweighted_field_fluctuations_insitu(particles_data);
  dn.fourier_transform();

  std::complex<double> sn_amp = double(particles_data.ntotal);

  FieldStats stats2pt(params);
  stats2pt.compute_ylm_wgtd_2pt_stats_in_config(
    dn, dn, sn_amp, params.ELL, 0, rbins
  );

  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    npairs_save[ibin] = stats2pt.npairs[ibin];
    r_save[ibin] = stats2pt.r[ibin];
    xi_save[ibin] += double(2*params.ELL + 1) * stats2pt.xi[ibin];
  }

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Two-point correlation function terms computed.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// Normalisation factor should redue to
  /// ``params.volume / particles_data.ntotal / particles_data.ntotal``
  /// in the case of uniform particle weights.

  /* * Output ************************************************************** */

  /// Fill in output struct.
  trv::TwoPCFMeasurements corrfunc_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    corrfunc_out.rbin.push_back(rbins.bin_centres[ibin]);
    corrfunc_out.reff.push_back(r_save[ibin]);
    corrfunc_out.npairs.push_back(npairs_save[ibin]);
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
      save_fileptr, params, particles_data, norm, norm_alt, 0., "config"
    );
    for (int ibin = 0; ibin < params.num_bins; ibin++) {
      std::fprintf(
        save_fileptr, "%.9e \t %.9e \t %d \t %.9e \t %.9e\n",
        corrfunc_out.rbin[ibin], corrfunc_out.reff[ibin],
        corrfunc_out.npairs[ibin],
        corrfunc_out.xi[ibin].real(), corrfunc_out.xi[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::sys::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  return corrfunc_out;
}

#ifdef TRV_USE_DISABLED_CODE
/**
 * Compute power spectrum window from a random catalogue and
 * optionally save the results.
 *
 * @param particles_rand (Random-source) particle catalogue.
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
trv::PowspecWindowMeasurements compute_powspec_window(
  ParticleCatalogue& particles_rand,
  LineOfSight* los_rand,
  trv::ParameterSet& params,
  std::vector<double> kbin,
  double alpha,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: power spectrum window "
      "from random catalogue.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up output.
  int* nmodes_save = new int[params.num_bins];
  double* k_save = new double[params.num_bins];
  std::complex<double>* pk_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute power spectrum window.
  MeshField dn_00(params);
  dn_00.compute_ylm_wgtd_field(particles_rand, los_rand, alpha, 0, 0);
  dn_00.fourier_transform();

  FieldStats stats2pt(params);
  std::complex<double> sn_amp = stats2pt.calc_ylm_wgtd_shotnoise_amp_for_powspec(
    particles_rand, los_rand, alpha, params.ELL, 0
  );

  stats2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
    dn_00, dn_00, sn_amp, params.ELL, 0, kbin
  );

  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    nmodes_save[ibin] = stats2pt.nmodes[ibin];
    k_save[ibin] = stats2pt.k[ibin];
    pk_save[ibin] += stats2pt.pk[ibin];
  }

  /// Renormalise the normalisation factor to be dimensionless.
  /// CAVEAT: Discretionary choice to eliminate the physical volume
  /// dimension by physical volume.
  norm /= params.volume;

  /* * Output ************************************************************** */

  /// Fill in output struct.
  trv::PowspecWindowMeasurements powwin_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    powwin_out.kbin.push_back(kbin[ibin]);
    powwin_out.keff.push_back(k_save[ibin]);
    powwin_out.nmodes.push_back(nmodes_save[ibin]);
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
    for (int ibin = 0; ibin < params.num_bins; ibin++) {
      std::fprintf(
        save_fileptr, "%.9e \t %.9e \t %d \t %.9e \t %.9e\n",
        powwin_out.kbin[ibin], powwin_out.keff[ibin],
        powwin_out.nmodes[ibin],
        powwin_out.pk[ibin].real(), powwin_out.pk[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::sys::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] nmodes_save; delete[] k_save; delete[] pk_save;

  return powwin_out;
}
#endif  // TRV_USE_DISABLED_CODE

/**
 * Compute two-point correlation function window from a random catalogue
 * and optionally save the results.
 *
 * @param particles_rand (Random-source) particle catalogue.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbins Separation bins.
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
trv::TwoPCFWindowMeasurements compute_corrfunc_window(
  trv::ParticleCatalogue& particles_rand,
  trv::LineOfSight* los_rand,
  trv::ParameterSet& params,
  trv::Binning rbins,
  double alpha,
  double norm,
  double norm_alt=0.,
  bool save=false
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: two-point correlation function window "
      "from random catalogue.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up input.
  int ell1 = params.ELL;

  /// Set up output.
  int* npairs_save = new int[params.num_bins];
  double* r_save = new double[params.num_bins];
  std::complex<double>* xi_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }

  /* * Measurement ********************************************************* */

  /// Compute 2PCF window.
  MeshField dn_00(params);
  dn_00.compute_ylm_wgtd_field(particles_rand, los_rand, alpha, 0, 0);
  dn_00.fourier_transform();

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params);
    dn_LM.compute_ylm_wgtd_field(
      particles_rand, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    FieldStats stats2pt(params);
    std::complex<double> sn_amp = stats2pt.calc_ylm_wgtd_shotnoise_amp_for_powspec(
      particles_rand, los_rand, alpha, params.ELL, M_
    );

    /// Calculate equivalent to (-1)^m_1 δ_{m_1, -M} which, after being
    /// summed over m_1, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = (2*params.ELL + 1) * (2*ell1 + 1)
          * wigner_3j(ell1, 0, params.ELL, 0, 0, 0)
          * wigner_3j(ell1, 0, params.ELL, m1, 0, M_);

      if (std::fabs(coupling) < EPS_COUPLING_2PT) {continue;}

      stats2pt.compute_ylm_wgtd_2pt_stats_in_config(
        dn_LM, dn_00, sn_amp, ell1, m1, rbins
      );

      for (int ibin = 0; ibin < params.num_bins; ibin++) {
        xi_save[ibin] += coupling * stats2pt.xi[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_bins; ibin++) {
          npairs_save[ibin] = stats2pt.npairs[ibin];
          r_save[ibin] = stats2pt.r[ibin];
        }
      }
    }

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Two-point correlation function window term "
        "at order M = %d computed.\n",
        trv::sys::show_timestamp().c_str(),
        M_
      );
    }
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  trv::TwoPCFWindowMeasurements corrfwin_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    corrfwin_out.rbin.push_back(rbins.bin_centres[ibin]);
    corrfwin_out.reff.push_back(r_save[ibin]);
    corrfwin_out.npairs.push_back(npairs_save[ibin]);
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
      save_fileptr, params, particles_rand, norm, norm_alt, 0., "config"
    );
    for (int ibin = 0; ibin < params.num_bins; ibin++) {
      std::fprintf(
        save_fileptr, "%.9e \t %.9e \t %d \t %.9e \t %.9e\n",
        corrfwin_out.rbin[ibin], corrfwin_out.reff[ibin],
        corrfwin_out.npairs[ibin],
        corrfwin_out.xi[ibin].real(), corrfwin_out.xi[ibin].imag()
      );
    }
    std::fclose(save_fileptr);

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s STAT] Measurements saved to %s.\n",
        trv::sys::show_timestamp().c_str(),
        save_filepath
      );
    }
  }

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  return corrfwin_out;
}

}  // trv::algo::
}  // trv::

#endif  // TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_
