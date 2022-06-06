/**
 * @file threept.hpp
 * @brief Three-point correlator computations.
 *
 * @note Unless otherwise specified, the Paper hereafter in this method's
 *       comments refers to Sugiyama et al. (2018)
 *       [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_

#include "monitor.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "particles.hpp"
#include "field.hpp"
#include "twopt.hpp"

using namespace trv::maths;
using namespace trv::obj;

const double EPS_COUPLING_3PT = 1.e-10;  /**< zero-tolerance for three-point
                                              coupling coefficients */
                                         // CAVEAT: discretionary choice

namespace trv {
namespace algo {

/**
 * Bispectrum measurements.
 *
 */
struct BispecMeasurements {
  std::vector<double> kbin1;  ///< first central wavenumber in bins
  std::vector<double> kbin2;  ///< second central wavenumber in bins
  std::vector<double> keff1;  ///< first effective wavenumber in bins
  std::vector<double> keff2;  ///< second effective wavenumber in bins
  std::vector<int> ntriag;  /**< number of contributing
                                 triangular configurations */
  std::vector< std::complex<double> > bk_raw;  /**< bispectrum
                                                    raw measurements */
  std::vector< std::complex<double> > bk_shot;  ///< bispectrum shot noise
};

/**
 * Three-point correlation function measurements.
 *
 */
struct ThreePCFMeasurements {
  std::vector<double> rbin1;  ///< first central separation in bins
  std::vector<double> rbin2;  ///< second central separation in bins
  std::vector<double> reff1;  ///< first effective separation in bins
  std::vector<double> reff2;  ///< second effective separation in bins
  std::vector<int> ntriag;  /**< number of contributing
                                 triangular configurations */
  std::vector< std::complex<double> > zeta_raw;  /**< three-point correlation
                                                      function raw
                                                      measurements */
  std::vector< std::complex<double> > zeta_shot;  /**< three-point correlation
                                                       function shot noise */
};

/**
 * Three-point correlation function window measurements.
 *
 */
struct ThreePCFWindowMeasurements {
  std::vector<double> rbin1;  ///< first central separation in bins
  std::vector<double> rbin2;  ///< second central separation in bins
  std::vector<double> reff1;  ///< first effective separation in bins
  std::vector<double> reff2;  ///< second effective separation in bins
  std::vector<int> ntriag;  /**< number of contributing
                                 triangular configurations */
  std::vector< std::complex<double> > zeta_raw;  /**< three-point correlation
                                                      function window raw
                                                      measurements */
  std::vector< std::complex<double> > zeta_shot;  /**< three-point correlation
                                                       function window
                                                       shot noise */
};

/**
 * Calculate mesh-based bispectrum normalisation.
 *
 * @param catalogue Particle catalogue.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @returns norm_factor Bispectrum normalisation factor.
 */
double calc_bispec_normalisation_from_mesh(
  ParticleCatalogue& catalogue,
  trv::scheme::ParameterSet& params,
  double alpha=1.
) {
  PseudoDensityField<ParticleCatalogue> catalogue_mesh(params);

  double norm_factor = catalogue_mesh._calc_wgt_cu_volume_norm(catalogue)
    / std::pow(alpha, 3);

  catalogue_mesh.finalise_density_field();

  return norm_factor;
}

/**
 * Calculate particle-based bispectrum normalisation.
 *
 * @param catalogue Particle catalogue.
 * @param alpha Alpha ratio.
 * @returns norm_factor Bispectrum normalisation factor.
 */
double calc_bispec_normalisation_from_particles(
  ParticleCatalogue& catalogue, double alpha=1.
) {
  double norm_factor = catalogue._calc_bispec_normalisation() / alpha;

  return norm_factor;
}

/// NOBUG: Standard naming convention is not always followed for
/// intermediary quantities in the functions below.

/// TODO: Shuffle functions.

/**
 * Compute bispectrum from paired catalogues and
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
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns bispec_out Output bispectrum measurements.
 */
BispecMeasurements compute_bispec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::scheme::ParameterSet& params,
  double* kbin,
  double alpha,
  double norm,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measuring bispectrum from data and random catalogues.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    std::fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (trv::runtime::currTask == 0) {
      throw trv::runtime::InvalidParameter(
        "[%s ERRO] Specified bispectrum multipole vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        trv::runtime::show_timestamp().c_str()
      );
    }
  }

  /// Set up output.
  int* ntriag_save = new int[params.num_kbin];
  double* k1_save = new double[params.num_kbin];
  double* k2_save = new double[params.num_kbin];
  std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
  std::complex<double>* sn_save = new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    ntriag_save[ibin] = 0;
    k1_save[ibin] = 0.;
    k2_save[ibin] = 0.;
    bk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

	/// Set up intermediary quantities.
  double vol_cell = params.volume / double(params.nmesh);

  double dr[3];
  dr[0] = params.boxsize[0] / double(params.ngrid[0]);
  dr[1] = params.boxsize[1] / double(params.ngrid[1]);
  dr[2] = params.boxsize[2] / double(params.ngrid[2]);

  /* * Measurement ********************************************************* */

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing shot noise terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute three shot noise terms out of four.
  PseudoDensityField<ParticleCatalogue> dn_00_for_sn(params);  // dn_00
  dn_00_for_sn.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00_for_sn.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
  			Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        /// Calculate \bar{S}_LM in eq. (46) in the Paper.
        /// QUEST: This is possibly redundant as `L` and `M`
        /// both need to be zero.
        std::complex<double> barS_LM =
          stats_sn.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );  // barS_LM

        /// Compute N_LM in eq. (46) in the Paper.
        PseudoDensityField<ParticleCatalogue> N_LM(params);  // N_LM
        N_LM.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
          particles_data, particles_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        N_LM.fourier_transform();

        /// Compute S_{ell1 ell2 L; i = j = k} in eq. (45) in the Paper.
        /// QUEST: This is possibly redundant and thus factorisable
        /// as `L` and `M` both need to be zero.
        if (params.ell1 == 0 && params.ell2 == 0) {
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            sn_save[ibin] += coupling * barS_LM;
          }
        }

        /// Compute S_{ell1 ell2 L; i ≠ j = k} in eq. (45) in the Paper.
        /// QUEST: This is possibly redundant and thus factorisable
        /// as `L` and `M` both need to be zero.
        if (params.ell2 == 0) {
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_sn, N_LM, barS_LM, kbin, params.ell1, m1_
          );
          if (params.form == "diag") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              sn_save[ibin] += coupling * (
                stats_sn.pk[ibin] - stats_sn.sn[ibin]
              );
            }
          } else
          if (params.form == "full") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              sn_save[ibin] += coupling * (
                stats_sn.pk[params.ith_kbin] - stats_sn.sn[params.ith_kbin]
              );
            }
          }
        }

        /// Compute S_{ell1 ell2 L; i = k ≠ j} in eq. (45) in the Paper.
        /// QUEST: This is possibly redundant and thus factorisable
        /// as `L` and `M` both need to be zero.
        if (params.ell1 == 0) {
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_sn, N_LM, barS_LM, kbin, params.ell2, m2_
          );
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            sn_save[ibin] += coupling * (stats_sn.pk[ibin] - stats_sn.sn[ibin]);
          }
        }

        if (trv::runtime::currTask == 0) {
          std::printf(
            "[%s STAT] Shot noise terms (3 out of 4) at order "
            "(m1, m2, M) = (%d, %d, %d) computed.\n",
            trv::runtime::show_timestamp().c_str(),
            m1_, m2_, M_
          );
        }
      }
    }
  }

  dn_00_for_sn.finalise_density_field();  // ~dn_00

  /// Compute the final shot noise term out of four.
  PseudoDensityField<ParticleCatalogue> N_00(params);  // N_00
  N_00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  N_00.fourier_transform();

  SphericalBesselCalculator sj1(params.ell1), sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      /// Check for vanishing cases where all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      if (flag_vanishing == "false") {
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_config_space(
            params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
          );
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_config_space(
            params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
          );
      }

      /// Compute S_{ell_1 ell_2 L; i = j ≠ k} in eq. (45) in the Paper.
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
  			Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        PseudoDensityField<ParticleCatalogue> dn_LM_for_sn(params);  // dn_LM
        dn_LM_for_sn.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM_for_sn.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in the Paper.
        std::complex<double> barS_LM =
          stats_sn.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );  // barS_LM

        /// Compute on mesh grids.
        fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
        trv::runtime::gbytesMem += double(params.nmesh)
          * sizeof(fftw_complex) / BYTES_PER_GBYTES;
        for (int gid = 0; gid < params.nmesh; gid++) {
          three_pt_holder[gid][0] = 0.;
          three_pt_holder[gid][1] = 0.;
        }

        stats_sn.compute_2pt_self_shotnoise_for_bispec_meshgrid(
          dn_LM_for_sn, N_00, barS_LM, three_pt_holder
          // params.ELL, M_
        );

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_a;
          double kmag_b = kbin[i_kbin];
          if (params.form == "diag") {
            kmag_a = kmag_b;
          } else
          if (params.form == "full") {
            kmag_a = kbin[params.ith_kbin];
          }

          double rvec[3];
          std::complex<double> S_ij_k = 0.;
          for (int i = 0; i < params.ngrid[0]; i++) {
            for (int j = 0; j < params.ngrid[1]; j++) {
              for (int k = 0; k < params.ngrid[2]; k++) {
                long long idx_grid =
                  (i * params.ngrid[1] + j) * params.ngrid[2] + k;

                /// This conforms to the absurd FFT array ordering
                /// convention that negative wavenumbers/frequencies come
                /// after zero and positive ones.
                rvec[0] = (i < params.ngrid[0]/2) ?
                  i * dr[0] : (i - params.ngrid[0]) * dr[0];
                rvec[1] = (j < params.ngrid[1]/2) ?
                  j * dr[1] : (j - params.ngrid[1]) * dr[1];
                rvec[2] = (k < params.ngrid[2]/2) ?
                  k * dr[2] : (k - params.ngrid[2]) * dr[2];
                double rmag = std::sqrt(
                  rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
                );

                double ja = sj1.eval(kmag_a * rmag);
                double jb = sj2.eval(kmag_b * rmag);
                std::complex<double> S_ij_k_gridpt(
                  three_pt_holder[idx_grid][0], three_pt_holder[idx_grid][1]
                );

                S_ij_k += ja * jb * ylm_a[idx_grid] * ylm_b[idx_grid]
                  * S_ij_k_gridpt;
              }
            }
          }

          S_ij_k *= vol_cell * std::pow(M_I, params.ell1 + params.ell2);

          sn_save[i_kbin] += coupling * S_ij_k;

          if (trv::runtime::currTask == 0) {
            std::printf(
              "[%s STAT] Shot noise term (last out of 4) at "
              "wavenumber k2 = %.4f and order (m1, m2, M) = (%d, %d, %d) "
              "computed.\n",
              trv::runtime::show_timestamp().c_str(),
              kmag_b, m1_, m2_, M_
            );
          }
        }

        fftw_free(three_pt_holder); three_pt_holder = NULL;
        trv::runtime::gbytesMem -= double(params.nmesh)
          * sizeof(fftw_complex) / BYTES_PER_GBYTES;
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed shot noise terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing bispectrum terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute bispectrum terms.
  PseudoDensityField<ParticleCatalogue> dn_00(params);  // dn_00
  dn_00.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      /// Check for vanishing cases where all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// QUEST: Escape early?
      // if (flag_vanishing == "true") {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      /// QUEST: This if-statement is redundant if early escape?
      if (flag_vanishing == "false") {
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_fourier_space(
            params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
          );
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_fourier_space(
            params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
          );
      }

      /// Compute a single term.
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        /// Compute G_LM in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> dn_LM(params);  // dn_LM
        dn_LM.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM.fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.inv_fourier_transform();

        /// Compute F_lm's in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> F_lm_a(params);  // F_lm_a
        double kmag_a, kmag_b;
        double dk = kbin[1] - kbin[0];  // HACK: only for linear binning
                                        // TODO: implement custom binning

        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_a, kmag_a, dk
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          kmag_b = kbin[i_kbin];

          PseudoDensityField<ParticleCatalogue> F_lm_b(params);  // F_lm_b
          F_lm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_b, kmag_b, dk
          );

          if (params.form == "diag") {
            kmag_a = kmag_b;
            F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
              dn_00, ylm_a, kmag_a, dk
            );
          }

          /// Add grid contribution.
          std::complex<double> bk_sum = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(dn_LM[gid][0], dn_LM[gid][1]);
            bk_sum += F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;
          }

          bk_save[i_kbin] += coupling * vol_cell * bk_sum;

          if (trv::runtime::currTask == 0) {
            std::printf(
              "[%s STAT] Bispectrum term at wavenumber k2 = %.4f and order "
              "(m1, m2, M) = (%d, %d, %d) computed.\n",
              trv::runtime::show_timestamp().c_str(),
              kmag_b, m1_, m2_, M_
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  dn_00.finalise_density_field();  // ~dn_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed bispectrum terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add bispectrum effective binning.
  BispecMeasurements bispec_out;
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    if (params.form == "diag") {
      bispec_out.kbin1.push_back(kbin[ibin]);
      bispec_out.keff1.push_back(kbin[ibin]);
    } else
    if (params.form == "full") {
      bispec_out.kbin1.push_back(kbin[params.ith_kbin]);
      bispec_out.keff1.push_back(kbin[params.ith_kbin]);
    }
    bispec_out.kbin2.push_back(kbin[ibin]);
    bispec_out.keff2.push_back(kbin[ibin]);
    bispec_out.ntriag.push_back(0);
    bispec_out.bk_raw.push_back(norm * bk_save[ibin]);
    bispec_out.bk_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      std::sprintf(
        save_filepath, "%s/bk%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      std::sprintf(
        save_filepath, "%s/bk%d%d%d_kbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_kbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        bispec_out.kbin1[ibin], bispec_out.keff1[ibin],
        bispec_out.kbin2[ibin], bispec_out.keff2[ibin],
        bispec_out.ntriag[ibin],
        bispec_out.bk_raw[ibin].real(), bispec_out.bk_raw[ibin].imag(),
        bispec_out.bk_shot[ibin].real(), bispec_out.bk_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] k1_save; delete[] k2_save;
  delete[] bk_save; delete[] sn_save;

  return bispec_out;
}

/**
 * Compute bispectrum in a periodic box and optionally save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param params Parameter set.
 * @param kbin Wavenumber bins.
 * @param norm Normalisation factor.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns bispec_out Output bispectrum measurements.
 */
BispecMeasurements compute_bispec_in_box(
  ParticleCatalogue& particles_data,
  trv::scheme::ParameterSet& params,
  double* kbin,
  double norm,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: bispectrum in a periodic box.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    std::fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (trv::runtime::currTask == 0) {
      throw trv::runtime::InvalidParameter(
        "[%s ERRO] Specified bispectrum multipole vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        trv::runtime::show_timestamp().c_str()
      );
    }
  }

  /// Set up output.
  int* ntriag_save = new int[params.num_kbin];
  double* k1_save = new double[params.num_kbin];
  double* k2_save = new double[params.num_kbin];
  std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
  std::complex<double>* sn_save = new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    ntriag_save[ibin] = 0;
    k1_save[ibin] = 0.;
    k2_save[ibin] = 0.;
    bk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

	/// Set up intermediary quantities.
  double vol_cell = params.volume / double(params.nmesh);

  double dr[3];
  dr[0] = params.boxsize[0] / double(params.ngrid[0]);
  dr[1] = params.boxsize[1] / double(params.ngrid[1]);
  dr[2] = params.boxsize[2] / double(params.ngrid[2]);

  /* * Measurement ********************************************************* */

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing shot noise terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute three shot noise terms out of four.
  /// QUEST: Why does this not include any weights?
  PseudoDensityField<ParticleCatalogue> dn_00_for_sn(params);  // dn_00
  dn_00_for_sn.compute_unweighted_fluctuation_insitu(
    particles_data, params.volume
  );
  dn_00_for_sn.fourier_transform();

  /// QUEST: Why is this not the sum of quadratic weights?
  PseudoDensityField<ParticleCatalogue> N_L0_for_sn(params);  // N_L0
  N_L0_for_sn.compute_unweighted_density(particles_data);
  N_L0_for_sn.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;
      Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

      /// Calculate Wigner-3j coupling coefficient.
      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Calculate \bar{S}_LM in eq. (46) in the Paper.
      /// QUEST: This is redundant and thus factorisable.
      std::complex<double> barS_LM = double(particles_data.ntotal);  // barS_LM

      /// Compute S_{ell1 ell2 L; i = j = k} in eq. (45) in the Paper.
      /// QUEST: This is possibly redundant and thus factorisable
      /// as `L` and `M` both need to be zero.
      /// QUEST: Why is this not the sum of cubic weights?
      if (params.ell1 == 0 && params.ell2 == 0) {
        for (int ibin = 0; ibin < params.num_kbin; ibin++) {
          sn_save[ibin] += coupling * barS_LM;
        }
      }

      /// Compute S_{ell1 ell2 L; i ≠ j = k} in eq. (45) in the Paper.
      /// QUEST: This is possibly redundant and thus factorisable
      /// as `L` and `M` both need to be zero.
      if (params.ell2 == 0) {
        stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
          dn_00_for_sn, N_L0_for_sn, barS_LM, kbin, params.ell1, m1_
        );
        if (params.form == "diag") {
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            sn_save[ibin] += coupling * (stats_sn.pk[ibin] - stats_sn.sn[ibin]);
          }
        } else
        if (params.form == "full") {
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            sn_save[ibin] += coupling * (
              stats_sn.pk[params.ith_kbin] - stats_sn.sn[params.ith_kbin]
            );
          }
        }
      }

      /// Compute S_{ell1 ell2 L; i = k ≠ j} in eq. (45) in the Paper.
      /// QUEST: This is possibly redundant and thus factorisable
      /// as `L` and `M` both need to be zero.
      if (params.ell1 == 0) {
        stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
          dn_00_for_sn, N_L0_for_sn, barS_LM, kbin, params.ell2, m2_
        );
        for (int ibin = 0; ibin < params.num_kbin; ibin++) {
          sn_save[ibin] += coupling * (stats_sn.pk[ibin] - stats_sn.sn[ibin]);
        }
      }

      if (trv::runtime::currTask == 0) {
        std::printf(
          "[%s STAT] Shot noise terms (3 out of 4) at order "
          "(m1, m2, M) = (%d, %d, %d) computed.\n",
          trv::runtime::show_timestamp().c_str(),
          m1_, m2_, M_
        );
      }
    }
  }

  dn_00_for_sn.finalise_density_field();  // ~dn_00
  N_L0_for_sn.finalise_density_field();  // ~N_L0

  /// Compute the final shot noise term out of four.
  /// QUEST: Why is this not the sum of quadratic weights?
  PseudoDensityField<ParticleCatalogue> N_00(params);  // N_00
	N_00.compute_unweighted_density(particles_data);
	N_00.fourier_transform();

  SphericalBesselCalculator sj1(params.ell1), sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;
      Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

      /// Calculate Wigner-3j coupling coefficient.
      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
      SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );

      PseudoDensityField<ParticleCatalogue> dn_L0_for_sn(params);  // dn_L0
      dn_L0_for_sn.compute_unweighted_fluctuation_insitu(
        particles_data, params.volume
      );
      dn_L0_for_sn.fourier_transform();

      /// Calculate \bar{S}_LM in eq. (46) in the Paper.
      /// QUEST: This is redundant and thus factorisable.
      std::complex<double> barS_L0 = double(particles_data.ntotal);  // barS_L0

      /// Compute on mesh grids.
      fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
        trv::runtime::gbytesMem += double(params.nmesh)
          * sizeof(fftw_complex) / BYTES_PER_GBYTES;
      for (int gid = 0; gid < params.nmesh; gid++) {
        three_pt_holder[gid][0] = 0.;
        three_pt_holder[gid][1] = 0.;
      }

      stats_sn.compute_2pt_self_shotnoise_for_bispec_meshgrid(
        dn_L0_for_sn, N_00, barS_L0, three_pt_holder
        // params.ELL, M_
      );

      for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
        double kmag_a;
        double kmag_b = kbin[i_kbin];
        if (params.form == "diag") {
          kmag_a = kmag_b;
        } else
        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
        }

        double rvec[3];
        std::complex<double> S_ij_k = 0.;
        for (int i = 0; i < params.ngrid[0]; i++) {
          for (int j = 0; j < params.ngrid[1]; j++) {
            for (int k = 0; k < params.ngrid[2]; k++) {
              long long idx_grid =
                (i * params.ngrid[1] + j) * params.ngrid[2] + k;

              /// This conforms to the absurd FFT array ordering
              /// convention that negative wavenumbers/frequencies come
              /// after zero and positive wavenumbers/frequencies.
              rvec[0] = (i < params.ngrid[0]/2) ?
                i * dr[0] : (i - params.ngrid[0]) * dr[0];
              rvec[1] = (j < params.ngrid[1]/2) ?
                j * dr[1] : (j - params.ngrid[1]) * dr[1];
              rvec[2] = (k < params.ngrid[2]/2) ?
                k * dr[2] : (k - params.ngrid[2]) * dr[2];
              double rmag = std::sqrt(
                rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
              );

              double ja = sj1.eval(kmag_a * rmag);
              double jb = sj2.eval(kmag_b * rmag);
              std::complex<double> S_ij_k_gridpt(
                three_pt_holder[idx_grid][0], three_pt_holder[idx_grid][1]
              );

              S_ij_k += ja * jb * ylm_a[idx_grid] * ylm_b[idx_grid]
                * S_ij_k_gridpt;
            }
          }
        }

        S_ij_k *= vol_cell * std::pow(M_I, params.ell1 + params.ell2);

        sn_save[i_kbin] += coupling * S_ij_k;

        if (trv::runtime::currTask == 0) {
          std::printf(
            "[%s STAT] Shot noise term (last out of 4) at "
            "wavenumber k2 = %.4f and order (m1, m2, M) = (%d, %d, %d) "
            "computed.\n",
            trv::runtime::show_timestamp().c_str(),
            kmag_b, m1_, m2_, M_
          );
        }
      }

      fftw_free(three_pt_holder); three_pt_holder = NULL;
        trv::runtime::gbytesMem -= double(params.nmesh)
          * sizeof(fftw_complex) / BYTES_PER_GBYTES;

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed shot noise terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing bispectrum terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute bispectrum terms.
  /// QUEST: Why does this not include any weights?
  PseudoDensityField<ParticleCatalogue> dn_00(params);  // dn_00
  dn_00.compute_unweighted_fluctuation_insitu(
    particles_data, params.volume
  );
  dn_00.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;

      /// Calculate Wigner-3j coupling coefficient.
      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
      SphericalHarmonicCalculator
        ::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );

      /// Compute G_00 in eq. (42) in the Paper (L, M = 0 in the global
      /// plane-parallel picture).
      /// QUEST: Why does this not include any weights?
      PseudoDensityField<ParticleCatalogue> dn_00_(params);
      dn_00_.compute_unweighted_fluctuation_insitu(
        particles_data, params.volume
      );
      dn_00_.fourier_transform();
      dn_00_.apply_assignment_compensation();
      dn_00_.inv_fourier_transform();

      /// Compute F_lm's in eq. (42) in the Paper.
      PseudoDensityField<ParticleCatalogue> F_lm_a(params);  // F_lm_a
      double kmag_a, kmag_b;
      double dk = kbin[1] - kbin[0];  // HACK: only for linear binning
                                      // TODO: implement custom binning

      if (params.form == "full") {
        kmag_a = kbin[params.ith_kbin];
        F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
          dn_00, ylm_a, kmag_a, dk
        );
      }

      for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
        kmag_b = kbin[i_kbin];

        PseudoDensityField<ParticleCatalogue> F_lm_b(params);  // F_lm_b
        F_lm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
          dn_00, ylm_b, kmag_b, dk
        );

        if (params.form == "diag") {
          kmag_a = kmag_b;
          F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_a, kmag_a, dk
          );
        }

        /// Add grid contribution.
        std::complex<double> bk_sum = 0.;
        for (int gid = 0; gid < params.nmesh; gid++) {
          std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
          std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
          std::complex<double> G_00_gridpt(dn_00_[gid][0], dn_00_[gid][1]);
          bk_sum += F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;
        }

        bk_save[i_kbin] += coupling * vol_cell * bk_sum;

        if (trv::runtime::currTask == 0) {
          std::printf(
            "[%s STAT] Bispectrum term at wavenumber k2 = %.4f and order "
            "(m1, m2, M) = (%d, %d, %d) computed.\n",
            trv::runtime::show_timestamp().c_str(),
            kmag_b, m1_, m2_, M_
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  dn_00.finalise_density_field();  // ~dn_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed bispectrum terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add bispectrum effective binning.
  BispecMeasurements bispec_out;
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    if (params.form == "diag") {
      bispec_out.kbin1.push_back(kbin[ibin]);
      bispec_out.keff1.push_back(kbin[ibin]);
    } else
    if (params.form == "full") {
      bispec_out.kbin1.push_back(kbin[params.ith_kbin]);
      bispec_out.keff1.push_back(kbin[params.ith_kbin]);
    }
    bispec_out.kbin2.push_back(kbin[ibin]);
    bispec_out.keff2.push_back(kbin[ibin]);
    bispec_out.ntriag.push_back(0);
    bispec_out.bk_raw.push_back(norm * bk_save[ibin]);
    bispec_out.bk_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      std::sprintf(
        save_filepath, "%s/bk%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      std::sprintf(
        save_filepath, "%s/bk%d%d%d_kbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_kbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t "
        "%.9e \t %.9e \t %.9e \t %.9e\n",
        bispec_out.kbin1[ibin], bispec_out.keff1[ibin],
        bispec_out.kbin2[ibin], bispec_out.keff2[ibin],
        bispec_out.ntriag[ibin],
        bispec_out.bk_raw[ibin].real(), bispec_out.bk_raw[ibin].imag(),
        bispec_out.bk_shot[ibin].real(), bispec_out.bk_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] k1_save; delete[] k2_save;
  delete[] bk_save; delete[] sn_save;

  return bispec_out;
}

/**
 * Compute three-point correlation function from paired catalogues and
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
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns threepcf_out Output three-point correlation function
 *                       measurements.
 */
ThreePCFMeasurements compute_3pcf(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::scheme::ParameterSet& params,
  double* rbin,
  double alpha,
  double norm,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measuring three-point correlation function "
      "from data and random catalogues.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    std::fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (trv::runtime::currTask == 0) {
      throw trv::runtime::InvalidParameter(
        "[%s ERRO] Specified three-point correlation function multipole"
        "vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        trv::runtime::show_timestamp().c_str()
      );
    }
  }

  /// Set up output.
  int* ntriag_save = new int[params.num_rbin];
  double* r1_save = new double[params.num_rbin];
  double* r2_save = new double[params.num_rbin];
  std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
  std::complex<double>* sn_save = new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    ntriag_save[ibin] = 0;
    r1_save[ibin] = 0.;
    r2_save[ibin] = 0.;
    zeta_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

	/// Set up intermediary quantities.
  double vol_cell = params.volume / double(params.nmesh);

  /* * Measurement ********************************************************* */

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing shot noise terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute the only shot noise term S_{ell1 ell2 L; i = j ≠ k}
  /// in eq. (51) in the Paper.
  PseudoDensityField<ParticleCatalogue> N_00(params);  // N_00
  N_00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  N_00.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      /// Check for vanishing cases where all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      if (flag_vanishing == "false") {
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_config_space(
            params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
          );
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_config_space(
            params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
          );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        PseudoDensityField<ParticleCatalogue> dn_LM_for_sn(params);  // dn_LM
        dn_LM_for_sn.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM_for_sn.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in the Paper.
        std::complex<double> barS_LM =
          stats_sn.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );  // barS_LM

        stats_sn.compute_uncoupled_shotnoise_for_3pcf(
          dn_LM_for_sn, N_00, ylm_a, ylm_b, barS_LM, rbin
          // params.ell1, m1_
        );

        for (int ibin = 0; ibin < params.num_rbin; ibin++) {
          if (params.form == "diag") {
            sn_save[ibin] += coupling * stats_sn.xi[ibin];
          } else
          if (params.form == "full") {
            /// Calculate shot noise contribution equivalent to the
            /// Kronecker delta in eq. (51) in the Paper.
            if (ibin == params.ith_rbin) {
              sn_save[ibin] += coupling * stats_sn.xi[ibin];
            } else {
              sn_save[ibin] += 0.;
            }
          }
        }

        if (trv::runtime::currTask == 0) {
          std::printf(
            "[%s STAT] Shot noise term at order "
            "(m1, m2, M) = (%d, %d, %d) computed.\n",
            trv::runtime::show_timestamp().c_str(),
            m1_, m2_, M_
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed shot noise terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing three-point correlation function terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute 3PCF terms.
  PseudoDensityField<ParticleCatalogue> dn_00(params);  // dn_00
  dn_00.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  SphericalBesselCalculator sj1(params.ell1), sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      /// Check for vanishing cases where all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// QUEST: Escape early?
      // if (flag_vanishing == "true") {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      /// QUEST: This if-statement is redundant if early escape?
      if (flag_vanishing == "false") {
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_fourier_space(
            params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
          );
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_fourier_space(
            params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
          );
      }

      /// Compute a single term.
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        /// Compute G_LM in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> dn_LM(params);  // dn_LM
        dn_LM.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM.fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.inv_fourier_transform();

        /// Compute F_lm's in eq. (49) in the Paper.
        PseudoDensityField<ParticleCatalogue> F_lm_a(params);  // F_lm_a
        double rmag_a, rmag_b;

        if (params.form == "full") {
          rmag_a = rbin[params.ith_rbin];
          F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            dn_00, ylm_a, sj1, rmag_a
          );
        }

        for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
          rmag_b = rbin[i_rbin];

          PseudoDensityField<ParticleCatalogue> F_lm_b(params);
          F_lm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            dn_00, ylm_b, sj2, rmag_b
          );

          if (params.form == "diag") {
            rmag_a = rmag_b;
            F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
              dn_00, ylm_a, sj1, rmag_a
            );
          }

          /// Add grid contribution.
          std::complex<double> zeta_sum = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(dn_LM[gid][0], dn_LM[gid][1]);
            zeta_sum += F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;
          }

          zeta_save[i_rbin] += std::pow(M_I, params.ell1 + params.ell2)
            * coupling * vol_cell * zeta_sum;

          if (trv::runtime::currTask == 0) {
            std::printf(
              "[%s STAT] Three-point correlation function term at "
              "separation r2 = %.3f and order (m1, m2, M) = (%d, %d, %d) "
              "computed.\n",
              trv::runtime::show_timestamp().c_str(),
              rmag_b, m1_, m2_, M_
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  dn_00.finalise_density_field();  // ~dn_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed three-point correlation function terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add 3PCF effective binning.
  ThreePCFMeasurements threepcf_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    if (params.form == "diag") {
      threepcf_out.rbin1.push_back(rbin[ibin]);
      threepcf_out.reff1.push_back(rbin[ibin]);
    } else
    if (params.form == "full") {
      threepcf_out.rbin1.push_back(rbin[params.ith_rbin]);
      threepcf_out.reff1.push_back(rbin[params.ith_rbin]);
    }
    threepcf_out.rbin2.push_back(rbin[ibin]);
    threepcf_out.reff2.push_back(rbin[ibin]);
    threepcf_out.ntriag.push_back(0);
    threepcf_out.zeta_raw.push_back(norm * zeta_save[ibin]);
    threepcf_out.zeta_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      std::sprintf(
        save_filepath, "%s/zeta%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      std::sprintf(
        save_filepath, "%s/zeta%d%d%d_rbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_rbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        threepcf_out.rbin1[ibin], threepcf_out.reff1[ibin],
        threepcf_out.rbin2[ibin], threepcf_out.reff2[ibin],
        threepcf_out.ntriag[ibin],
        threepcf_out.zeta_raw[ibin].real(), threepcf_out.zeta_raw[ibin].imag(),
        threepcf_out.zeta_shot[ibin].real(), threepcf_out.zeta_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  return threepcf_out;
}

/**
 * Compute three-point correlation function in a periodic box
 * and save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param params Parameter set.
 * @param rbin Separation bins.
 * @param norm Normalisation factor.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns threepcf_out Output three-point correlation
 *                       function measurements.
 */
ThreePCFMeasurements compute_3pcf_in_box(
  ParticleCatalogue& particles_data,
  trv::scheme::ParameterSet& params,
  double* rbin,
  double norm,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: three-point correlation function "
      "in a periodic box.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    std::fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (trv::runtime::currTask == 0) {
      throw trv::runtime::InvalidParameter(
        "[%s ERRO] Specified three-point correlation function multipole "
        "vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        trv::runtime::show_timestamp().c_str()
      );
    }
  }

  /// Set up output.
  int* ntriag_save = new int[params.num_rbin];
  double* r1_save = new double[params.num_rbin];
  double* r2_save = new double[params.num_rbin];
  std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
  std::complex<double>* sn_save = new std::complex<double>[params.num_rbin];
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    ntriag_save[ibin] = 0;
    r1_save[ibin] = 0.;
    r2_save[ibin] = 0.;
    zeta_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

	/// Set up intermediary quantities.
  double vol_cell = params.volume / double(params.nmesh);

  /* * Measurement ********************************************************* */

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing shot noise terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute the only shot noise term S_{ell1 ell2 L; i = j ≠ k}
  /// in eq. (51) in the Paper.
  /// QUEST: Why is this not the sum of quadratic weights?
  PseudoDensityField<ParticleCatalogue> N_00(params);
  N_00.compute_unweighted_density(particles_data);
  N_00.fourier_transform();

  /// Compute shot noise terms.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;
      Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

      /// Calculate Wigner-3j coupling coefficient.
      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
      SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );

      PseudoDensityField<ParticleCatalogue> dn_L0_for_sn(params);  // dn_L0
      dn_L0_for_sn.compute_unweighted_fluctuation_insitu(
        particles_data, params.volume
      );
      dn_L0_for_sn.fourier_transform();

      /// Calculate \bar{S}_LM in eq. (46) in the Paper.
      /// QUEST: This is redundant and thus factorisable.
      std::complex<double> barS_L0 = double(particles_data.ntotal);  // barS_L0

      stats_sn.compute_uncoupled_shotnoise_for_3pcf(
        dn_L0_for_sn, N_00, ylm_a, ylm_b, barS_L0, rbin
        // params.ell1, m1_
      );

      for (int ibin = 0; ibin < params.num_rbin; ibin++) {
        if (params.form == "diag") {
          sn_save[ibin] += coupling * stats_sn.xi[ibin];
        } else if (params.form == "full") {
          /// Calculate shot noise contribution equivalent to the
          /// Kronecker delta in eq. (51) in the Paper.
          if (ibin == params.ith_rbin) {
            sn_save[ibin] += coupling * stats_sn.xi[ibin];
          } else {
            sn_save[ibin] += 0.;
          }
        }
      }

      if (trv::runtime::currTask == 0) {
        std::printf(
          "[%s STAT] Shot noise term at order "
          "(m1, m2, M) = (%d, %d, %d) computed.\n",
          trv::runtime::show_timestamp().c_str(),
          m1_, m2_, M_
        );
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed shot noise terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing bispectrum terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute 3PCF terms.
  /// QUEST: Why does this not include any weights?
  PseudoDensityField<ParticleCatalogue> dn_00(params);  // dn_00
  dn_00.compute_unweighted_fluctuation_insitu(
    particles_data, params.volume
  );
  dn_00.fourier_transform();

  SphericalBesselCalculator sj1(params.ell1), sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;

      /// Calculate Wigner-3j coupling coefficient.
      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
      SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );

      /// Compute G_00 in eq. (42) in the Paper (L, M = 0 in the global
      /// plane-parallel picture).
      /// QUEST: Why does this not include any weights?
      PseudoDensityField<ParticleCatalogue> dn_00_(params);
      dn_00_.compute_unweighted_fluctuation_insitu(
        particles_data, params.volume
      );
      dn_00_.fourier_transform();
      dn_00_.apply_assignment_compensation();
      dn_00_.inv_fourier_transform();

      /// Compute F_lm's in eq. (42) in the Paper.
      PseudoDensityField<ParticleCatalogue> F_lm_a(params);
      double rmag_a, rmag_b;

      if (params.form == "full") {
        rmag_a = rbin[params.ith_rbin];
        F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
          dn_00, ylm_a, sj1, rmag_a
        );
      }

      for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
        rmag_b = rbin[i_rbin];

        PseudoDensityField<ParticleCatalogue> F_lm_b(params);  // F_lm_b
        F_lm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
          dn_00, ylm_b, sj2, rmag_b
        );

        if (params.form == "diag") {
          rmag_a = rmag_b;
          F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            dn_00, ylm_a, sj1, rmag_a
          );
        }

        /// Add grid contribution.
        std::complex<double> zeta_sum = 0.;
        for (int gid = 0; gid < params.nmesh; gid++) {
          std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
          std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
          std::complex<double> G_00_gridpt(dn_00_[gid][0], dn_00_[gid][1]);
          zeta_sum += F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;
        }

        zeta_save[i_rbin] += std::pow(M_I, params.ell1 + params.ell2)
          * coupling * vol_cell * zeta_sum;

        if (trv::runtime::currTask == 0) {
          std::printf(
            "[%s STAT] Three-point correlation function term at "
            "separation r2 = %.3f and order (m1, m2, M) = (%d, %d, %d) "
            "computed.\n",
            trv::runtime::show_timestamp().c_str(),
            rmag_b, m1_, m2_, M_
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  dn_00.finalise_density_field();  // ~dn_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed three-point correlation function terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add 3PCF effective binning.
  ThreePCFMeasurements threepcf_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    if (params.form == "diag") {
      threepcf_out.rbin1.push_back(rbin[ibin]);
      threepcf_out.reff1.push_back(rbin[ibin]);
    } else
    if (params.form == "full") {
      threepcf_out.rbin1.push_back(rbin[params.ith_rbin]);
      threepcf_out.reff1.push_back(rbin[params.ith_rbin]);
    }
    threepcf_out.rbin2.push_back(rbin[ibin]);
    threepcf_out.reff2.push_back(rbin[ibin]);
    threepcf_out.ntriag.push_back(0);
    threepcf_out.zeta_raw.push_back(norm * zeta_save[ibin]);
    threepcf_out.zeta_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      std::sprintf(
        save_filepath, "%s/zeta%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      std::sprintf(
        save_filepath, "%s/zeta%d%d%d_rbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_rbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        threepcf_out.rbin1[ibin], threepcf_out.reff1[ibin],
        threepcf_out.rbin2[ibin], threepcf_out.reff2[ibin],
        threepcf_out.ntriag[ibin],
        threepcf_out.zeta_raw[ibin].real(), threepcf_out.zeta_raw[ibin].imag(),
        threepcf_out.zeta_shot[ibin].real(), threepcf_out.zeta_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  return threepcf_out;
}

/**
 * Compute three-point correlation function window from a random
 * catalogue and optionally save the results.
 *
 * @param particles_rand (Random-source) particle container.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbin Separation bins.
 * @param alpha Alpha ratio.
 * @param norm Normalisation factor.
 * @param wide_angle If `true` (default is `false`), compute wide-angle
 *                   correction terms as set by `params`.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns threepcfwin_out Output three-point correlation function
 *                          window measurements.
 */
ThreePCFWindowMeasurements compute_3pcf_window(
  ParticleCatalogue& particles_rand,
  LineOfSight* los_rand,
  trv::scheme::ParameterSet& params,
  double* rbin,
  double alpha,
  double norm,
  bool wide_angle=false,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measurement: three-point correlation function window "
      "from random catalogue.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    std::fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (trv::runtime::currTask == 0) {
      throw trv::runtime::InvalidParameter(
        "[%s ERRO] Specified three-point correlation function window "
        "multipole vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        trv::runtime::show_timestamp().c_str()
      );
    }
  }

  /// Set up output.
  int* ntriag_save = new int[params.num_rbin];
  double* r1_save = new double[params.num_rbin];
  double* r2_save = new double[params.num_rbin];
  std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
  std::complex<double>* sn_save = new std::complex<double>[params.num_rbin];
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    ntriag_save[ibin] = 0;
    r1_save[ibin] = 0.;
    r2_save[ibin] = 0.;
    zeta_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

	/// Set up intermediary quantities.
  double vol_cell = params.volume / double(params.nmesh);

  /* * Measurement ********************************************************* */

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing shot noise terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute the only shot noise term S_{ell1 ell2 L; i = j ≠ k}
  /// in eq. (51) in the Paper.
  PseudoDensityField<ParticleCatalogue> N_00(params);  // N_00
  N_00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
    particles_rand, los_rand, alpha, 0, 0
  );
  N_00.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      /// Check for vanishing cases where all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      if (flag_vanishing == "false") {
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_config_space(
            params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
          );
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_config_space(
            params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
          );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        PseudoDensityField<ParticleCatalogue> n_LM_for_sn(params);  // n_LM
        n_LM_for_sn.compute_ylm_wgtd_density(
          particles_rand, los_rand, alpha, params.ELL, M_
        );
        n_LM_for_sn.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in the Paper.
        /// QUEST: Check the power index in isotropic shot noise.
        std::complex<double> barS_LM =
          stats_sn.calc_ylm_wgtd_shotnoise_for_powspec(
            particles_rand, los_rand, alpha, params.ELL, M_
          );  // barS_LM

        stats_sn.compute_uncoupled_shotnoise_for_3pcf(
          n_LM_for_sn, N_00, ylm_a, ylm_b, barS_LM, rbin
          // params.ell1, m1_
        );

        for (int ibin = 0; ibin < params.num_rbin; ibin++) {
          if (params.form == "diag") {
            sn_save[ibin] += coupling * stats_sn.xi[ibin];
          } else if (params.form == "full") {
            /// Calculate shot noise contribution equivalent to
            /// the Kronecker delta in eq. (51) in the Paper.
            if (ibin == params.ith_rbin) {
              sn_save[ibin] += coupling * stats_sn.xi[ibin];
            } else {
              sn_save[ibin] += 0.;
            }
          }
        }

        if (trv::runtime::currTask == 0) {
          std::printf(
            "[%s STAT] Shot noise term at order "
            "(m1, m2, M) = (%d, %d, %d) computed.\n",
            trv::runtime::show_timestamp().c_str(),
            m1_, m2_, M_
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed shot noise terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing three-point correlation function window terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute 3PCF window terms.
  PseudoDensityField<ParticleCatalogue> n_00(params);  // n_00
  n_00.compute_ylm_wgtd_density(particles_rand, los_rand, alpha, 0, 0);
  n_00.fourier_transform();

  SphericalBesselCalculator sj1(params.ell1), sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      /// Check for vanishing cases where all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// QUEST: Escape early?
      // if (flag_vanishing == "true") {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      /// QUEST: This if-statement is redundant if early escape?
      if (flag_vanishing == "false") {
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_fourier_space(
            params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
          );
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_fourier_space(
            params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
          );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        /// Compute G_LM in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> n_LM(params);  // n_LM
        n_LM.compute_ylm_wgtd_density(
          particles_rand, los_rand, alpha, params.ELL, M_
        );
        n_LM.fourier_transform();
        n_LM.apply_assignment_compensation();
        n_LM.inv_fourier_transform();

        /// Perform wide-angle corrections if required.
        if (wide_angle) {
          n_LM.apply_power_law_weight_for_wide_angle();
        }

        /// Compute F_lm's in eq. (49) in the Paper.
        PseudoDensityField<ParticleCatalogue> F_lm_a(params);  // F_lm_a
        double rmag_a, rmag_b;

        if (params.form == "full") {
          rmag_a = rbin[params.ith_rbin];
          F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            n_00, ylm_a, sj1, rmag_a
          );
        }

        for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
          rmag_b = rbin[i_rbin];

          PseudoDensityField<ParticleCatalogue> F_lm_b(params);
          F_lm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            n_00, ylm_b, sj2, rmag_b
          );

          if (params.form == "diag") {
            rmag_a = rmag_b;
            F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
              n_00, ylm_a, sj1, rmag_a
            );
          }

          /// Add grid contribution.
          std::complex<double> zeta_sum = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(n_LM[gid][0], n_LM[gid][1]);
            zeta_sum += F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;
          }

          zeta_save[i_rbin] += std::pow(M_I, params.ell1 + params.ell2)
            * coupling * vol_cell * zeta_sum;

          if (trv::runtime::currTask == 0) {
            std::printf(
              "[%s STAT] Three-point correlation function window term at "
              "separation r2 = %.3f and order (m1, m2, M) = (%d, %d, %d) "
              "computed.\n",
              trv::runtime::show_timestamp().c_str(),
              rmag_b, m1_, m2_, M_
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  n_00.finalise_density_field();  // ~n_00

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed three-point correlation function window terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add 3PCF effective binning.
  ThreePCFWindowMeasurements threepcfwin_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    if (params.form == "diag") {
      threepcfwin_out.rbin1.push_back(rbin[ibin]);
      threepcfwin_out.reff1.push_back(rbin[ibin]);
    } else
    if (params.form == "full") {
      threepcfwin_out.rbin1.push_back(rbin[params.ith_rbin]);
      threepcfwin_out.reff1.push_back(rbin[params.ith_rbin]);
    }
    threepcfwin_out.rbin2.push_back(rbin[ibin]);
    threepcfwin_out.reff2.push_back(rbin[ibin]);
    threepcfwin_out.ntriag.push_back(0);
    threepcfwin_out.zeta_raw.push_back(norm * zeta_save[ibin]);
    threepcfwin_out.zeta_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  std::string window_tag;
  if (wide_angle) {
    char window_tag_[32];
    std::sprintf(window_tag_, "window_wa%d%d", params.i_wa, params.j_wa);
    window_tag = std::string(window_tag_);
  } else {
    window_tag = "window";
  }

  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      std::sprintf(
        save_filepath, "%s/zeta%d%d%d_%s%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        window_tag.c_str(), params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      std::sprintf(
        save_filepath, "%s/zeta%d%d%d_%s_rbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        window_tag.c_str(),
        params.ith_rbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        threepcfwin_out.rbin1[ibin], threepcfwin_out.reff1[ibin],
        threepcfwin_out.rbin2[ibin], threepcfwin_out.reff2[ibin],
        threepcfwin_out.ntriag[ibin],
        threepcfwin_out.zeta_raw[ibin].real(),
        threepcfwin_out.zeta_raw[ibin].imag(),
        threepcfwin_out.zeta_shot[ibin].real(),
        threepcfwin_out.zeta_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  return threepcfwin_out;
}

#ifdef TRIUMVIRATE_USE_DISABLED_CODE
/**
 * Compute bispectrum from paired catalogues with respect to a choice of
 * line of sight and optionally save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param particles_rand (Random-source) particle container.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param kbin Wavenumber bins.
 * @param alpha Alpha ratio.
 * @param norm Survey volume normalisation constant.
 * @param los_choice Choice of line of sight {0, 1 or 2}.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns bispec_out Output bispectrum measurements.
 */
BispecMeasurements compute_bispec_for_los_choice(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::scheme::ParameterSet& params,
  double* kbin,
  double alpha,
  double norm,
  int los_choice,
  bool save=false
) {
  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Measuring bispectrum from data and random catalogues "
      "for the specified choice of line of sight.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    std::fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (trv::runtime::currTask == 0) {
      throw trv::runtime::InvalidParameter(
        "[%s ERRO] Specified bispectrum multipole vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        trv::runtime::show_timestamp().c_str()
      );
    }
  }

  /// Set up output.
  int* ntriag_save = new int[params.num_kbin];
  double* k1_save = new double[params.num_kbin];
  double* k2_save = new double[params.num_kbin];
  std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
  std::complex<double>* sn_save = new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    ntriag_save[ibin] = 0;
    k1_save[ibin] = 0.;
    k2_save[ibin] = 0.;
    bk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }

	/// Set up intermediary quantities.
  double vol_cell = params.volume / double(params.nmesh);

  double dr[3];
  dr[0] = params.boxsize[0] / double(params.ngrid[0]);
  dr[1] = params.boxsize[1] / double(params.ngrid[1]);
  dr[2] = params.boxsize[2] / double(params.ngrid[2]);

  /* * Measurement ********************************************************* */

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing shot noise terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute three shot noise terms out of four.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        PseudoDensityField<ParticleCatalogue> dn_LM00_for_sn(params);
          // dn_LM or dn_00
        if (los_choice == 0) {
          dn_LM00_for_sn.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM00_for_sn.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM00_for_sn.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in the Paper.
        std::complex<double> barS_LM =
          stats_sn.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );

        /// Compute N_LM in eq. (46) in the Paper.
        PseudoDensityField<ParticleCatalogue> N_LM00(params);  // N_LM or N_00
        if (los_choice == 0) {
          N_LM00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand,
            alpha, 0, 0
          );
        } else {
          N_LM00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand,
            alpha, params.ELL, M_
          );
        }
        N_LM00.fourier_transform();

        /// Compute S_{ell1 ell2 L; i = j = k} in eq. (45) in the Paper.
        if (params.ell1 == 0 && params.ell2 == 0) {
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            sn_save[ibin] += coupling * barS_LM;
          }
        }

        /// Compute S_{ell1 ell2 L; i ≠ j = k} in eq. (45) in the Paper.
        if (params.ell2 == 0) {
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_LM00_for_sn, N_LM00, barS_LM, kbin,
            params.ell1, m1_
          );
          if (params.form == "diag") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              sn_save[ibin] += coupling * (
                stats_sn.pk[ibin] - stats_sn.sn[ibin]
              );
            }
          } else if (params.form == "full") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              sn_save[ibin] += coupling * (
                stats_sn.pk[params.ith_kbin] - stats_sn.sn[params.ith_kbin]
              );
            }
          }
        }

        /// WARNING: This was inherited from Sugiyama et al. without a
        /// matching equation.  S_{ell1 ell2 L; i = k ≠ j}
        /// (i.e. ell1 == 0 case) has been shuffled (see below).
      }
    }
  }

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        PseudoDensityField<ParticleCatalogue> dn_LM00_for_sn(params);
          // dn_LM or dn_00
        if (los_choice == 1) {
          dn_LM00_for_sn.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM00_for_sn.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM00_for_sn.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in the Paper.
        std::complex<double> barS_LM =
          stats_sn.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );  // barS_LM

        /// Compute N_LM in eq. (46) in the Paper.
        PseudoDensityField<ParticleCatalogue> N_LM00(params);  // N_LM or N_00
        if (los_choice == 1) {
          N_LM00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            0, 0
          );
        } else {
          N_LM00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        }
        N_LM00.fourier_transform();

        /// Compute S_{ell1 ell2 L; i = k ≠ j} in eq. (45) in the Paper.
        if (params.ell1 == 0) {
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_LM00_for_sn, N_LM00, barS_LM, kbin,
            params.ell2, m2_
          );
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            sn_save[ibin] += coupling * (stats_sn.pk[ibin] - stats_sn.sn[ibin]);
          }
        }

        /// WARNING: This was inherited from Sugiyama et al. without a
        /// matching equation.  S_{ell1 ell2 L; i = k ≠ j}
        /// (i.e. ell1 == 0 case) has been shuffled (see above).

        if (trv::runtime::currTask == 0) {
          std::printf(
            "[%s STAT] Shot noise terms (3 out of 4) at order "
            "(m1, m2, M) = (%d, %d, %d) computed.\n",
            trv::runtime::show_timestamp().c_str(),
            m1_, m2_, M_
          );
        }
      }
    }
  }

  /// Compute the final shot noise term out of four.
  SphericalBesselCalculator sj1(params.ell1), sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      /// Check for vanishing cases where all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      if (flag_vanishing == "false") {
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_config_space(
            params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
          );
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_config_space(
            params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
          );
      }

      /// Compute S_{ell_1 ell_2 L; i = j ≠ k} in eq. (45) in the Paper.
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        Pseudo2ptStats<ParticleCatalogue> stats(params);

        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        PseudoDensityField<ParticleCatalogue> N_LM00(params);  // N_LM or N_00
        if (los_choice == 2) {
          N_LM00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            0, 0
          );
        } else {
          N_LM00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        }
        N_LM00.fourier_transform();

        PseudoDensityField<ParticleCatalogue> dn_LM00_for_sn(params);
          // dn_LM or dn_00
        if (los_choice == 2) {
          dn_LM00_for_sn.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM00_for_sn.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM00_for_sn.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in the Paper.
        std::complex<double> barS_LM =
          stats.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );  // barS_LM

        /// Compute on mesh grids.
        fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
        trv::runtime::gbytesMem += double(params.nmesh)
          * sizeof(fftw_complex) / BYTES_PER_GBYTES;
        for (int gid = 0; gid < params.nmesh; gid++) {
          three_pt_holder[gid][0] = 0.;
          three_pt_holder[gid][1] = 0.;
        }

        stats.compute_2pt_self_shotnoise_for_bispec_meshgrid(
          dn_LM00_for_sn, N_LM00, barS_LM, three_pt_holder
          // params.ELL, M_
        );

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_a;
          double kmag_b = kbin[i_kbin];
          if (params.form == "diag") {
            kmag_a = kmag_b;
          } else
          if (params.form == "full") {
            kmag_a = kbin[params.ith_kbin];
          }

          double rvec[3];
          std::complex<double> S_ij_k = 0.;
          for (int i = 0; i < params.ngrid[0]; i++) {
            for (int j = 0; j < params.ngrid[1]; j++) {
              for (int k = 0; k < params.ngrid[2]; k++) {
                long long idx_grid =
                  (i * params.ngrid[1] + j) * params.ngrid[2] + k;

                /// This conforms to the absurd FFT array ordering
                /// convention that negative wavenumbers/frequencies come
                /// after zero and positive wavenumbers/frequencies.
                rvec[0] = (i < params.ngrid[0]/2) ?
                  i * dr[0] : (i - params.ngrid[0]) * dr[0];
                rvec[1] = (j < params.ngrid[1]/2) ?
                  j * dr[1] : (j - params.ngrid[1]) * dr[1];
                rvec[2] = (k < params.ngrid[2]/2) ?
                  k * dr[2] : (k - params.ngrid[2]) * dr[2];
                double rmag = std::sqrt(
                  rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
                );

                double j1 = sj1.eval(kmag_a * rmag);
                double j2 = sj2.eval(kmag_b * rmag);
                std::complex<double> S_ij_k_gridpt(
                  three_pt_holder[idx_grid][0],
                  three_pt_holder[idx_grid][1]
                );

                S_ij_k += j1 * j2 * ylm_a[idx_grid] * ylm_b[idx_grid]
                  * S_ij_k_gridpt;
              }
            }
          }

          S_ij_k *= vol_cell * std::pow(M_I, params.ell1 + params.ell2);

          sn_save[i_kbin] += coupling * S_ij_k;

          if (trv::runtime::currTask == 0) {
            std::printf(
              "[%s STAT] Shot noise term (last out of 4) at "
              "wavenumber k2 = %.4f and order (m1, m2, M) = (%d, %d, %d) "
              "computed.\n",
              trv::runtime::show_timestamp().c_str(),
              kmag_b, m1_, m2_, M_
            );
          }
        }

        fftw_free(three_pt_holder); three_pt_holder = NULL;
        trv::runtime::gbytesMem -= double(params.nmesh)
          * sizeof(fftw_complex) / BYTES_PER_GBYTES;
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed shot noise terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] Computing bispectrum terms...\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /// Compute bispectrum terms.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      /// Check for vanishing cases where all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// QUEST: Escape early?
      // if (flag_vanishing == "true") {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      trv::runtime::gbytesMem += 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;

      /// QUEST: This if-statement is redundant if early escape?
      if (flag_vanishing == "false") {
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_fourier_space(
            params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
          );
        SphericalHarmonicCalculator::
          store_reduced_spherical_harmonic_in_fourier_space(
            params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
          );
      }

      /// Compute a single term.
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (std::fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        /// Compute G_LM in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> dn_LM00_for_los1(params);
          // dn_LM or dn_00
        if (los_choice == 0) {
          dn_LM00_for_los1.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM00_for_los1.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM00_for_los1.fourier_transform();

        PseudoDensityField<ParticleCatalogue> dn_LM00_for_los2(params);
          // dn_LM or dn_00
        if (los_choice == 1) {
          dn_LM00_for_los2.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM00_for_los2.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM00_for_los2.fourier_transform();

        PseudoDensityField<ParticleCatalogue> dn_LM00_for_los3(params);
          // dn_LM or dn_00
        if (los_choice == 2) {
          dn_LM00_for_los3.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM00_for_los3.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM00_for_los3.fourier_transform();
        dn_LM00_for_los3.apply_assignment_compensation();
        dn_LM00_for_los3.inv_fourier_transform();

        /// Compute F_lm's in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> F_lm_a(params);  // F_lm_a
        double kmag_a, kmag_b;
        double dk = kbin[1] - kbin[0];  // HACK: only for linear binning
                                        // TODO: implement custom binning

        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_LM00_for_los1, ylm_a, kmag_a, dk
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          kmag_b = kbin[i_kbin];

          PseudoDensityField<ParticleCatalogue> F_lm_b(params);  // F_lm_b
          F_lm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_LM00_for_los2, ylm_b, kmag_b, dk
          );

          if (params.form == "diag") {
            kmag_a = kmag_b;
            F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
              dn_LM00_for_los1, ylm_a, kmag_a, dk
            );
          }

          /// Add grid contribution.
          std::complex<double> bk_sum = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(
              dn_LM00_for_los3[gid][0], dn_LM00_for_los3[gid][1]
            );
            bk_sum += F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;
          }

          bk_save[i_kbin] += coupling * vol_cell * bk_sum;

          if (trv::runtime::currTask == 0) {
            std::printf(
              "[%s STAT] Bispectrum term at wavenumber k2 = %.4f and order "
              "(m1, m2, M) = (%d, %d, %d) computed.\n",
              trv::runtime::show_timestamp().c_str(),
              kmag_b, m1_, m2_, M_
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      trv::runtime::gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  if (trv::runtime::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed bispectrum terms.\n",
      trv::runtime::show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add bispectrum effective binning.
  BispecMeasurements bispec_out;
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    if (params.form == "diag") {
      bispec_out.kbin1.push_back(kbin[ibin]);
      bispec_out.keff1.push_back(kbin[ibin]);
    } else
    if (params.form == "full") {
      bispec_out.kbin1.push_back(kbin[params.ith_kbin]);
      bispec_out.keff1.push_back(kbin[params.ith_kbin]);
    }
    bispec_out.kbin2.push_back(kbin[ibin]);
    bispec_out.keff2.push_back(kbin[ibin]);
    bispec_out.ntriag.push_back(0);
    bispec_out.bk_raw.push_back(norm * bk_save[ibin]);
    bispec_out.bk_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      std::sprintf(
        save_filepath, "%s/bk%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      std::sprintf(
        save_filepath, "%s/bk%d%d%d_kbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_kbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      std::fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        bispec_out.kbin1[ibin], bispec_out.keff1[ibin],
        bispec_out.kbin2[ibin], bispec_out.keff2[ibin],
        bispec_out.ntriag[ibin],
        bispec_out.bk_raw[ibin].real(), bispec_out.bk_raw[ibin].imag(),
        bispec_out.bk_shot[ibin].real(), bispec_out.bk_shot[ibin].imag()
      );
    }
    std::fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] k1_save; delete[] k2_save;
  delete[] bk_save; delete[] sn_save;

  return bispec_out;
}
#endif  // TRIUMVIRATE_USE_DISABLED_CODE

}  // trv::algo::
}  // trv::

#endif  // TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_
