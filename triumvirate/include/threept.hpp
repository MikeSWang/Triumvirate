/**
 * @file threept.hpp
 * @brief Three-point correlator computations.
 *
 * @note Unless otherwise specified, the Paper hereafter in this method's
 *       comments refers to Sugiyama et al. (2018) [arXiv:1803.02132].
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_

#include "common.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "particles.hpp"
#include "field.hpp"
#include "twopt.hpp"

const std::complex<double> M_I(0., 1.);  ///< imaginary unit
const double EPS_COUPLING_3PT = 1.e-10;  /**< zero-tolerance for three-point
                                              coupling coefficients */
                                         // CAVEAT: discretionary choice

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
  ParticleCatalogue& catalogue, ParameterSet& params, double alpha=1.
) {
  PseudoDensityField<ParticleCatalogue> catalogue_mesh(params);

  double norm_factor = catalogue_mesh._calc_wgt_cu_volume_norm(catalogue)
    / pow(alpha, 3);

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
 * Calculate bispectrum from paired catalogues and
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
BispecMeasurements calc_bispec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  ParameterSet& params,
  double* kbin,
  double alpha,
  double norm,
  bool save=false
) {
  if (currTask == 0) {
    printf(
      "[%s STAT] Measuring bispectrum from data and random catalogues.\n",
      show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[%s ERRO] Specified bispectrum multipole vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        show_timestamp().c_str()
      );
    }
    exit(1);
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

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing shot noise terms...\n",
      show_timestamp().c_str()
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
        if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

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
              sn_save[ibin] += coupling * stats_sn.pk[ibin];
            }
          } else
          if (params.form == "full") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              sn_save[ibin] += coupling * stats_sn.pk[params.ith_kbin];
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
            sn_save[ibin] += coupling * stats_sn.pk[ibin];
          }
        }

        if (currTask == 0) {
          printf(
            "[%s STAT] Shot noise terms (3 out of 4) at order "
            "(m1, m2, M) = (%d, %d, %d)` computed.\n",
            show_timestamp().c_str(),
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
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * double(params.nmesh)
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

      /// Compute S_{\ell_1 \ell_2 L; i = j ≠ k} in eq. (45) in the Paper.
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
  			Pseudo2ptStats<ParticleCatalogue> stats_sn(params);

        /// Calculate Wigner-3j coupling coefficient.
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

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
        gbytesMem += double(params.nmesh)
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
                double rmag = sqrt(
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

          S_ij_k *= vol_cell * pow(M_I, params.ell1 + params.ell2);

          sn_save[i_kbin] += coupling * S_ij_k;
        }

        if (currTask == 0) {
          printf(
            "[%s STAT] Shot noise terms (last out of 4) at order "
            "(m1, m2, M) = (%d, %d, %d)` computed.\n",
            show_timestamp().c_str(),
            m1_, m2_, M_
          );
        }

        fftw_free(three_pt_holder); three_pt_holder = NULL;
        gbytesMem -= double(params.nmesh)
          * sizeof(fftw_complex) / BYTES_PER_GBYTES;
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed shot noise terms.\n",
      show_timestamp().c_str()
    );
  }

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing bispectrum terms...\n",
      show_timestamp().c_str()
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
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// QUEST: Escape early?
      // if (flag_vanishing == "true") {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * double(params.nmesh)
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
        if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        /// Compute G_LM in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> dn_LM(params);  // dn_LM
        dn_LM.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM.fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.inv_fourier_transform();

        /// Calculate F_lm's in eq. (42) in the Paper.
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

          if (currTask == 0) {
            printf(
              "[%s STAT] Bispectrum term at wavenumber k2 = %.4f and order "
              "(m1, m2, M) = (%d, %d, %d) computed.\n",
              show_timestamp().c_str(),
              kmag_b, m1_, m2_, M_
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  dn_00.finalise_density_field();  // ~dn_00

  // /// Recalculate normalisation.
  /// double norm = pow(norm, 2) / pow(particles_data.wtotal, 3);

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed bispectrum terms.\n",
      show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add bispectrum effective binning.
  BispecMeasurements bispec_out;
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    if (params.form == "diag") {
      bispec_out.kbin1.push_back(kbin[ibin]);
    } else
    if (params.form == "full") {
      bispec_out.kbin1.push_back(kbin[params.ith_kbin]);
    }
    // bispec_out.keff1.push_back(k1_save[ibin]);
    bispec_out.kbin2.push_back(kbin[ibin]);
    // bispec_out.keff2.push_back(k2_save[ibin]);
    // bispec_out.ntriag.push_back(ntriag_save[ibin]);
    bispec_out.bk_raw.push_back(norm * bk_save[ibin]);
    bispec_out.bk_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      sprintf(
        save_filepath, "%s/bk%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      sprintf(
        save_filepath, "%s/bk%d%d%d_kbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_kbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    FILE* save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        bispec_out.kbin1[ibin], bispec_out.keff1[ibin],
        bispec_out.kbin2[ibin], bispec_out.keff2[ibin],
        bispec_out.ntriag[ibin],
        bispec_out.bk_raw[ibin].real(), bispec_out.bk_raw[ibin].imag(),
        bispec_out.bk_shot[ibin].real(), bispec_out.bk_shot[ibin].imag()
      );
    }
    fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] k1_save; delete[] k2_save;
  delete[] bk_save; delete[] sn_save;

  return bispec_out;
}

/**
 * Calculate bispectrum in a periodic box and optionally save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param params Parameter set.
 * @param kbin Wavenumber bins.
 * @param norm Normalisation factor.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns bispec_out Output bispectrum measurements.
 */
BispecMeasurements calc_bispec_in_box(
  ParticleCatalogue& particles_data,
  ParameterSet& params,
  double* kbin, double norm,
  bool save=false
) {
  if (currTask == 0) {
    printf(
      "[%s STAT] Measurement: bispectrum in a periodic box.\n",
      show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[%s ERRO] Specified bispectrum multipole vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        show_timestamp().c_str()
      );
    }
    exit(1);
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

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing shot noise terms...\n",
      show_timestamp().c_str()
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
      if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

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
            sn_save[ibin] += coupling * stats_sn.pk[ibin];
          }
        } else
        if (params.form == "full") {
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            sn_save[ibin] += coupling * stats_sn.pk[params.ith_kbin];
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
          sn_save[ibin] += coupling * stats_sn.pk[ibin];
        }
      }

      if (currTask == 0) {
        printf(
          "[%s STAT] Shot noise terms (3 out of 4) at order "
          "(m1, m2, M) = (%d, %d, %d)` computed.\n",
          show_timestamp().c_str(),
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
      if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * double(params.nmesh)
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
        gbytesMem += double(params.nmesh)
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
              double rmag = sqrt(
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

        S_ij_k *= vol_cell * pow(M_I, params.ell1 + params.ell2);

        sn_save[i_kbin] += coupling * S_ij_k;
      }

      if (currTask == 0) {
        printf(
          "[%s STAT] Shot noise terms (last out of 4) at order "
          "(m1, m2, M) = (%d, %d, %d)` computed.\n",
          show_timestamp().c_str(),
          m1_, m2_, M_
        );
      }

      fftw_free(three_pt_holder); three_pt_holder = NULL;
        gbytesMem -= double(params.nmesh)
          * sizeof(fftw_complex) / BYTES_PER_GBYTES;

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed shot noise terms.\n",
      show_timestamp().c_str()
    );
  }

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing bispectrum terms...\n",
      show_timestamp().c_str()
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
      if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * double(params.nmesh)
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
      PseudoDensityField<ParticleCatalogue> dn_00(params);
      dn_00.compute_unweighted_fluctuation_insitu(
        particles_data, params.volume
      );
      dn_00.fourier_transform();
      dn_00.apply_assignment_compensation();
      dn_00.inv_fourier_transform();

      /// Calculate F_lm's in eq. (42) in the Paper.
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
          std::complex<double> G_00_gridpt(dn_00[gid][0], dn_00[gid][1]);
          bk_sum += F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;
        }

        bk_save[i_kbin] += coupling * vol_cell * bk_sum;

        if (currTask == 0) {
          printf(
            "[%s STAT] Bispectrum term at wavenumber k2 = %.4f and order "
            "(m1, m2, M) = (%d, %d, %d) computed.\n",
            show_timestamp().c_str(),
            kmag_b, m1_, m2_, M_
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  dn_00.finalise_density_field();  // ~dn_00

  // /// Recalculate normalisation.
  // double norm = pow(params.volume, 2) / pow(double(particles_data.ntotal), 3);

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed bispectrum terms.\n",
      show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add bispectrum effective binning.
  BispecMeasurements bispec_out;
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    if (params.form == "diag") {
      bispec_out.kbin1.push_back(kbin[ibin]);
    } else
    if (params.form == "full") {
      bispec_out.kbin1.push_back(kbin[params.ith_kbin]);
    }
    // bispec_out.keff1.push_back(k1_save[ibin]);
    bispec_out.kbin2.push_back(kbin[ibin]);
    // bispec_out.keff2.push_back(k2_save[ibin]);
    // bispec_out.ntriag.push_back(ntriag_save[ibin]);
    bispec_out.bk_raw.push_back(norm * bk_save[ibin]);
    bispec_out.bk_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      sprintf(
        save_filepath, "%s/bk%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      sprintf(
        save_filepath, "%s/bk%d%d%d_kbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_kbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    FILE* save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      fprintf(
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
    fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] k1_save; delete[] k2_save;
  delete[] bk_save; delete[] sn_save;

  return bispec_out;
}

/**
 * Calculate three-point correlation function from paired catalogues and
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
ThreePCFMeasurements calc_3pcf(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  ParameterSet& params,
  double* rbin,
  double alpha,
  double norm,
  bool save=false
) {
  if (currTask == 0) {
    printf(
      "[%s STAT] Measuring three-point correlation function "
      "from data and random catalogues.\n",
      show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[%s ERRO] Specified three-point correlation function multipole"
        "vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        show_timestamp().c_str()
      );
    }
    exit(1);
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

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing shot noise terms...\n",
      show_timestamp().c_str()
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
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * double(params.nmesh)
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
        if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

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

        if (currTask == 0) {
          printf(
            "[%s STAT] Shot noise term at order "
            "(m1, m2, M) = (%d, %d, %d)` computed.\n",
            show_timestamp().c_str(),
            m1_, m2_, M_
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed shot noise terms.\n",
      show_timestamp().c_str()
    );
  }

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing three-point correlation function terms...\n",
      show_timestamp().c_str()
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
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// QUEST: Escape early?
      // if (flag_vanishing == "true") {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * double(params.nmesh)
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
        if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        /// Compute G_LM in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> dn_LM(params);  // dn_LM
        dn_LM.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM.fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.inv_fourier_transform();

        /// Calculate F_lm's in eq. (49) in the Paper.
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

          zeta_save[i_rbin] += pow(M_I, params.ell1 + params.ell2)
            * coupling * vol_cell * zeta_sum;

          if (currTask == 0) {
            printf(
              "[%s STAT] Three-point correlation function term at "
              "separation r2 = %.3f and order (m1, m2, M) = (%d, %d, %d) "
              "computed.\n",
              show_timestamp().c_str(),
              rmag_b, m1_, m2_, M_
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  dn_00.finalise_density_field();  // ~dn_00

  // /// Recalculate normalisation.
  /// double norm = pow(norm, 2) / pow(particles_data.ntotal, 3);

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed three-point correlation function terms.\n",
      show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add 3PCF effective binning.
  ThreePCFMeasurements threepcf_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    if (params.form == "diag") {
      threepcf_out.rbin1.push_back(rbin[ibin]);
    } else
    if (params.form == "full") {
      threepcf_out.rbin1.push_back(rbin[params.ith_rbin]);
    }
    // bispec_out.reff1.push_back(r1_save[ibin]);
    threepcf_out.rbin2.push_back(rbin[ibin]);
    // bispec_out.reff2.push_back(r2_save[ibin]);
    // bispec_out.ntriag.push_back(ntriag_save[ibin]);
    threepcf_out.zeta_raw.push_back(norm * zeta_save[ibin]);
    threepcf_out.zeta_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      sprintf(
        save_filepath, "%s/zeta%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      sprintf(
        save_filepath, "%s/zeta%d%d%d_rbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_rbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    FILE* save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        threepcf_out.rbin1[ibin], threepcf_out.reff1[ibin],
        threepcf_out.rbin2[ibin], threepcf_out.reff2[ibin],
        threepcf_out.ntriag[ibin],
        threepcf_out.zeta_raw[ibin].real(), threepcf_out.zeta_raw[ibin].imag(),
        threepcf_out.zeta_shot[ibin].real(), threepcf_out.zeta_shot[ibin].imag()
      );
    }
    fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  return threepcf_out;
}

/**
 * Calculate three-point correlation function in a periodic box
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
ThreePCFMeasurements calc_3pcf_in_box(
  ParticleCatalogue& particles_data,
  ParameterSet& params,
  double* rbin, double norm,
  bool save=false
) {
  if (currTask == 0) {
    printf(
      "[%s STAT] Measurement: three-point correlation function "
      "in a periodic box.\n",
      show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[%s ERRO] Specified three-point correlation function multipole "
        "vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        show_timestamp().c_str()
      );
    }
    exit(1);
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

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing shot noise terms...\n",
      show_timestamp().c_str()
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
      if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * double(params.nmesh)
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

      if (currTask == 0) {
        printf(
          "[%s STAT] Shot noise term at order "
          "(m1, m2, M) = (%d, %d, %d)` computed.\n",
          show_timestamp().c_str(),
          m1_, m2_, M_
        );
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed shot noise terms.\n",
      show_timestamp().c_str()
    );
  }

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing bispectrum terms...\n",
      show_timestamp().c_str()
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
      if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * double(params.nmesh)
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
      PseudoDensityField<ParticleCatalogue> dn_00(params);
      dn_00.compute_unweighted_fluctuation_insitu(
        particles_data, params.volume
      );
      dn_00.fourier_transform();
      dn_00.apply_assignment_compensation();
      dn_00.inv_fourier_transform();

      /// Calculate F_lm's in eq. (42) in the Paper.
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
          std::complex<double> G_00_gridpt(dn_00[gid][0], dn_00[gid][1]);
          zeta_sum += F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;
        }

        zeta_save[i_rbin] += pow(M_I, params.ell1 + params.ell2)
          * coupling * vol_cell * zeta_sum;

        if (currTask == 0) {
          printf(
            "[%s STAT] Three-point correlation function term at "
            "separation r2 = %.3f and order (m1, m2, M) = (%d, %d, %d) "
            "computed.\n",
            show_timestamp().c_str(),
            rmag_b, m1_, m2_, M_
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  dn_00.finalise_density_field();  // ~dn_00

  // /// Recalculate normalisation.
  /// double norm = pow(norm, 2) / pow(particles_data.ntotal, 3);

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed three-point correlation function terms.\n",
      show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add 3PCF effective binning.
  ThreePCFMeasurements threepcf_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    if (params.form == "diag") {
      threepcf_out.rbin1.push_back(rbin[ibin]);
    } else
    if (params.form == "full") {
      threepcf_out.rbin1.push_back(rbin[params.ith_rbin]);
    }
    // bispec_out.reff1.push_back(r1_save[ibin]);
    threepcf_out.rbin2.push_back(rbin[ibin]);
    // bispec_out.reff2.push_back(r2_save[ibin]);
    // bispec_out.ntriag.push_back(ntriag_save[ibin]);
    threepcf_out.zeta_raw.push_back(norm * zeta_save[ibin]);
    threepcf_out.zeta_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      sprintf(
        save_filepath, "%s/zeta%d%d%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      sprintf(
        save_filepath, "%s/zeta%d%d%d_rbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_rbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    FILE* save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      fprintf(
        save_fileptr,
        "%.9e \t %.9e \t %.9e \t %.9e \t %d \t %.9e \t %.9e \t %.9e \t %.9e\n",
        threepcf_out.rbin1[ibin], threepcf_out.reff1[ibin],
        threepcf_out.rbin2[ibin], threepcf_out.reff2[ibin],
        threepcf_out.ntriag[ibin],
        threepcf_out.zeta_raw[ibin].real(), threepcf_out.zeta_raw[ibin].imag(),
        threepcf_out.zeta_shot[ibin].real(), threepcf_out.zeta_shot[ibin].imag()
      );
    }
    fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  return threepcf_out;
}

/**
 * Calculate three-point correlation function window from a random
 * catalogue and optionally save the results.
 *
 * @param particles_rand (Random-source) particle container.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbin Separation bins.
 * @param alpha Alpha ratio.
 * @param norm Normalisation factor.
 * @param save If `true` (default is `false`), write computed results
 *             to the measurement output file set by `params`.
 * @returns threepcfwin_out Output three-point correlation function
 *                          window measurements.
 */
ThreePCFWindowMeasurements calc_3pcf_window(
  ParticleCatalogue& particles_rand,
  LineOfSight* los_rand,
  ParameterSet& params,
  double* rbin,
  double alpha,
  double norm,
  bool save=false
) {
  if (currTask == 0) {
    printf(
      "[%s STAT] Measurement: three-point correlation function window "
      "from random catalogue.\n",
      show_timestamp().c_str()
    );
  }

  /* * Set-up ************************************************************** */

  /// Set up/check input.
  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0))
    < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[%s ERRO] Specified three-point correlation function window "
        "multipole vanishes identically. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n",
        show_timestamp().c_str()
      );
    }
    exit(1);
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

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing shot noise terms...\n",
      show_timestamp().c_str()
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
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
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
        if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

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

        if (currTask == 0) {
          printf(
            "[%s STAT] Shot noise term at order "
            "(m1, m2, M) = (%d, %d, %d)` computed.\n",
            show_timestamp().c_str(),
            m1_, m2_, M_
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  N_00.finalise_density_field();  // ~N_00

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed shot noise terms.\n",
      show_timestamp().c_str()
    );
  }

  if (currTask == 0) {
    printf(
      "[%s STAT] Computing three-point correlation function window terms...\n",
      show_timestamp().c_str()
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
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_vanishing = "false";
        }
      }

      /// QUEST: Escape early?
      // if (flag_vanishing == "true") {continue;}

      /// Initialise/reset spherical harmonic grids.
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
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
        if (fabs(coupling) < EPS_COUPLING_3PT) {continue;}

        /// Compute G_LM in eq. (42) in the Paper.
        PseudoDensityField<ParticleCatalogue> n_LM(params);  // n_LM
        n_LM.compute_ylm_wgtd_density(
          particles_rand, los_rand, alpha, params.ELL, M_
        );
        n_LM.fourier_transform();
        n_LM.apply_assignment_compensation();
        n_LM.inv_fourier_transform();

        /// Calculate F_lm's in eq. (49) in the Paper.
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

          zeta_save[i_rbin] += pow(M_I, params.ell1 + params.ell2)
            * coupling * vol_cell * zeta_sum;

          if (currTask == 0) {
            printf(
              "[%s STAT] Three-point correlation function window term at "
              "separation r2 = %.3f and order (m1, m2, M) = (%d, %d, %d) "
              "computed.\n",
              show_timestamp().c_str(),
              rmag_b, m1_, m2_, M_
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * double(params.nmesh)
        * sizeof(std::complex<double>) / BYTES_PER_GBYTES;
    }
  }

  n_00.finalise_density_field();  // ~n_00

  // /// Recalculate normalisation.
  /// double norm = pow(norm, 2) / pow(particles_data.wtotal, 3);

  if (currTask == 0) {
    printf(
      "[%s STAT] ... computed three-point correlation function window terms.\n",
      show_timestamp().c_str()
    );
  }

  /* * Output ************************************************************** */

  /// Fill in output struct.
  /// FIXME: Add 3PCF effective binning.
  ThreePCFWindowMeasurements threepcfwin_out;
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    if (params.form == "diag") {
      threepcfwin_out.rbin1.push_back(rbin[ibin]);
    } else
    if (params.form == "full") {
      threepcfwin_out.rbin1.push_back(rbin[params.ith_rbin]);
    }
    // bispec_out.reff1.push_back(r1_save[ibin]);
    threepcfwin_out.rbin2.push_back(rbin[ibin]);
    // bispec_out.reff2.push_back(r2_save[ibin]);
    // bispec_out.ntriag.push_back(ntriag_save[ibin]);
    threepcfwin_out.zeta_raw.push_back(norm * zeta_save[ibin]);
    threepcfwin_out.zeta_shot.push_back(norm * sn_save[ibin]);
  }

  /// Save (optionally) to output file.
  if (save) {
    /// Set output path.
    char save_filepath[1024];
    if (params.form == "diag") {
      sprintf(
        save_filepath, "%s/zeta%d%d%d_window%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "full") {
      sprintf(
        save_filepath, "%s/zeta%d%d%d_window_rbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.ith_rbin,
        params.output_tag.c_str()
      );
    }

    /// Write output.
    FILE* save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      fprintf(
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
    fclose(save_fileptr);
  }

  delete[] ntriag_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  return threepcfwin_out;
}

/**
 * Calculate three-point correlation function window from random
 * catalogues and save the results.  This function parallelises the
 * tasks into different processes.
 *
 * @param particles_rand (Random-source) particle container.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pcf_window_mpi(
  ParticleCatalogue& particles_rand,
  LineOfSight* los_rand,
  ParameterSet& params,
  double* rbin,
  double alpha,
  double norm
) {  // WARNING: inherited from Sugiyama et al. without matching equation
  if (currTask == 0) {
    printf(
      "[STAT] :: Measuring three-point correlation function window "
      "from random catalogues.\n"
    );
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[ERRO] :: Disallowed multipole degree combination "
        "for three-point correlation function window measurements. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
      );
    }
    exit(1);
  }

  int n_temp = 10;  // CAVEAT: discretionary choice
  int NR = 3;  // CAVEAT: discretionary choice

  params.ith_rbin = currTask;

  /// Set up binning.
  /// CAVEAT: Discretionary choices.
  rbin[0] = 0.;
  rbin[1] = 1.;
  rbin[2] = 10.;
  rbin[3] = 20.;
  rbin[4] = 30.;
  rbin[5] = 40.;
  rbin[6] = 50.;
  rbin[7] = 60.;
  double rmin = 70.;
  double dlnr = (log(params.rmax) - log(rmin))
    / double((params.num_rbin - 8) - 1);

  for (int ibin = 8; ibin < params.num_rbin; ibin++) {
    rbin[ibin] = rmin * exp(dlnr * (ibin - 8));
  }

  /// Initialise output shot noise terms.
  if (currTask == 0) {
    printf(
      "[STAT] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* shotnoise_save = new std::complex<double>[n_temp];
  for (int ibin = 0; ibin < n_temp; ibin++) {
    shotnoise_save[ibin] = 0.;
  }

  /// Compute shot noise terms.
  PseudoDensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
  shotnoise_quadratic_00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
    particles_rand, los_rand, alpha, 0, 0
  );
  shotnoise_quadratic_00.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          *  wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        PseudoDensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
        dn_LM_for_shotnoise.compute_ylm_wgtd_density(
          particles_rand, los_rand, alpha, params.ELL, M_
        );
        dn_LM_for_shotnoise.fourier_transform();

        /// QUEST: Check the power index in isotropic shot noise.
        Pseudo2ptStats<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise =
          stats.calc_ylm_wgtd_shotnoise_for_powspec(
            particles_rand, los_rand, alpha, params.ELL, M_
          );

        stats.compute_uncoupled_shotnoise_for_3pcf(
          dn_LM_for_shotnoise, shotnoise_quadratic_00,
          ylm_a, ylm_b,
          shotnoise,
          rbin
          // params.ell1, m1_
        );

        for (int idx = 0; idx < n_temp; idx++) {
          if (params.form == "diag") {
            shotnoise_save[idx] += coupling * stats.xi[idx + NR * n_temp];
          } else if (params.form == "full") {
            /// Calculate shot noise contribution equivalent to
            /// the Kronecker delta.  See `calc_3pcf`.
            if (idx + NR * n_temp == params.ith_rbin) {
              shotnoise_save[idx] += coupling * stats.xi[idx + NR * n_temp];
            } else {
              shotnoise_save[idx] += 0.;
            }
          }
        }

        if (currTask == 0) {
          printf(
            "[STAT] :: Computed shot noise term for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
            m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;
    }
  }

  shotnoise_quadratic_00.finalise_density_field();

  if (currTask == 0) {
    printf(
      "[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Initialise output three-point correlation function window.
  if (currTask == 0) {
    printf(
      "[STAT] :: Computing three-point correlation function window terms "
      "(%.3f seconds elapsed...).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* zeta_save = new std::complex<double>[n_temp];
  for (int idx = 0; idx < n_temp; idx++) {
    zeta_save[idx] = 0.;
  }

  /// Compute three-point correlation function window.
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_density(particles_rand, los_rand, alpha, 0, 0);
  dn_00.fourier_transform();

  SphericalBesselCalculator sj1(params.ell1);
  SphericalBesselCalculator sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        /// Calculate G_{LM}.  See `calc_3pcf`.
        PseudoDensityField<ParticleCatalogue> dn_LM(params);
        dn_LM.compute_ylm_wgtd_density(
          particles_rand, los_rand, alpha, params.ELL, M_
        );
        dn_LM.fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.inv_fourier_transform();

        /// Calculate F_{\ell m}.  See `calc_3pcf`.
        PseudoDensityField<ParticleCatalogue> F_lm_a(params);
        double rmag_a;
        if (params.form == "full") {
          rmag_a = rbin[params.ith_rbin];
          F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            dn_00, ylm_a, sj1, rmag_a
          );
        }

        for (int i_rbin = 0; i_rbin < n_temp; i_rbin++) {
          double rmag_b = rbin[i_rbin + NR * n_temp];

          if (params.form == "diag") {
            rmag_a = rmag_b;
            F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
              dn_00, ylm_a, sj1, rmag_a
            );
          }

          PseudoDensityField<ParticleCatalogue> F_lm_b(params);
          F_lm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            dn_00, ylm_b, sj2, rmag_b
          );

          std::complex<double> zeta_sum = 0.;
          double factor = params.volume / double(params.nmesh);
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_ellm_1(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_ellm_2(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM(dn_LM[gid][0], dn_LM[gid][1]);
            zeta_sum += pow(M_I, params.ell1 + params.ell2) * factor
              * F_ellm_1 * F_ellm_2 * G_LM;
          }

          zeta_save[i_rbin] += coupling * zeta_sum;

          if (currTask == 0) {
            printf(
              "[STAT] :: Computed three-point correlation function window "
              "term for separation and orders "
              "`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
              "(%.3f seconds elapsed).\n",
              rmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;
    }
  }

  if (currTask == 0) {
    printf(
      "[STAT] :: Computed three-point correlation function window terms "
      "(... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  FILE* save_fileptr;
  char save_filepath[1024];
  if (params.form == "diag") {
    sprintf(
      save_filepath, "%s/zeta_window_%d%d%d_%d%s",
      params.measurement_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      NR,
      params.output_tag.c_str()
    );
    save_fileptr = fopen(save_filepath, "w");
    for (int idx = 0; idx < n_temp; idx++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
        rbin[idx + NR * n_temp], rbin[idx + NR * n_temp],
        norm * (zeta_save[idx].real() - shotnoise_save[idx].real()),
        norm * shotnoise_save[idx].real()
      );
    }
  } else if (params.form == "full") {
    sprintf(
      save_filepath, "%s/zeta_window_%d%d%d_rbin%02d_%d%s",
      params.measurement_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      params.ith_rbin, NR,
      params.output_tag.c_str()
    );
    save_fileptr = fopen(save_filepath, "w");
    for (int idx = 0; idx < n_temp; idx++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
        rbin[params.ith_rbin], rbin[idx + NR * n_temp],
        norm * (zeta_save[idx].real() - shotnoise_save[idx].real()),
        norm * shotnoise_save[idx].real()
      );
    }
  }
  fclose(save_fileptr);

  delete[] shotnoise_save;
  delete[] zeta_save;

  return 0;
}

/**
 * Calculate three-point correlation function window wide-angle correction
 * terms for three-point correlation function from random catalogues and
 * save the results.  This function uses logarithmic binning except on
 * the smallest scales, where discretionary bins are used.
 *
 * @param particles_rand (Random-source) particle container.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pcf_window_for_wide_angle(
  ParticleCatalogue& particles_rand,
  LineOfSight* los_rand,
  ParameterSet& params,
  double* rbin,
  double alpha,
  double norm
) {
  if (currTask == 0) {
    printf(
      "[STAT] :: Measuring three-point correlation function window wide-angle"
      "correction terms for three-point correlation function "
      "from random catalogues.\n"
    );
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[ERRO] :: Disallowed multipole degree combination "
        "for three-point correlation function window measurements. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
      );
    }
    exit(1);
  }

  /// Initialise output shot noise terms.
  if (currTask == 0) {
    printf(
      "[STAT] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* shotnoise_save =
    new std::complex<double>[params.num_rbin];
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    shotnoise_save[ibin] = 0.;
  }

  /// Compute shot noise terms.
  PseudoDensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
  shotnoise_quadratic_00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
    particles_rand, los_rand, alpha, 0, 0
  );
  shotnoise_quadratic_00.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          *  wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        PseudoDensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
        dn_LM_for_shotnoise.compute_ylm_wgtd_density(
          particles_rand, los_rand, alpha, params.ELL, M_
        );
        dn_LM_for_shotnoise.fourier_transform();

        /// QUEST: Check the power index in isotropic shot noise.
        Pseudo2ptStats<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise =
          stats.calc_ylm_wgtd_shotnoise_for_powspec(
            particles_rand, los_rand, alpha, params.ELL, M_
          );

        stats.compute_uncoupled_shotnoise_for_3pcf(
          dn_LM_for_shotnoise, shotnoise_quadratic_00,
          ylm_a, ylm_b,
          shotnoise,
          rbin
          // params.ell1, m1_
        );

        for (int ibin = 0; ibin < params.num_rbin; ibin++) {
          if (params.form == "diag") {
            shotnoise_save[ibin] += coupling * stats.xi[ibin];
          } else if (params.form == "full") {
            /// Calculate shot noise contribution equivalent to
            /// the Kronecker delta.  See `calc_3pcf`.
            if (ibin == params.ith_rbin) {
              shotnoise_save[ibin] += coupling * stats.xi[ibin];
            } else {
              shotnoise_save[ibin] += 0.;
            }
          }
        }

        if (currTask == 0) {
          printf(
            "[STAT] :: Computed shot noise term for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
            m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
          );
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;
    }
  }

  shotnoise_quadratic_00.finalise_density_field();

  if (currTask == 0) {
    printf(
      "[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Initialise output three-point correlation function window.
  if (currTask == 0) {
    printf(
      "[STAT] :: Computing three-point correlation function window "
      "wide-angle correction terms (%.3f seconds elapsed...).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
  for (int ibin = 0; ibin < params.num_rbin; ibin++) {
    zeta_save[ibin] = 0.;
  }

  /// Compute three-point correlation function window.
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_density(particles_rand, los_rand, alpha, 0, 0);
  dn_00.fourier_transform();

  SphericalBesselCalculator sj1(params.ell1);
  SphericalBesselCalculator sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        /// Calculate x^{-i-j} G_{LM}.  See `calc_3pcf`.s
        PseudoDensityField<ParticleCatalogue> dn_LM(params);
        dn_LM.compute_ylm_wgtd_density(
          particles_rand, los_rand, alpha, params.ELL, M_
        );
        dn_LM.fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.inv_fourier_transform();
        dn_LM.apply_power_law_weight_for_wide_angle();

        /// Calculate F_{\ell m}.  See `calc_3pcf`.
        PseudoDensityField<ParticleCatalogue> F_lm_a(params);
        double rmag_a;
        if (params.form == "full") {
          rmag_a = rbin[params.ith_rbin];
          F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            dn_00, ylm_a, sj1, rmag_a
          );
        }

        for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
          double rmag_b = rbin[i_rbin];
          if (params.form == "diag") {
            rmag_a = rmag_b;
            F_lm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
              dn_00, ylm_a, sj1, rmag_a
            );
          }

          PseudoDensityField<ParticleCatalogue> F_lm_b(params);
          F_lm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
            dn_00, ylm_b, sj2, rmag_b
          );

          std::complex<double> zeta_sum = 0.;
          double factor = params.volume / double(params.nmesh);
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_ellm_1(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_ellm_2(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> x_G_LM(dn_LM[gid][0], dn_LM[gid][1]);
            zeta_sum += pow(M_I, params.ell1 + params.ell2) * factor
              * F_ellm_1 * F_ellm_2 * x_G_LM;
          }

          zeta_save[i_rbin] += coupling * zeta_sum;

          if (currTask == 0) {
            printf(
              "[STAT] :: Computed three-point correlation function window "
              "term for separation and orders "
              "`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
              "(%.3f seconds elapsed).\n",
              rmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * sizeof(std::complex<double>)
          * double(params.nmesh) / BYTES_PER_GBYTES;
    }
  }

  if (currTask == 0) {
    printf(
      "[STAT] :: Computed three-point correlation function window terms "
      "(... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Normalise and then save the output.
  double norm = pow(norm, 2) / pow(particles_rand.wtotal, 3);
  norm /= alpha * alpha * alpha;

  FILE* save_fileptr;
  char save_filepath[1024];
  if (params.form == "diag") {
    sprintf(
      save_filepath, "%s/zeta%d%d%d_window-wa%s",
      params.measurement_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      params.output_tag.c_str()
    );
    save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
        rbin[ibin], rbin[ibin],
        norm * (zeta_save[ibin].real() - shotnoise_save[ibin].real()),
        norm * shotnoise_save[ibin].real()
      );
    }
  } else if (params.form == "full") {
    sprintf(
      save_filepath, "%s/zeta%d%d%d_window-wa%d%d_rbin%02d%s",
      params.measurement_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      params.i_wa, params.j_wa,
      params.ith_rbin,
      params.output_tag.c_str()
    );
    save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
        rbin[params.ith_rbin], rbin[ibin],
        norm * (zeta_save[ibin].real() - shotnoise_save[ibin].real()),
        norm * shotnoise_save[ibin].real()
      );
    }
  }
  fclose(save_fileptr);

  delete[] shotnoise_save;
  delete[] zeta_save;

  return 0;
}

// #ifdef TRIUMVIRATE_USE_DISABLED_CODE_
/**
 * Calculate bispectrum from catalogues with respect to a choice of
 * line of sight and save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param particles_rand (Random-source) particle container.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param los Choice of line of sight.
 * @param norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec_for_los_choice(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  ParameterSet& params,
  double alpha,
  double* kbin,
  int los,
  double norm
) {  // WARNING: inherited from Sugiyama et al. without matching equation
  if (currTask == 0) {
    printf(
      "[STAT] :: Measuring bispectrum for the choice of line of sight "
      "from data and random catalogues.\n"
    );
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[ERRO] :: Disallowed multipole degree combination "
        "for bispectrum measurements. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
      );
    }
    exit(1);
  }

  /// Initialise output shot noise terms.
  if (currTask == 0) {
    printf(
      "[STAT] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* shotnoise_save =
    new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    shotnoise_save[ibin] = 0.;
  }

  /// Compute shot noise terms.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        PseudoDensityField<ParticleCatalogue> dn_for_shotnoise(params);
        if (los == 0) {
          dn_for_shotnoise.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        } else {
          dn_for_shotnoise.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        }
        dn_for_shotnoise.fourier_transform();

        /// Calculate N_LM in eq. (46) in arXiv:1803.02132.
        PseudoDensityField<ParticleCatalogue> shotnoise_quadratic(params);
        if (los == 0) {
          shotnoise_quadratic.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        } else {
          shotnoise_quadratic.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        }
        shotnoise_quadratic.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
        Pseudo2ptStats<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
          stats.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );

        /// Calculate S_{\ell_1 \ell_2 L; i = j = k} in eq. (45)
        /// in arXiv:1803.02132.
        if (params.ell1 == 0 && params.ell2 == 0) {
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            shotnoise_save[ibin] += coupling * shotnoise_cubic_LM;
          }
        }

        /// Calculate S_{\ell_1 \ell_2 L; i ≠ j = k} in eq. (45)
        /// in arXiv:1803.02132.
        if (params.ell2 == 0) {
          stats.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_for_shotnoise, shotnoise_quadratic,
            shotnoise_cubic_LM,
            kbin,
            params.ell1, m1_
          );
          if (params.form == "diag") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              shotnoise_save[ibin] += coupling * stats.pk[ibin];
            }
          } else if (params.form == "full") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              shotnoise_save[ibin] += coupling * stats.pk[params.ith_kbin];
            }
          }
        }

        /// WARNING: This was inherited from Sugiyama et al. without
        /// matching equation.  S_{\ell_1 \ell_2 L; i = k ≠ j}
        /// (i.e. ell1 == 0 case) may have been shuffled.

        if (currTask == 0) {
          printf(
            "[STAT] :: Computed shot noise term for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
            m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC);
        }
      }
    }
  }

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        PseudoDensityField<ParticleCatalogue> dn_for_shotnoise(params);
        if (los == 1) {
          dn_for_shotnoise.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        } else {
          dn_for_shotnoise.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        }
        dn_for_shotnoise.fourier_transform();

        /// Calculate N_LM in eq. (46) in arXiv:1803.02132.
        PseudoDensityField<ParticleCatalogue> shotnoise_quadratic(params);
        if (los == 1) {
          shotnoise_quadratic.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        } else {
          shotnoise_quadratic.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        }
        shotnoise_quadratic.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
        Pseudo2ptStats<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
          stats.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );

        /// WARNING: This was inherited from Sugiyama et al. without
        /// matching equation.  S_{\ell_1 \ell_2 L; i = k ≠ j}
        /// (i.e. ell1 == 0 case) may have been shuffled.

        /// Calculate S_{\ell_1 \ell_2 L; i = k ≠ j} in eq. (45)
        /// in arXiv:1803.02132.
        if (params.ell1 == 0) {
          stats.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_for_shotnoise, shotnoise_quadratic,
            shotnoise_cubic_LM,
            kbin,
            params.ell2, m2_
          );
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            shotnoise_save[ibin] += coupling * stats.pk[ibin];
          }
        }

        if (currTask == 0) {
          printf(
            "[STAT] :: Computed for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` "
            "(... %.3f seconds elapsed in total).\n",
            m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC);
        }
      }
    }
  }

  SphericalBesselCalculator sj1(params.ell1);
  SphericalBesselCalculator sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          *  wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        PseudoDensityField<ParticleCatalogue> shotnoise_quadratic(params);
        if (los == 2) {
          shotnoise_quadratic.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        } else {
          shotnoise_quadratic.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        }
        shotnoise_quadratic.fourier_transform();

        PseudoDensityField<ParticleCatalogue> dn_shotnoise(params);
        if (los == 2) {
          dn_shotnoise.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        } else {
          dn_shotnoise.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        }
        dn_shotnoise.fourier_transform();

        /// Calculate S_{\ell_1 \ell_2 L; i = j ≠ k} in eq. (45)
        /// in arXiv:1803.02132.
        Pseudo2ptStats<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
          stats.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );

        fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
        gbytesMem += sizeof(fftw_complex)
          * double(params.nmesh) / BYTES_PER_GBYTES;
        for (int ibin = 0; ibin < params.nmesh; ibin++) {
          three_pt_holder[ibin][0] = 0.;
          three_pt_holder[ibin][1] = 0.;
        }

        stats.compute_2pt_self_shotnoise_for_bispec_meshgrid(
          dn_shotnoise, shotnoise_quadratic, shotnoise_cubic_LM, three_pt_holder
          // params.ELL, M_
        );

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];
          double kmag_a;
          if (params.form == "diag") {
            kmag_a = kmag_b;
          } else if (params.form == "full") {
            kmag_a = kbin[params.ith_kbin];
          }

          double dr[3];
          dr[0] = params.boxsize[0] / double(params.ngrid[0]);
          dr[1] = params.boxsize[1] / double(params.ngrid[1]);
          dr[2] = params.boxsize[2] / double(params.ngrid[2]);

          double rvec[3];
          std::complex<double> shotnoise_sum = 0.;
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
                double rmag = sqrt(
                  rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
                );

                std::complex<double> shotnoise_grid(
                  three_pt_holder[idx_grid][0],
                  three_pt_holder[idx_grid][1]
                );
                double j1 = sj1.eval(kmag_a * rmag);
                double j2 = sj2.eval(kmag_b * rmag);

                shotnoise_sum +=
                  j1 * j2 * ylm_a[idx_grid] * ylm_b[idx_grid]
                  * shotnoise_grid;
              }
            }
          }

          double factor = params.volume / double(params.nmesh);
          shotnoise_sum *= factor * pow(M_I, params.ell1 + params.ell2);

          shotnoise_save[i_kbin] += coupling * shotnoise_sum;

          if (currTask == 0) {
            printf(
              "[STAT] :: Computed shot noise term for wavenumber and orders "
              "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
              "(%.3f seconds elapsed).\n",
              kmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
            );
          }
        }

        fftw_free(three_pt_holder); three_pt_holder = NULL;
        gbytesMem -= sizeof(fftw_complex)
          * double(params.nmesh) / BYTES_PER_GBYTES;
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;
    }
  }

  if (currTask == 0) {
    printf(
      "[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Initialise output bispectrum.
  if (currTask == 0) {
    printf(
      "[STAT] :: Computing bispectrum terms (%.3f seconds elapsed...).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    bk_save[ibin] = 0.;
  }

  /// Compute bispectrum.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          *  wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        /// Calculate G_{LM} in eq. (42) in arXiv:1803.02132.
        PseudoDensityField<ParticleCatalogue> dn_los1(params);
        if (los == 0) {
          dn_los1.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        } else {
          dn_los1.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        }
        dn_los1.fourier_transform();

        PseudoDensityField<ParticleCatalogue> dn_los2(params);
        if (los == 1) {
          dn_los2.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        } else {
          dn_los2.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        }
        dn_los2.fourier_transform();

        PseudoDensityField<ParticleCatalogue> dn_los3(params);
        if (los == 2) {
          dn_los3.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );
        } else {
          dn_los3.compute_ylm_wgtd_fluctuation(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            0, 0
          );
        }
        dn_los3.fourier_transform();
        dn_los3.apply_assignment_compensation();
        dn_los3.inv_fourier_transform();

        /// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
        PseudoDensityField<ParticleCatalogue> F_lm_a(params);
        double kmag_a;
        double dk = kbin[1] - kbin[0];
        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_los1, ylm_a, kmag_a, dk
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];

          if (params.form == "diag") {
            kmag_a = kmag_b;
            F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
              dn_los1, ylm_a, kmag_a, dk
            );
          }

          PseudoDensityField<ParticleCatalogue> F_lm_b(params);
          F_lm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_los2, ylm_b, kmag_b, dk
          );

          double factor = params.volume / double(params.nmesh);
          std::complex<double> bk_sum = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_ellm_1(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_ellm_2(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM(dn_los3[gid][0], dn_los3[gid][1]);
            bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
          }

          bk_save[i_kbin] += coupling * bk_sum;

          if (currTask == 0) {
            printf(
              "[STAT] :: Computed bispectrum term for wavenumber and orders "
              "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
              "(%.3f seconds elapsed).\n",
              kmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;
    }
  }

  if (currTask == 0) {
    printf(
      "[STAT] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Normalise and then save the output.
  double norm = pow(norm, 2) / pow(particles_data.wtotal, 3);

  FILE* save_fileptr;
  char save_filepath[1024];
  if (params.form == "diag") {
    sprintf(
      save_filepath, "%s/bk%d%d%d%s",
      params.measurement_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      params.output_tag.c_str()
    );
    save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
        kbin[ibin], kbin[ibin],
        norm * (bk_save[ibin].real() - shotnoise_save[ibin].real()),
        norm * shotnoise_save[ibin].real()
      );
    }
  } else if (params.form == "full") {
    sprintf(
      save_filepath, "%s/bk%d%d%d_kbin%02d%s",
      params.measurement_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      params.ith_kbin,
      params.output_tag.c_str()
    );
    save_fileptr = fopen(save_filepath, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
        kbin[params.ith_kbin], kbin[ibin],
        norm * (bk_save[ibin].real() - shotnoise_save[ibin].real()),
        norm * shotnoise_save[ibin].real()
      );
    }
  }
  fclose(save_fileptr);

  delete[] shotnoise_save;
  delete[] bk_save;

  return 0;
}

/**
 * Calculate bispectrum from catalogues for modes with individual
 * orders @f$M@f$ and save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param particles_rand (Random-source) particle container.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec_for_M_mode(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  ParameterSet& params,
  double alpha,
  double* kbin,
  double norm
) {
  if (currTask == 0) {
    printf(
      "[STAT] :: Measuring bispectrum for individual modes of order `m` "
      "from data and random catalogues.\n"
    );
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < EPS_COUPLING_3PT
  ) {
    if (currTask == 0) {
      printf(
        "[ERRO] :: Disallowed multipole degree combination "
        "for bispectrum measurements. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
      );
    }
    exit(1);
  }

  /// Initialise output shot noise terms.
  if (currTask == 0) {
    printf(
      "[STAT] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  std::vector< std::vector< std::complex<double> > > shotnoise_save;

  shotnoise_save.resize(2*params.ELL + 1);
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    shotnoise_save[ibin].resize(params.num_kbin);
  }

  for (int m_idx = 0; m_idx < 2*params.ELL + 1; m_idx++) {
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      shotnoise_save[m_idx][ibin] = 0.;
    }
  }

  /// Compute shot noise terms.
  PseudoDensityField<ParticleCatalogue> dn_00_for_shotnoise(params);
  dn_00_for_shotnoise.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  dn_00_for_shotnoise.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        /// Calculate N_LM in eq. (46) in arXiv:1803.02132.
        PseudoDensityField<ParticleCatalogue> shotnoise_quadratic_LM(params);
        shotnoise_quadratic_LM.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        shotnoise_quadratic_LM.fourier_transform();

        /// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
        Pseudo2ptStats<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
          stats.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );

        /// Calculate S_{\ell_1 \ell_2 L; i = j = k} in eq. (45)
        /// in arXiv:1803.02132.
        if (params.ell1 == 0 && params.ell2 == 0) {
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            shotnoise_save[M_ + params.ELL][ibin] += shotnoise_cubic_LM;
          }
        }

        /// Calculate S_{\ell_1 \ell_2 L; i ≠ j = k} in eq. (45)
        /// in arXiv:1803.02132.
        if (params.ell2 == 0) {
          stats.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_shotnoise, shotnoise_quadratic_LM,
            shotnoise_cubic_LM,
            kbin,
            params.ell1, m1_
          );
          if (params.form == "diag") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              shotnoise_save[M_ + params.ELL][ibin] += stats.pk[ibin];
            }
          } else if (params.form == "full") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              shotnoise_save[M_ + params.ELL][ibin] += stats.pk[params.ith_kbin];
            }
          }
        }

        /// Calculate S_{\ell_1 \ell_2 L; i = k ≠ j} in eq. (45)
        /// in arXiv:1803.02132.
        if (params.ell1 == 0) {
          stats.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_shotnoise, shotnoise_quadratic_LM,
            shotnoise_cubic_LM,
            kbin,
            params.ell2, m2_
          );
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            shotnoise_save[M_ + params.ELL][ibin] += stats.pk[ibin];
          }
        }

        if (currTask == 0) {
          printf(
            "[STAT] :: Computed shot noise term for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
            m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
          );
        }
      }
    }
  }

  dn_00_for_shotnoise.finalise_density_field();

  /// Calculate N_00 in eq. (45) in arXiv:1803.02132.
  PseudoDensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
  shotnoise_quadratic_00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  shotnoise_quadratic_00.fourier_transform();

  SphericalBesselCalculator sj1(params.ell1);
  SphericalBesselCalculator sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          *  wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        };

        PseudoDensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
        dn_LM_for_shotnoise.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM_for_shotnoise.fourier_transform();

        /// Calculate S_{\ell_1 \ell_2 L; i = j ≠ k} in eq. (45)
        /// in arXiv:1803.02132.
        Pseudo2ptStats<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
          stats.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
            particles_data, particles_rand,
            los_data, los_rand,
            alpha,
            params.ELL, M_
          );

        fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
        gbytesMem += sizeof(fftw_complex) *
          double(params.nmesh) / BYTES_PER_GBYTES;
        for (int gid = 0; gid < params.nmesh; gid++) {
          three_pt_holder[gid][0] = 0.;
          three_pt_holder[gid][1] = 0.;
        }

        stats.compute_2pt_self_shotnoise_for_bispec_meshgrid(
          dn_LM_for_shotnoise, shotnoise_quadratic_00, shotnoise_cubic_LM, three_pt_holder
          // params.ELL, M_
        );

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];
          double kmag_a;
          if (params.form == "diag") {
            kmag_a = kmag_b;
          } else if (params.form == "full") {
            kmag_a = kbin[params.ith_kbin];
          }

          double dr[3];
          dr[0] = params.boxsize[0] / double(params.ngrid[0]);
          dr[1] = params.boxsize[1] / double(params.ngrid[1]);
          dr[2] = params.boxsize[2] / double(params.ngrid[2]);

          double rvec[3];
          std::complex<double> shotnoise_sum = 0.;
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
                double rmag = sqrt(
                  rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
                );

                std::complex<double> shotnoise_grid(
                  three_pt_holder[idx_grid][0],
                  three_pt_holder[idx_grid][1]
                );
                double j1 = sj1.eval(kmag_a * rmag);
                double j2 = sj2.eval(kmag_b * rmag);

                shotnoise_sum +=
                  j1 * j2 * ylm_a[idx_grid] * ylm_b[idx_grid]
                  * shotnoise_grid;
              }
            }
          }

          double factor = params.volume / double(params.nmesh);
          shotnoise_sum *= factor * pow(M_I, params.ell1 + params.ell2);

          shotnoise_save[M_ + params.ELL][i_kbin] += shotnoise_sum;

          if (currTask == 0) {
            printf(
              "[STAT] :: Computed shot noise term for wavenumber and orders "
              "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
              "(%.3f seconds elapsed).\n",
              kmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
            );
          }
        }

        fftw_free(three_pt_holder); three_pt_holder = NULL;
        gbytesMem -= sizeof(fftw_complex)
          * double(params.nmesh) / BYTES_PER_GBYTES;
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;
    }
  }

  shotnoise_quadratic_00.finalise_density_field();

  if (currTask == 0) {
    printf(
      "[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Initialise output bispectrum.
  if (currTask == 0) {
    printf(
      "[STAT] :: Computing bispectrum terms (%.3f seconds elapsed...).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  std::vector< std::vector< std::complex<double> > > bk_save;

  /// Compute bispectrum.
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  dn_00.fourier_transform();

  bk_save.resize(2*params.ELL + 1);
  for (int ibin = 0; ibin < params.num_kbin; ibin++) {
    bk_save[ibin].resize(params.num_kbin);
  }

  for (int m = 0; m < 2*params.ELL + 1; m++) {
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      bk_save[m][ibin] = 0.;
    }
  }

  /// Compute bispectrum.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      gbytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > EPS_COUPLING_3PT) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
        );
        SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
        );
      }

      /// Calculate G_{LM} in eq. (42) in arXiv:1803.02132.
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < EPS_COUPLING_3PT) {
          continue;
        }

        PseudoDensityField<ParticleCatalogue> dn_LM(params);
        dn_LM.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM.fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.inv_fourier_transform();

        /// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
        PseudoDensityField<ParticleCatalogue> F_lm_a(params);
        double kmag_a;
        double dk = kbin[1] - kbin[0];
        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_a, kmag_a, dk
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];

          if (params.form == "diag") {
            kmag_a = kmag_b;
            F_lm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
              dn_00, ylm_a, kmag_a, dk
            );
          }

          PseudoDensityField<ParticleCatalogue> F_lm_b(params);
          F_lm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_b, kmag_b, dk
          );

          double factor = params.volume / double(params.nmesh);
          std::complex<double> bk_sum = 0.;
          for (int ibin = 0; ibin < params.nmesh; ibin++) {
            std::complex<double> F_ellm_1(F_lm_a[ibin][0], F_lm_a[ibin][1]);
            std::complex<double> F_ellm_2(F_lm_b[ibin][0], F_lm_b[ibin][1]);
            std::complex<double> G_LM(dn_LM[ibin][0], dn_LM[ibin][1]);
            bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
          }

          /// No coupling multiplication.
          bk_save[M_ + params.ELL][i_kbin] += bk_sum;

          if (currTask == 0) {
            printf(
              "[STAT] :: Computed bispectrum term for wavenumber and orders "
              "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
              "(%.3f seconds elapsed).\n",
              kmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      gbytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / BYTES_PER_GBYTES;
    }
  }

  if (currTask == 0) {
    printf(
      "[STAT] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Normalise and then save the output.
  double norm = pow(norm, 2) / pow(particles_data.wtotal, 3);

  for (int M_ = 0; M_<2*params.ELL + 1; M_++) {
    FILE* save_fileptr;
    char save_filepath[1024];
    if (params.form == "diag") {
      sprintf(
        save_filepath, "%s/bk%d%d%d_M%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2,
        params.ELL, M_,
        params.output_tag.c_str()
      );
      save_fileptr = fopen(save_filepath, "w");
      for (int ibin = 0; ibin < params.num_kbin; ibin++) {
        fprintf(
          save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
          kbin[ibin], kbin[ibin],
          norm * (bk_save[M_][ibin].real() - shotnoise_save[M_][ibin].real()),
          norm * shotnoise_save[M_][ibin].real()
        );
      }
    } else if (params.form == "full") {
      sprintf(
        save_filepath, "%s/bk%d%d%d_M%d_kbin%02d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2,
        params.ELL, M_,
        params.ith_kbin,
        params.output_tag.c_str()
      );
      save_fileptr = fopen(save_filepath, "w");
      for (int ibin = 0; ibin < params.num_kbin; ibin++) {
        fprintf(
          save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
          kbin[params.ith_kbin], kbin[ibin],
          norm * (bk_save[M_][ibin].real() - shotnoise_save[M_][ibin].real()),
          norm * shotnoise_save[M_][ibin].real()
        );
      }
    }
    fclose(save_fileptr);
  }

  return 0;
}
// #endif  // TRIUMVIRATE_USE_DISABLED_CODE_

#endif  // TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_
