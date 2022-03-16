#ifndef TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>

#include <fftw3.h>

#include "common.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "harmonic.hpp"
#include "bessel.hpp"
#include "particles.hpp"
#include "field.hpp"
#include "twopt.hpp"

/**
 * Calculate bispectrum from catalogues and save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param particles_rand (Random-source) particle container.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  ParameterSet& params,
  double alpha,
  double* kbin,
  double survey_vol_norm
) {
  if (currTask == 0) {
    printf(
			"[STAT] :: Measuring bispectrum from data and random catalogues.\n"
		);
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
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
  clockElapsed = double(clock() - clockStart);
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
        if (fabs(coupling) < 1.e-10) {
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
            shotnoise_save[ibin] += coupling * shotnoise_cubic_LM;
          }
        }

    		/// Calculate S_{\ell_1 \ell_2 L; i != j = k} in eq. (45)
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
              shotnoise_save[ibin] += coupling * stats.pk[ibin];
            }
          } else if (params.form == "full") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              shotnoise_save[ibin] += coupling * stats.pk[params.ith_kbin];
            }
          }
        }

    		/// Calculate S_{\ell_1 \ell_2 L; i = k != j} in eq. (45)
				/// in arXiv:1803.02132.
        if (params.ell1 == 0) {
          stats.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_shotnoise, shotnoise_quadratic_LM,
						shotnoise_cubic_LM,
						kbin,
						params.ell2, m2_
          );
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            shotnoise_save[ibin] += coupling * stats.pk[ibin];
          }
        }

        clockElapsed = double(clock() - clockStart);
        if (currTask == 0) {
          printf(
            "[STAT] :: Computed shot noise term for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
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
      bytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > 1.e-10) {
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
        if (fabs(coupling) < 1.e-10) {
          continue;
        }

        PseudoDensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
        dn_LM_for_shotnoise.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM_for_shotnoise.fourier_transform();

    		/// Calculate S_{\ell_1 \ell_2 L; i = j != k} in eq. (45)
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
        bytesMem += sizeof(fftw_complex) *
          double(params.nmesh) / 1024. / 1024. / 1024.;
        for (int gid = 0; gid < params.nmesh; gid++) {
          three_pt_holder[gid][0] = 0.;
          three_pt_holder[gid][1] = 0.;
        }

        stats.compute_2pt_self_shotnoise_for_bispec_meshgrid(
          dn_LM_for_shotnoise, shotnoise_quadratic_00,
					shotnoise_cubic_LM,
					params.ELL, M_,
					three_pt_holder
        );

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_a;
          double kmag_b = kbin[i_kbin];
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

								/// NOTE: This conforms to the absurd FFT array ordering
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

          std::complex<double> I_(0., 1.);
          double factor = params.volume / double(params.nmesh);
          shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

          shotnoise_save[i_kbin] += coupling * shotnoise_sum;

          clockElapsed = double(clock() - clockStart);
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
        bytesMem -= sizeof(fftw_complex)
          * double(params.nmesh) / 1024. / 1024. / 1024.;
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      bytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;
    }
  }

  shotnoise_quadratic_00.finalise_density_field();

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Initialise output bispectrum.
	double clockElapsed = double(clock() - clockStart);
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
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
  );
  dn_00.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      bytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > 1.e-10) {
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
        if (fabs(coupling) < 1.e-10) {
          continue;
        }

				/// Calculate G_{LM} in eq. (42) in arXiv:1803.02132.
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

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
        PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
        double kmag_a;
        double dk = kbin[1] - kbin[0];
        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_a, kmag_a, dk
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];

          if (params.form == "diag") {
            kmag_a = kmag_b;
            F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
              dn_00, ylm_a, kmag_a, dk
            );
          }

          PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
          F_ellm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_b, kmag_b, dk
          );

          double factor = params.volume / double(params.nmesh);
          std::complex<double> bk_sum = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
            std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
            std::complex<double> G_LM(dn_LM[gid][0], dn_LM[gid][1]);
            bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
          }

          bk_save[i_kbin] += coupling * bk_sum;

          double clockElapsed = double(clock() - clockStart);
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
      bytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;
    }
  }

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_data.wtotal, 3);

  FILE* save_fileptr;
  char buf[1024];
  if (params.form == "diag") {
    sprintf(
      buf, "%s/bk%d%d%d%s",
      params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
    );
    save_fileptr = fopen(buf, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[ibin], kbin[ibin],
        1. * (bk_save[ibin].real() - shotnoise_save[ibin].real()),
        1. * (bk_save[ibin].imag() - shotnoise_save[ibin].imag()),
        1. * shotnoise_save[ibin].real(),
        1. * shotnoise_save[ibin].imag()
      );
    }
  } else if (params.form == "full") {
    sprintf(
      buf, "%s/bk%d%d%d_kbin%02d%s",
      params.measurement_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      params.ith_kbin,
			params.output_tag.c_str()
    );
    save_fileptr = fopen(buf, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[params.ith_kbin], kbin[ibin],
        norm * (bk_save[ibin].real() - shotnoise_save[ibin].real()),
        norm * (bk_save[ibin].imag() - shotnoise_save[ibin].imag()),
        norm * shotnoise_save[ibin].real(),
        norm * shotnoise_save[ibin].imag()
      );
    }
  }
  fclose(save_fileptr);

  delete[] shotnoise_save;
  delete[] bk_save;

  return 0;
}

/**
 * Calculate bispectrum in a periodic box and save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param params Parameter set.
 * @param kbin Wavenumber bins.
 * @returns Exit status.
 */
int calc_bispec_in_box(
  ParticleCatalogue& particles_data,
  ParameterSet& params,
  double* kbin
) {
  if (currTask == 0) {
    printf("[STAT] :: Measuring bispectrum in a periodic box.\n");
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
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
  clockElapsed = double(clock() - clockStart);
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

	/// Calculate normal density field in the global plane-parallel
	/// picture in contrast with `shotnoise_cubic_LM`.
	PseudoDensityField<ParticleCatalogue> density(params);
	density.compute_unweighted_density(particles_data);
	density.fourier_transform();

  /// Compute shot noise terms.
	/// WARNING: This was inherited from Sugiyama et al. without
	/// matching equation.
  PseudoDensityField<ParticleCatalogue> dn_00_for_shotnoise(params);
  dn_00_for_shotnoise.compute_unweighted_fluctuation_insitu(
		particles_data, params.volume
	);
  dn_00_for_shotnoise.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;

      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (fabs(coupling) < 1.e-10) {
        continue;
      }

      Pseudo2ptStats<ParticleCatalogue> stats(params);
      std::complex<double> shotnoise = double(particles_data.ntotal);
        // NOTE: ``double`` conversion essential here

      if (params.ell1 == 0 && params.ell2 == 0) {
        for (int ibin = 0; ibin < params.num_kbin; ibin++) {
          shotnoise_save[ibin] += coupling * shotnoise;
        }
      }

      if (params.ell2 == 0) {
        stats.compute_ylm_wgtd_2pt_stats_in_fourier(
          dn_00_for_shotnoise, density,
          shotnoise,
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

      if (params.ell1 == 0) {
        stats.compute_ylm_wgtd_2pt_stats_in_fourier(
          dn_00_for_shotnoise, density,
          shotnoise,
          kbin,
          params.ell2, m2_
        );
        for (int ibin = 0; ibin < params.num_kbin; ibin++) {
          shotnoise_save[ibin] += coupling * stats.pk[ibin];
        }
      }

      clockElapsed = double(clock() - clockStart);
      if (currTask == 0) {
        printf(
          "[STAT] :: Computed shot noise term for orders "
          "`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
          m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
        );
      }
    }
  }

  dn_00_for_shotnoise.finalise_density_field();

  SphericalBesselCalculator sj1(params.ell1);
  SphericalBesselCalculator sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;

      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (fabs(coupling) < 1.e-10) {
        continue;
      }

      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      bytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;

      SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
        params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
      );
      SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
        params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
      );

      PseudoDensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
      dn_LM_for_shotnoise.compute_unweighted_fluctuation_insitu(
				particles_data, params.volume
			);
      dn_LM_for_shotnoise.fourier_transform();

      Pseudo2ptStats<ParticleCatalogue> stats(params);
      std::complex<double> shotnoise = double(particles_data.ntotal);

      fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
      bytesMem += sizeof(fftw_complex)
        * double(params.nmesh) / 1024. / 1024. / 1024.;
      for (int gid = 0; gid < params.nmesh; gid++) {
        three_pt_holder[gid][0] = 0.;
        three_pt_holder[gid][1] = 0.;
      }

      stats.compute_2pt_self_shotnoise_for_bispec_meshgrid(
        dn_LM_for_shotnoise, density,
				shotnoise,
        params.ELL, M_,
        three_pt_holder
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

							/// NOTE: This conforms to the absurd FFT array ordering
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
								three_pt_holder[idx_grid][0], three_pt_holder[idx_grid][1]
							);
              double j1 = sj1.eval(kmag_a * rmag);
              double j2 = sj2.eval(kmag_b * rmag);

              shotnoise_sum +=
								j1 * j2 * ylm_a[idx_grid] * ylm_b[idx_grid]
								* shotnoise_grid;
            }
          }
        }

        std::complex<double> I_(0., 1.);
        double factor = params.volume / double(params.nmesh);
        shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

        shotnoise_save[i_kbin] += coupling * shotnoise_sum;

        clockElapsed = double(clock() - clockStart);
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
      bytesMem -= sizeof(fftw_complex)
        * double(params.nmesh) / 1024. / 1024. / 1024.;

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
  		bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	density.finalise_density_field();

	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	double clockElapsed = double(clock() - clockStart);
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
	PseudoDensityField<ParticleCatalogue> dn_00(params);
	dn_00.compute_unweighted_fluctuation_insitu(
		particles_data, params.volume
	);
	dn_00.fourier_transform();

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			int M_ = 0;

			double coupling = double(2*params.ELL + 1)
				* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
				* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
			if (fabs(coupling) < 1.e-10) {
				continue;
			}

			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
				params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
			);
			SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
				params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
			);

			/// Calculate G_{LM} in eq. (42) in arXiv:1803.02132 (equivalent to
			/// normal density fluctuation with L, M = 0, 0 in the global
			/// plane-parallel picture).
			PseudoDensityField<ParticleCatalogue> dn_LM(params);
			dn_LM.compute_unweighted_fluctuation_insitu(
				particles_data, params.volume
			);
			dn_LM.fourier_transform();
			dn_LM.apply_assignment_compensation();
			dn_LM.inv_fourier_transform();

			/// NOTE: Standard naming convention is overriden below.

			/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
			PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
			double kmag_a;
			double dk = kbin[1] - kbin[0];
			if (params.form == "full") {
				kmag_a = kbin[params.ith_kbin];
				F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
					dn_00, ylm_a, kmag_a, dk
				);
			}

			for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
				double kmag_b = kbin[i_kbin];

				if (params.form == "diag") {
					kmag_a = kmag_b;
					F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
						dn_00, ylm_a, kmag_a, dk
					);
				}

				PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
				F_ellm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
					dn_00, ylm_b, kmag_b, dk
				);

				double factor = params.volume / double(params.nmesh);
				std::complex<double> bk_sum = 0.;
				for (int gid = 0; gid < params.nmesh; gid++) {
					std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
					std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
					std::complex<double> G_LM(dn_LM[gid][0], dn_LM[gid][1]);
					bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
				}

				bk_save[i_kbin] += coupling * bk_sum;

				double clockElapsed = double(clock() - clockStart);
				if (currTask == 0) {
					printf(
						"[STAT] :: Computed bispectrum term for wavenumber and orders "
						"`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
						kmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
					);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
	/// NOTE: Save the real parts only.
	double norm = params.volume
		/ double(particles_data.ntotal) / double(particles_data.ntotal);
	norm *= params.volume / double(particles_data.ntotal);

	FILE* save_fileptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/bk%d%d%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
		for (int ibin = 0; ibin < params.num_kbin; ibin++) {
			fprintf(
				save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[ibin], kbin[ibin],
				1. * (bk_save[ibin].real() - shotnoise_save[ibin].real()),
				1. * shotnoise_save[ibin].real()
			);
		}
	} else if (params.form == "full") {
		sprintf(
			buf, "%s/bk%d%d%d_kbin%02d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_kbin,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
		for (int ibin = 0; ibin < params.num_kbin; ibin++) {
			fprintf(
				save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[params.ith_kbin], kbin[ibin],
				1. * (bk_save[ibin].real() - shotnoise_save[ibin].real()),
				1. * shotnoise_save[ibin].real()
			);
		}
	}
	fclose(save_fileptr);

	delete[] shotnoise_save;
	delete[] bk_save;

	return 0;
}

/**
 * Calculate three-point correlation function from catalogues
 * and save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param particles_rand (Random-source) particle container.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pcf(
	ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (currTask == 0) {
		printf(
			"[STAT] :: Measuring three-point correlation function "
			"from data and random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[ERRO] :: Disallowed multipole degree combination "
				"for three-point correlation function measurements. "
				"Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	clockElapsed = double(clock() - clockStart);
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

	/// Compute shot noise terms, including only
	/// S_{\ell_1 \ell_2 L; i = j != k} in eq. (45) in arXiv:1803.02132
	/// (see eq. 51).
	PseudoDensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
	shotnoise_quadratic_00.compute_ylm_wgtd_2pt_self_component_for_shotnoise(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	shotnoise_quadratic_00.fourier_transform();

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				PseudoDensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
				dn_LM_for_shotnoise.compute_ylm_wgtd_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM_for_shotnoise.fourier_transform();

				Pseudo2ptStats<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise_cubic_LM =
					stats.calc_ylm_wgtd_3pt_self_component_for_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

				stats.compute_uncoupled_shotnoise_for_3pcf(
					dn_LM_for_shotnoise, shotnoise_quadratic_00,
					ylm_a, ylm_b,
					shotnoise_cubic_LM,
					rbin
					// params.ell1, m1_
				);

				for (int ibin = 0; ibin < params.num_rbin; ibin++) {
					if (params.form == "diag") {
						shotnoise_save[ibin] += coupling * stats.xi[ibin];
					} else if (params.form == "full") {
						/// Calculate shot noise contribution equivalent to the
						/// Kronecker delta in eq. (51) in arXiv:1803.02132.
						if (ibin == params.ith_rbin) {
							shotnoise_save[ibin] += coupling * stats.xi[ibin];
						} else {
							shotnoise_save[ibin] += 0.;
						}
					}
				}

				clockElapsed = double(clock() - clockStart);
				if (currTask == 0) {
					printf(
						"[STAT] :: Computed shot noise term for orders "
						"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
						m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_quadratic_00.finalise_density_field();

	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function.
	double clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computing three-point correlation terms "
			"(%.3f seconds elapsed...).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int ibin = 0; ibin < params.num_rbin; ibin++) {
		zeta_save[ibin] = 0.;
	}

  /// Compute three-point correlation function.
	PseudoDensityField<ParticleCatalogue> dn_00(params);
	dn_00.compute_ylm_wgtd_fluctuation(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	dn_00.fourier_transform();

	SphericalBesselCalculator sj1(params.ell1);
	SphericalBesselCalculator sj2(params.ell2);
	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				/// Calculate G_{LM} in eq. (42) in arXiv:1803.02132.
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

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (49) in arXiv:1803.02132.
				PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_a, sj1, rmag_a
					);
				}

				for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
					double rmag_b = rbin[i_rbin];

					if (params.form == "diag") {
						rmag_a = rmag_b;
						F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
							dn_00, ylm_a, sj1, rmag_a
						);
					}

					PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_b, sj2, rmag_b
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh);
					for (int gid = 0; gid < params.nmesh; gid++) {
						std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
						std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
						std::complex<double> G_LM(dn_LM[gid][0], dn_LM[gid][1]);
						zeta_sum += pow(I_, params.ell1 + params.ell2) * factor
							* F_ellm_1 * F_ellm_2 * G_LM;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double clockElapsed = double(clock() - clockStart);
					if (currTask == 0) {
						printf(
							"[STAT] :: Computed three-point correlation function term "
							"for separation and orders "
							"`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							rmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
						);
					}
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed three-point correlation function terms "
      "(... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_data.wtotal, 3);

	FILE* save_fileptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
			buf, "%s/zeta%d%d%d_rbin%02d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_rbin,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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

/**
 * Calculate three-point correlation function in a periodic box
 * and save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param params Parameter set.
 * @param rbin Separation bins.
 * @returns Exit status.
 */
int calc_3pcf_in_box(
	ParticleCatalogue& particles_data,
	ParameterSet& params,
	double* rbin
) {
	if (currTask == 0) {
		printf(
			"[STAT] :: Measuring three-point correlation function "
      "in a periodic box.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[ERRO] :: Disallowed multipole degree combination "
				"for three-point correlation function measurements. "
				"Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	clockElapsed = double(clock() - clockStart);
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

	/// Calculate normal density field in the global plane-parallel
	/// picture in contrast with `shotnoise_cubic_LM`.
	/// WARNING: This was inherited from Sugiyama et al. without
	/// matching equation.
	PseudoDensityField<ParticleCatalogue> density(params);
	density.compute_unweighted_density(particles_data);
	density.fourier_transform();

	/// Compute shot noise terms.
	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			int M_ = 0;

			double coupling = double(2*params.ELL + 1)
				* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
				* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
			if (fabs(coupling) < 1.e-10) {
				continue;
			}

			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
				params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
			);
			SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
				params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
			);

			PseudoDensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
			dn_LM_for_shotnoise.compute_unweighted_fluctuation_insitu(
				particles_data, params.volume
			);
			dn_LM_for_shotnoise.fourier_transform();

			Pseudo2ptStats<ParticleCatalogue> stats(params);
			std::complex<double> shotnoise = double(particles_data.ntotal);

			stats.compute_uncoupled_shotnoise_for_3pcf(
				dn_LM_for_shotnoise, density,
				ylm_a, ylm_b,
				shotnoise,
				rbin
				// params.ell1, m1_
			);

			for (int ibin = 0; ibin < params.num_rbin; ibin++) {
				if (params.form == "diag") {
					shotnoise_save[ibin] += coupling * stats.xi[ibin];
				} else if (params.form == "full") {
					/// Calculate shot noise contribution equivalent to the
					/// Kronecker delta in eq. (51) in arXiv:1803.02132.
					if (ibin == params.ith_rbin) {
						shotnoise_save[ibin] += coupling * stats.xi[ibin];
					} else {
						shotnoise_save[ibin] += 0.;
					}
				}
			}

			clockElapsed = double(clock() - clockStart);
			if (currTask == 0) {
				printf(
					"[STAT] :: Computed shot noise term for wavenumber and orders "
					"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
					m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
				);
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
					* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	density.finalise_density_field();

	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function.
	double clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computing three-point correlation function terms "
			"(%.3f seconds elapsed...).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int ibin = 0; ibin < params.num_rbin; ibin++) {
		zeta_save[ibin] = 0.;
	}

	/// Compute three-point correlation function.
	PseudoDensityField<ParticleCatalogue> dn_00(params);
	dn_00.compute_unweighted_fluctuation_insitu(
		particles_data, params.volume
	);
	dn_00.fourier_transform();

	SphericalBesselCalculator sj1(params.ell1);
	SphericalBesselCalculator sj2(params.ell2);
	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			int M_ = 0;

			double coupling = double(2*params.ELL + 1)
				* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
				* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
			if (fabs(coupling) < 1.e-10) {
				continue;
			}

			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
				params.ell1, m1_, params.boxsize, params.ngrid, ylm_a
			);
			SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
				params.ell2, m2_, params.boxsize, params.ngrid, ylm_b
			);

			/// Calculate G_{LM} in eq. (42) in arXiv:1803.02132 (equivalent to
			/// normal density fluctuation with L, M = 0, 0 in the global
			/// plane-parallel picture).
			PseudoDensityField<ParticleCatalogue> dn_LM(params);
			dn_LM.compute_unweighted_fluctuation_insitu(
				particles_data, params.volume
			);
			dn_LM.fourier_transform();
			dn_LM.apply_assignment_compensation();
			dn_LM.inv_fourier_transform();

			/// NOTE: Standard naming convention is overriden below.

			/// Calculate F_{\ell m} in eq. (49) in arXiv:1803.02132.
			PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
			double rmag_a;
			if (params.form == "full") {
				rmag_a = rbin[params.ith_rbin];
				F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
					dn_00, ylm_a, sj1, rmag_a
				);
			}

			for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
				double rmag_b = rbin[i_rbin];

				if (params.form == "diag") {
					rmag_a = rmag_b;
					F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_a, sj1, rmag_a
					);
				}

				PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
				F_ellm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
					dn_00, ylm_b, sj2, rmag_b
				);

				double factor = params.volume / double(params.nmesh);
				std::complex<double> I_(0., 1.);
				std::complex<double> zeta_sum = 0.;
				for (int gid = 0; gid < params.nmesh; gid++) {
					std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
					std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
					std::complex<double> G_LM(dn_LM[gid][0], dn_LM[gid][1]);
					zeta_sum += factor * pow(I_, params.ell1 + params.ell2)
						* F_ellm_1 * F_ellm_2 * G_LM;
				}

				zeta_save[i_rbin] += coupling * zeta_sum;

				double clockElapsed = double(clock() - clockStart);
				if (currTask == 0) {
					printf(
						"[STAT] :: Computed three-point correlation function term "
						"for separation and orders "
						"`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
						rmag_b, m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
					);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed three-point correlation function terms "
      "(... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
	double norm = params.volume
		/ double(particles_data.ntotal) / double(particles_data.ntotal);
	norm *= params.volume / double(particles_data.ntotal);

	FILE* save_fileptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
			buf, "%s/zeta%d%d%d_rbin%02d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_rbin,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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

/**
 * Calculate three-point correlation function window for three-point
 * correlation function from random catalogues and save the results.
 *
 * @param particles_rand (Random-source) particle container.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pcf_window(
	ParticleCatalogue& particles_rand,
	LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (currTask == 0) {
		printf(
			"[STAT] :: Measuring three-point correlation function window "
			"for three-point correlation function from random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
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
	clockElapsed = double(clock() - clockStart);
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

				clockElapsed = double(clock() - clockStart);
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
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_quadratic_00.finalise_density_field();

	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function window.
	double clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computing three-point correlation function window terms "
			"(%.3f seconds elapsed...).\n",
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				/// Calculate G_{LM}.  See `calc_3pcf`.s
				PseudoDensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.compute_ylm_wgtd_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM.fourier_transform();
				dn_LM.apply_assignment_compensation();
				dn_LM.inv_fourier_transform();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m}.  See `calc_3pcf`.
				PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_a, sj1, rmag_a
					);
				}

				for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
					double rmag_b = rbin[i_rbin];
					if (params.form == "diag") {
						rmag_a = rmag_b;
						F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
							dn_00, ylm_a, sj1, rmag_a
						);
					}

					PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_b, sj2, rmag_b
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh);
					for (int gid = 0; gid < params.nmesh; gid++) {
						std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
						std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
						std::complex<double> G_LM(dn_LM[gid][0], dn_LM[gid][1]);
						zeta_sum += pow(I_, params.ell1 + params.ell2) * factor
							* F_ellm_1 * F_ellm_2 * G_LM;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double clockElapsed = double(clock() - clockStart);
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
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed three-point correlation function window terms "
      "(... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_rand.wtotal, 3);
	norm /= alpha * alpha * alpha;

	FILE* save_fileptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d_window%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
			buf, "%s/zeta%d%d%d_window_rbin%02d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_rbin,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pcf_window_mpi(
	ParticleCatalogue& particles_rand,
	LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {  // WARNING: inherited from Sugiyama et al. without matching equation
	if (currTask == 0) {
		printf(
			"[STAT] :: Measuring three-point correlation function window "
			"from random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
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

	int n_temp = 10;  // NOTE: discretionary choice
	int NR = 3;  // NOTE: discretionary choice

	params.ith_rbin = currTask;

	/// Set up binning.
	/// NOTE: The setup binning is discretionary.
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
	clockElapsed = double(clock() - clockStart);
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

				clockElapsed = double(clock() - clockStart);
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
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_quadratic_00.finalise_density_field();

	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function window.
	double clockElapsed = double(clock() - clockStart);
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m}.  See `calc_3pcf`.
				PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_a, sj1, rmag_a
					);
				}

				for (int i_rbin = 0; i_rbin < n_temp; i_rbin++) {
					double rmag_b = rbin[i_rbin + NR * n_temp];

					if (params.form == "diag") {
						rmag_a = rmag_b;
						F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
							dn_00, ylm_a, sj1, rmag_a
						);
					}

					PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_b, sj2, rmag_b
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh);
					for (int gid = 0; gid < params.nmesh; gid++) {
						std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
						std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
						std::complex<double> G_LM(dn_LM[gid][0], dn_LM[gid][1]);
						zeta_sum += pow(I_, params.ell1 + params.ell2) * factor
							* F_ellm_1 * F_ellm_2 * G_LM;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double clockElapsed = double(clock() - clockStart);
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
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed three-point correlation function window terms "
      "(... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_rand.wtotal, 3);
	norm /= alpha * alpha * alpha;

	FILE* save_fileptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta_window_%d%d%d_%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			NR,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
			buf, "%s/zeta_window_%d%d%d_rbin%02d_%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_rbin, NR,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pcf_window_for_wide_angle(
	ParticleCatalogue& particles_rand,
	LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (currTask == 0) {
		printf(
			"[STAT] :: Measuring three-point correlation function window wide-angle"
			"correction terms for three-point correlation function "
			"from random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
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
	clockElapsed = double(clock() - clockStart);
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

				clockElapsed = double(clock() - clockStart);
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
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_quadratic_00.finalise_density_field();

	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function window.
	double clockElapsed = double(clock() - clockStart);
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m}.  See `calc_3pcf`.
				PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_a, sj1, rmag_a
					);
				}

				for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
					double rmag_b = rbin[i_rbin];
					if (params.form == "diag") {
						rmag_a = rmag_b;
						F_ellm_a.inv_fourier_transform_for_sjl_ylm_wgtd_field(
							dn_00, ylm_a, sj1, rmag_a
						);
					}

					PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.inv_fourier_transform_for_sjl_ylm_wgtd_field(
						dn_00, ylm_b, sj2, rmag_b
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh);
					for (int gid = 0; gid < params.nmesh; gid++) {
						std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
						std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
						std::complex<double> x_G_LM(dn_LM[gid][0], dn_LM[gid][1]);
						zeta_sum += pow(I_, params.ell1 + params.ell2) * factor
							* F_ellm_1 * F_ellm_2 * x_G_LM;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double clockElapsed = double(clock() - clockStart);
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
			bytesMem -= 2 * sizeof(std::complex<double>)
					* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed three-point correlation function window terms "
      "(... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_rand.wtotal, 3);
	norm /= alpha * alpha * alpha;

	FILE* save_fileptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d_window-wa%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
			buf, "%s/zeta%d%d%d_window-wa%d%d_rbin%02d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.i_wa, params.j_wa,
			params.ith_rbin,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec_for_los_choice(
	ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* kbin,
	int los,
	double survey_vol_norm
) {  // WARNING: inherited from Sugiyama et al. without matching equation
	if (currTask == 0) {
		printf(
			"[STAT] :: Measuring bispectrum for the choice of line of sight "
			"from data and random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
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
	clockElapsed = double(clock() - clockStart);
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
				if (fabs(coupling) < 1.e-10) {
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

    		/// Calculate S_{\ell_1 \ell_2 L; i != j = k} in eq. (45)
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
				/// matching equation.  S_{\ell_1 \ell_2 L; i = k != j}
				/// (i.e. ell1 == 0 case) may have been shuffled.

				clockElapsed = double(clock() - clockStart);
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
				if (fabs(coupling) < 1.e-10) {
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
				/// matching equation.  S_{\ell_1 \ell_2 L; i = k != j}
				/// (i.e. ell1 == 0 case) may have been shuffled.

    		/// Calculate S_{\ell_1 \ell_2 L; i = k != j} in eq. (45)
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

				clockElapsed = double(clock() - clockStart);
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

    		/// Calculate S_{\ell_1 \ell_2 L; i = j != k} in eq. (45)
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
				bytesMem += sizeof(fftw_complex)
					* double(params.nmesh) / 1024. / 1024. / 1024.;
				for (int ibin = 0; ibin < params.nmesh; ibin++) {
					three_pt_holder[ibin][0] = 0.;
					three_pt_holder[ibin][1] = 0.;
				}

				stats.compute_2pt_self_shotnoise_for_bispec_meshgrid(
					dn_shotnoise, shotnoise_quadratic,
					shotnoise_cubic_LM,
					params.ELL, M_,
					three_pt_holder
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

								/// NOTE: This conforms to the absurd FFT array ordering
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

					std::complex<double> I_(0., 1.);
					double factor = params.volume / double(params.nmesh);
					shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

					shotnoise_save[i_kbin] += coupling * shotnoise_sum;

					clockElapsed = double(clock() - clockStart);
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
				bytesMem -= sizeof(fftw_complex)
					* double(params.nmesh) / 1024. / 1024. / 1024.;
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	double clockElapsed = double(clock() - clockStart);
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				  *	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
				PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
				double kmag_a;
				double dk = kbin[1] - kbin[0];
				if (params.form == "full") {
					kmag_a = kbin[params.ith_kbin];
					F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
						dn_los1, ylm_a, kmag_a, dk
					);
				}

				for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
					double kmag_b = kbin[i_kbin];

					if (params.form == "diag") {
						kmag_a = kmag_b;
						F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
							dn_los1, ylm_a, kmag_a, dk
						);
					}

					PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
						dn_los2, ylm_b, kmag_b, dk
					);

					double factor = params.volume / double(params.nmesh);
					std::complex<double> bk_sum = 0.;
					for (int gid = 0; gid < params.nmesh; gid++) {
						std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
						std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
						std::complex<double> G_LM(dn_los3[gid][0], dn_los3[gid][1]);
						bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
					}

					bk_save[i_kbin] += coupling * bk_sum;

					double clockElapsed = double(clock() - clockStart);
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
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
	/// NOTE: Save the real parts only.
  double norm = pow(survey_vol_norm, 2) / pow(particles_data.wtotal, 3);

	FILE* save_fileptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/bk%d%d%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
			buf, "%s/bk%d%d%d_kbin%02d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_kbin,
			params.output_tag.c_str()
		);
		save_fileptr = fopen(buf, "w");
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
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec_for_M_mode(
	ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* kbin,
	double survey_vol_norm
) {
	if (currTask == 0) {
		printf(
			"[STAT] :: Measuring bispectrum for individual modes of order `m` "
			"from data and random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
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
	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// NOTE: Additional spaces are needed below to avoid clash with
	/// `>>` opertor.
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
				if (fabs(coupling) < 1.e-10) {
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

    		/// Calculate S_{\ell_1 \ell_2 L; i != j = k} in eq. (45)
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

    		/// Calculate S_{\ell_1 \ell_2 L; i = k != j} in eq. (45)
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

				clockElapsed = double(clock() - clockStart);
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				  *	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

    		/// Calculate S_{\ell_1 \ell_2 L; i = j != k} in eq. (45)
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
				bytesMem += sizeof(fftw_complex) *
					double(params.nmesh) / 1024. / 1024. / 1024.;
				for (int gid = 0; gid < params.nmesh; gid++) {
					three_pt_holder[gid][0] = 0.;
					three_pt_holder[gid][1] = 0.;
				}

				stats.compute_2pt_self_shotnoise_for_bispec_meshgrid(
					dn_LM_for_shotnoise, shotnoise_quadratic_00,
					shotnoise_cubic_LM,
					params.ELL, M_,
					three_pt_holder
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

								/// NOTE: This conforms to the absurd FFT array ordering
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

					std::complex<double> I_(0., 1.);
					double factor = params.volume / double(params.nmesh);
					shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

					shotnoise_save[M_ + params.ELL][i_kbin] += shotnoise_sum;

					clockElapsed = double(clock() - clockStart);
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
				bytesMem -= sizeof(fftw_complex)
					* double(params.nmesh) / 1024. / 1024. / 1024.;
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_quadratic_00.finalise_density_field();

	clockElapsed = double(clock() - clockStart);
	if (currTask == 0) {
		printf(
			"[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	double clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
			"[STAT] :: Computing bispectrum terms (%.3f seconds elapsed...).\n",
			clockElapsed / CLOCKS_PER_SEC
		);
  }

	/// NOTE: Additional spaces are needed below to avoid clash with
	/// `>>` opertor.
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
			bytesMem += 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
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
				if (fabs(coupling) < 1.e-10) {
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

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
				PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
				double kmag_a;
				double dk = kbin[1] - kbin[0];
				if (params.form == "full") {
					kmag_a = kbin[params.ith_kbin];
					F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
						dn_00, ylm_a, kmag_a, dk
					);
				}

				for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
					double kmag_b = kbin[i_kbin];

					if (params.form == "diag") {
						kmag_a = kmag_b;
						F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
							dn_00, ylm_a, kmag_a, dk
						);
					}

					PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
						dn_00, ylm_b, kmag_b, dk
					);

					double factor = params.volume / double(params.nmesh);
					std::complex<double> bk_sum = 0.;
					for (int ibin = 0; ibin < params.nmesh; ibin++) {
						std::complex<double> F_ellm_1(F_ellm_a[ibin][0], F_ellm_a[ibin][1]);
						std::complex<double> F_ellm_2(F_ellm_b[ibin][0], F_ellm_b[ibin][1]);
						std::complex<double> G_LM(dn_LM[ibin][0], dn_LM[ibin][1]);
						bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
					}

					/// No coupling multiplication.
					bk_save[M_ + params.ELL][i_kbin] += bk_sum;

					double clockElapsed = double(clock() - clockStart);
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
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
	/// NOTE: Save the real parts only.
  double norm = pow(survey_vol_norm, 2) / pow(particles_data.wtotal, 3);

	for (int M_ = 0; M_<2*params.ELL + 1; M_++) {
		FILE* save_fileptr;
		char buf[1024];
		if (params.form == "diag") {
			sprintf(
				buf, "%s/bk%d%d%d_M%d%s",
				params.measurement_dir.c_str(),
				params.ell1, params.ell2,
				params.ELL, M_,
				params.output_tag.c_str()
			);
			save_fileptr = fopen(buf, "w");
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
				buf, "%s/bk%d%d%d_M%d_kbin%02d%s",
				params.measurement_dir.c_str(),
				params.ell1, params.ell2,
				params.ELL, M_,
				params.ith_kbin,
				params.output_tag.c_str()
			);
			save_fileptr = fopen(buf, "w");
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

/**
 * Calculate bispectrum from catalogues and save the results.
 *
 * @param particles_data (Data-source) particle container.
 * @param particles_rand (Random-source) particle container.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec_(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  std::vector< std::vector<double> > los_data_arr,
	std::vector< std::vector<double> > los_rand_arr,
  ParameterSet& params,
  double alpha,
  double* kbin
) {
  if (currTask == 0) {
    printf(
			"[STAT] :: Measuring bispectrum from data and random catalogues.\n"
		);
  }

	LineOfSight* los_data = new LineOfSight[los_data_arr.size()];
	LineOfSight* los_rand = new LineOfSight[los_rand_arr.size()];
  for (int id = 0; id < particles_data.ntotal; id++) {
      los_data[id].pos[0] = los_data_arr[id][0];
      los_data[id].pos[1] = los_data_arr[id][1];
      los_data[id].pos[2] = los_data_arr[id][2];
  }
  for (int id = 0; id < particles_rand.ntotal; id++) {
      los_rand[id].pos[0] = los_rand_arr[id][0];
      los_rand[id].pos[1] = los_rand_arr[id][1];
      los_rand[id].pos[2] = los_rand_arr[id][2];
  }

	ParticleCatalogue::boxify_catalogues_for_fft(
		particles_data, particles_rand, params.boxsize, params.ngrid
	);

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
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
  clockElapsed = double(clock() - clockStart);
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
        if (fabs(coupling) < 1.e-10) {
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
            shotnoise_save[ibin] += coupling * shotnoise_cubic_LM;
          }
        }

    		/// Calculate S_{\ell_1 \ell_2 L; i != j = k} in eq. (45)
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
              shotnoise_save[ibin] += coupling * stats.pk[ibin];
            }
          } else if (params.form == "full") {
            for (int ibin = 0; ibin < params.num_kbin; ibin++) {
              shotnoise_save[ibin] += coupling * stats.pk[params.ith_kbin];
            }
          }
        }

    		/// Calculate S_{\ell_1 \ell_2 L; i = k != j} in eq. (45)
				/// in arXiv:1803.02132.
        if (params.ell1 == 0) {
          stats.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_shotnoise, shotnoise_quadratic_LM,
						shotnoise_cubic_LM,
						kbin,
						params.ell2, m2_
          );
          for (int ibin = 0; ibin < params.num_kbin; ibin++) {
            shotnoise_save[ibin] += coupling * stats.pk[ibin];
          }
        }

        clockElapsed = double(clock() - clockStart);
        if (currTask == 0) {
          printf(
            "[STAT] :: Computed shot noise term for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
            m1_, m2_, M_, clockElapsed / CLOCKS_PER_SEC
          );
        }
      }
    }
  }

	// MARKER
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
      bytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > 1.e-10) {
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

        if (fabs(coupling) < 1.e-10) {
          continue;
        }

        PseudoDensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
        dn_LM_for_shotnoise.compute_ylm_wgtd_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM_for_shotnoise.fourier_transform();

    		/// Calculate S_{\ell_1 \ell_2 L; i = j != k} in eq. (45)
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
        bytesMem += sizeof(fftw_complex) *
          double(params.nmesh) / 1024. / 1024. / 1024.;
        for (int gid = 0; gid < params.nmesh; gid++) {
          three_pt_holder[gid][0] = 0.;
          three_pt_holder[gid][1] = 0.;
        }

        stats.compute_2pt_self_shotnoise_for_bispec_meshgrid(
          dn_LM_for_shotnoise, shotnoise_quadratic_00,
					shotnoise_cubic_LM,
					params.ELL, M_,
					three_pt_holder
        );

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_a;
          double kmag_b = kbin[i_kbin];
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

								/// NOTE: This conforms to the absurd FFT array ordering
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

          std::complex<double> I_(0., 1.);
          double factor = params.volume / double(params.nmesh);
          shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

          shotnoise_save[i_kbin] += coupling * shotnoise_sum;

          clockElapsed = double(clock() - clockStart);
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
        bytesMem -= sizeof(fftw_complex)
          * double(params.nmesh) / 1024. / 1024. / 1024.;
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      bytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;
    }
  }

  shotnoise_quadratic_00.finalise_density_field();

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Initialise output bispectrum.
	double clockElapsed = double(clock() - clockStart);
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
  PseudoDensityField<ParticleCatalogue> dn_00(params);
  dn_00.compute_ylm_wgtd_fluctuation(
    particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
  );
  dn_00.fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh];
      bytesMem += 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > 1.e-10) {
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
        if (fabs(coupling) < 1.e-10) {
          continue;
        }

				/// Calculate G_{LM} in eq. (42) in arXiv:1803.02132.
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

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
        PseudoDensityField<ParticleCatalogue> F_ellm_a(params);
        double kmag_a;
        double dk = kbin[1] - kbin[0];
        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_a, kmag_a, dk
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];

          if (params.form == "diag") {
            kmag_a = kmag_b;
            F_ellm_a.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
              dn_00, ylm_a, kmag_a, dk
            );
          }

          PseudoDensityField<ParticleCatalogue> F_ellm_b(params);
          F_ellm_b.inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
            dn_00, ylm_b, kmag_b, dk
          );

          double factor = params.volume / double(params.nmesh);
          std::complex<double> bk_sum = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_ellm_1(F_ellm_a[gid][0], F_ellm_a[gid][1]);
            std::complex<double> F_ellm_2(F_ellm_b[gid][0], F_ellm_b[gid][1]);
            std::complex<double> G_LM(dn_LM[gid][0], dn_LM[gid][1]);
            bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
          }

          bk_save[i_kbin] += coupling * bk_sum;

          double clockElapsed = double(clock() - clockStart);
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
      bytesMem -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh) / 1024. / 1024. / 1024.;
    }
  }

  clockElapsed = double(clock() - clockStart);
  if (currTask == 0) {
    printf(
      "[STAT] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      clockElapsed / CLOCKS_PER_SEC
    );
  }

  /// Normalise and then save the output.
  FILE* save_fileptr;

  char buf[1024];

  if (params.form == "diag") {
    sprintf(
      buf, "%s/bk%d%d%d%s",
      params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
    );

    save_fileptr = fopen(buf, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[ibin], kbin[ibin],
        (bk_save[ibin].real() - shotnoise_save[ibin].real()),
        (bk_save[ibin].imag() - shotnoise_save[ibin].imag()),
        shotnoise_save[ibin].real(),
        shotnoise_save[ibin].imag()
      );
    }
  } else if (params.form == "full") {
    sprintf(
      buf, "%s/bk%d%d%d_kbin%02d%s",
      params.measurement_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      params.ith_kbin,
			params.output_tag.c_str()
    );

    save_fileptr = fopen(buf, "w");
    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      fprintf(
        save_fileptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[params.ith_kbin], kbin[ibin],
        (bk_save[ibin].real() - shotnoise_save[ibin].real()),
        (bk_save[ibin].imag() - shotnoise_save[ibin].imag()),
        shotnoise_save[ibin].real(),
        shotnoise_save[ibin].imag()
      );
    }
  }
  fclose(save_fileptr);

  delete[] shotnoise_save;
  delete[] bk_save;

  return 0;
}

#endif
