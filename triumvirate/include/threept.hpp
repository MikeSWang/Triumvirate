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
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
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
			"[Info] :: Measuring bispectrum from data and random catalogues.\n"
		);
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
  ) {
    if (currTask == 0) {
      printf(
        "[Error] :: Disallowed multipole degree combination "
        "for bispectrum measurements. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
      );
    }
    exit(1);
  }

  /// Initialise output shot noise terms.
  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* shotnoise_save =
    new std::complex<double>[params.num_kbin];
  for (int i = 0; i < params.num_kbin; i++) {
    shotnoise_save[i] = 0.;
  }

  /// Compute shot noise terms.
  DensityField<ParticleCatalogue> dn_00_for_shotnoise(params);
  dn_00_for_shotnoise.calc_ylm_weighted_fluctuation(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  dn_00_for_shotnoise.calc_fourier_transform();

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
        DensityField<ParticleCatalogue> shotnoise_quadratic_LM(params);
        shotnoise_quadratic_LM.calc_ylm_weighted_field_for_bispec_shotnoise(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        shotnoise_quadratic_LM.calc_fourier_transform();

    		/// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
        TwoPointStatistics<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

    		/// Calculate S_{\ell_1 \ell_2 L; i = j = k} in eq. (45)
				/// in arXiv:1803.02132.
        if (params.ell1 == 0 && params.ell2 == 0) {
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * shotnoise_cubic_LM;
          }
        }

    		/// Calculate S_{\ell_1 \ell_2 L; i != j = k} in eq. (45)
				/// in arXiv:1803.02132.
        if (params.ell2 == 0) {
          stats.calc_2pt_func_in_fourier(
            dn_00_for_shotnoise, shotnoise_quadratic_LM,
						kbin,
						shotnoise_cubic_LM,
						params.ell1, m1_
          );
          if (params.form == "diag") {
            for (int i = 0; i < params.num_kbin; i++) {
              shotnoise_save[i] += coupling * stats.pk[i];
            }
          } else if (params.form == "full") {
            for (int i = 0; i < params.num_kbin; i++) {
              shotnoise_save[i] += coupling * stats.pk[params.ith_kbin];
            }
          }
        }

    		/// Calculate S_{\ell_1 \ell_2 L; i = k != j} in eq. (45)
				/// in arXiv:1803.02132.
        if (params.ell1 == 0) {
          stats.calc_2pt_func_in_fourier(
            dn_00_for_shotnoise, shotnoise_quadratic_LM,
						kbin,
						shotnoise_cubic_LM,
						params.ell2, m2_
          );
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * stats.pk[i];
          }
        }

        durationInSec = double(clock() - timeStart);
        if (currTask == 0) {
          printf(
            "[Status] :: Computed shot noise term for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
            m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
          );
        }
      }
    }
  }

  dn_00_for_shotnoise.finalise_density_field();

	/// Calculate N_00 in eq. (45) in arXiv:1803.02132.
  DensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
  shotnoise_quadratic_00.calc_ylm_weighted_field_for_bispec_shotnoise(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  shotnoise_quadratic_00.calc_fourier_transform();

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

        DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
        dn_LM_for_shotnoise.calc_ylm_weighted_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM_for_shotnoise.calc_fourier_transform();

    		/// Calculate S_{\ell_1 \ell_2 L; i = j != k} in eq. (45)
				/// in arXiv:1803.02132.
        TwoPointStatistics<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

        fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
        bytesMem += sizeof(fftw_complex) *
          double(params.nmesh) / 1024. / 1024. / 1024.;
        for (int i = 0; i < params.nmesh; i++) {
          three_pt_holder[i][0] = 0.;
          three_pt_holder[i][1] = 0.;
        }

        stats.calc_shotnoise_for_bispec_on_grid(
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

          durationInSec = double(clock() - timeStart);
          if (currTask == 0) {
            printf(
              "[Status] :: Computed shot noise term for wavenumber and orders "
              "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
              kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  /// Initialise output bispectrum.
	double durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
			"[Status] :: Computing bispectrum terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
  }

  std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
  for (int i = 0; i < params.num_kbin; i++) {
    bk_save[i] = 0.;
  }

  /// Compute bispectrum.
  DensityField<ParticleCatalogue> dn_00(params);
  dn_00.calc_ylm_weighted_fluctuation(
    particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
  );
  dn_00.calc_fourier_transform();

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
        DensityField<ParticleCatalogue> dn_LM(params);
        dn_LM.calc_ylm_weighted_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM.calc_fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.calc_inverse_fourier_transform();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
        DensityField<ParticleCatalogue> F_ellm_a(params);
        double kmag_a;
        double dk = kbin[1] - kbin[0];
        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          F_ellm_a.calc_inverse_fourier_transform_for_bispec(
            dn_00, kmag_a, dk, ylm_a
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];

          if (params.form == "diag") {
            kmag_a = kmag_b;
            F_ellm_a.calc_inverse_fourier_transform_for_bispec(
              dn_00, kmag_a, dk, ylm_a
            );
          }

          DensityField<ParticleCatalogue> F_ellm_b(params);
          F_ellm_b.calc_inverse_fourier_transform_for_bispec(
            dn_00, kmag_b, dk, ylm_b
          );

          double factor = params.volume / double(params.nmesh);
          std::complex<double> bk_sum = 0.;
          for (int i = 0; i < params.nmesh; i++) {
            std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
            std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
            std::complex<double> G_LM(dn_LM[i][0], dn_LM[i][1]);
            bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
          }

          bk_save[i_kbin] += coupling * bk_sum;

          double durationInSec = double(clock() - timeStart);
          if (currTask == 0) {
            printf(
              "[Status] :: Computed bispectrum term for wavenumber and orders "
              "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
              kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  /// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_data.wtotal, 3);

  FILE* saved_file_ptr;
  char buf[1024];
  if (params.form == "diag") {
    sprintf(
      buf, "%s/bk%d%d%d%s",
      params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
    );
    saved_file_ptr = fopen(buf, "w");
    for (int i = 0; i < params.num_kbin; i++) {
      fprintf(
        saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[i], kbin[i],
        1. * (bk_save[i].real() - shotnoise_save[i].real()),
        1. * (bk_save[i].imag() - shotnoise_save[i].imag()),
        1. * shotnoise_save[i].real(),
        1. * shotnoise_save[i].imag()
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
    saved_file_ptr = fopen(buf, "w");
    for (int i = 0; i < params.num_kbin; i++) {
      fprintf(
        saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[params.ith_kbin], kbin[i],
        norm * (bk_save[i].real() - shotnoise_save[i].real()),
        norm * (bk_save[i].imag() - shotnoise_save[i].imag()),
        norm * shotnoise_save[i].real(),
        norm * shotnoise_save[i].imag()
      );
    }
  }
  fclose(saved_file_ptr);

  delete[] shotnoise_save;
  delete[] bk_save;

  return 0;
}

/**
 * Calculate bispectrum in a periodic box and save the results.
 *
 * @param particles_data (Reference to) the data-source particle container.
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
    printf("[Info] :: Measuring bispectrum in a periodic box.\n");
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
  ) {
    if (currTask == 0) {
      printf(
        "[Error] :: Disallowed multipole degree combination "
        "for bispectrum measurements. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
      );
    }
    exit(1);
  }

  /// Initialise output shot noise terms.
  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* shotnoise_save =
    new std::complex<double>[params.num_kbin];
  for (int i = 0; i < params.num_kbin; i++) {
    shotnoise_save[i] = 0.;
  }

	/// Calculate normal density field in the global plane-parallel
	/// picture in contrast with `shotnoise_cubic_LM`.
	DensityField<ParticleCatalogue> density(params);
	density.calc_density_field_in_box_for_bispec(particles_data);
	density.calc_fourier_transform();

  /// Compute shot noise terms.
	/// WARNING: This was inherited from Sugiyama et al. without
	/// matching equation.
  DensityField<ParticleCatalogue> dn_00_for_shotnoise(params);
  dn_00_for_shotnoise.calc_fluctuation_in_box(particles_data, params);
  dn_00_for_shotnoise.calc_fourier_transform();

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

      TwoPointStatistics<ParticleCatalogue> stats(params);
      std::complex<double> shotnoise = double(particles_data.ntotal);
        // NOTE: ``double`` conversion essential here

      if (params.ell1 == 0 && params.ell2 == 0) {
        for (int i = 0; i < params.num_kbin; i++) {
          shotnoise_save[i] += coupling * shotnoise;
        }
      }

      if (params.ell2 == 0) {
        stats.calc_2pt_func_in_fourier(
          dn_00_for_shotnoise, density,
          kbin,
          shotnoise,
          params.ell1, m1_
        );
        if (params.form == "diag") {
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * stats.pk[i];
          }
        } else if (params.form == "full") {
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * stats.pk[params.ith_kbin];
          }
        }
      }

      if (params.ell1 == 0) {
        stats.calc_2pt_func_in_fourier(
          dn_00_for_shotnoise, density,
          kbin,
          shotnoise,
          params.ell2, m2_
        );
        for (int i = 0; i < params.num_kbin; i++) {
          shotnoise_save[i] += coupling * stats.pk[i];
        }
      }

      durationInSec = double(clock() - timeStart);
      if (currTask == 0) {
        printf(
          "[Status] :: Computed shot noise term for orders "
          "`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
          m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

      DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
      dn_LM_for_shotnoise.calc_fluctuation_in_box(particles_data, params);
      dn_LM_for_shotnoise.calc_fourier_transform();

      TwoPointStatistics<ParticleCatalogue> stats(params);
      std::complex<double> shotnoise = double(particles_data.ntotal);

      fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
      bytesMem += sizeof(fftw_complex)
        * double(params.nmesh) / 1024. / 1024. / 1024.;
      for (int i = 0; i < params.nmesh; i++) {
        three_pt_holder[i][0] = 0.;
        three_pt_holder[i][1] = 0.;
      }

      stats.calc_shotnoise_for_bispec_on_grid(
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

        durationInSec = double(clock() - timeStart);
        if (currTask == 0) {
          printf(
            "[Status] :: Computed shot noise term for wavenumber and orders "
            "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
            kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	double durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
    printf(
			"[Status] :: Computing bispectrum terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
	for (int i = 0; i < params.num_kbin; i++) {
		bk_save[i] = 0.;
	}

	/// Compute bispectrum.
	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_fluctuation_in_box(particles_data, params);
	dn_00.calc_fourier_transform();

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
			DensityField<ParticleCatalogue> dn_LM(params);
			dn_LM.calc_fluctuation_in_box(particles_data, params);
			dn_LM.calc_fourier_transform();
			dn_LM.apply_assignment_compensation();
			dn_LM.calc_inverse_fourier_transform();

			/// NOTE: Standard naming convention is overriden below.

			/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
			DensityField<ParticleCatalogue> F_ellm_a(params);
			double kmag_a;
			double dk = kbin[1] - kbin[0];
			if (params.form == "full") {
				kmag_a = kbin[params.ith_kbin];
				F_ellm_a.calc_inverse_fourier_transform_for_bispec(
					dn_00, kmag_a, dk, ylm_a
				);
			}

			for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
				double kmag_b = kbin[i_kbin];

				if (params.form == "diag") {
					kmag_a = kmag_b;
					F_ellm_a.calc_inverse_fourier_transform_for_bispec(
						dn_00, kmag_a, dk, ylm_a
					);
				}

				DensityField<ParticleCatalogue> F_ellm_b(params);
				F_ellm_b.calc_inverse_fourier_transform_for_bispec(
					dn_00, kmag_b, dk, ylm_b
				);

				double factor = params.volume / double(params.nmesh);
				std::complex<double> bk_sum = 0.;
				for (int i = 0; i < params.nmesh; i++) {
					std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
					std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
					std::complex<double> G_LM(dn_LM[i][0], dn_LM[i][1]);
					bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
				}

				bk_save[i_kbin] += coupling * bk_sum;

				double durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed bispectrum term for wavenumber and orders "
						"`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
						kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
	/// NOTE: Save the real parts only.
	double norm = params.volume
		/ double(particles_data.ntotal) / double(particles_data.ntotal);
	norm *= params.volume / double(particles_data.ntotal);

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/bk%d%d%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_kbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[i], kbin[i],
				1. * (bk_save[i].real() - shotnoise_save[i].real()),
				1. * shotnoise_save[i].real()
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
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_kbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[params.ith_kbin], kbin[i],
				1. * (bk_save[i].real() - shotnoise_save[i].real()),
				1. * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] bk_save;

	return 0;
}

/**
 * Calculate three-point correlation function from catalogues
 * and save the results.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pt_corrfunc(
	ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (currTask == 0) {
		printf(
			"[Info] :: Measuring three-point correlation function "
			"from data and random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for three-point correlation function measurements. "
				"Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms, including only
	/// S_{\ell_1 \ell_2 L; i = j != k} in eq. (45) in arXiv:1803.02132
	/// (see eq. 51).
	DensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
	shotnoise_quadratic_00.calc_ylm_weighted_field_for_bispec_shotnoise(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	shotnoise_quadratic_00.calc_fourier_transform();

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

				DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
				dn_LM_for_shotnoise.calc_ylm_weighted_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM_for_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

				stats.calc_2pt_func_for_3pt_corrfunc(
					dn_LM_for_shotnoise, shotnoise_quadratic_00,
					rbin,
					shotnoise_cubic_LM,
					params.ell1, m1_,
					ylm_a, ylm_b
				);

				for (int i = 0; i < params.num_rbin; i++) {
					if (params.form == "diag") {
						shotnoise_save[i] += coupling * stats.xi[i];
					} else if (params.form == "full") {
						/// Calculate shot noise contribution equivalent to the
						/// Kronecker delta in eq. (51) in arXiv:1803.02132.
						if (i == params.ith_rbin) {
							shotnoise_save[i] += coupling * stats.xi[i];
						} else {
							shotnoise_save[i] += 0.;
						}
					}
				}

				durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed shot noise term for orders "
						"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_quadratic_00.finalise_density_field();

	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function.
	double durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing three-point correlation terms "
			"(%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		zeta_save[i] = 0.;
	}

  /// Compute three-point correlation function.
	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_fluctuation(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	dn_00.calc_fourier_transform();

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
				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM.calc_fourier_transform();
				dn_LM.apply_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (49) in arXiv:1803.02132.
				DensityField<ParticleCatalogue> F_ellm_a(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_a, ylm_a, sj1
					);
				}

				for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
					double rmag_b = rbin[i_rbin];

					if (params.form == "diag") {
						rmag_a = rmag_b;
						F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
							dn_00, rmag_a, ylm_a, sj1
						);
					}

					DensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_b, ylm_b, sj2
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh);
					for (int i = 0; i < params.nmesh; i++) {
						std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
						std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
						std::complex<double> G_LM(dn_LM[i][0], dn_LM[i][1]);
						zeta_sum += pow(I_, params.ell1 + params.ell2) * factor
							* F_ellm_1 * F_ellm_2 * G_LM;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double durationInSec = double(clock() - timeStart);
					if (currTask == 0) {
						printf(
							"[Status] :: Computed three-point correlation function term "
							"for separation and orders "
							"`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed three-point correlation function terms "
      "(... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_data.wtotal, 3);

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
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
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] zeta_save;

	return 0;
}

/**
 * Calculate three-point correlation function in a periodic box
 * and save the results.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param params Parameter set.
 * @param rbin Separation bins.
 * @returns Exit status.
 */
int calc_3pt_corrfunc_in_box(
	ParticleCatalogue& particles_data,
	ParameterSet& params,
	double* rbin
) {
	if (currTask == 0) {
		printf(
			"[Info] :: Measuring three-point correlation function "
      "in a periodic box.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for three-point correlation function measurements. "
				"Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Calculate normal density field in the global plane-parallel
	/// picture in contrast with `shotnoise_cubic_LM`.
	/// WARNING: This was inherited from Sugiyama et al. without
	/// matching equation.
	DensityField<ParticleCatalogue> density(params);
	density.calc_density_field_in_box_for_bispec(particles_data);
	density.calc_fourier_transform();

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

			DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
			dn_LM_for_shotnoise.calc_fluctuation_in_box(particles_data, params);
			dn_LM_for_shotnoise.calc_fourier_transform();

			TwoPointStatistics<ParticleCatalogue> stats(params);
			std::complex<double> shotnoise = double(particles_data.ntotal);

			stats.calc_2pt_func_for_3pt_corrfunc(
				dn_LM_for_shotnoise, density,
				rbin,
				shotnoise,
				params.ell1, m1_,
				ylm_a, ylm_b
			);

			for (int i = 0; i < params.num_rbin; i++) {
				if (params.form == "diag") {
					shotnoise_save[i] += coupling * stats.xi[i];
				} else if (params.form == "full") {
					/// Calculate shot noise contribution equivalent to the
					/// Kronecker delta in eq. (51) in arXiv:1803.02132.
					if (i == params.ith_rbin) {
						shotnoise_save[i] += coupling * stats.xi[i];
					} else {
						shotnoise_save[i] += 0.;
					}
				}
			}

			durationInSec = double(clock() - timeStart);
			if (currTask == 0) {
				printf(
					"[Status] :: Computed shot noise term for wavenumber and orders "
					"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
					m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
				);
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
					* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

	density.finalise_density_field();

	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function.
	double durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing three-point correlation function terms "
			"(%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		zeta_save[i] = 0.;
	}

	/// Compute three-point correlation function.
	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_fluctuation_in_box(particles_data, params);
	dn_00.calc_fourier_transform();

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
			DensityField<ParticleCatalogue> dn_LM(params);
			dn_LM.calc_fluctuation_in_box(particles_data, params);
			dn_LM.calc_fourier_transform();
			dn_LM.apply_assignment_compensation();
			dn_LM.calc_inverse_fourier_transform();

			/// NOTE: Standard naming convention is overriden below.

			/// Calculate F_{\ell m} in eq. (49) in arXiv:1803.02132.
			DensityField<ParticleCatalogue> F_ellm_a(params);
			double rmag_a;
			if (params.form == "full") {
				rmag_a = rbin[params.ith_rbin];
				F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
					dn_00, rmag_a, ylm_a, sj1
				);
			}

			for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
				double rmag_b = rbin[i_rbin];

				if (params.form == "diag") {
					rmag_a = rmag_b;
					F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_a, ylm_a, sj1
					);
				}

				DensityField<ParticleCatalogue> F_ellm_b(params);
				F_ellm_b.calc_inverse_fourier_transform_for_3pt_corrfunc(
					dn_00, rmag_b, ylm_b, sj2
				);

				double factor = params.volume / double(params.nmesh);
				std::complex<double> I_(0., 1.);
				std::complex<double> zeta_sum = 0.;
				for (int i = 0; i < params.nmesh; i++) {
					std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
					std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
					std::complex<double> G_LM(dn_LM[i][0], dn_LM[i][1]);
					zeta_sum += factor * pow(I_, params.ell1 + params.ell2)
						* F_ellm_1 * F_ellm_2 * G_LM;
				}

				zeta_save[i_rbin] += coupling * zeta_sum;

				double durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed three-point correlation function term "
						"for separation and orders "
						"`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
						rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytesMem -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh) / 1024. / 1024. / 1024.;
		}
	}

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed three-point correlation function terms "
      "(... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
	double norm = params.volume
		/ double(particles_data.ntotal) / double(particles_data.ntotal);
	norm *= params.volume / double(particles_data.ntotal);

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
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
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] zeta_save;

	return 0;
}

/**
 * Calculate three-point correlation function window for three-point
 * correlation function from random catalogues and save the results.
 *
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_rand Random-source particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pt_corrfunc_window(
	ParticleCatalogue& particles_rand,
	LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (currTask == 0) {
		printf(
			"[Info] :: Measuring three-point correlation function window "
			"for three-point correlation function from random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for three-point correlation function window measurements. "
				"Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
	shotnoise_quadratic_00.calc_ylm_weighted_mean_density_for_3pcf_window_shotnoise(
		particles_rand, los_rand, alpha, 0, 0
	);
	shotnoise_quadratic_00.calc_fourier_transform();

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

				DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
				dn_LM_for_shotnoise.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM_for_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise =
					stats.calc_shotnoise_for_corrfunc_window(
						particles_rand, los_rand, alpha, params.ELL, M_
					);

				stats.calc_2pt_func_for_3pt_corrfunc(
					dn_LM_for_shotnoise, shotnoise_quadratic_00,
					rbin,
					shotnoise,
					params.ell1, m1_,
					ylm_a, ylm_b
				);

				for (int i = 0; i < params.num_rbin; i++) {
					if (params.form == "diag") {
						shotnoise_save[i] += coupling * stats.xi[i];
					} else if (params.form == "full") {
						/// Calculate shot noise contribution equivalent to
						/// the Kronecker delta.  See `calc_3pt_corrfunc`.
						if (i == params.ith_rbin) {
							shotnoise_save[i] += coupling * stats.xi[i];
						} else {
							shotnoise_save[i] += 0.;
						}
					}
				}

				durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed shot noise term for orders "
						"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function window.
	double durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing three-point correlation function window terms "
			"(%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		zeta_save[i] = 0.;
	}

  /// Compute three-point correlation function window.
	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, 0, 0);
	dn_00.calc_fourier_transform();

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

				/// Calculate G_{LM}.  See `calc_3pt_corrfunc`.s
				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM.calc_fourier_transform();
				dn_LM.apply_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m}.  See `calc_3pt_corrfunc`.
				DensityField<ParticleCatalogue> F_ellm_a(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_a, ylm_a, sj1
					);
				}

				for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
					double rmag_b = rbin[i_rbin];
					if (params.form == "diag") {
						rmag_a = rmag_b;
						F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
							dn_00, rmag_a, ylm_a, sj1
						);
					}

					DensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_b, ylm_b, sj2
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh);
					for (int i = 0; i < params.nmesh; i++) {
						std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
						std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
						std::complex<double> G_LM(dn_LM[i][0], dn_LM[i][1]);
						zeta_sum += pow(I_, params.ell1 + params.ell2) * factor
							* F_ellm_1 * F_ellm_2 * G_LM;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double durationInSec = double(clock() - timeStart);
					if (currTask == 0) {
						printf(
							"[Status] :: Computed three-point correlation function window "
							"term for separation and orders "
							"`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed three-point correlation function window terms "
      "(... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_rand.wtotal, 3);
	norm /= alpha * alpha * alpha;

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d_window%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
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
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] zeta_save;

	return 0;
}

/**
 * Calculate three-point correlation function window from random
 * catalogues and save the results.  This function parallelises the
 * tasks into different processes.
 *
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_rand Random-source particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pt_corrfunc_window_mpi(
	ParticleCatalogue& particles_rand,
	LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {  // WARNING: inherited from Sugiyama et al. without matching equation
	if (currTask == 0) {
		printf(
			"[Info] :: Measuring three-point correlation function window "
			"from random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
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

	for (int i = 8; i < params.num_rbin; i++) {
		rbin[i] = rmin * exp(dlnr * (i - 8));
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save = new std::complex<double>[n_temp];
	for (int i = 0; i < n_temp; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
	shotnoise_quadratic_00.calc_ylm_weighted_mean_density_for_3pcf_window_shotnoise(
		particles_rand, los_rand, alpha, 0, 0
	);
	shotnoise_quadratic_00.calc_fourier_transform();

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

				DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
				dn_LM_for_shotnoise.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM_for_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise =
					stats.calc_shotnoise_for_corrfunc_window(
						particles_rand, los_rand, alpha, params.ELL, M_
					);

				stats.calc_2pt_func_for_3pt_corrfunc(
					dn_LM_for_shotnoise, shotnoise_quadratic_00,
					rbin,
					shotnoise,
					params.ell1, m1_,
					ylm_a, ylm_b
				);

				for (int i = 0; i < n_temp; i++) {
					if (params.form == "diag") {
						shotnoise_save[i] += coupling * stats.xi[i + NR * n_temp];
					} else if (params.form == "full") {
						/// Calculate shot noise contribution equivalent to
						/// the Kronecker delta.  See `calc_3pt_corrfunc`.
						if (i + NR * n_temp == params.ith_rbin) {
							shotnoise_save[i] += coupling * stats.xi[i + NR * n_temp];
						} else {
							shotnoise_save[i] += 0.;
						}
					}
				}

				durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed shot noise term for orders "
						"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function window.
	double durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing three-point correlation function window terms "
			"(%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* zeta_save = new std::complex<double>[n_temp];
	for (int i = 0; i < n_temp; i++) {
		zeta_save[i] = 0.;
	}

  /// Compute three-point correlation function window.
	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, 0, 0);
	dn_00.calc_fourier_transform();

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

				/// Calculate G_{LM}.  See `calc_3pt_corrfunc`.
				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM.calc_fourier_transform();
				dn_LM.apply_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m}.  See `calc_3pt_corrfunc`.
				DensityField<ParticleCatalogue> F_ellm_a(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_a, ylm_a, sj1
					);
				}

				for (int i_rbin = 0; i_rbin < n_temp; i_rbin++) {
					double rmag_b = rbin[i_rbin + NR * n_temp];

					if (params.form == "diag") {
						rmag_a = rmag_b;
						F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
							dn_00, rmag_a, ylm_a, sj1
						);
					}

					DensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_b, ylm_b, sj2
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh);
					for (int i = 0; i < params.nmesh; i++) {
						std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
						std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
						std::complex<double> G_LM(dn_LM[i][0], dn_LM[i][1]);
						zeta_sum += pow(I_, params.ell1 + params.ell2) * factor
							* F_ellm_1 * F_ellm_2 * G_LM;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double durationInSec = double(clock() - timeStart);
					if (currTask == 0) {
						printf(
							"[Status] :: Computed three-point correlation function window "
							"term for separation and orders "
							"`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed three-point correlation function window terms "
      "(... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_rand.wtotal, 3);
	norm /= alpha * alpha * alpha;

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta_window_%d%d%d_%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			NR,
			params.output_tag.c_str()
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < n_temp; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i + NR * n_temp], rbin[i + NR * n_temp],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
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
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < n_temp; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i + NR * n_temp],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

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
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_rand Random-source particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pt_corrfunc_window_for_wide_angle(
	ParticleCatalogue& particles_rand,
	LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (currTask == 0) {
		printf(
			"[Info] :: Measuring three-point correlation function window wide-angle"
			"correction terms for three-point correlation function "
			"from random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for three-point correlation function window measurements. "
				"Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
	shotnoise_quadratic_00.calc_ylm_weighted_mean_density_for_3pcf_window_shotnoise(
		particles_rand, los_rand, alpha, 0, 0
	);
	shotnoise_quadratic_00.calc_fourier_transform();

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

				DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
				dn_LM_for_shotnoise.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM_for_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise =
					stats.calc_shotnoise_for_corrfunc_window(
						particles_rand, los_rand, alpha, params.ELL, M_
					);

				stats.calc_2pt_func_for_3pt_corrfunc(
					dn_LM_for_shotnoise, shotnoise_quadratic_00,
					rbin,
					shotnoise,
					params.ell1, m1_,
					ylm_a, ylm_b
				);

				for (int i = 0; i < params.num_rbin; i++) {
					if (params.form == "diag") {
						shotnoise_save[i] += coupling * stats.xi[i];
					} else if (params.form == "full") {
						/// Calculate shot noise contribution equivalent to
						/// the Kronecker delta.  See `calc_3pt_corrfunc`.
						if (i == params.ith_rbin) {
							shotnoise_save[i] += coupling * stats.xi[i];
						} else {
							shotnoise_save[i] += 0.;
						}
					}
				}

				durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed shot noise term for orders "
						"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point correlation function window.
	double durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing three-point correlation function window "
			"wide-angle correction terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		zeta_save[i] = 0.;
	}

  /// Compute three-point correlation function window.
	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, 0, 0);
	dn_00.calc_fourier_transform();

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

				/// Calculate x^{-i-j} G_{LM}.  See `calc_3pt_corrfunc`.s
				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM.calc_fourier_transform();
				dn_LM.apply_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();
				dn_LM.apply_separation_weight_for_wide_angle();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m}.  See `calc_3pt_corrfunc`.
				DensityField<ParticleCatalogue> F_ellm_a(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_a, ylm_a, sj1
					);
				}

				for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
					double rmag_b = rbin[i_rbin];
					if (params.form == "diag") {
						rmag_a = rmag_b;
						F_ellm_a.calc_inverse_fourier_transform_for_3pt_corrfunc(
							dn_00, rmag_a, ylm_a, sj1
						);
					}

					DensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.calc_inverse_fourier_transform_for_3pt_corrfunc(
						dn_00, rmag_b, ylm_b, sj2
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh);
					for (int i = 0; i < params.nmesh; i++) {
						std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
						std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
						std::complex<double> x_G_LM(dn_LM[i][0], dn_LM[i][1]);
						zeta_sum += pow(I_, params.ell1 + params.ell2) * factor
							* F_ellm_1 * F_ellm_2 * x_G_LM;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double durationInSec = double(clock() - timeStart);
					if (currTask == 0) {
						printf(
							"[Status] :: Computed three-point correlation function window "
							"term for separation and orders "
							"`r2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed three-point correlation function window terms "
      "(... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
  double norm = pow(survey_vol_norm, 2) / pow(particles_rand.wtotal, 3);
	norm /= alpha * alpha * alpha;

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d_window-wa%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
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
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] zeta_save;

	return 0;
}

/**
 * Calculate bispectrum from catalogues with respect to a choice of
 * line of sight and save the results.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
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
			"[Info] :: Measuring bispectrum for the choice of line of sight "
			"from data and random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for bispectrum measurements. "
				"Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_kbin];
	for (int i = 0; i < params.num_kbin; i++) {
		shotnoise_save[i] = 0.;
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

				DensityField<ParticleCatalogue> dn_for_shotnoise(params);
				if (los == 0) {
					dn_for_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_for_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_for_shotnoise.calc_fourier_transform();

				/// Calculate N_LM in eq. (46) in arXiv:1803.02132.
				DensityField<ParticleCatalogue> shotnoise_quadratic(params);
				if (los == 0) {
					shotnoise_quadratic.calc_ylm_weighted_field_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				} else {
					shotnoise_quadratic.calc_ylm_weighted_field_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				}
				shotnoise_quadratic.calc_fourier_transform();

    		/// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

    		/// Calculate S_{\ell_1 \ell_2 L; i = j = k} in eq. (45)
				/// in arXiv:1803.02132.
				if (params.ell1 == 0 && params.ell2 == 0) {
					for (int i = 0; i < params.num_kbin; i++) {
						shotnoise_save[i] += coupling * shotnoise_cubic_LM;
					}
				}

    		/// Calculate S_{\ell_1 \ell_2 L; i != j = k} in eq. (45)
				/// in arXiv:1803.02132.
				if (params.ell2 == 0) {
					stats.calc_2pt_func_in_fourier(
						dn_for_shotnoise, shotnoise_quadratic,
						kbin,
						shotnoise_cubic_LM,
						params.ell1, m1_
					);
					if (params.form == "diag") {
						for (int i = 0; i < params.num_kbin; i++) {
							shotnoise_save[i] += coupling * stats.pk[i];
						}
					} else if (params.form == "full") {
						for (int i = 0; i < params.num_kbin; i++) {
							shotnoise_save[i] += coupling * stats.pk[params.ith_kbin];
						}
					}
				}

				/// WARNING: This was inherited from Sugiyama et al. without
				/// matching equation.  S_{\ell_1 \ell_2 L; i = k != j}
				/// (i.e. ell1 == 0 case) may have been shuffled.

				durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed shot noise term for orders "
						"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC);
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

				DensityField<ParticleCatalogue> dn_for_shotnoise(params);
				if (los == 1) {
					dn_for_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_for_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_for_shotnoise.calc_fourier_transform();

				/// Calculate N_LM in eq. (46) in arXiv:1803.02132.
				DensityField<ParticleCatalogue> shotnoise_quadratic(params);
				if (los == 1) {
					shotnoise_quadratic.calc_ylm_weighted_field_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				} else {
					shotnoise_quadratic.calc_ylm_weighted_field_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				}
				shotnoise_quadratic.calc_fourier_transform();

    		/// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
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
					stats.calc_2pt_func_in_fourier(
						dn_for_shotnoise, shotnoise_quadratic,
						kbin,
						shotnoise_cubic_LM,
						params.ell2, m2_
					);
					for (int i = 0; i < params.num_kbin; i++) {
						shotnoise_save[i] += coupling * stats.pk[i];
					}
				}

				durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed for orders "
						"`m1 = %d`, `m2 = %d`, `M = %d` "
						"(... %.3f seconds elapsed in total).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC);
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

				DensityField<ParticleCatalogue> shotnoise_quadratic(params);
				if (los == 2) {
					shotnoise_quadratic.calc_ylm_weighted_field_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				} else {
					shotnoise_quadratic.calc_ylm_weighted_field_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				}
				shotnoise_quadratic.calc_fourier_transform();

				DensityField<ParticleCatalogue> dn_shotnoise(params);
				if (los == 2) {
					dn_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_shotnoise.calc_fourier_transform();

    		/// Calculate S_{\ell_1 \ell_2 L; i = j != k} in eq. (45)
				/// in arXiv:1803.02132.
				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

				fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
				bytesMem += sizeof(fftw_complex)
					* double(params.nmesh) / 1024. / 1024. / 1024.;
				for (int i = 0; i < params.nmesh; i++) {
					three_pt_holder[i][0] = 0.;
					three_pt_holder[i][1] = 0.;
				}

				stats.calc_shotnoise_for_bispec_on_grid(
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

					durationInSec = double(clock() - timeStart);
					if (currTask == 0) {
						printf(
							"[Status] :: Computed shot noise term for wavenumber and orders "
							"`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	double durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
			"[Status] :: Computing bispectrum terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
  }

	std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
	for (int i = 0; i < params.num_kbin; i++) {
		bk_save[i] = 0.;
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
				DensityField<ParticleCatalogue> dn_los1(params);
				if (los == 0) {
					dn_los1.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
					  los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_los1.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_los1.calc_fourier_transform();

				DensityField<ParticleCatalogue> dn_los2(params);
				if (los == 1) {
					dn_los2.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_los2.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_los2.calc_fourier_transform();

				DensityField<ParticleCatalogue> dn_los3(params);
				if (los == 2) {
					dn_los3.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_los3.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_los3.calc_fourier_transform();
				dn_los3.apply_assignment_compensation();
				dn_los3.calc_inverse_fourier_transform();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
				DensityField<ParticleCatalogue> F_ellm_a(params);
				double kmag_a;
				double dk = kbin[1] - kbin[0];
				if (params.form == "full") {
					kmag_a = kbin[params.ith_kbin];
					F_ellm_a.calc_inverse_fourier_transform_for_bispec(
						dn_los1, kmag_a, dk, ylm_a
					);
				}

				for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
					double kmag_b = kbin[i_kbin];

					if (params.form == "diag") {
						kmag_a = kmag_b;
						F_ellm_a.calc_inverse_fourier_transform_for_bispec(
							dn_los1, kmag_a, dk, ylm_a
						);
					}

					DensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.calc_inverse_fourier_transform_for_bispec(
						dn_los2, kmag_b, dk, ylm_b
					);

					double factor = params.volume / double(params.nmesh);
					std::complex<double> bk_sum = 0.;
					for (int i = 0; i < params.nmesh; i++) {
						std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
						std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
						std::complex<double> G_LM(dn_los3[i][0], dn_los3[i][1]);
						bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
					}

					bk_save[i_kbin] += coupling * bk_sum;

					double durationInSec = double(clock() - timeStart);
					if (currTask == 0) {
						printf(
							"[Status] :: Computed bispectrum term for wavenumber and orders "
							"`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
	/// NOTE: Save the real parts only.
  double norm = pow(survey_vol_norm, 2) / pow(particles_data.wtotal, 3);

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/bk%d%d%d%s",
			params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_kbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[i], kbin[i],
				norm * (bk_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
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
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_kbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[params.ith_kbin], kbin[i],
				norm * (bk_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] bk_save;

	return 0;
}

/**
 * Calculate bispectrum from catalogues for modes with individual
 * orders @f$M@f$ and save the results.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
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
			"[Info] :: Measuring bispectrum for individual modes of order `m` "
			"from data and random catalogues.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (currTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for bispectrum measurements. "
				"Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// NOTE: Additional spaces are needed below to avoid clash with
	/// `>>` opertor.
	std::vector< std::vector< std::complex<double> > > shotnoise_save;

	shotnoise_save.resize(2*params.ELL + 1);
	for (int i = 0; i < params.num_kbin; i++) {
		shotnoise_save[i].resize(params.num_kbin);
	}

	for (int m_idx = 0; m_idx < 2*params.ELL + 1; m_idx++) {
		for (int i = 0; i < params.num_kbin; i++) {
			shotnoise_save[m_idx][i] = 0.;
		}
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> dn_00_for_shotnoise(params);
	dn_00_for_shotnoise.calc_ylm_weighted_fluctuation(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	dn_00_for_shotnoise.calc_fourier_transform();

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
				DensityField<ParticleCatalogue> shotnoise_quadratic_LM(params);
				shotnoise_quadratic_LM.calc_ylm_weighted_field_for_bispec_shotnoise(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				shotnoise_quadratic_LM.calc_fourier_transform();

    		/// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

    		/// Calculate S_{\ell_1 \ell_2 L; i = j = k} in eq. (45)
				/// in arXiv:1803.02132.
				if (params.ell1 == 0 && params.ell2 == 0) {
					for (int i = 0; i < params.num_kbin; i++) {
						shotnoise_save[M_ + params.ELL][i] += shotnoise_cubic_LM;
					}
				}

    		/// Calculate S_{\ell_1 \ell_2 L; i != j = k} in eq. (45)
				/// in arXiv:1803.02132.
				if (params.ell2 == 0) {
					stats.calc_2pt_func_in_fourier(
						dn_00_for_shotnoise, shotnoise_quadratic_LM,
						kbin,
						shotnoise_cubic_LM,
						params.ell1, m1_
					);
					if (params.form == "diag") {
						for (int i = 0; i < params.num_kbin; i++) {
							shotnoise_save[M_ + params.ELL][i] += stats.pk[i];
						}
					} else if (params.form == "full") {
						for (int i = 0; i < params.num_kbin; i++) {
							shotnoise_save[M_ + params.ELL][i] += stats.pk[params.ith_kbin];
						}
					}
				}

    		/// Calculate S_{\ell_1 \ell_2 L; i = k != j} in eq. (45)
				/// in arXiv:1803.02132.
				if (params.ell1 == 0) {
					stats.calc_2pt_func_in_fourier(
						dn_00_for_shotnoise, shotnoise_quadratic_LM,
						kbin,
						shotnoise_cubic_LM,
						params.ell2, m2_
					);
					for (int i = 0; i < params.num_kbin; i++) {
						shotnoise_save[M_ + params.ELL][i] += stats.pk[i];
					}
				}

				durationInSec = double(clock() - timeStart);
				if (currTask == 0) {
					printf(
						"[Status] :: Computed shot noise term for orders "
						"`m1 = %d`, `m2 = %d`, `M = %d` (%.3f seconds elapsed).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}
		}
	}

	dn_00_for_shotnoise.finalise_density_field();

	/// Calculate N_00 in eq. (45) in arXiv:1803.02132.
	DensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
	shotnoise_quadratic_00.calc_ylm_weighted_field_for_bispec_shotnoise(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	shotnoise_quadratic_00.calc_fourier_transform();

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

				DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
				dn_LM_for_shotnoise.calc_ylm_weighted_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM_for_shotnoise.calc_fourier_transform();

    		/// Calculate S_{\ell_1 \ell_2 L; i = j != k} in eq. (45)
				/// in arXiv:1803.02132.
				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

				fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
				bytesMem += sizeof(fftw_complex) *
					double(params.nmesh) / 1024. / 1024. / 1024.;
				for (int i = 0; i < params.nmesh; i++) {
					three_pt_holder[i][0] = 0.;
					three_pt_holder[i][1] = 0.;
				}

				stats.calc_shotnoise_for_bispec_on_grid(
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

					durationInSec = double(clock() - timeStart);
					if (currTask == 0) {
						printf(
							"[Status] :: Computed shot noise term for wavenumber and orders "
							"`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

	durationInSec = double(clock() - timeStart);
	if (currTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	double durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
			"[Status] :: Computing bispectrum terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
  }

	/// NOTE: Additional spaces are needed below to avoid clash with
	/// `>>` opertor.
	std::vector< std::vector< std::complex<double> > > bk_save;

  /// Compute bispectrum.
	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_fluctuation(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	dn_00.calc_fourier_transform();

	bk_save.resize(2*params.ELL + 1);
	for (int i = 0; i < params.num_kbin; i++) {
		bk_save[i].resize(params.num_kbin);
	}

	for (int m = 0; m < 2*params.ELL + 1; m++) {
		for (int i = 0; i < params.num_kbin; i++) {
			bk_save[m][i] = 0.;
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

				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM.calc_fourier_transform();
				dn_LM.apply_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
				DensityField<ParticleCatalogue> F_ellm_a(params);
				double kmag_a;
				double dk = kbin[1] - kbin[0];
				if (params.form == "full") {
					kmag_a = kbin[params.ith_kbin];
					F_ellm_a.calc_inverse_fourier_transform_for_bispec(
						dn_00, kmag_a, dk, ylm_a
					);
				}

				for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
					double kmag_b = kbin[i_kbin];

					if (params.form == "diag") {
						kmag_a = kmag_b;
						F_ellm_a.calc_inverse_fourier_transform_for_bispec(
							dn_00, kmag_a, dk, ylm_a
						);
					}

					DensityField<ParticleCatalogue> F_ellm_b(params);
					F_ellm_b.calc_inverse_fourier_transform_for_bispec(
						dn_00, kmag_b, dk, ylm_b
					);

					double factor = params.volume / double(params.nmesh);
					std::complex<double> bk_sum = 0.;
					for (int i = 0; i < params.nmesh; i++) {
						std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
						std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
						std::complex<double> G_LM(dn_LM[i][0], dn_LM[i][1]);
						bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
					}

					/// No coupling multiplication.
					bk_save[M_ + params.ELL][i_kbin] += bk_sum;

					double durationInSec = double(clock() - timeStart);
					if (currTask == 0) {
						printf(
							"[Status] :: Computed bispectrum term for wavenumber and orders "
							"`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
							kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

	/// Normalise and then save the output.
	/// NOTE: Save the real parts only.
  double norm = pow(survey_vol_norm, 2) / pow(particles_data.wtotal, 3);

	for (int M_ = 0; M_<2*params.ELL + 1; M_++) {
		FILE* saved_file_ptr;
		char buf[1024];
		if (params.form == "diag") {
			sprintf(
				buf, "%s/bk%d%d%d_M%d%s",
				params.measurement_dir.c_str(),
				params.ell1, params.ell2,
				params.ELL, M_,
				params.output_tag.c_str()
			);
			saved_file_ptr = fopen(buf, "w");
			for (int i = 0; i < params.num_kbin; i++) {
				fprintf(
					saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
					kbin[i], kbin[i],
					norm * (bk_save[M_][i].real() - shotnoise_save[M_][i].real()),
					norm * shotnoise_save[M_][i].real()
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
			saved_file_ptr = fopen(buf, "w");
			for (int i = 0; i < params.num_kbin; i++) {
				fprintf(
					saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
					kbin[params.ith_kbin], kbin[i],
					norm * (bk_save[M_][i].real() - shotnoise_save[M_][i].real()),
					norm * shotnoise_save[M_][i].real()
				);
			}
		}
		fclose(saved_file_ptr);
	}

	return 0;
}

/**
 * Calculate bispectrum from catalogues and save the results.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
 * @param params Parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec_(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  std::vector<std::vector<double> > los_data_arr,
	std::vector<std::vector<double> > los_rand_arr,
  ParameterSet& params,
  double alpha,
  double* kbin
) {
  if (currTask == 0) {
    printf(
			"[Info] :: Measuring bispectrum from data and random catalogues.\n"
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

	ParticleCatalogue::align_catalogues_for_fft(
		particles_data, particles_rand, params.boxsize, params.ngrid
	);

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
  ) {
    if (currTask == 0) {
      printf(
        "[Error] :: Disallowed multipole degree combination "
        "for bispectrum measurements. "
        "Please ensure `wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0`.\n"
      );
    }
    exit(1);
  }

  /// Initialise output shot noise terms.
  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computing shot noise terms (%.3f seconds elapsed...).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* shotnoise_save =
    new std::complex<double>[params.num_kbin];
  for (int i = 0; i < params.num_kbin; i++) {
    shotnoise_save[i] = 0.;
  }

  /// Compute shot noise terms.
  DensityField<ParticleCatalogue> dn_00_for_shotnoise(params);
  dn_00_for_shotnoise.calc_ylm_weighted_fluctuation(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  dn_00_for_shotnoise.calc_fourier_transform();

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
        DensityField<ParticleCatalogue> shotnoise_quadratic_LM(params);
        shotnoise_quadratic_LM.calc_ylm_weighted_field_for_bispec_shotnoise(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        shotnoise_quadratic_LM.calc_fourier_transform();

    		/// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
        TwoPointStatistics<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

    		/// Calculate S_{\ell_1 \ell_2 L; i = j = k} in eq. (45)
				/// in arXiv:1803.02132.
        if (params.ell1 == 0 && params.ell2 == 0) {
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * shotnoise_cubic_LM;
          }
        }

    		/// Calculate S_{\ell_1 \ell_2 L; i != j = k} in eq. (45)
				/// in arXiv:1803.02132.
        if (params.ell2 == 0) {
          stats.calc_2pt_func_in_fourier(
            dn_00_for_shotnoise, shotnoise_quadratic_LM,
						kbin,
						shotnoise_cubic_LM,
						params.ell1, m1_
          );
          if (params.form == "diag") {
            for (int i = 0; i < params.num_kbin; i++) {
              shotnoise_save[i] += coupling * stats.pk[i];
            }
          } else if (params.form == "full") {
            for (int i = 0; i < params.num_kbin; i++) {
              shotnoise_save[i] += coupling * stats.pk[params.ith_kbin];
            }
          }
        }

    		/// Calculate S_{\ell_1 \ell_2 L; i = k != j} in eq. (45)
				/// in arXiv:1803.02132.
        if (params.ell1 == 0) {
          stats.calc_2pt_func_in_fourier(
            dn_00_for_shotnoise, shotnoise_quadratic_LM,
						kbin,
						shotnoise_cubic_LM,
						params.ell2, m2_
          );
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * stats.pk[i];
          }
        }

        durationInSec = double(clock() - timeStart);
        if (currTask == 0) {
          printf(
            "[Status] :: Computed shot noise term for orders "
            "`m1 = %d`, `m2 = %d`, `M = %d` "
						"(%.3f seconds elapsed).\n",
            m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
          );
        }
      }
    }
  }

	// MARKER
  dn_00_for_shotnoise.finalise_density_field();

	/// Calculate N_00 in eq. (45) in arXiv:1803.02132.
  DensityField<ParticleCatalogue> shotnoise_quadratic_00(params);
  shotnoise_quadratic_00.calc_ylm_weighted_field_for_bispec_shotnoise(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  shotnoise_quadratic_00.calc_fourier_transform();

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

        DensityField<ParticleCatalogue> dn_LM_for_shotnoise(params);
        dn_LM_for_shotnoise.calc_ylm_weighted_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM_for_shotnoise.calc_fourier_transform();

    		/// Calculate S_{\ell_1 \ell_2 L; i = j != k} in eq. (45)
				/// in arXiv:1803.02132.
        TwoPointStatistics<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise_cubic_LM =
					stats.calc_shotnoise_for_bispec_from_self(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);

        fftw_complex* three_pt_holder = fftw_alloc_complex(params.nmesh);
        bytesMem += sizeof(fftw_complex) *
          double(params.nmesh) / 1024. / 1024. / 1024.;
        for (int i = 0; i < params.nmesh; i++) {
          three_pt_holder[i][0] = 0.;
          three_pt_holder[i][1] = 0.;
        }

        stats.calc_shotnoise_for_bispec_on_grid(
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

          durationInSec = double(clock() - timeStart);
          if (currTask == 0) {
            printf(
              "[Status] :: Computed shot noise term for wavenumber and orders "
              "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
              kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed shot noise terms (... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  /// Initialise output bispectrum.
	double durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
			"[Status] :: Computing bispectrum terms (%.3f seconds elapsed...).\n",
			durationInSec / CLOCKS_PER_SEC
		);
  }

  std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
  for (int i = 0; i < params.num_kbin; i++) {
    bk_save[i] = 0.;
  }

  /// Compute bispectrum.
  DensityField<ParticleCatalogue> dn_00(params);
  dn_00.calc_ylm_weighted_fluctuation(
    particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
  );
  dn_00.calc_fourier_transform();

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
        DensityField<ParticleCatalogue> dn_LM(params);
        dn_LM.calc_ylm_weighted_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM.calc_fourier_transform();
        dn_LM.apply_assignment_compensation();
        dn_LM.calc_inverse_fourier_transform();

				/// NOTE: Standard naming convention is overriden below.

				/// Calculate F_{\ell m} in eq. (42) in arXiv:1803.02132.
        DensityField<ParticleCatalogue> F_ellm_a(params);
        double kmag_a;
        double dk = kbin[1] - kbin[0];
        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          F_ellm_a.calc_inverse_fourier_transform_for_bispec(
            dn_00, kmag_a, dk, ylm_a
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];

          if (params.form == "diag") {
            kmag_a = kmag_b;
            F_ellm_a.calc_inverse_fourier_transform_for_bispec(
              dn_00, kmag_a, dk, ylm_a
            );
          }

          DensityField<ParticleCatalogue> F_ellm_b(params);
          F_ellm_b.calc_inverse_fourier_transform_for_bispec(
            dn_00, kmag_b, dk, ylm_b
          );

          double factor = params.volume / double(params.nmesh);
          std::complex<double> bk_sum = 0.;
          for (int i = 0; i < params.nmesh; i++) {
            std::complex<double> F_ellm_1(F_ellm_a[i][0], F_ellm_a[i][1]);
            std::complex<double> F_ellm_2(F_ellm_b[i][0], F_ellm_b[i][1]);
            std::complex<double> G_LM(dn_LM[i][0], dn_LM[i][1]);
            bk_sum += factor * F_ellm_1 * F_ellm_2 * G_LM;
          }

          bk_save[i_kbin] += coupling * bk_sum;

          double durationInSec = double(clock() - timeStart);
          if (currTask == 0) {
            printf(
              "[Status] :: Computed bispectrum term for wavenumber and orders "
              "`k2 = %.3f`, `m1 = %d`, `m2 = %d`, `M = %d` "
							"(%.3f seconds elapsed).\n",
              kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
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

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: Computed bispectrum terms (... %.3f seconds elapsed).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  /// Normalise and then save the output.
  FILE* saved_file_ptr;

  char buf[1024];

  if (params.form == "diag") {
    sprintf(
      buf, "%s/bk%d%d%d%s",
      params.measurement_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.output_tag.c_str()
    );

    saved_file_ptr = fopen(buf, "w");
    for (int i = 0; i < params.num_kbin; i++) {
      fprintf(
        saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[i], kbin[i],
        (bk_save[i].real() - shotnoise_save[i].real()),
        (bk_save[i].imag() - shotnoise_save[i].imag()),
        shotnoise_save[i].real(),
        shotnoise_save[i].imag()
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

    saved_file_ptr = fopen(buf, "w");
    for (int i = 0; i < params.num_kbin; i++) {
      fprintf(
        saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[params.ith_kbin], kbin[i],
        (bk_save[i].real() - shotnoise_save[i].real()),
        (bk_save[i].imag() - shotnoise_save[i].imag()),
        shotnoise_save[i].real(),
        shotnoise_save[i].imag()
      );
    }
  }
  fclose(saved_file_ptr);

  delete[] shotnoise_save;
  delete[] bk_save;

  return 0;
}

#endif
