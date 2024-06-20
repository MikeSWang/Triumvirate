// Copyright (C) [GPLv3 Licence]
//
// This file is part of the Triumvirate program. See the COPYRIGHT
// and LICENCE files at the top-level directory of this distribution
// for details of copyright and licensing.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

/**
 * @file threept.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "threept.hpp"

namespace trvs = trv::sys;
namespace trvm = trv::maths;

namespace trv {

// ***********************************************************************
// Coupling coefficients
// ***********************************************************************

double calc_coupling_coeff_3pt(
  int ell1, int ell2, int ELL, int m1, int m2, int M
) {
  return double(2*ell1 + 1) * double(2*ell2 + 1) * double(2*ELL + 1)
    * trvm::wigner_3j(ell1, ell2, ELL, 0, 0, 0)
    * trvm::wigner_3j(ell1, ell2, ELL, m1, m2, M);
}

void validate_multipole_coupling(trv::ParameterSet& params) {
  double coupling_ = trvm::wigner_3j(
    params.ell1, params.ell2, params.ELL, 0, 0, 0
  );
  if (std::fabs(coupling_) < trvm::eps_coupling) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Specified three-point correlator multipole "
        "vanishes identically owing to zero-valued Wigner 3-j symbol."
      );
    }
    throw trvs::InvalidParameterError(
      "Specified three-point correlator multipole "
      "vanishes identically owing to zero-valued Wigner 3-j symbol.\n"
    );
  }
}


// ***********************************************************************
// Normalisation
// ***********************************************************************

double calc_bispec_normalisation_from_particles(
  ParticleCatalogue& particles, double alpha
) {
  if (particles.pdata == nullptr) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Particle data are uninitialised.");
    }
    throw trvs::InvalidDataError("Particle data are uninitialised.\n");
  }

  double norm = 0.;  // I₃

#ifdef TRV_USE_OMP
#pragma omp parallel for simd reduction(+:norm)
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    norm += particles[pid].ws
      * std::pow(particles[pid].nz, 2) * std::pow(particles[pid].wc, 3);
  }

  if (norm == 0.) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Particle 'nz' values appear to be all zeros. "
        "Check the input catalogue contains valid 'nz' field."
      );
    }
    throw trvs::InvalidDataError(
      "Particle 'nz' values appear to be all zeros. "
      "Check the input catalogue contains valid 'nz' field.\n"
    );
  }

  double norm_factor = 1. / (alpha * norm);  // 1/I₃

  return norm_factor;
}

double calc_bispec_normalisation_from_mesh(
  ParticleCatalogue& particles, trv::ParameterSet& params, double alpha
) {
  MeshField catalogue_mesh(params, false, "`catalogue_mesh`");

  double norm_factor =
    catalogue_mesh.calc_grid_based_powlaw_norm(particles, 3);

  norm_factor /= std::pow(alpha, 3);

  return norm_factor;
}


// ***********************************************************************
// Shot noise
// ***********************************************************************

std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_bispec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  double alpha, int ell, int m
) {
  double sn_data_real = 0., sn_data_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:sn_data_real, sn_data_imag)
#endif
  for (int pid = 0; pid < particles_data.ntotal; pid++) {
    double los_[3] = {
      los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    std::complex<double> sn_part = ylm * std::pow(particles_data[pid].w, 3);
    double sn_part_real = sn_part.real();
    double sn_part_imag = sn_part.imag();

    sn_data_real += sn_part_real;
    sn_data_imag += sn_part_imag;
  }

  std::complex<double> sn_data(sn_data_real, sn_data_imag);

  double sn_rand_real = 0., sn_rand_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:sn_rand_real, sn_rand_imag)
#endif
  for (int pid = 0; pid < particles_rand.ntotal; pid++) {
    double los_[3] = {
      los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    std::complex<double> sn_part = ylm * std::pow(particles_rand[pid].w, 3);
    double sn_part_real = sn_part.real();
    double sn_part_imag = sn_part.imag();

    sn_rand_real += sn_part_real;
    sn_rand_imag += sn_part_imag;
  }

  std::complex<double> sn_rand(sn_rand_real, sn_rand_imag);

  return sn_data + std::pow(alpha, 3) * sn_rand;
}

std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_bispec(
  ParticleCatalogue& particles, LineOfSight* los,
  double alpha, int ell, int m
) {
  double sn_real = 0., sn_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:sn_real, sn_imag)
#endif
  for (int pid = 0; pid < particles.ntotal; pid++) {
    double los_[3] = {los[pid].pos[0], los[pid].pos[1], los[pid].pos[2]};

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    std::complex<double> sn_part = ylm * std::pow(particles[pid].w, 3);
    double sn_part_real = sn_part.real();
    double sn_part_imag = sn_part.imag();

    sn_real += sn_part_real;
    sn_imag += sn_part_imag;
  }

  std::complex<double> sn(sn_real, sn_imag);

  return std::pow(alpha, 3) * sn;
}


// ***********************************************************************
// Full statistics
// ***********************************************************************

// STYLE: Standard naming convention is not always followed for
// intermediary quantities in the functions below.

// Hereafter 'the Paper' refers to Sugiyama et al. (2019) [1803.02132].

trv::BispecMeasurements compute_bispec(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing bispectrum from paired survey-type catalogues..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  double alpha = catalogue_data.wstotal / catalogue_rand.wstotal;

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int dv_dim = 0;  // data vector dimension
  if (params.shape == "diag" || params.shape == "row" ) {
    dv_dim = kbinning.num_bins;
  } else
  if (params.shape == "off-diag") {
    dv_dim = kbinning.num_bins - std::abs(params.idx_bin);
  } else
  if (params.shape == "full") {
    dv_dim = kbinning.num_bins * kbinning.num_bins;
  } else
  if (params.shape == "triu") {
    dv_dim = kbinning.num_bins * (kbinning.num_bins + 1) / 2;
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Three-point statistic form is not recognised: `form` = '%s'.",
        params.form.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Three-point statistic form is not recognised: `form` = '%s'.\n",
      params.form.c_str()
    );
  }

  int* nmodes1_dv = new int[dv_dim];
  int* nmodes2_dv = new int[dv_dim];
  double* k1bin_dv = new double[dv_dim];
  double* k2bin_dv = new double[dv_dim];
  double* k1eff_dv = new double[dv_dim];
  double* k2eff_dv = new double[dv_dim];
  std::complex<double>* bk_dv = new std::complex<double>[dv_dim];
  std::complex<double>* sn_dv = new std::complex<double>[dv_dim];
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    nmodes1_dv[idx_dv] = 0;
    nmodes2_dv[idx_dv] = 0;
    k1bin_dv[idx_dv] = 0.;
    k2bin_dv[idx_dv] = 0.;
    k1eff_dv[idx_dv] = 0.;
    k2eff_dv[idx_dv] = 0.;
    bk_dv[idx_dv] = 0.;
    sn_dv[idx_dv] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params, true, "`dn_00`");  // δn_00(k)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  MeshField& dn_00_for_sn = dn_00;  // δn_00(k) (for shot noise)

  double vol_cell = dn_00.vol_cell;

  MeshField N_00(params, true, "`N_00`");  // N_00(k)
  N_00.compute_ylm_wgtd_quad_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  N_00.fourier_transform();

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  FieldStats stats_sn(params);

  // Initialise reduced-spherical-harmonic weights on mesh grids.
  std::vector< std::complex<double> > ylm_k_a(params.nmesh);
  std::vector< std::complex<double> > ylm_k_b(params.nmesh);
  std::vector< std::complex<double> > ylm_r_a(params.nmesh);
  std::vector< std::complex<double> > ylm_r_b(params.nmesh);
  trvs::count_cgrid += 4;
  trvs::count_grid += 4;
  trvs::update_maxcntgrid();
  trvs::gbytesMem +=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);
  trvs::update_maxmem();

  // Compute bispectrum terms including shot noise.
  int count_terms = 0;
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Check for if all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );
        if (std::fabs(coupling) > trvm::eps_coupling) {
          flag_vanishing = "false";
          break;
        }
      }
      if (flag_vanishing == "true") {continue;}

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        // Calculate the coupling coefficient.
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );  // Wigner 3-j's
        if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

        // ·······························································
        // Raw bispectrum
        // ·······························································

        // Compute bispectrum components in eqs. (41) & (42) in the Paper.
        MeshField G_LM(params, true, "`G_LM`");  // G_LM
        G_LM.compute_ylm_wgtd_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        G_LM.fourier_transform();
        G_LM.apply_assignment_compensation();
        G_LM.inv_fourier_transform();

        MeshField F_lm_a(params, true, "`F_lm_a`");  // F_lm_a
        MeshField F_lm_b(params, true, "`F_lm_b`");  // F_lm_b

        if (params.shape == "diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;

            double k_lower = kbinning.bin_edges[ibin];
            double k_upper = kbinning.bin_edges[ibin + 1];

            double k_eff_a_, k_eff_b_;
            int nmodes_a_, nmodes_b_;

            F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
            );
            F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_b, k_lower, k_upper, k_eff_b_, nmodes_b_
            );

            if (count_terms == 0) {
              k1bin_dv[idx_dv] = kbinning.bin_centres[ibin];
              k2bin_dv[idx_dv] = kbinning.bin_centres[ibin];
              k1eff_dv[idx_dv] = k_eff_a_;
              k2eff_dv[idx_dv] = k_eff_b_;
              nmodes1_dv[idx_dv] = nmodes_a_;
              nmodes2_dv[idx_dv] = nmodes_b_;
            }

            // B_{l₁ l₂ L}^{m₁ m₂ M}
            double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
            for (long long gid = 0; gid < params.nmesh; gid++) {
              std::complex<double> F_lm_a_gridpt(
                F_lm_a[gid][0], F_lm_a[gid][1]
              );
              std::complex<double> F_lm_b_gridpt(
                F_lm_b[gid][0], F_lm_b[gid][1]
              );
              std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
              std::complex<double> bk_gridpt =
                F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

              bk_comp_real += bk_gridpt.real();
              bk_comp_imag += bk_gridpt.imag();
            }

            std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

            bk_dv[idx_dv] += coupling * vol_cell * bk_component;
          }
        }

        if (params.shape == "off-diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_row, ibin_col;
            if (params.idx_bin >= 0) {
              ibin_row = idx_dv;
              ibin_col = idx_dv + std::abs(params.idx_bin);
            } else {
              ibin_row = idx_dv + std::abs(params.idx_bin);
              ibin_col = idx_dv;
            }

            double k_lower_a = kbinning.bin_edges[ibin_row];
            double k_upper_a = kbinning.bin_edges[ibin_row + 1];
            double k_lower_b = kbinning.bin_edges[ibin_col];
            double k_upper_b = kbinning.bin_edges[ibin_col + 1];

            double k_eff_a_, k_eff_b_;
            int nmodes_a_, nmodes_b_;

            F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
            );
            F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
            );

            if (count_terms == 0) {
              k1bin_dv[idx_dv] = kbinning.bin_centres[ibin_row];
              k2bin_dv[idx_dv] = kbinning.bin_centres[ibin_col];
              k1eff_dv[idx_dv] = k_eff_a_;
              k2eff_dv[idx_dv] = k_eff_b_;
              nmodes1_dv[idx_dv] = nmodes_a_;
              nmodes2_dv[idx_dv] = nmodes_b_;
            }

            // B_{l₁ l₂ L}^{m₁ m₂ M}
            double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
            for (long long gid = 0; gid < params.nmesh; gid++) {
              std::complex<double> F_lm_a_gridpt(
                F_lm_a[gid][0], F_lm_a[gid][1]
              );
              std::complex<double> F_lm_b_gridpt(
                F_lm_b[gid][0], F_lm_b[gid][1]
              );
              std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
              std::complex<double> bk_gridpt =
                F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

              bk_comp_real += bk_gridpt.real();
              bk_comp_imag += bk_gridpt.imag();
            }

            std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

            bk_dv[idx_dv] += coupling * vol_cell * bk_component;
          }
        }

        if (params.shape == "row") {
          int ibin_row = params.idx_bin;

          double k_lower_a = kbinning.bin_edges[ibin_row];
          double k_upper_a = kbinning.bin_edges[ibin_row + 1];

          double k_eff_a_;
          int nmodes_a_;

          F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
          );

          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_col = idx_dv;

            double k_lower_b = kbinning.bin_edges[ibin_col];
            double k_upper_b = kbinning.bin_edges[ibin_col + 1];

            double k_eff_b_;
            int nmodes_b_;

            F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
            );

            if (count_terms == 0) {
              k1bin_dv[idx_dv] = kbinning.bin_centres[ibin_row];
              k2bin_dv[idx_dv] = kbinning.bin_centres[ibin_col];
              k1eff_dv[idx_dv] = k_eff_a_;
              k2eff_dv[idx_dv] = k_eff_b_;
              nmodes1_dv[idx_dv] = nmodes_a_;
              nmodes2_dv[idx_dv] = nmodes_b_;
            }

            // B_{l₁ l₂ L}^{m₁ m₂ M}
            double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
            for (long long gid = 0; gid < params.nmesh; gid++) {
              std::complex<double> F_lm_a_gridpt(
                F_lm_a[gid][0], F_lm_a[gid][1]
              );
              std::complex<double> F_lm_b_gridpt(
                F_lm_b[gid][0], F_lm_b[gid][1]
              );
              std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
              std::complex<double> bk_gridpt =
                F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

              bk_comp_real += bk_gridpt.real();
              bk_comp_imag += bk_gridpt.imag();
            }

            std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

            bk_dv[idx_dv] += coupling * vol_cell * bk_component;
          }
        }

        if (params.shape == "full") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
              int idx_dv = idx_row * params.num_bins + idx_col;

              double k_lower_a = kbinning.bin_edges[idx_row];
              double k_upper_a = kbinning.bin_edges[idx_row + 1];
              double k_lower_b = kbinning.bin_edges[idx_col];
              double k_upper_b = kbinning.bin_edges[idx_col + 1];

              double k_eff_a_, k_eff_b_;
              int nmodes_a_, nmodes_b_;

              F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
                dn_00, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
              );
              F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
                dn_00, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
              );

              if (count_terms == 0) {
                k1bin_dv[idx_dv] = kbinning.bin_centres[idx_row];
                k2bin_dv[idx_dv] = kbinning.bin_centres[idx_col];
                k1eff_dv[idx_dv] = k_eff_a_;
                k2eff_dv[idx_dv] = k_eff_b_;
                nmodes1_dv[idx_dv] = nmodes_a_;
                nmodes2_dv[idx_dv] = nmodes_b_;
              }

              // B_{l₁ l₂ L}^{m₁ m₂ M}
              double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
              for (long long gid = 0; gid < params.nmesh; gid++) {
                std::complex<double> F_lm_a_gridpt(
                  F_lm_a[gid][0], F_lm_a[gid][1]
                );
                std::complex<double> F_lm_b_gridpt(
                  F_lm_b[gid][0], F_lm_b[gid][1]
                );
                std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
                std::complex<double> bk_gridpt =
                  F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

                bk_comp_real += bk_gridpt.real();
                bk_comp_imag += bk_gridpt.imag();
              }

              std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

              bk_dv[idx_dv] += coupling * vol_cell * bk_component;
            }
          }
        }

        if (params.shape == "triu") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
              int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                + (idx_col - idx_row);

              double k_lower_a = kbinning.bin_edges[idx_row];
              double k_upper_a = kbinning.bin_edges[idx_row + 1];
              double k_lower_b = kbinning.bin_edges[idx_col];
              double k_upper_b = kbinning.bin_edges[idx_col + 1];

              double k_eff_a_, k_eff_b_;
              int nmodes_a_, nmodes_b_;

              F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
                dn_00, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
              );
              F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
                dn_00, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
              );

              if (count_terms == 0) {
                k1bin_dv[idx_dv] = kbinning.bin_centres[idx_row];
                k2bin_dv[idx_dv] = kbinning.bin_centres[idx_col];
                k1eff_dv[idx_dv] = k_eff_a_;
                k2eff_dv[idx_dv] = k_eff_b_;
                nmodes1_dv[idx_dv] = nmodes_a_;
                nmodes2_dv[idx_dv] = nmodes_b_;
              }

              // B_{l₁ l₂ L}^{m₁ m₂ M}
              double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
              for (long long gid = 0; gid < params.nmesh; gid++) {
                std::complex<double> F_lm_a_gridpt(
                  F_lm_a[gid][0], F_lm_a[gid][1]
                );
                std::complex<double> F_lm_b_gridpt(
                  F_lm_b[gid][0], F_lm_b[gid][1]
                );
                std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
                std::complex<double> bk_gridpt =
                  F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

                bk_comp_real += bk_gridpt.real();
                bk_comp_imag += bk_gridpt.imag();
              }

              std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

              bk_dv[idx_dv] += coupling * vol_cell * bk_component;
            }
          }
        }

        // ·······························································
        // Shot noise
        // ·······························································

        // Compute shot noise components in eqs. (45) & (46) in the Paper.
        MeshField dn_LM_for_sn(params, true, "`dn_LM_for_sn`");  // δn_LM(k)
                                                                 // (for
                                                                 // shot noise)
        dn_LM_for_sn.compute_ylm_wgtd_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM_for_sn.fourier_transform();

        MeshField N_LM(params, true, "`N_LM`");  // N_LM(k)
        N_LM.compute_ylm_wgtd_quad_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        N_LM.fourier_transform();

        std::complex<double> Sbar_LM = calc_ylm_wgtd_shotnoise_amp_for_bispec(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );  // \bar{S}_LM

        if (params.ell1 == 0 && params.ell2 == 0) {
          // When l₁ = l₂ = 0, the Wigner 3-j symbol enforces L = 0
          // and the pre-factors involving degrees and orders become 1.
          std::complex<double> S_ijk = coupling * Sbar_LM;  // S|{i = j = k}
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            sn_dv[idx_dv] += S_ijk;
          }
        }

        if (params.ell2 == 0) {  // S|{i ≠ j = k}
          // When l₂ = 0, the Wigner 3-j symbol enforces L = l₁.
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_sn, N_LM, Sbar_LM, params.ell1, m1_, kbinning
          );

          if (params.shape == "diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin = idx_dv;
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin] - stats_sn.sn[ibin]
              );
            }
          }

          if (params.shape == "off-diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_row;
              if (params.idx_bin >= 0) {
                ibin_row = idx_dv;
              } else {
                ibin_row = idx_dv + std::abs(params.idx_bin);
              }
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin_row] - stats_sn.sn[ibin_row]
              );
            }
          }

          if (params.shape == "row") {
            std::complex<double> sn_row_ = coupling * (
              stats_sn.pk[params.idx_bin] - stats_sn.sn[params.idx_bin]
            );
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              sn_dv[idx_dv] += sn_row_;
            }
          }

          if (params.shape == "full") {
            for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
              std::complex<double> sn_row_ = coupling * (
                stats_sn.pk[idx_row] - stats_sn.sn[idx_row]
              );
              for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
                int idx_dv = idx_row * params.num_bins + idx_col;
                sn_dv[idx_dv] += sn_row_;
              }
            }
          }

          if (params.shape == "triu") {
            for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
              std::complex<double> sn_row_ = coupling * (
                stats_sn.pk[idx_row] - stats_sn.sn[idx_row]
              );
              for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
                int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                  + (idx_col - idx_row);
                sn_dv[idx_dv] += sn_row_;
              }
            }
          }
        }

        if (params.ell1 == 0) {  // S|{j ≠ i = k}
          // When l₁ = 0, the Wigner 3-j symbol enforces L = l₂.
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_sn, N_LM, Sbar_LM, params.ell2, m2_, kbinning
          );

          if (params.shape == "diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin = idx_dv;
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin] - stats_sn.sn[ibin]
              );
            }
          }

          if (params.shape == "off-diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_col;
              if (params.idx_bin >= 0) {
                ibin_col = idx_dv + params.idx_bin;
              } else {
                ibin_col = idx_dv;
              }
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin_col] - stats_sn.sn[ibin_col]
              );
            }
          }

          if (params.shape == "row") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_col = idx_dv;
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin_col] - stats_sn.sn[ibin_col]
              );
            }
          }

          if (params.shape == "full") {
            for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
              std::complex<double> sn_col_ = coupling * (
                stats_sn.pk[idx_col] - stats_sn.sn[idx_col]
              );
              for (int idx_row = 0; idx_row <= params.num_bins; idx_row++) {
                int idx_dv = idx_row * params.num_bins + idx_col;
                sn_dv[idx_dv] += sn_col_;
              }
            }
          }

          if (params.shape == "triu") {
            for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
              std::complex<double> sn_col_ = coupling * (
                stats_sn.pk[idx_col] - stats_sn.sn[idx_col]
              );
              for (int idx_row = 0; idx_row <= idx_col; idx_row++) {
                int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                  + (idx_col - idx_row);
                sn_dv[idx_dv] += sn_col_;
              }
            }
          }
        }

        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          double k_a = k1eff_dv[idx_dv];
          double k_b = k2eff_dv[idx_dv];

          std::complex<double> S_ij_k = parity *
            stats_sn.compute_uncoupled_shotnoise_for_bispec_per_bin(
              dn_LM_for_sn, N_00, ylm_r_a, ylm_r_b, sj_a, sj_b,
              Sbar_LM, k_a, k_b
            );  // S|{i = j ≠ k}

          sn_dv[idx_dv] += coupling * S_ij_k;
        }

        count_terms++;
        if (trvs::currTask == 0) {
          trvs::logger.stat(
            "Bispectrum term computed at orders (m1, m2, M) = (%d, %d, %d).",
            m1_, m2_, M_
          );
        }
      }
    }
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::BispecMeasurements bispec_out;
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    bispec_out.k1_bin.push_back(k1bin_dv[idx_dv]);
    bispec_out.k1_eff.push_back(k1eff_dv[idx_dv]);
    bispec_out.nmodes_1.push_back(nmodes1_dv[idx_dv]);
    bispec_out.k2_bin.push_back(k2bin_dv[idx_dv]);
    bispec_out.k2_eff.push_back(k2eff_dv[idx_dv]);
    bispec_out.nmodes_2.push_back(nmodes2_dv[idx_dv]);
    bispec_out.bk_raw.push_back(norm_factor * bk_dv[idx_dv]);
    bispec_out.bk_shot.push_back(norm_factor * sn_dv[idx_dv]);
  }
  bispec_out.dim = dv_dim;

  delete[] nmodes1_dv; delete[] nmodes2_dv;
  delete[] k1bin_dv; delete[] k2bin_dv;
  delete[] k1eff_dv; delete[] k2eff_dv;
  delete[] bk_dv; delete[] sn_dv;

  trvs::count_cgrid -= 4;
  trvs::count_grid -= 4;
  trvs::gbytesMem -=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed bispectrum from paired survey-type catalogues."
    );
  }

  return bispec_out;
}

trv::ThreePCFMeasurements compute_3pcf(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing three-point correlation function "
      "from paired survey-type catalogues..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  double alpha = catalogue_data.wstotal / catalogue_rand.wstotal;

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int dv_dim = 0;  // data vector dimension
  if (params.shape == "diag" || params.shape == "row" ) {
    dv_dim = rbinning.num_bins;
  } else
  if (params.shape == "off-diag") {
    dv_dim = rbinning.num_bins - std::abs(params.idx_bin);
  } else
  if (params.shape == "full") {
    dv_dim = rbinning.num_bins * rbinning.num_bins;
  } else
  if (params.shape == "triu") {
    dv_dim = rbinning.num_bins * (rbinning.num_bins + 1) / 2;
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Three-point statistic form is not recognised: `form` = '%s'.",
        params.form.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Three-point statistic form is not recognised: `form` = '%s'.\n",
      params.form.c_str()
    );
  }

  int* npairs1_dv = new int[dv_dim];
  int* npairs2_dv = new int[dv_dim];
  double* r1bin_dv = new double[dv_dim];
  double* r2bin_dv = new double[dv_dim];
  double* r1eff_dv = new double[dv_dim];
  double* r2eff_dv = new double[dv_dim];
  std::complex<double>* zeta_dv = new std::complex<double>[dv_dim];
  std::complex<double>* sn_dv = new std::complex<double>[dv_dim];
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    npairs1_dv[idx_dv] = 0;
    npairs2_dv[idx_dv] = 0;
    r1bin_dv[idx_dv] = 0.;
    r2bin_dv[idx_dv] = 0.;
    r1eff_dv[idx_dv] = 0.;
    r2eff_dv[idx_dv] = 0.;
    zeta_dv[idx_dv] = 0.;
    sn_dv[idx_dv] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params, true, "`dn_00`");  // δn_00(k)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  double vol_cell = dn_00.vol_cell;

  MeshField N_00(params, true, "`N_00`");  // N_00(k)
  N_00.compute_ylm_wgtd_quad_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  N_00.fourier_transform();

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  FieldStats stats_sn(params);

  // Initialise reduced-spherical-harmonic weights on mesh grids.
  std::vector< std::complex<double> > ylm_r_a(params.nmesh);
  std::vector< std::complex<double> > ylm_r_b(params.nmesh);
  std::vector< std::complex<double> > ylm_k_a(params.nmesh);
  std::vector< std::complex<double> > ylm_k_b(params.nmesh);
  trvs::count_cgrid += 4;
  trvs::count_grid += 4;
  trvs::update_maxcntgrid();
  trvs::gbytesMem +=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);
  trvs::update_maxmem();

  // Compute 3PCF terms including shot noise.
  int count_terms = 0;
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Check for vanishing cases where all Wigner-3j symbols are zero.
      // Check for if all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );
        if (std::fabs(coupling) > trvm::eps_coupling) {
          flag_vanishing = "false";
          break;
        }
      }
      if (flag_vanishing == "true") {continue;}

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        // Calculate the coupling coefficient.
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );  // Wigner 3-j's
        if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

        // ·······························································
        // Shot noise
        // ·······························································

        // Compute shot noise components in eq. (51) in the Paper.
        MeshField dn_LM_for_sn(params, true, "`dn_LM_for_sn`");  // δn_LM(k)
                                                                 // (for
                                                                 // shot noise)
        dn_LM_for_sn.compute_ylm_wgtd_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM_for_sn.fourier_transform();

        std::complex<double> Sbar_LM = calc_ylm_wgtd_shotnoise_amp_for_bispec(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );  // \bar{S}_LM

        stats_sn.compute_uncoupled_shotnoise_for_3pcf(
          dn_LM_for_sn, N_00, ylm_r_a, ylm_r_b, Sbar_LM, rbinning
        );  // S|{i = j ≠ k}

        // Enforce the Kronecker delta in eq. (51) in the Paper.
        if (params.shape == "diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;
            sn_dv[idx_dv] += coupling * stats_sn.xi[ibin];
          }
        }

        if (params.shape == "off-diag" && params.idx_bin == 0) {
          // Note that ``idx_col == idx_row``.
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;
            sn_dv[idx_dv] += coupling * stats_sn.xi[ibin];
          }
        }

        if (params.shape == "row") {
          sn_dv[params.idx_bin] += coupling * stats_sn.xi[params.idx_bin];
        }

        if (params.shape == "full") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            // Note that ``idx_col == idx_row``.
            int idx_dv = idx_row * params.num_bins + idx_row;
            sn_dv[idx_dv] += coupling * stats_sn.xi[idx_row];
          }
        }

        if (params.shape == "triu") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            // Note that ``idx_col == idx_row``.
            int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2;
            sn_dv[idx_dv] += coupling * stats_sn.xi[idx_row];
          }
        }

        // Only record the binned coordinates and counts once.
        if (count_terms == 0) {
          if (params.shape == "diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin = idx_dv;

              r1bin_dv[idx_dv] = rbinning.bin_centres[ibin];
              r2bin_dv[idx_dv] = rbinning.bin_centres[ibin];
              r1eff_dv[idx_dv] = stats_sn.r[ibin];
              r2eff_dv[idx_dv] = stats_sn.r[ibin];
              npairs1_dv[idx_dv] = stats_sn.npairs[ibin];
              npairs2_dv[idx_dv] = stats_sn.npairs[ibin];
            }
          }

          if (params.shape == "off-diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_row, ibin_col;
              if (params.idx_bin >= 0) {
                ibin_row = idx_dv;
                ibin_col = idx_dv + std::abs(params.idx_bin);
              } else {
                ibin_row = idx_dv + std::abs(params.idx_bin);
                ibin_col = idx_dv;
              }

              r1bin_dv[idx_dv] = rbinning.bin_centres[ibin_row];
              r2bin_dv[idx_dv] = rbinning.bin_centres[ibin_col];
              r1eff_dv[idx_dv] = stats_sn.r[ibin_row];
              r2eff_dv[idx_dv] = stats_sn.r[ibin_col];
              npairs1_dv[idx_dv] = stats_sn.npairs[ibin_row];
              npairs2_dv[idx_dv] = stats_sn.npairs[ibin_col];
            }
          }

          if (params.shape == "row") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_row = params.idx_bin;
              int ibin_col = idx_dv;

              r1bin_dv[idx_dv] = rbinning.bin_centres[ibin_row];
              r2bin_dv[idx_dv] = rbinning.bin_centres[ibin_col];
              r1eff_dv[idx_dv] = stats_sn.r[ibin_row];
              r2eff_dv[idx_dv] = stats_sn.r[ibin_col];
              npairs1_dv[idx_dv] = stats_sn.npairs[ibin_row];
              npairs2_dv[idx_dv] = stats_sn.npairs[ibin_col];
            }
          }

          if (params.shape == "full") {
            for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
              for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
                int idx_dv = idx_row * params.num_bins + idx_col;

                r1bin_dv[idx_dv] = rbinning.bin_centres[idx_row];
                r2bin_dv[idx_dv] = rbinning.bin_centres[idx_col];
                r1eff_dv[idx_dv] = stats_sn.r[idx_row];
                r2eff_dv[idx_dv] = stats_sn.r[idx_col];
                npairs1_dv[idx_dv] = stats_sn.npairs[idx_row];
                npairs2_dv[idx_dv] = stats_sn.npairs[idx_col];
              }
            }
          }

          if (params.shape == "triu") {
            for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
              for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
                int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                + (idx_col - idx_row);

                r1bin_dv[idx_dv] = rbinning.bin_centres[idx_row];
                r2bin_dv[idx_dv] = rbinning.bin_centres[idx_col];
                r1eff_dv[idx_dv] = stats_sn.r[idx_row];
                r2eff_dv[idx_dv] = stats_sn.r[idx_col];
                npairs1_dv[idx_dv] = stats_sn.npairs[idx_row];
                npairs2_dv[idx_dv] = stats_sn.npairs[idx_col];
              }
            }
          }
        }

        // ·······························································
        // Raw 3PCF
        // ·······························································

        // Compute 3PCF components in eqs. (42), (48) & (49) in the Paper.
        MeshField G_LM(params, true, "`G_LM`");  // G_LM
        G_LM.compute_ylm_wgtd_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        G_LM.fourier_transform();
        G_LM.apply_assignment_compensation();
        G_LM.inv_fourier_transform();

        MeshField F_lm_a(params, true, "`F_lm_a`");  // F_lm_a
        MeshField F_lm_b(params, true, "`F_lm_b`");  // F_lm_b

        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          double r_a = r1eff_dv[idx_dv];
          F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
            dn_00, ylm_k_a, sj_a, r_a
          );

          double r_b = r2eff_dv[idx_dv];
          F_lm_b.inv_fourier_transform_sjl_ylm_wgtd_field(
            dn_00, ylm_k_b, sj_b, r_b
          );

          // ζ_{l₁ l₂ L}^{m₁ m₂ M}
          double zeta_comp_real = 0., zeta_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:zeta_comp_real, zeta_comp_imag)
#endif  // TRV_USE_OMP
          for (long long gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
            std::complex<double> zeta_gridpt =
              F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

            zeta_comp_real += zeta_gridpt.real();
            zeta_comp_imag += zeta_gridpt.imag();
          }

          std::complex<double> zeta_component(zeta_comp_real, zeta_comp_imag);

          zeta_dv[idx_dv] += parity * coupling * vol_cell * zeta_component;
        }

        count_terms++;
        if (trvs::currTask == 0) {
          trvs::logger.stat(
            "Three-point correlation function term computed at orders "
            "(m1, m2, M) = (%d, %d, %d).",
            m1_, m2_, M_
          );
        }
      }
    }
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::ThreePCFMeasurements threepcf_out;
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    threepcf_out.r1_bin.push_back(r1bin_dv[idx_dv]);
    threepcf_out.r1_eff.push_back(r1eff_dv[idx_dv]);
    threepcf_out.npairs_1.push_back(npairs1_dv[idx_dv]);
    threepcf_out.r2_bin.push_back(r2bin_dv[idx_dv]);
    threepcf_out.r2_eff.push_back(r2eff_dv[idx_dv]);
    threepcf_out.npairs_2.push_back(npairs2_dv[idx_dv]);
    threepcf_out.zeta_raw.push_back(norm_factor * zeta_dv[idx_dv]);
    threepcf_out.zeta_shot.push_back(norm_factor * sn_dv[idx_dv]);
  }
  threepcf_out.dim = dv_dim;

  delete[] npairs1_dv; delete[] npairs2_dv;
  delete[] r1bin_dv; delete[] r2bin_dv;
  delete[] r1eff_dv; delete[] r2eff_dv;
  delete[] zeta_dv; delete[] sn_dv;

  trvs::count_cgrid -= 4;
  trvs::count_grid -= 4;
  trvs::gbytesMem -=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed three-point correlation function "
      "from paired survey-type catalogues."
    );
  }

  return threepcf_out;
}

trv::BispecMeasurements compute_bispec_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning kbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing bispectrum from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int dv_dim = 0;  // data vector dimension
  if (params.shape == "diag" || params.shape == "row" ) {
    dv_dim = kbinning.num_bins;
  } else
  if (params.shape == "off-diag") {
    dv_dim = kbinning.num_bins - std::abs(params.idx_bin);
  } else
  if (params.shape == "full") {
    dv_dim = kbinning.num_bins * kbinning.num_bins;
  } else
  if (params.shape == "triu") {
    dv_dim = kbinning.num_bins * (kbinning.num_bins + 1) / 2;
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Three-point statistic form is not recognised: `form` = '%s'.",
        params.form.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Three-point statistic form is not recognised: `form` = '%s'.\n",
      params.form.c_str()
    );
  }

  int* nmodes1_dv = new int[dv_dim];
  int* nmodes2_dv = new int[dv_dim];
  double* k1bin_dv = new double[dv_dim];
  double* k2bin_dv = new double[dv_dim];
  double* k1eff_dv = new double[dv_dim];
  double* k2eff_dv = new double[dv_dim];
  std::complex<double>* bk_dv = new std::complex<double>[dv_dim];
  std::complex<double>* sn_dv = new std::complex<double>[dv_dim];
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    nmodes1_dv[idx_dv] = 0;
    nmodes2_dv[idx_dv] = 0;
    k1bin_dv[idx_dv] = 0.;
    k2bin_dv[idx_dv] = 0.;
    k1eff_dv[idx_dv] = 0.;
    k2eff_dv[idx_dv] = 0.;
    bk_dv[idx_dv] = 0.;
    sn_dv[idx_dv] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params, true, "`dn_00`");  // δn_00(k)
  dn_00.compute_unweighted_field_fluctuations_insitu(catalogue_data);
  dn_00.fourier_transform();

  MeshField& dn_00_for_sn = dn_00;  // δn_00(k) (for shot noise)
  MeshField& dn_L0_for_sn = dn_00;  // δn_L0(k) (for shot noise)

  double vol_cell = dn_00.vol_cell;

  // Under the global plane-parallel approximation, y_{LM} = δᴰ_{M0}
  // (L-invariant) for the line-of-sight spherical harmonic.
  MeshField N_L0(params, true, "`N_L0`");  // N_L0(k)
  N_L0.compute_unweighted_field(catalogue_data);
  N_L0.fourier_transform();

  MeshField& N_00 = N_L0;  // N_00(k)

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  FieldStats stats_sn(params);

  // Initialise/reset spherical harmonic mesh grids.
  std::vector< std::complex<double> > ylm_k_a(params.nmesh);
  std::vector< std::complex<double> > ylm_k_b(params.nmesh);
  std::vector< std::complex<double> > ylm_r_a(params.nmesh);
  std::vector< std::complex<double> > ylm_r_b(params.nmesh);
  trvs::count_cgrid += 4;
  trvs::count_grid += 4;
  trvs::update_maxcntgrid();
  trvs::gbytesMem +=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);
  trvs::update_maxmem();

  // Compute bispectrum terms including shot noise.
  int count_terms = 0;
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Under the global plane-parallel approximation, δᴰ_{M0} enforces
      // M = 0 for any spherical-harmonic-weighted field fluctuations.
      int M_ = 0;

      // Calculate the coupling coefficient.
      double coupling = trv::calc_coupling_coeff_3pt(
        params.ell1, params.ell2, params.ELL, m1_, m2_, M_
      );  // Wigner 3-j's
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator
        ::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );

      // ·································································
      // Raw bispectrum
      // ·································································

      MeshField G_00(params, true, "`G_00`");  // G_00
      G_00.compute_unweighted_field_fluctuations_insitu(catalogue_data);
      G_00.fourier_transform();
      G_00.apply_assignment_compensation();
      G_00.inv_fourier_transform();

      MeshField F_lm_a(params, true, "`F_lm_a`");  // F_lm_a
      MeshField F_lm_b(params, true, "`F_lm_b`");  // F_lm_b

      if (params.shape == "diag") {
        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          int ibin = idx_dv;

          double k_lower = kbinning.bin_edges[ibin];
          double k_upper = kbinning.bin_edges[ibin + 1];

          double k_eff_a_, k_eff_b_;
          int nmodes_a_, nmodes_b_;

          F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
          );
          F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_b, k_lower, k_upper, k_eff_b_, nmodes_b_
          );

          if (count_terms == 0) {
            k1bin_dv[idx_dv] = kbinning.bin_centres[ibin];
            k2bin_dv[idx_dv] = kbinning.bin_centres[ibin];
            k1eff_dv[idx_dv] = k_eff_a_;
            k2eff_dv[idx_dv] = k_eff_b_;
            nmodes1_dv[idx_dv] = nmodes_a_;
            nmodes2_dv[idx_dv] = nmodes_b_;
          }

          // B_{l₁ l₂ L}^{m₁ m₂ M}
          double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
          for (long long gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_00_gridpt(G_00[gid][0], G_00[gid][1]);
            std::complex<double> bk_gridpt =
              F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;

            bk_comp_real += bk_gridpt.real();
            bk_comp_imag += bk_gridpt.imag();
          }

          std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

          bk_dv[idx_dv] += coupling * vol_cell * bk_component;
        }
      }

      if (params.shape == "off-diag") {
        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          int ibin_row, ibin_col;
          if (params.idx_bin >= 0) {
            ibin_row = idx_dv;
            ibin_col = idx_dv + std::abs(params.idx_bin);
          } else {
            ibin_row = idx_dv + std::abs(params.idx_bin);
            ibin_col = idx_dv;
          }

          double k_lower_a = kbinning.bin_edges[ibin_row];
          double k_upper_a = kbinning.bin_edges[ibin_row + 1];
          double k_lower_b = kbinning.bin_edges[ibin_col];
          double k_upper_b = kbinning.bin_edges[ibin_col + 1];

          double k_eff_a_, k_eff_b_;
          int nmodes_a_, nmodes_b_;

          F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
          );
          F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
          );

          if (count_terms == 0) {
            k1bin_dv[idx_dv] = kbinning.bin_centres[ibin_row];
            k2bin_dv[idx_dv] = kbinning.bin_centres[ibin_col];
            k1eff_dv[idx_dv] = k_eff_a_;
            k2eff_dv[idx_dv] = k_eff_b_;
            nmodes1_dv[idx_dv] = nmodes_a_;
            nmodes2_dv[idx_dv] = nmodes_b_;
          }

          // B_{l₁ l₂ L}^{m₁ m₂ M}
          double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
          for (long long gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_00_gridpt(G_00[gid][0], G_00[gid][1]);
            std::complex<double> bk_gridpt =
              F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;

            bk_comp_real += bk_gridpt.real();
            bk_comp_imag += bk_gridpt.imag();
          }

          std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

          bk_dv[idx_dv] += coupling * vol_cell * bk_component;
        }
      }

      if (params.shape == "row") {
        int ibin_row = params.idx_bin;

        double k_lower_a = kbinning.bin_edges[ibin_row];
        double k_upper_a = kbinning.bin_edges[ibin_row + 1];

        double k_eff_a_;
        int nmodes_a_;

        F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
          dn_00, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
        );

        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          int ibin_col = idx_dv;

          double k_lower_b = kbinning.bin_edges[ibin_col];
          double k_upper_b = kbinning.bin_edges[ibin_col + 1];

          double k_eff_b_;
          int nmodes_b_;

          F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
          );

          if (count_terms == 0) {
            k1bin_dv[idx_dv] = kbinning.bin_centres[ibin_row];
            k2bin_dv[idx_dv] = kbinning.bin_centres[ibin_col];
            k1eff_dv[idx_dv] = k_eff_a_;
            k2eff_dv[idx_dv] = k_eff_b_;
            nmodes1_dv[idx_dv] = nmodes_a_;
            nmodes2_dv[idx_dv] = nmodes_b_;
          }

          // B_{l₁ l₂ L}^{m₁ m₂ M}
          double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
          for (long long gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_00_gridpt(G_00[gid][0], G_00[gid][1]);
            std::complex<double> bk_gridpt =
              F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;

            bk_comp_real += bk_gridpt.real();
            bk_comp_imag += bk_gridpt.imag();
          }

          std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

          bk_dv[idx_dv] += coupling * vol_cell * bk_component;
        }
      }

      if (params.shape == "full") {
        for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
          for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
            int idx_dv = idx_row * params.num_bins + idx_col;

            double k_lower_a = kbinning.bin_edges[idx_row];
            double k_upper_a = kbinning.bin_edges[idx_row + 1];
            double k_lower_b = kbinning.bin_edges[idx_col];
            double k_upper_b = kbinning.bin_edges[idx_col + 1];

            double k_eff_a_, k_eff_b_;
            int nmodes_a_, nmodes_b_;

            F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
            );
            F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
            );

            if (count_terms == 0) {
              k1bin_dv[idx_dv] = kbinning.bin_centres[idx_row];
              k2bin_dv[idx_dv] = kbinning.bin_centres[idx_col];
              k1eff_dv[idx_dv] = k_eff_a_;
              k2eff_dv[idx_dv] = k_eff_b_;
              nmodes1_dv[idx_dv] = nmodes_a_;
              nmodes2_dv[idx_dv] = nmodes_b_;
            }

            // B_{l₁ l₂ L}^{m₁ m₂ M}
            double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
            for (long long gid = 0; gid < params.nmesh; gid++) {
              std::complex<double> F_lm_a_gridpt(
                F_lm_a[gid][0], F_lm_a[gid][1]
              );
              std::complex<double> F_lm_b_gridpt(
                F_lm_b[gid][0], F_lm_b[gid][1]
              );
              std::complex<double> G_00_gridpt(G_00[gid][0], G_00[gid][1]);
              std::complex<double> bk_gridpt =
                F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;

              bk_comp_real += bk_gridpt.real();
              bk_comp_imag += bk_gridpt.imag();
            }

            std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

            bk_dv[idx_dv] += coupling * vol_cell * bk_component;
          }
        }
      }

      if (params.shape == "triu") {
        for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
          for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
            int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
              + (idx_col - idx_row);

            double k_lower_a = kbinning.bin_edges[idx_row];
            double k_upper_a = kbinning.bin_edges[idx_row + 1];
            double k_lower_b = kbinning.bin_edges[idx_col];
            double k_upper_b = kbinning.bin_edges[idx_col + 1];

            double k_eff_a_, k_eff_b_;
            int nmodes_a_, nmodes_b_;

            F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
            );
            F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
            );

            if (count_terms == 0) {
              k1bin_dv[idx_dv] = kbinning.bin_centres[idx_row];
              k2bin_dv[idx_dv] = kbinning.bin_centres[idx_col];
              k1eff_dv[idx_dv] = k_eff_a_;
              k2eff_dv[idx_dv] = k_eff_b_;
              nmodes1_dv[idx_dv] = nmodes_a_;
              nmodes2_dv[idx_dv] = nmodes_b_;
            }

            // B_{l₁ l₂ L}^{m₁ m₂ M}
            double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
            for (long long gid = 0; gid < params.nmesh; gid++) {
              std::complex<double> F_lm_a_gridpt(
                F_lm_a[gid][0], F_lm_a[gid][1]
              );
              std::complex<double> F_lm_b_gridpt(
                F_lm_b[gid][0], F_lm_b[gid][1]
              );
              std::complex<double> G_00_gridpt(G_00[gid][0], G_00[gid][1]);
              std::complex<double> bk_gridpt =
                F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;

              bk_comp_real += bk_gridpt.real();
              bk_comp_imag += bk_gridpt.imag();
            }

            std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

            bk_dv[idx_dv] += coupling * vol_cell * bk_component;
          }
        }
      }

      // ·································································
      // Shot noise
      // ·································································

      // Under the global plane-parallel approximation, y_{LM} = δᴰ_{M0}
      // (L-invariant) for the line-of-sight spherical harmonic.
      // Also note the field is unweighted from simulation sources.
      std::complex<double> Sbar_LM =
        double(catalogue_data.ntotal);  // \bar{S}_LM
      std::complex<double> Sbar_L0 = Sbar_LM;  // \bar{S}_L0

      if (params.ell1 == 0 && params.ell2 == 0) {
        std::complex<double> S_ijk = coupling * Sbar_LM;  // S|{i = j = k}
        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          sn_dv[idx_dv] += S_ijk;
        }
      }

      if (params.ell2 == 0) {  // S|{i ≠ j = k}
        // When l₂ = 0, the Wigner 3-j symbol enforces L = l₁.
        stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
          dn_00_for_sn, N_L0, Sbar_LM, params.ell1, m1_, kbinning
        );

        if (params.shape == "diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;
            sn_dv[idx_dv] += coupling * (
              stats_sn.pk[ibin] - stats_sn.sn[ibin]
            );
          }
        }

        if (params.shape == "off-diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_row;
            if (params.idx_bin >= 0) {
              ibin_row = idx_dv;
            } else {
              ibin_row = idx_dv + std::abs(params.idx_bin);
            }
            sn_dv[idx_dv] += coupling * (
              stats_sn.pk[ibin_row] - stats_sn.sn[ibin_row]
            );
          }
        }

        if (params.shape == "row") {
          std::complex<double> sn_row_ = coupling * (
            stats_sn.pk[params.idx_bin] - stats_sn.sn[params.idx_bin]
          );
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            sn_dv[idx_dv] += sn_row_;
          }
        }

        if (params.shape == "full") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            std::complex<double> sn_row_ = coupling * (
              stats_sn.pk[idx_row] - stats_sn.sn[idx_row]
            );
            for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
              int idx_dv = idx_row * params.num_bins + idx_col;
              sn_dv[idx_dv] += sn_row_;
            }
          }
        }

        if (params.shape == "triu") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            std::complex<double> sn_row_ = coupling * (
              stats_sn.pk[idx_row] - stats_sn.sn[idx_row]
            );
            for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
              int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                + (idx_col - idx_row);
              sn_dv[idx_dv] += sn_row_;
            }
          }
        }
      }

      if (params.ell1 == 0) {  // S|{j ≠ i = k}
        // When l₁ = 0, the Wigner 3-j symbol enforces L = l₂.
        stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
          dn_00_for_sn, N_L0, Sbar_LM, params.ell2, m2_, kbinning
        );

        if (params.shape == "diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;
            sn_dv[idx_dv] += coupling * (
              stats_sn.pk[ibin] - stats_sn.sn[ibin]
            );
          }
        }

        if (params.shape == "off-diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_col;
            if (params.idx_bin >= 0) {
              ibin_col = idx_dv + params.idx_bin;
            } else {
              ibin_col = idx_dv;
            }
            sn_dv[idx_dv] += coupling * (
              stats_sn.pk[ibin_col] - stats_sn.sn[ibin_col]
            );
          }
        }

        if (params.shape == "row") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_col = idx_dv;
            sn_dv[idx_dv] += coupling * (
              stats_sn.pk[ibin_col] - stats_sn.sn[ibin_col]
            );
          }
        }

        if (params.shape == "full") {
          for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
            std::complex<double> sn_col_ = coupling * (
              stats_sn.pk[idx_col] - stats_sn.sn[idx_col]
            );
            for (int idx_row = 0; idx_row <= params.num_bins; idx_row++) {
              int idx_dv = idx_row * params.num_bins + idx_col;
              sn_dv[idx_dv] += sn_col_;
            }
          }
        }

        if (params.shape == "triu") {
          for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
            std::complex<double> sn_col_ = coupling * (
              stats_sn.pk[idx_col] - stats_sn.sn[idx_col]
            );
            for (int idx_row = 0; idx_row <= idx_col; idx_row++) {
              int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                + (idx_col - idx_row);
              sn_dv[idx_dv] += sn_col_;
            }
          }
        }
      }

      for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
        double k_a = k1eff_dv[idx_dv];
        double k_b = k2eff_dv[idx_dv];

        std::complex<double> S_ij_k = parity *
          stats_sn.compute_uncoupled_shotnoise_for_bispec_per_bin(
            dn_L0_for_sn, N_00, ylm_r_a, ylm_r_b, sj_a, sj_b,
            Sbar_L0, k_a, k_b
          );  // S|{i = j ≠ k}

        sn_dv[idx_dv] += coupling * S_ij_k;
      }

      count_terms++;
      if (trvs::currTask == 0) {
        trvs::logger.stat(
          "Bispectrum term computed at orders (m1, m2, M) = (%d, %d, 0).",
          m1_, m2_
        );
      }
    }
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::BispecMeasurements bispec_out;
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    bispec_out.k1_bin.push_back(k1bin_dv[idx_dv]);
    bispec_out.k1_eff.push_back(k1eff_dv[idx_dv]);
    bispec_out.nmodes_1.push_back(nmodes1_dv[idx_dv]);
    bispec_out.k2_bin.push_back(k2bin_dv[idx_dv]);
    bispec_out.k2_eff.push_back(k2eff_dv[idx_dv]);
    bispec_out.nmodes_2.push_back(nmodes2_dv[idx_dv]);
    bispec_out.bk_raw.push_back(norm_factor * bk_dv[idx_dv]);
    bispec_out.bk_shot.push_back(norm_factor * sn_dv[idx_dv]);
  }
  bispec_out.dim = dv_dim;

  delete[] nmodes1_dv; delete[] nmodes2_dv;
  delete[] k1bin_dv; delete[] k2bin_dv;
  delete[] k1eff_dv; delete[] k2eff_dv;
  delete[] bk_dv; delete[] sn_dv;

  trvs::count_cgrid -= 4;
  trvs::count_grid -= 4;
  trvs::gbytesMem -=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed bispectrum from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation."
    );
  }

  return bispec_out;
}

trv::ThreePCFMeasurements compute_3pcf_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing three-point correlation function "
      "from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int dv_dim = 0;  // data vector dimension
  if (params.shape == "diag" || params.shape == "row" ) {
    dv_dim = rbinning.num_bins;
  } else
  if (params.shape == "off-diag") {
    dv_dim = rbinning.num_bins - std::abs(params.idx_bin);
  } else
  if (params.shape == "full") {
    dv_dim = rbinning.num_bins * rbinning.num_bins;
  } else
  if (params.shape == "triu") {
    dv_dim = rbinning.num_bins * (rbinning.num_bins + 1) / 2;
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Three-point statistic form is not recognised: `form` = '%s'.",
        params.form.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Three-point statistic form is not recognised: `form` = '%s'.\n",
      params.form.c_str()
    );
  }

  int* npairs1_dv = new int[dv_dim];
  int* npairs2_dv = new int[dv_dim];
  double* r1bin_dv = new double[dv_dim];
  double* r2bin_dv = new double[dv_dim];
  double* r1eff_dv = new double[dv_dim];
  double* r2eff_dv = new double[dv_dim];
  std::complex<double>* zeta_dv = new std::complex<double>[dv_dim];
  std::complex<double>* sn_dv = new std::complex<double>[dv_dim];
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    npairs1_dv[idx_dv] = 0;
    npairs2_dv[idx_dv] = 0;
    r1bin_dv[idx_dv] = 0.;
    r2bin_dv[idx_dv] = 0.;
    r1eff_dv[idx_dv] = 0.;
    r2eff_dv[idx_dv] = 0.;
    zeta_dv[idx_dv] = 0.;
    sn_dv[idx_dv] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params, true, "`dn_00`");  // δn_00(k)
  dn_00.compute_unweighted_field_fluctuations_insitu(catalogue_data);
  dn_00.fourier_transform();

  MeshField& dn_L0_for_sn = dn_00;  // δn_L0(k)

  double vol_cell = dn_00.vol_cell;

  MeshField N_00(params, true, "`N_00`");  // N_00(k)
  N_00.compute_unweighted_field(catalogue_data);
  N_00.fourier_transform();

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  FieldStats stats_sn(params);

  // Initialise/reset spherical harmonic mesh grids.
  std::vector< std::complex<double> > ylm_r_a(params.nmesh);
  std::vector< std::complex<double> > ylm_r_b(params.nmesh);
  std::vector< std::complex<double> > ylm_k_a(params.nmesh);
  std::vector< std::complex<double> > ylm_k_b(params.nmesh);
  trvs::count_cgrid += 4;
  trvs::count_grid += 4;
  trvs::update_maxcntgrid();
  trvs::gbytesMem +=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);
  trvs::update_maxmem();

  // Compute 3PCF terms including shot noise.
  int count_terms = 0;
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Under the global plane-parallel approximation, δᴰ_{M0} enforces
      // M = 0 for any spherical-harmonic-weighted field fluctuations.
      int M_ = 0;

      // Calculate the coupling coefficient.
      double coupling = trv::calc_coupling_coeff_3pt(
        params.ell1, params.ell2, params.ELL, m1_, m2_, M_
      );  // Wigner 3-j's
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );

      // ·································································
      // Shot noise
      // ·································································

      // Under the global plane-parallel approximation, y_{LM} = δᴰ_{M0}
      // (L-invariant) for the line-of-sight spherical harmonic.
      // Also note the field is unweighted from simulation sources.
      std::complex<double> Sbar_L0 =
        double(catalogue_data.ntotal);  // \bar{S}_L0

      stats_sn.compute_uncoupled_shotnoise_for_3pcf(
        dn_L0_for_sn, N_00, ylm_r_a, ylm_r_b, Sbar_L0, rbinning
      );  // S|{i = j ≠ k}

      // Enforce the Kronecker delta in eq. (51) in the Paper.
      if (params.shape == "diag") {
        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          int ibin = idx_dv;
          sn_dv[idx_dv] += coupling * stats_sn.xi[ibin];
        }
      }

      if (params.shape == "off-diag" && params.idx_bin == 0) {
        // Note that ``idx_col == idx_row``.
        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          int ibin = idx_dv;
          sn_dv[idx_dv] += coupling * stats_sn.xi[ibin];
        }
      }

      if (params.shape == "row") {
        sn_dv[params.idx_bin] += coupling * stats_sn.xi[params.idx_bin];
      }

      if (params.shape == "full") {
        for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
          // Note that ``idx_col == idx_row``.
          int idx_dv = idx_row * params.num_bins + idx_row;
          sn_dv[idx_dv] += coupling * stats_sn.xi[idx_row];
        }
      }

      if (params.shape == "triu") {
        for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
          // Note that ``idx_col == idx_row``.
          int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2;
          sn_dv[idx_dv] += coupling * stats_sn.xi[idx_row];
        }
      }

      // Only record the binned coordinates and counts once.
      if (count_terms == 0) {
        if (params.shape == "diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;

            r1bin_dv[idx_dv] = rbinning.bin_centres[ibin];
            r2bin_dv[idx_dv] = rbinning.bin_centres[ibin];
            r1eff_dv[idx_dv] = stats_sn.r[ibin];
            r2eff_dv[idx_dv] = stats_sn.r[ibin];
            npairs1_dv[idx_dv] = stats_sn.npairs[ibin];
            npairs2_dv[idx_dv] = stats_sn.npairs[ibin];
          }
        }

        if (params.shape == "off-diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_row, ibin_col;
            if (params.idx_bin >= 0) {
              ibin_row = idx_dv;
              ibin_col = idx_dv + std::abs(params.idx_bin);
            } else {
              ibin_row = idx_dv + std::abs(params.idx_bin);
              ibin_col = idx_dv;
            }

            r1bin_dv[idx_dv] = rbinning.bin_centres[ibin_row];
            r2bin_dv[idx_dv] = rbinning.bin_centres[ibin_col];
            r1eff_dv[idx_dv] = stats_sn.r[ibin_row];
            r2eff_dv[idx_dv] = stats_sn.r[ibin_col];
            npairs1_dv[idx_dv] = stats_sn.npairs[ibin_row];
            npairs2_dv[idx_dv] = stats_sn.npairs[ibin_col];
          }
        }

        if (params.shape == "row") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_row = params.idx_bin;
            int ibin_col = idx_dv;

            r1bin_dv[idx_dv] = rbinning.bin_centres[ibin_row];
            r2bin_dv[idx_dv] = rbinning.bin_centres[ibin_col];
            r1eff_dv[idx_dv] = stats_sn.r[ibin_row];
            r2eff_dv[idx_dv] = stats_sn.r[ibin_col];
            npairs1_dv[idx_dv] = stats_sn.npairs[ibin_row];
            npairs2_dv[idx_dv] = stats_sn.npairs[ibin_col];
          }
        }

        if (params.shape == "full") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
              int idx_dv = idx_row * params.num_bins + idx_col;

              r1bin_dv[idx_dv] = rbinning.bin_centres[idx_row];
              r2bin_dv[idx_dv] = rbinning.bin_centres[idx_col];
              r1eff_dv[idx_dv] = stats_sn.r[idx_row];
              r2eff_dv[idx_dv] = stats_sn.r[idx_col];
              npairs1_dv[idx_dv] = stats_sn.npairs[idx_row];
              npairs2_dv[idx_dv] = stats_sn.npairs[idx_col];
            }
          }
        }

        if (params.shape == "triu") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
              int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
              + (idx_col - idx_row);

              r1bin_dv[idx_dv] = rbinning.bin_centres[idx_row];
              r2bin_dv[idx_dv] = rbinning.bin_centres[idx_col];
              r1eff_dv[idx_dv] = stats_sn.r[idx_row];
              r2eff_dv[idx_dv] = stats_sn.r[idx_col];
              npairs1_dv[idx_dv] = stats_sn.npairs[idx_row];
              npairs2_dv[idx_dv] = stats_sn.npairs[idx_col];
            }
          }
        }
      }

      // ·································································
      // Raw 3PCF
      // ·································································

      // Compute 3PCF components in eqs. (42), (48) & (49) in the Paper.
      MeshField G_00(params, true, "`G_00`");  // G_00
      G_00.compute_unweighted_field_fluctuations_insitu(catalogue_data);
      G_00.fourier_transform();
      G_00.apply_assignment_compensation();
      G_00.inv_fourier_transform();

      MeshField F_lm_a(params, true, "`F_lm_a`");  // F_lm_a
      MeshField F_lm_b(params, true, "`F_lm_b`");  // F_lm_b

      for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
        double r_a = r1eff_dv[idx_dv];
        F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
          dn_00, ylm_k_a, sj_a, r_a
        );

        double r_b = r2eff_dv[idx_dv];
        F_lm_b.inv_fourier_transform_sjl_ylm_wgtd_field(
          dn_00, ylm_k_b, sj_b, r_b
        );

        // ζ_{l₁ l₂ L}^{m₁ m₂ M}
        double zeta_comp_real = 0., zeta_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:zeta_comp_real, zeta_comp_imag)
#endif  // TRV_USE_OMP
        for (long long gid = 0; gid < params.nmesh; gid++) {
          std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
          std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
          std::complex<double> G_00_gridpt(G_00[gid][0], G_00[gid][1]);
          std::complex<double> zeta_gridpt =
            F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;

          zeta_comp_real += zeta_gridpt.real();
          zeta_comp_imag += zeta_gridpt.imag();
        }

        std::complex<double> zeta_component(zeta_comp_real, zeta_comp_imag);

        zeta_dv[idx_dv] += parity * coupling * vol_cell * zeta_component;
      }

      count_terms++;
      if (trvs::currTask == 0) {
        trvs::logger.stat(
          "Three-point correlation function term computed at orders "
          "(m1, m2, M) = (%d, %d, 0).",
          m1_, m2_
        );
      }
    }
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::ThreePCFMeasurements threepcf_out;
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    threepcf_out.r1_bin.push_back(r1bin_dv[idx_dv]);
    threepcf_out.r1_eff.push_back(r1eff_dv[idx_dv]);
    threepcf_out.npairs_1.push_back(npairs1_dv[idx_dv]);
    threepcf_out.r2_bin.push_back(r2bin_dv[idx_dv]);
    threepcf_out.r2_eff.push_back(r2eff_dv[idx_dv]);
    threepcf_out.npairs_2.push_back(npairs2_dv[idx_dv]);
    threepcf_out.zeta_raw.push_back(norm_factor * zeta_dv[idx_dv]);
    threepcf_out.zeta_shot.push_back(norm_factor * sn_dv[idx_dv]);
  }
  threepcf_out.dim = dv_dim;

  delete[] npairs1_dv; delete[] npairs2_dv;
  delete[] r1bin_dv; delete[] r2bin_dv;
  delete[] r1eff_dv; delete[] r2eff_dv;
  delete[] zeta_dv; delete[] sn_dv;

  trvs::count_cgrid -= 4;
  trvs::count_grid -= 4;
  trvs::gbytesMem -=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed three-point correlation function "
      "from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation."
    );
  }

  return threepcf_out;
}

trv::ThreePCFWindowMeasurements compute_3pcf_window(
  ParticleCatalogue& catalogue_rand, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double alpha, double norm_factor, bool wide_angle
) {
  std::string msg_tag = wide_angle ? "wide-angle corrections " : "";

  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing three-point correlation function window %s"
      "from random catalogue...",
      msg_tag.c_str()
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int dv_dim = 0;  // data vector dimension
  if (params.shape == "diag" || params.shape == "row" ) {
    dv_dim = rbinning.num_bins;
  } else
  if (params.shape == "off-diag") {
    dv_dim = rbinning.num_bins - std::abs(params.idx_bin);
  } else
  if (params.shape == "full") {
    dv_dim = rbinning.num_bins * rbinning.num_bins;
  } else
  if (params.shape == "triu") {
    dv_dim = rbinning.num_bins * (rbinning.num_bins + 1) / 2;
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Three-point statistic form is not recognised: `form` = '%s'.",
        params.form.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Three-point statistic form is not recognised: `form` = '%s'.\n",
      params.form.c_str()
    );
  }

  int* npairs1_dv = new int[dv_dim];
  int* npairs2_dv = new int[dv_dim];
  double* r1bin_dv = new double[dv_dim];
  double* r2bin_dv = new double[dv_dim];
  double* r1eff_dv = new double[dv_dim];
  double* r2eff_dv = new double[dv_dim];
  std::complex<double>* zeta_dv = new std::complex<double>[dv_dim];
  std::complex<double>* sn_dv = new std::complex<double>[dv_dim];
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    npairs1_dv[idx_dv] = 0;
    npairs2_dv[idx_dv] = 0;
    r1bin_dv[idx_dv] = 0.;
    r2bin_dv[idx_dv] = 0.;
    r1eff_dv[idx_dv] = 0.;
    r2eff_dv[idx_dv] = 0.;
    zeta_dv[idx_dv] = 0.;
    sn_dv[idx_dv] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField n_00(params, true, "`n_00`");  // n_00(k)
  n_00.compute_ylm_wgtd_field(catalogue_rand, los_rand, alpha, 0, 0);
  n_00.fourier_transform();

  double vol_cell = n_00.vol_cell;

  MeshField N_00(params, true, "`N_00`");  // N_00(k)
  N_00.compute_ylm_wgtd_quad_field(catalogue_rand, los_rand, alpha, 0, 0);
  N_00.fourier_transform();

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  FieldStats stats_sn(params);

  // Initialise reduced-spherical-harmonic weights on mesh grids.
  std::vector< std::complex<double> > ylm_r_a(params.nmesh);
  std::vector< std::complex<double> > ylm_r_b(params.nmesh);
  std::vector< std::complex<double> > ylm_k_a(params.nmesh);
  std::vector< std::complex<double> > ylm_k_b(params.nmesh);
  trvs::count_cgrid += 4;
  trvs::count_grid += 4;
  trvs::update_maxcntgrid();
  trvs::gbytesMem +=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);
  trvs::update_maxmem();

  // Compute 3PCF window terms including shot noise.
  int count_terms = 0;
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Check for vanishing cases where all Wigner-3j symbols are zero.
      // Check for if all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );
        if (std::fabs(coupling) > trvm::eps_coupling) {
          flag_vanishing = "false";
          break;
        }
      }
      if (flag_vanishing == "true") {continue;}

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        // Calculate the coupling coefficient.
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );  // Wigner 3-j's
        if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

        // ·······························································
        // Shot noise
        // ·······························································

        // Compute shot noise components in eq. (51) in the Paper.
        MeshField n_LM_for_sn(params, true, "`n_LM_for_sn`");  // δn_LM(k)
                                                               // (for
                                                               // shot noise)
        n_LM_for_sn.compute_ylm_wgtd_field(
          catalogue_rand, los_rand, alpha, params.ELL, M_
        );
        n_LM_for_sn.fourier_transform();

        std::complex<double> Sbar_LM = calc_ylm_wgtd_shotnoise_amp_for_bispec(
          catalogue_rand, los_rand, alpha, params.ELL, M_
        );  // \bar{S}_LM

        stats_sn.compute_uncoupled_shotnoise_for_3pcf(
          n_LM_for_sn, N_00, ylm_r_a, ylm_r_b, Sbar_LM, rbinning
        );  // S|{i = j ≠ k}

        // Enforce the Kronecker delta in eq. (51) in the Paper.
        if (params.shape == "diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;
            sn_dv[idx_dv] += coupling * stats_sn.xi[ibin];
          }
        }

        if (params.shape == "off-diag" && params.idx_bin == 0) {
          // Note that ``idx_col == idx_row``.
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;
            sn_dv[idx_dv] += coupling * stats_sn.xi[ibin];
          }
        }

        if (params.shape == "row") {
          sn_dv[params.idx_bin] += coupling * stats_sn.xi[params.idx_bin];
        }

        if (params.shape == "full") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            // Note that ``idx_col == idx_row``.
            int idx_dv = idx_row * params.num_bins + idx_row;
            sn_dv[idx_dv] += coupling * stats_sn.xi[idx_row];
          }
        }

        if (params.shape == "triu") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            // Note that ``idx_col == idx_row``.
            int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2;
            sn_dv[idx_dv] += coupling * stats_sn.xi[idx_row];
          }
        }

        if (count_terms == 0) {
          if (params.shape == "diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin = idx_dv;

              r1bin_dv[idx_dv] = rbinning.bin_centres[ibin];
              r2bin_dv[idx_dv] = rbinning.bin_centres[ibin];
              r1eff_dv[idx_dv] = stats_sn.r[ibin];
              r2eff_dv[idx_dv] = stats_sn.r[ibin];
              npairs1_dv[idx_dv] = stats_sn.npairs[ibin];
              npairs2_dv[idx_dv] = stats_sn.npairs[ibin];
            }
          }

          if (params.shape == "off-diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_row, ibin_col;
              if (params.idx_bin >= 0) {
                ibin_row = idx_dv;
                ibin_col = idx_dv + std::abs(params.idx_bin);
              } else {
                ibin_row = idx_dv + std::abs(params.idx_bin);
                ibin_col = idx_dv;
              }

              r1bin_dv[idx_dv] = rbinning.bin_centres[ibin_row];
              r2bin_dv[idx_dv] = rbinning.bin_centres[ibin_col];
              r1eff_dv[idx_dv] = stats_sn.r[ibin_row];
              r2eff_dv[idx_dv] = stats_sn.r[ibin_col];
              npairs1_dv[idx_dv] = stats_sn.npairs[ibin_row];
              npairs2_dv[idx_dv] = stats_sn.npairs[ibin_col];
            }
          }

          if (params.shape == "row") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_row = params.idx_bin;
              int ibin_col = idx_dv;

              r1bin_dv[idx_dv] = rbinning.bin_centres[ibin_row];
              r2bin_dv[idx_dv] = rbinning.bin_centres[ibin_col];
              r1eff_dv[idx_dv] = stats_sn.r[ibin_row];
              r2eff_dv[idx_dv] = stats_sn.r[ibin_col];
              npairs1_dv[idx_dv] = stats_sn.npairs[ibin_row];
              npairs2_dv[idx_dv] = stats_sn.npairs[ibin_col];
            }
          }

          if (params.shape == "full") {
            for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
              for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
                int idx_dv = idx_row * params.num_bins + idx_col;

                r1bin_dv[idx_dv] = rbinning.bin_centres[idx_row];
                r2bin_dv[idx_dv] = rbinning.bin_centres[idx_col];
                r1eff_dv[idx_dv] = stats_sn.r[idx_row];
                r2eff_dv[idx_dv] = stats_sn.r[idx_col];
                npairs1_dv[idx_dv] = stats_sn.npairs[idx_row];
                npairs2_dv[idx_dv] = stats_sn.npairs[idx_col];
              }
            }
          }

          if (params.shape == "triu") {
            for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
              for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
                int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                + (idx_col - idx_row);

                r1bin_dv[idx_dv] = rbinning.bin_centres[idx_row];
                r2bin_dv[idx_dv] = rbinning.bin_centres[idx_col];
                r1eff_dv[idx_dv] = stats_sn.r[idx_row];
                r2eff_dv[idx_dv] = stats_sn.r[idx_col];
                npairs1_dv[idx_dv] = stats_sn.npairs[idx_row];
                npairs2_dv[idx_dv] = stats_sn.npairs[idx_col];
              }
            }
          }
        }

        // ·······························································
        // Raw 3PCF
        // ·······························································

        // Compute 3PCF components in eqs. (42), (48) & (49) in the Paper.
        MeshField G_LM(params, true, "`G_LM`");  // G_LM
        G_LM.compute_ylm_wgtd_field(
          catalogue_rand, los_rand, alpha, params.ELL, M_
        );
        G_LM.fourier_transform();
        G_LM.apply_assignment_compensation();
        G_LM.inv_fourier_transform();

        // Perform wide-angle corrections if required.
        if (wide_angle) {
          G_LM.apply_wide_angle_pow_law_kernel();
        }

        MeshField F_lm_a(params, true, "`F_lm_a`");  // F_lm_a
        MeshField F_lm_b(params, true, "`F_lm_b`");  // F_lm_b

        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          double r_a = r1eff_dv[idx_dv];
          F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
            n_00, ylm_k_a, sj_a, r_a
          );

          double r_b = r2eff_dv[idx_dv];
          F_lm_b.inv_fourier_transform_sjl_ylm_wgtd_field(
            n_00, ylm_k_b, sj_b, r_b
          );

          // ζ_{l₁ l₂ L}^{m₁ m₂ M}
          double zeta_comp_real = 0., zeta_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:zeta_comp_real, zeta_comp_imag)
#endif  // TRV_USE_OMP
          for (long long gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
            std::complex<double> zeta_gridpt =
              F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

            zeta_comp_real += zeta_gridpt.real();
            zeta_comp_imag += zeta_gridpt.imag();
          }

          std::complex<double> zeta_component(zeta_comp_real, zeta_comp_imag);

          zeta_dv[idx_dv] += parity * coupling * vol_cell * zeta_component;
        }

        count_terms++;
        if (trvs::currTask == 0) {
          trvs::logger.stat(
            "Three-point correlation function window term computed at orders "
            "(m1, m2, M) = (%d, %d, %d).",
            m1_, m2_, M_
          );
        }
      }
    }
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::ThreePCFWindowMeasurements threepcfwin_out;
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    threepcfwin_out.r1_bin.push_back(r1bin_dv[idx_dv]);
    threepcfwin_out.r1_eff.push_back(r1eff_dv[idx_dv]);
    threepcfwin_out.npairs_1.push_back(npairs1_dv[idx_dv]);
    threepcfwin_out.r2_bin.push_back(r2bin_dv[idx_dv]);
    threepcfwin_out.r2_eff.push_back(r2eff_dv[idx_dv]);
    threepcfwin_out.npairs_2.push_back(npairs2_dv[idx_dv]);
    threepcfwin_out.zeta_raw.push_back(norm_factor * zeta_dv[idx_dv]);
    threepcfwin_out.zeta_shot.push_back(norm_factor * sn_dv[idx_dv]);
  }
  threepcfwin_out.dim = dv_dim;

  delete[] npairs1_dv; delete[] npairs2_dv;
  delete[] r1bin_dv; delete[] r2bin_dv;
  delete[] r1eff_dv; delete[] r2eff_dv;
  delete[] zeta_dv; delete[] sn_dv;

  trvs::count_cgrid -= 4;
  trvs::count_grid -= 4;
  trvs::gbytesMem -=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed three-point correlation function window %s"
      "from random catalogue.",
      msg_tag.c_str()
    );
  }

  return threepcfwin_out;
}

#ifdef TRV_USE_LEGACY_CODE
trv::BispecMeasurements compute_bispec_for_los_choice(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  int los_choice,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing bispectrum from paired survey-type catalogues "
      "for line-of-sight choice %d...",
      los_choice
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  double alpha = catalogue_data.wstotal / catalogue_rand.wstotal;

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int dv_dim = 0;  // data vector dimension
  if (params.shape == "diag" || params.shape == "row" ) {
    dv_dim = kbinning.num_bins;
  } else
  if (params.shape == "off-diag") {
    dv_dim = kbinning.num_bins - std::abs(params.idx_bin);
  } else
  if (params.shape == "full") {
    dv_dim = kbinning.num_bins * kbinning.num_bins;
  } else
  if (params.shape == "triu") {
    dv_dim = kbinning.num_bins * (kbinning.num_bins + 1) / 2;
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Three-point statistic form is not recognised: `form` = '%s'.",
        params.form.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Three-point statistic form is not recognised: `form` = '%s'.\n",
      params.form.c_str()
    );
  }

  int* nmodes1_dv = new int[dv_dim];
  int* nmodes2_dv = new int[dv_dim];
  double* k1bin_dv = new double[dv_dim];
  double* k2bin_dv = new double[dv_dim];
  double* k1eff_dv = new double[dv_dim];
  double* k2eff_dv = new double[dv_dim];
  std::complex<double>* bk_dv = new std::complex<double>[dv_dim];
  std::complex<double>* sn_dv = new std::complex<double>[dv_dim];
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    nmodes1_dv[idx_dv] = 0;
    nmodes2_dv[idx_dv] = 0;
    k1bin_dv[idx_dv] = 0.;
    k2bin_dv[idx_dv] = 0.;
    k1eff_dv[idx_dv] = 0.;
    k2eff_dv[idx_dv] = 0.;
    bk_dv[idx_dv] = 0.;
    sn_dv[idx_dv] = 0.;
  }  // likely redundant but safe

  // ---->
  // // Set up FFTW master plans.
  // fftw_complex* array_holder = fftw_alloc_complex(params.nmesh);
  // fftw_plan fwd_master_plan = fftw_plan_dft_3d(
  //   params.ngrid[0], params.ngrid[1], params.ngrid[2],
  //   array_holder, array_holder,
  //   FFTW_FORWARD, FFTW_MEASURE
  // );
  // fftw_plan bwd_master_plan = fftw_plan_dft_3d(
  //   params.ngrid[0], params.ngrid[1], params.ngrid[2],
  //   array_holder, array_holder,
  //   FFTW_BACKWARD, FFTW_MEASURE
  // );
  // ----<

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // RFE: Not adopted until copy-assignment constructor is checked.
  // ---->
  // // Compute common field quantities.
  // MeshField dn_00(params, true, "`dn_00`");  // δn_00(r)
  // dn_00.compute_ylm_wgtd_field(
  //   catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  // );

  // double vol_cell = dn_00.vol_cell;

  // MeshField N_00(params, true, "`N_00`");  // N_00(r)
  // N_00.compute_ylm_wgtd_quad_field(
  //   catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  // );
  // ----<

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  FieldStats stats_sn(params);

  // Initialise reduced-spherical-harmonic weights on mesh grids.
  std::vector< std::complex<double> > ylm_k_a(params.nmesh);
  std::vector< std::complex<double> > ylm_k_b(params.nmesh);
  std::vector< std::complex<double> > ylm_r_a(params.nmesh);
  std::vector< std::complex<double> > ylm_r_b(params.nmesh);
  trvs::count_cgrid += 4;
  trvs::count_grid += 4;
  trvs::update_maxcntgrid();
  trvs::gbytesMem +=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);
  trvs::update_maxmem();

  // Compute bispectrum terms.
  int count_terms = 0;
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Check for if all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );
        if (std::fabs(coupling) > trvm::eps_coupling) {
          flag_vanishing = "false";
          break;
        }
      }
      if (flag_vanishing == "true") {continue;}

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        // Calculate the coupling coefficient.
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );  // Wigner 3-j's
        if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

        // ·······························································
        // Raw bispectrum
        // ·······························································

        // Compute bispectrum components in eqs. (41) & (42) in the Paper.
        MeshField dn_LM_a(params, true, "`dn_LM_a`");  // δn_LM_a
        if (los_choice == 0) {
          dn_LM_a.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_a.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM_a.fourier_transform();

        MeshField dn_LM_b(params, true, "`dn_LM_b`");  // δn_LM_b
        if (los_choice == 1) {
          dn_LM_b.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_b.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM_b.fourier_transform();

        MeshField G_LM(params, true, "`G_LM`");  // G_LM
        if (los_choice == 2) {
          G_LM.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          G_LM.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        G_LM.fourier_transform();
        G_LM.apply_assignment_compensation();
        G_LM.inv_fourier_transform();

        double vol_cell = G_LM.vol_cell;

        MeshField F_lm_a(params, true, "`F_lm_a`");  // F_lm_a
        MeshField F_lm_b(params, true, "`F_lm_b`");  // F_lm_b

        if (params.shape == "diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin = idx_dv;

            double k_lower = kbinning.bin_edges[ibin];
            double k_upper = kbinning.bin_edges[ibin + 1];

            double k_eff_a_, k_eff_b_;
            int nmodes_a_, nmodes_b_;

            F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_LM_a, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
            );
            F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_LM_b, ylm_k_b, k_lower, k_upper, k_eff_b_, nmodes_b_
            );

            if (count_terms == 0) {
              k1bin_dv[idx_dv] = kbinning.bin_centres[ibin];
              k2bin_dv[idx_dv] = kbinning.bin_centres[ibin];
              k1eff_dv[idx_dv] = k_eff_a_;
              k2eff_dv[idx_dv] = k_eff_b_;
              nmodes1_dv[idx_dv] = nmodes_a_;
              nmodes2_dv[idx_dv] = nmodes_b_;
            }

            // B_{l₁ l₂ L}^{m₁ m₂ M}
            double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
            for (long long gid = 0; gid < params.nmesh; gid++) {
              std::complex<double> F_lm_a_gridpt(
                F_lm_a[gid][0], F_lm_a[gid][1]
              );
              std::complex<double> F_lm_b_gridpt(
                F_lm_b[gid][0], F_lm_b[gid][1]
              );
              std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
              std::complex<double> bk_gridpt =
                F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

              bk_comp_real += bk_gridpt.real();
              bk_comp_imag += bk_gridpt.imag();
            }

            std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

            bk_dv[idx_dv] += coupling * vol_cell * bk_component;
          }
        }

        if (params.shape == "off-diag") {
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_row, ibin_col;
            if (params.idx_bin >= 0) {
              ibin_row = idx_dv;
              ibin_col = idx_dv + std::abs(params.idx_bin);
            } else {
              ibin_row = idx_dv + std::abs(params.idx_bin);
              ibin_col = idx_dv;
            }

            double k_lower_a = kbinning.bin_edges[ibin_row];
            double k_upper_a = kbinning.bin_edges[ibin_row + 1];
            double k_lower_b = kbinning.bin_edges[ibin_col];
            double k_upper_b = kbinning.bin_edges[ibin_col + 1];

            double k_eff_a_, k_eff_b_;
            int nmodes_a_, nmodes_b_;

            F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_LM_a, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
            );
            F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_LM_b, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
            );

            if (count_terms == 0) {
              k1bin_dv[idx_dv] = kbinning.bin_centres[ibin_row];
              k2bin_dv[idx_dv] = kbinning.bin_centres[ibin_col];
              k1eff_dv[idx_dv] = k_eff_a_;
              k2eff_dv[idx_dv] = k_eff_b_;
              nmodes1_dv[idx_dv] = nmodes_a_;
              nmodes2_dv[idx_dv] = nmodes_b_;
            }

            // B_{l₁ l₂ L}^{m₁ m₂ M}
            double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
            for (long long gid = 0; gid < params.nmesh; gid++) {
              std::complex<double> F_lm_a_gridpt(
                F_lm_a[gid][0], F_lm_a[gid][1]
              );
              std::complex<double> F_lm_b_gridpt(
                F_lm_b[gid][0], F_lm_b[gid][1]
              );
              std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
              std::complex<double> bk_gridpt =
                F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

              bk_comp_real += bk_gridpt.real();
              bk_comp_imag += bk_gridpt.imag();
            }

            std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

            bk_dv[idx_dv] += coupling * vol_cell * bk_component;
          }
        }

        if (params.shape == "row") {
          int ibin_row = params.idx_bin;

          double k_lower_a = kbinning.bin_edges[params.idx_bin];
          double k_upper_a = kbinning.bin_edges[params.idx_bin + 1];

          double k_eff_a_;
          int nmodes_a_;

          F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_LM_a, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
          );

          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            int ibin_col = idx_dv;

            double k_lower_b = kbinning.bin_edges[ibin_col];
            double k_upper_b = kbinning.bin_edges[ibin_col + 1];

            double k_eff_b_;
            int nmodes_b_;

            F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_LM_b, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
            );

            if (count_terms == 0) {
              k1bin_dv[idx_dv] = kbinning.bin_centres[ibin_row];
              k2bin_dv[idx_dv] = kbinning.bin_centres[ibin_col];
              k1eff_dv[idx_dv] = k_eff_a_;
              k2eff_dv[idx_dv] = k_eff_b_;
              nmodes1_dv[idx_dv] = nmodes_a_;
              nmodes2_dv[idx_dv] = nmodes_b_;
            }

            // B_{l₁ l₂ L}^{m₁ m₂ M}
            double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
            for (long long gid = 0; gid < params.nmesh; gid++) {
              std::complex<double> F_lm_a_gridpt(
                F_lm_a[gid][0], F_lm_a[gid][1]
              );
              std::complex<double> F_lm_b_gridpt(
                F_lm_b[gid][0], F_lm_b[gid][1]
              );
              std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
              std::complex<double> bk_gridpt =
                F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

              bk_comp_real += bk_gridpt.real();
              bk_comp_imag += bk_gridpt.imag();
            }

            std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

            bk_dv[idx_dv] += coupling * vol_cell * bk_component;
          }
        }

        if (params.shape == "full") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
              int idx_dv = idx_row * params.num_bins + idx_col;

              double k_lower_a = kbinning.bin_edges[idx_row];
              double k_upper_a = kbinning.bin_edges[idx_row + 1];
              double k_lower_b = kbinning.bin_edges[idx_col];
              double k_upper_b = kbinning.bin_edges[idx_col + 1];

              double k_eff_a_, k_eff_b_;
              int nmodes_a_, nmodes_b_;

              F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
                dn_LM_a, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
              );
              F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
                dn_LM_b, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
              );

              if (count_terms == 0) {
                k1bin_dv[idx_dv] = kbinning.bin_centres[idx_row];
                k2bin_dv[idx_dv] = kbinning.bin_centres[idx_col];
                k1eff_dv[idx_dv] = k_eff_a_;
                k2eff_dv[idx_dv] = k_eff_b_;
                nmodes1_dv[idx_dv] = nmodes_a_;
                nmodes2_dv[idx_dv] = nmodes_b_;
              }

              // B_{l₁ l₂ L}^{m₁ m₂ M}
              double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
              for (long long gid = 0; gid < params.nmesh; gid++) {
                std::complex<double> F_lm_a_gridpt(
                  F_lm_a[gid][0], F_lm_a[gid][1]
                );
                std::complex<double> F_lm_b_gridpt(
                  F_lm_b[gid][0], F_lm_b[gid][1]
                );
                std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
                std::complex<double> bk_gridpt =
                  F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

                bk_comp_real += bk_gridpt.real();
                bk_comp_imag += bk_gridpt.imag();
              }

              std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

              bk_dv[idx_dv] += coupling * vol_cell * bk_component;
            }
          }
        }

        if (params.shape == "triu") {
          for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
            for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
              int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                + (idx_col - idx_row);

              double k_lower_a = kbinning.bin_edges[idx_row];
              double k_upper_a = kbinning.bin_edges[idx_row + 1];
              double k_lower_b = kbinning.bin_edges[idx_col];
              double k_upper_b = kbinning.bin_edges[idx_col + 1];

              double k_eff_a_, k_eff_b_;
              int nmodes_a_, nmodes_b_;

              F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
                dn_LM_a, ylm_k_a, k_lower_a, k_upper_a, k_eff_a_, nmodes_a_
              );
              F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
                dn_LM_b, ylm_k_b, k_lower_b, k_upper_b, k_eff_b_, nmodes_b_
              );

              if (count_terms == 0) {
                k1bin_dv[idx_dv] = kbinning.bin_centres[idx_row];
                k2bin_dv[idx_dv] = kbinning.bin_centres[idx_col];
                k1eff_dv[idx_dv] = k_eff_a_;
                k2eff_dv[idx_dv] = k_eff_b_;
                nmodes1_dv[idx_dv] = nmodes_a_;
                nmodes2_dv[idx_dv] = nmodes_b_;
              }

              // B_{l₁ l₂ L}^{m₁ m₂ M}
              double bk_comp_real = 0., bk_comp_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:bk_comp_real, bk_comp_imag)
#endif  // TRV_USE_OMP
              for (long long gid = 0; gid < params.nmesh; gid++) {
                std::complex<double> F_lm_a_gridpt(
                  F_lm_a[gid][0], F_lm_a[gid][1]
                );
                std::complex<double> F_lm_b_gridpt(
                  F_lm_b[gid][0], F_lm_b[gid][1]
                );
                std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
                std::complex<double> bk_gridpt =
                  F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;

                bk_comp_real += bk_gridpt.real();
                bk_comp_imag += bk_gridpt.imag();
              }

              std::complex<double> bk_component(bk_comp_real, bk_comp_imag);

              bk_dv[idx_dv] += coupling * vol_cell * bk_component;
            }
          }
        }

        // ·······························································
        // Shot noise
        // ·······························································

        // Compute shot noise components in eqs. (45) & (46) in the Paper.
        MeshField dn_LM_a_for_sn(params, true, "`dn_LM_a_for_sn`");  // δn_LM_a(k)
                                                               // (for shot
                                                               // noise)
        if (los_choice == 0) {
          dn_LM_a_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_a_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM_a_for_sn.fourier_transform();

        MeshField dn_LM_b_for_sn(params, true, "`dn_LM_b_for_sn`");  // δn_LM_b(k)
                                                               // (for shot
                                                               // noise)
        if (los_choice == 1) {
          dn_LM_b_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_b_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM_b_for_sn.fourier_transform();

        // δn_LM_c(k) (for shot noise)
        MeshField dn_LM_c_for_sn(params, true, "`dn_LM_c_for_sn`");
        if (los_choice == 2) {
          dn_LM_c_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_c_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        }
        dn_LM_c_for_sn.fourier_transform();

        MeshField N_LM_a(params, true, "`N_LM_a`");  // N_LM_a(k)
        if (los_choice == 0) {
          N_LM_a.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        } else {
          N_LM_a.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        }
        N_LM_a.fourier_transform();

        MeshField N_LM_b(params, true, "`N_LM_b`");  // N_LM_b(k)
        if (los_choice == 1) {
          N_LM_b.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        } else {
          N_LM_b.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        }
        N_LM_b.fourier_transform();

        MeshField N_LM_c(params, true, "`N_LM_c`");  // N_LM_c(k)
        if (los_choice == 2) {
          N_LM_c.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            0, 0
          );
        } else {
          N_LM_c.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        }
        N_LM_c.fourier_transform();

        std::complex<double> Sbar_LM = calc_ylm_wgtd_shotnoise_amp_for_bispec(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );  // \bar{S}_LM

        if (params.ell1 == 0 && params.ell2 == 0) {
          // When l₁ = l₂ = 0, the Wigner 3-j symbol enforces L = 0
          // and the pre-factors involving degrees and orders become 1.
          std::complex<double> S_ijk = coupling * Sbar_LM;  // S|{i = j = k}
          for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
            sn_dv[idx_dv] += S_ijk;
          }
        }

        if (params.ell2 == 0) {  // S|{i ≠ j = k}
          // When l₂ = 0, the Wigner 3-j symbol enforces L = l₁.
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_LM_a_for_sn, N_LM_a, Sbar_LM, params.ell1, m1_, kbinning
          );

          if (params.shape == "diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin = idx_dv;
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin] - stats_sn.sn[ibin]
              );
            }
          }

          if (params.shape == "off-diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_row;
              if (params.idx_bin >= 0) {
                ibin_row = idx_dv;
              } else {
                ibin_row = idx_dv + std::abs(params.idx_bin);
              }
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin_row] - stats_sn.sn[ibin_row]
              );
            }
          }

          if (params.shape == "row") {
            std::complex<double> sn_row_ = coupling * (
              stats_sn.pk[params.idx_bin] - stats_sn.sn[params.idx_bin]
            );
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              sn_dv[idx_dv] += sn_row_;
            }
          }

          if (params.shape == "full") {
            for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
              std::complex<double> sn_row_ = coupling * (
                stats_sn.pk[idx_row] - stats_sn.sn[idx_row]
              );
              for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
                int idx_dv = idx_row * params.num_bins + idx_col;
                sn_dv[idx_dv] += sn_row_;
              }
            }
          }

          if (params.shape == "triu") {
            for (int idx_row = 0; idx_row < params.num_bins; idx_row++) {
              std::complex<double> sn_row_ = coupling * (
                stats_sn.pk[idx_row] - stats_sn.sn[idx_row]
              );
              for (int idx_col = idx_row; idx_col < params.num_bins; idx_col++) {
                int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                  + (idx_col - idx_row);
                sn_dv[idx_dv] += sn_row_;
              }
            }
          }
        }

        if (params.ell1 == 0) {  // S|{j ≠ i = k}
          // When l₁ = 0, the Wigner 3-j symbol enforces L = l₂.
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_LM_b_for_sn, N_LM_b, Sbar_LM,
            params.ell2, m2_, kbinning
          );

          if (params.shape == "diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin = idx_dv;
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin] - stats_sn.sn[ibin]
              );
            }
          }

          if (params.shape == "off-diag") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_col;
              if (params.idx_bin >= 0) {
                ibin_col = idx_dv + params.idx_bin;
              } else {
                ibin_col = idx_dv;
              }
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin_col] - stats_sn.sn[ibin_col]
              );
            }
          }

          if (params.shape == "row") {
            for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
              int ibin_col = idx_dv;
              sn_dv[idx_dv] += coupling * (
                stats_sn.pk[ibin_col] - stats_sn.sn[ibin_col]
              );
            }
          }

          if (params.shape == "full") {
            for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
              std::complex<double> sn_col_ = coupling * (
                stats_sn.pk[idx_col] - stats_sn.sn[idx_col]
              );
              for (int idx_row = 0; idx_row <= params.num_bins; idx_row++) {
                int idx_dv = idx_row * params.num_bins + idx_col;
                sn_dv[idx_dv] += sn_col_;
              }
            }
          }

          if (params.shape == "triu") {
            for (int idx_col = 0; idx_col < params.num_bins; idx_col++) {
              std::complex<double> sn_col_ = coupling * (
                stats_sn.pk[idx_col] - stats_sn.sn[idx_col]
              );
              for (int idx_row = 0; idx_row <= idx_col; idx_row++) {
                int idx_dv = (2*params.num_bins - idx_row + 1) * idx_row / 2
                  + (idx_col - idx_row);
                sn_dv[idx_dv] += sn_col_;
              }
            }
          }
        }

        for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
          double k_a = k1eff_dv[idx_dv];
          double k_b = k2eff_dv[idx_dv];

          std::complex<double> S_ij_k = parity *
            stats_sn.compute_uncoupled_shotnoise_for_bispec_per_bin(
              dn_LM_c_for_sn, N_LM_c, ylm_r_a, ylm_r_b, sj_a, sj_b,
              Sbar_LM, k_a, k_b
            );  // S|{i = j ≠ k}

          sn_dv[idx_dv] += coupling * S_ij_k;
        }

        count_terms++;
        if (trvs::currTask == 0) {
          trvs::logger.stat(
            "Bispectrum term computed at orders (m1, m2, M) = (%d, %d, %d).",
            m1_, m2_, M_
          );
        }
      }
    }
  }

  // ---->
  // // Clean up FFTW master plans.
  // fftw_destroy_plan(fwd_master_plan);
  // fftw_destroy_plan(bwd_master_plan);
  // fftw_free(array_holder);
  // ----<

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::BispecMeasurements bispec_out;
  for (int idx_dv = 0; idx_dv < dv_dim; idx_dv++) {
    bispec_out.k1_bin.push_back(k1bin_dv[idx_dv]);
    bispec_out.k1_eff.push_back(k1eff_dv[idx_dv]);
    bispec_out.nmodes_1.push_back(nmodes1_dv[idx_dv]);
    bispec_out.k2_bin.push_back(k2bin_dv[idx_dv]);
    bispec_out.k2_eff.push_back(k2eff_dv[idx_dv]);
    bispec_out.nmodes_2.push_back(nmodes2_dv[idx_dv]);
    bispec_out.bk_raw.push_back(norm_factor * bk_dv[idx_dv]);
    bispec_out.bk_shot.push_back(norm_factor * sn_dv[idx_dv]);
  }
  bispec_out.dim = dv_dim;

  delete[] nmodes1_dv; delete[] nmodes2_dv;
  delete[] k1bin_dv; delete[] k2bin_dv;
  delete[] k1eff_dv; delete[] k2eff_dv;
  delete[] bk_dv; delete[] sn_dv;

  trvs::count_cgrid -= 4;
  trvs::count_grid -= 4;
  trvs::gbytesMem -=
    trvs::size_in_gb< std::complex<double> >(4*params.nmesh);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed bispectrum from paired survey-type catalogues "
      "for line-of-sight choice %d.",
      los_choice
    );
  }

  return bispec_out;
}
#endif  // TRV_USE_LEGACY_CODE

}  // namespace trv
