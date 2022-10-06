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
 * @file twopt.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang)
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "twopt.hpp"

/// CAVEAT: Discretionary choice.
const double eps_norm = 1.e-5;

namespace trvm = trv::maths;

namespace trv {

/// **********************************************************************
/// Coupling coefficients
/// **********************************************************************

double calc_coupling_coeff_2pt(int ell, int ELL, int m, int M) {
  return (2*ell + 1) * (2*ELL + 1)
    * trvm::wigner_3j(ell, 0, ELL, 0, 0, 0)
    * trvm::wigner_3j(ell, 0, ELL, m, 0, M);
}


/// **********************************************************************
/// Normalisation
/// **********************************************************************

double calc_powspec_normalisation_from_mesh(
  trv::ParticleCatalogue& particles, trv::ParameterSet& params, double alpha
) {
  trv::MeshField catalogue_mesh(params);

  double norm_factor = catalogue_mesh.calc_grid_based_powlaw_norm(particles, 2);

  catalogue_mesh.finalise_density_field();  // likely redundant but safe

  norm_factor /= std::pow(alpha, 2);

  return norm_factor;
}

double calc_powspec_normalisation_from_particles(
  ParticleCatalogue& particles, double alpha
) {
  if (particles.pdata == nullptr) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle data are uninitialised.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  double norm = 0.;  // I₂
  for (int pid = 0; pid < particles.ntotal; pid++) {
    norm += particles[pid].ws
      * particles[pid].nz * std::pow(particles[pid].wc, 2);
  }

  if (norm == 0.) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle 'nz' values appear to be all zeros. "
        "Check the input catalogue contains valid 'nz' field.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  double norm_factor = 1. / (alpha * norm);  // 1/I₂

  return norm_factor;
}


/// **********************************************************************
/// Shot noise
/// **********************************************************************

double calc_powspec_shotnoise_from_particles(
  ParticleCatalogue& particles, double alpha
) {
  if (particles.pdata == nullptr) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle data are uninitialised.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  double shotnoise = 0.;
  for (int pid = 0; pid < particles.ntotal; pid++) {
    shotnoise +=
      std::pow(particles[pid].ws, 2) * std::pow(particles[pid].wc, 2);
  }

  return shotnoise;
}

std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_powspec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  double alpha, int ell, int m
) {
  std::complex<double> sn_data = 0.;
  for (int pid = 0; pid < particles_data.ntotal; pid++) {
    std::vector<double> los_{
      los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    sn_data += ylm * std::pow(particles_data[pid].w, 2);
  }

  std::complex<double> sn_rand = 0.;
  for (int pid = 0; pid < particles_rand.ntotal; pid++) {
    std::vector<double> los_{
      los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    sn_rand += ylm * std::pow(particles_rand[pid].w, 2);
  }

  return sn_data + std::pow(alpha, 2) * sn_rand;
}

std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_powspec(
  ParticleCatalogue& particles, LineOfSight* los,
  double alpha, int ell, int m
) {
  std::complex<double> sn = 0.;
  for (int pid = 0; pid < particles.ntotal; pid++) {
    std::vector<double> los_{los[pid].pos[0], los[pid].pos[1], los[pid].pos[2]};

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    sn += ylm * std::pow(particles[pid].w, 2);
  }

  return std::pow(alpha, 2) * sn;
}


/// **********************************************************************
/// Full statistics
/// **********************************************************************

/// STYLE: Standard naming convention is not always followed for
/// intermediary quantities in the functions below.

trv::PowspecMeasurements compute_powspec(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measuring power spectrum from "
      "paired survey-type catalogues...\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// --------------------------------------------------------------------
  /// Set-up
  /// --------------------------------------------------------------------

  /// Set up input.
  double alpha = catalogue_data.wtotal / catalogue_rand.wtotal;
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
  }  // likely redundant but safe

  /// --------------------------------------------------------------------
  /// Measurement
  /// --------------------------------------------------------------------

  MeshField dn_00(params);  // δn_00(k)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params);  // δn_LM(k)
    dn_LM.compute_ylm_wgtd_field(
      catalogue_data, catalogue_rand, los_data, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    std::complex<double> sn_amp = trv::calc_ylm_wgtd_shotnoise_amp_for_powspec(
      catalogue_data, catalogue_rand, los_data, los_rand, alpha, params.ELL, M_
    );  // \bar{N}_LM(k)

    /// Compute quantity equivalent to (-1)^m₁ δᴰ_{m₁, -M} which, after
    /// being summed over m₁, agrees with Hand et al. (2017) [1704.02357].
    FieldStats stats_2pt(params);
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = calc_coupling_coeff_2pt(ell1, params.ELL, m1, M_);
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      stats_2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
        dn_LM, dn_00, sn_amp, ell1, m1, kbinning
      );

      for (int ibin = 0; ibin < params.num_bins; ibin++) {
        pk_save[ibin] += coupling * stats_2pt.pk[ibin];
        sn_save[ibin] += coupling * stats_2pt.sn[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_bins; ibin++) {
          nmodes_save[ibin] = stats_2pt.nmodes[ibin];
          k_save[ibin] = stats_2pt.k[ibin];
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

  /// --------------------------------------------------------------------
  /// Results
  /// --------------------------------------------------------------------

  trv::PowspecMeasurements powspec_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    powspec_out.kbin.push_back(kbinning.bin_centres[ibin]);
    powspec_out.keff.push_back(k_save[ibin]);
    powspec_out.nmodes.push_back(nmodes_save[ibin]);
    powspec_out.pk_raw.push_back(norm_factor * pk_save[ibin]);
    powspec_out.pk_shot.push_back(norm_factor * sn_save[ibin]);
  }

  delete[] nmodes_save; delete[] k_save; delete[] pk_save; delete[] sn_save;

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] ... measured power spectrum "
      "from paired survey-type catalogues.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  return powspec_out;
}

trv::TwoPCFMeasurements compute_corrfunc(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measuring two-point correlation function "
      "from paired survey-type catalogues...\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// --------------------------------------------------------------------
  /// Set-up
  /// --------------------------------------------------------------------

  /// Set up input.
  double alpha = catalogue_data.wtotal / catalogue_rand.wtotal;
  int ell1 = params.ELL;

  /// Set up output.
  int* npairs_save = new int[params.num_bins];
  double* r_save = new double[params.num_bins];
  std::complex<double>* xi_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }  // likely redundant but safe

  /// --------------------------------------------------------------------
  /// Measurement
  /// --------------------------------------------------------------------

  MeshField dn_00(params);  // δn_00(k)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params);  // δn_LM(k)
    dn_LM.compute_ylm_wgtd_field(
      catalogue_data, catalogue_rand, los_data, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    std::complex<double> sn_amp = trv::calc_ylm_wgtd_shotnoise_amp_for_powspec(
      catalogue_data, catalogue_rand, los_data, los_rand, alpha, params.ELL, M_
    );  // \bar{N}_LM(k)

    /// Compute quantity equivalent to (-1)^m₁ δᴰ_{m₁, -M} which, after
    /// being summed over m₁, agrees with Hand et al. (2017) [1704.02357].
    FieldStats stats_2pt(params);
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = calc_coupling_coeff_2pt(ell1, params.ELL, m1, M_);
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      stats_2pt.compute_ylm_wgtd_2pt_stats_in_config(
        dn_LM, dn_00, sn_amp, ell1, m1, rbinning
      );

      for (int ibin = 0; ibin < params.num_bins; ibin++) {
        xi_save[ibin] += coupling * stats_2pt.xi[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_bins; ibin++) {
          npairs_save[ibin] = stats_2pt.npairs[ibin];
          r_save[ibin] = stats_2pt.r[ibin];
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

  /// --------------------------------------------------------------------
  /// Results
  /// --------------------------------------------------------------------

  trv::TwoPCFMeasurements corrfunc_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    corrfunc_out.rbin.push_back(rbinning.bin_centres[ibin]);
    corrfunc_out.reff.push_back(r_save[ibin]);
    corrfunc_out.npairs.push_back(npairs_save[ibin]);
    corrfunc_out.xi.push_back(norm_factor * xi_save[ibin]);
  }

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] ... measured two-point correlation function "
      "from paired survey-type catalogues.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  return corrfunc_out;
}

trv::PowspecMeasurements compute_powspec_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning kbinning,
  double norm_factor
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measuring power spectrum from "
      "a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// --------------------------------------------------------------------
  /// Set-up
  /// --------------------------------------------------------------------

  int* nmodes_save = new int[params.num_bins];
  double* k_save = new double[params.num_bins];
  std::complex<double>* pk_save = new std::complex<double>[params.num_bins];
  std::complex<double>* sn_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  /// Check input normalisation matches expectation.
  double norm = catalogue_data.ntotal * catalogue_data.ntotal / params.volume;
  if (std::fabs(1 - norm * norm_factor) < eps_norm) {
    std::printf(
      "[%s WARN] Power spectrum normalisation input differs from "
      "expected value for an unweight field in a periodic box.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// --------------------------------------------------------------------
  /// Measurement
  /// --------------------------------------------------------------------

  /// Compute power spectrum.
  MeshField dn(params);  // δn(k)
  dn.compute_unweighted_field_fluctuations_insitu(catalogue_data);
  dn.fourier_transform();

  std::complex<double> sn_amp = double(catalogue_data.ntotal);  // \bar{N}

  /// Under the global plane-parallel approximation, δᴰ_{M0} enforces
  /// M = 0 for any spherical-harmonic-weighted field fluctuations.
  FieldStats stats_2pt(params);
  stats_2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
    dn, dn, sn_amp, params.ELL, 0, kbinning
  );

  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    nmodes_save[ibin] = stats_2pt.nmodes[ibin];
    k_save[ibin] = stats_2pt.k[ibin];
    pk_save[ibin] += double(2*params.ELL + 1) * stats_2pt.pk[ibin];
    sn_save[ibin] += double(2*params.ELL + 1) * stats_2pt.sn[ibin];
  }

  /// --------------------------------------------------------------------
  /// Results
  /// --------------------------------------------------------------------

  /// Fill in output struct.
  trv::PowspecMeasurements powspec_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    powspec_out.kbin.push_back(kbinning.bin_centres[ibin]);
    powspec_out.keff.push_back(k_save[ibin]);
    powspec_out.nmodes.push_back(nmodes_save[ibin]);
    powspec_out.pk_raw.push_back(norm_factor * pk_save[ibin]);
    powspec_out.pk_shot.push_back(norm_factor * sn_save[ibin]);
  }

  delete[] nmodes_save; delete[] k_save; delete[] pk_save; delete[] sn_save;

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] ... measured power spectrum from "
      "a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  return powspec_out;
}

trv::TwoPCFMeasurements compute_corrfunc_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measuring two-point correlation function "
      "a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// --------------------------------------------------------------------
  /// Set-up
  /// --------------------------------------------------------------------

  int* npairs_save = new int[params.num_bins];
  double* r_save = new double[params.num_bins];
  std::complex<double>* xi_save = new std::complex<double>[params.num_bins];
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }  // likely redundant but safe

  /// Check input normalisation matches expectation.
  double norm = catalogue_data.ntotal * catalogue_data.ntotal / params.volume;
  if (std::fabs(1 - norm * norm_factor) < eps_norm) {
    std::printf(
      "[%s WARN] Two-point correlation function normalisation input differs "
      "from expected value for an unweight field in a periodic box.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// --------------------------------------------------------------------
  /// Measurement
  /// --------------------------------------------------------------------

  /// Compute 2PCF.
  MeshField dn(params);  // δn(k)
  dn.compute_unweighted_field_fluctuations_insitu(catalogue_data);
  dn.fourier_transform();

  std::complex<double> sn_amp = double(catalogue_data.ntotal);  // \bar{N}

  /// Under the global plane-parallel approximation, δᴰ_{M0} enforces
  /// M = 0 for any spherical-harmonic-weighted field fluctuations.
  FieldStats stats_2pt(params);
  stats_2pt.compute_ylm_wgtd_2pt_stats_in_config(
    dn, dn, sn_amp, params.ELL, 0, rbinning
  );

  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    npairs_save[ibin] = stats_2pt.npairs[ibin];
    r_save[ibin] = stats_2pt.r[ibin];
    xi_save[ibin] += double(2*params.ELL + 1) * stats_2pt.xi[ibin];
  }

  /// --------------------------------------------------------------------
  /// Results
  /// --------------------------------------------------------------------

  /// Fill in output struct.
  trv::TwoPCFMeasurements corrfunc_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    corrfunc_out.rbin.push_back(rbinning.bin_centres[ibin]);
    corrfunc_out.reff.push_back(r_save[ibin]);
    corrfunc_out.npairs.push_back(npairs_save[ibin]);
    corrfunc_out.xi.push_back(norm_factor * xi_save[ibin]);
  }

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] ... measured two-point correlation function "
      "a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  return corrfunc_out;
}

trv::TwoPCFWindowMeasurements compute_corrfunc_window(
  trv::ParticleCatalogue& catalogue_rand, trv::LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning rbinning,
  double alpha, double norm_factor
) {
  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Measuring two-point correlation function window "
      "from a random catalogue...\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  /// --------------------------------------------------------------------
  /// Set-up
  /// --------------------------------------------------------------------

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
  }  // likely redundant but safe

  /// --------------------------------------------------------------------
  /// Measurement
  /// --------------------------------------------------------------------

  MeshField dn_00(params);
  dn_00.compute_ylm_wgtd_field(catalogue_rand, los_rand, alpha, 0, 0);
  dn_00.fourier_transform();  // δn_00(k)

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params);
    dn_LM.compute_ylm_wgtd_field(
      catalogue_rand, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();  // δn_LM(k)

    FieldStats stats_2pt(params);
    std::complex<double> sn_amp = trv::calc_ylm_wgtd_shotnoise_amp_for_powspec(
      catalogue_rand, los_rand, alpha, params.ELL, M_
    );  // \bar{N}_LM(k)

    /// Calculate equivalent to (-1)^m₁ δᴰ_{m₁, -M} which, after being
    /// summed over m₁, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = trv::calc_coupling_coeff_2pt(ell1, params.ELL, m1, M_);
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      stats_2pt.compute_ylm_wgtd_2pt_stats_in_config(
        dn_LM, dn_00, sn_amp, ell1, m1, rbinning
      );

      for (int ibin = 0; ibin < params.num_bins; ibin++) {
        xi_save[ibin] += coupling * stats_2pt.xi[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < params.num_bins; ibin++) {
          npairs_save[ibin] = stats_2pt.npairs[ibin];
          r_save[ibin] = stats_2pt.r[ibin];
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

  /// --------------------------------------------------------------------
  /// Results
  /// --------------------------------------------------------------------

  trv::TwoPCFWindowMeasurements corrfunc_win_out;
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    corrfunc_win_out.rbin.push_back(rbinning.bin_centres[ibin]);
    corrfunc_win_out.reff.push_back(r_save[ibin]);
    corrfunc_win_out.npairs.push_back(npairs_save[ibin]);
    corrfunc_win_out.xi.push_back(norm_factor * xi_save[ibin]);
  }

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] ... measured two-point correlation function window "
      "from a random catalogue.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  return corrfunc_win_out;
}

}  // namespace trv
