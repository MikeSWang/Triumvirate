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
 * @file threept.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Three-point statistic computations.
 *
 * This module provides the computation of three-point
 * statistics including:
 * - bispectrum normalisation (mesh- or particle-based);
 * - reduced-spherical-harmonic-weighted bispectrum shot noise;
 * - bispectrum and three-point correlation function for paired
 *   survey-type catalogues;
 * - bispectrum and three-point correlation function for periodic-box
 *   simulation-type catalogues in the global plane-parallel approximation.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_

#include <fftw3.h>

#include <cmath>
#include <complex>
#include <cstdio>

#include "monitor.hpp"
#include "parameters.hpp"
#include "maths.hpp"
#include "particles.hpp"
#include "dataobjs.hpp"
#include "field.hpp"
#include "twopt.hpp"

namespace trv {

// ***********************************************************************
// Coupling coefficients
// ***********************************************************************

/**
 * @brief Calculate the coupling coefficient for spherical-harmonic
 *        components of full three-point statistics.
 *
 * The calculated quantity is
 * @f[
 *   (2\ell_1 + 1) (2\ell_2 + 1) (2L + 1)
 *   W(\ell_1, \ell_2, L; 0, 0, 0) W(\ell_1, \ell_2, L; m, 0, M) \,.
 * @f]
 *
 * @param ell1, ell2, ELL, m1, m2, M Wigner 3-j components.
 * @returns Coupling coefficient for spherical-harmonic components in
 *         three-point statistics.
 */
double calc_coupling_coeff_3pt(
  int ell1, int ell2, int ELL, int m1, int m2, int M
);

/**
 * @brief Validate three-point correlator multipoles are non-vanishing.
 *
 * @param params Parameter set.
 */
void validate_multipole_coupling(trv::ParameterSet& params);


// ***********************************************************************
// Normalisation
// ***********************************************************************

/**
 * @brief Calculate particle-based bispectrum normalisation.
 *
 * @param particles Particle catalogue.
 * @param alpha Alpha contrast (default is 1.).
 * @returns Bispectrum normalisation factor.
 */
double calc_bispec_normalisation_from_particles(
  ParticleCatalogue& particles, double alpha = 1.
);

/**
 * @brief Calculate mesh-based bispectrum normalisation.
 *
 * @param particles Particle catalogue.
 * @param params Parameter set.
 * @param alpha Alpha contrast (default is 1.).
 * @returns Bispectrum normalisation factor.
 */
double calc_bispec_normalisation_from_mesh(
  ParticleCatalogue& particles, trv::ParameterSet& params, double alpha = 1.
);


// ***********************************************************************
// Shot noise
// ***********************************************************************

/**
 * @brief Calculate bispectrum shot noise amplitude weighted by
 *        reduced spherical harmonics.
 *
 * This calculates the quantity
 * @f[
 *   \bar{S}_{LM} = {\sum_i}{\vphantom{\sum}}'
 *     y_{LM}^*(\vec{x}_i) w(\vec{x}_i)^3 \,,
 * @f]
 * where if a pair of catalogues are provided,
 * @f[
 *   {\sum_i}{\vphantom{\sum}}' =
 *     \sum_{i \in \mathrm{data}} - \alpha^3 \sum_{i \in \mathrm{rand}} \,,
 * @f]
 * and otherwise
 * @f[
 *   {\sum_i}{\vphantom{\sum}}' =
 *     \alpha^3  \sum_{i \in \mathrm{data\ or\ rand}} \,.
 * @f]
 *
 * @see Eq. (46) in Sugiyama et al. (2019)
 *      [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
 *
 * @param particles_data (Data-source) particle catalogue.
 * @param particles_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param alpha Alpha contrast.
 * @param ell Degree of the spherical harmonic.
 * @param m Order of the spherical harmonic.
 * @returns Weighted shot-noise contribution for bispectrum.
 */
std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_bispec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  double alpha, int ell, int m
);

/**
 * @brief Calculate bispectrum shot noise amplitude weighted by
 *        reduced spherical harmonics.
 *
 * @param particles Particle catalogue.
 * @param los Particle lines of sight.
 * @param alpha Alpha contrast.
 * @param ell Degree of the spherical harmonic.
 * @param m Order of the spherical harmonic.
 * @returns Weighted shot-noise amplitude for bispectrum.
 *
 * @overload
 */
std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_bispec(
  ParticleCatalogue& particles, LineOfSight* los,
  double alpha, int ell, int m
);


// ***********************************************************************
// Full statistics
// ***********************************************************************

// STYLE: Standard naming convention is not always followed for
// intermediary quantities in the functions below.

// Hereafter 'the Paper' refers to Sugiyama et al. (2019) [1803.02132].

/**
 * @brief Compute bispectrum from paired survey-type catalogues.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param catalogue_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param kbinning Wavenumber binning.
 * @param norm_factor Normalisation factor.
 * @returns Bispectrum measurements.
 */
trv::BispecMeasurements compute_bispec(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
);

/**
 * @brief Compute three-point correlation function from paired
 *        survey-type catalogues.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param catalogue_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbinning Separation binning.
 * @param norm_factor Normalisation factor.
 * @returns Three-point correlation function measurements.
 */
trv::ThreePCFMeasurements compute_3pcf(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
);

/**
 * @brief Compute bispectrum in a periodic box in the global
 *        plane-parallel approximation.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param params Parameter set.
 * @param kbinning Wavenumber binning.
 * @param norm_factor Normalisation factor.
 * @returns Bispectrum measurements.
 */
trv::BispecMeasurements compute_bispec_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning kbinning,
  double norm_factor
);

/**
 * @brief Compute three-point correlation function in a periodic box
 *        in the global plane-parallel approximation.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param params Parameter set.
 * @param rbinning Separation binning.
 * @param norm_factor Normalisation factor.
 * @returns Three-point correlation function measurements.
 */
trv::ThreePCFMeasurements compute_3pcf_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
);

/**
 * @brief Compute three-point correlation function window from
 *        a random catalogue.
 *
 * @param catalogue_rand (Random-source) particle catalogue.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbinning Separation binning.
 * @param alpha Alpha contrast.
 * @param norm_factor Normalisation factor.
 * @param wide_angle Whether wide-angle corretions or not
 *                   (default is `false`).
 * @returns Three-point correlation function window measurements.
 */
trv::ThreePCFWindowMeasurements compute_3pcf_window(
  ParticleCatalogue& catalogue_rand, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double alpha, double norm_factor, bool wide_angle = false
);

#ifdef TRV_USE_LEGACY_CODE
/**
 * @brief Compute bispectrum from paired survey-type catalogues with
 *        a particular choice of line of sight.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param catalogue_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param los_choice Choice of line of sight in {0, 1, 2}.
 * @param params Parameter set.
 * @param kbin Wavenumber binning.
 * @param alpha Alpha contrast.
 * @param norm Normalisation factor.
 * @returns Bispectrum measurements.
 */
trv::BispecMeasurements compute_bispec_for_los_choice(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  int los_choice,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
);
#endif  // TRV_USE_LEGACY_CODE

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_THREEPT_HPP_INCLUDED_
