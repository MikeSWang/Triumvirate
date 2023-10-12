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
 * @file twopt.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Two-point statistic computations.
 *
 * This module provides the computation of two-point statistics including:
 * - power spectrum normalisation (mesh- or particle-based);
 * - power spectrum shot noise level and
 *   reduced-spherical-harmonic-weighted shot noise;
 * - power spectrum and two-point correlation function for paired
 *   survey-type catalogues;
 * - power spectrum and two-point correlation function for periodic-box
 *   simulation catalogues.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>

#include "monitor.hpp"
#include "maths.hpp"
#include "parameters.hpp"
#include "dataobjs.hpp"
#include "particles.hpp"
#include "field.hpp"

namespace trv {

// ***********************************************************************
// Coupling coefficients
// ***********************************************************************

/**
 * @brief Calculate the coupling coefficient for spherical-harmonic
 *        components of full two-point statistics.
 *
 * The calculated quantity is
 * @f[
 *   (2\ell + 1) (2L + 1) W(\ell, 0, L; m, 0, M) \,.
 * @f]
 *
 * @param ell, ELL, m, M Wigner 3-j components.
 * @returns Coupling coefficient for spherical-harmonic components in
 *         two-point statistics.
 */
double calc_coupling_coeff_2pt(int ell, int ELL, int m, int M);


// ***********************************************************************
// Normalisation
// ***********************************************************************

/**
 * @brief Calculate particle-based power spectrum normalisation.
 *
 * @param particles Particle catalogue.
 * @param alpha Alpha contrast (default is 1.).
 * @returns Power spectrum normalisation factor.
 */
double calc_powspec_normalisation_from_particles(
  ParticleCatalogue& particles, double alpha = 1.
);

/**
 * @brief Calculate mesh-based power spectrum normalisation.
 *
 * @param particles Particle catalogue.
 * @param params Parameter set.
 * @param alpha Alpha contrast (default is 1.).
 * @returns Power spectrum normalisation factor.
 */
double calc_powspec_normalisation_from_mesh(
  trv::ParticleCatalogue& particles,
  trv::ParameterSet& params, double alpha = 1.
);

/**
 * @brief Calculate power spectrum normalisation from mixed meshes.
 *
 * This normalisation convention is used for paired survey-type catalogues
 * only, where the two meshes come from separate data- and random-source
 * particle catalogues.
 *
 * @see <a href="https://pypower.readthedocs.io/en/latest/api/api.html
 * #pypower.fft_power.normalization">@c pypower documentation</a>.
 *
 * @param particles_data Data-source particle catalogue.
 * @param particles_rand Random-source particle catalogue.
 * @param params Parameter set.
 * @param alpha Alpha contrast.
 * @returns Power spectrum normalisation factor.
 */
double calc_powspec_normalisation_from_meshes(
  trv::ParticleCatalogue& particles_data,
  trv::ParticleCatalogue& particles_rand,
  trv::ParameterSet& params, double alpha
);

/**
 * @brief Calculate power spectrum normalisation from mixed meshes.
 *
 * The box size per dimension is internally determined by the particle
 * coordinate extents, with a padding factor applied. The mesh number in
 * each dimension is determined by the box size and the cell size.
 *
 * @param particles_data Data-source particle catalogue.
 * @param particles_rand Random-source particle catalogue.
 * @param params Parameter set.
 * @param alpha Alpha contrast.
 * @param padding Padding factor for the box size in each dimension.
 *                This overrides the box size set in the parameter set.
 * @param cellsize Cell size of the mesh grid in each dimension.  This
 *                 overrides the mesh number set in the parameter set.
 * @param assignment Assignment scheme for particles to mesh.  This
 *                   overrides the assignment scheme set in the
 *                   parameter set.
 * @returns Power spectrum normalisation factor.
 *
 * @overload
 */
double calc_powspec_normalisation_from_meshes(
  trv::ParticleCatalogue& particles_data,
  trv::ParticleCatalogue& particles_rand,
  trv::ParameterSet& params, double alpha,
  double padding, double cellsize, std::string assignment
);


// ***********************************************************************
// Shot noise
// ***********************************************************************

/**
 * @brief Calculate particle-based power spectrum shot noise level.
 *
 * @param particles Particle catalogue.
 * @param alpha Alpha contrast (default is 1.).
 * @returns Power spectrum normalisation factor.
 */
double calc_powspec_shotnoise_from_particles(
  ParticleCatalogue& particles, double alpha = 1.
);

/**
 * @brief Calculate power spectrum shot noise weighted by
 *        reduced spherical harmonics.
 *
 * This calculates the quantity
 * @f[
 *   \bar{N}_{LM}(\vec{x}) = {\sum_i}{\vphantom{\sum}}'
 *     y_{LM}^*(\hat{\vec{x}}) w(\vec{x})^2 \,,
 * @f]
 * where if a pair of catalogues are provided,
 * @f[
 *   {\sum_i}{\vphantom{\sum}}' =
 *     \sum_{i \in \mathrm{data}} + \alpha^2 \sum_{i \in \mathrm{rand}} \,,
 * @f]
 * and otherwise
 * @f[
 *   {\sum_i}{\vphantom{\sum}}' =
 *     \alpha^2 \sum_{i \in \mathrm{data\ or\ rand}} \,.
 * @f]
 *
 * @param particles_data (Data-source) particle catalogue.
 * @param particles_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param alpha Alpha contrast.
 * @param ell Degree of the spherical harmonic.
 * @param m Order of the spherical harmonic.
 * @returns Weighted shot noise for power spectrum.
 */
std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_powspec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  double alpha, int ell, int m
);

/**
 * @brief Calculate power spectrum shot noise weighted by
 *        reduced spherical harmonics.
 *
 * @param particles Particle catalogue.
 * @param los Particle lines of sight.
 * @param alpha Alpha contrast.
 * @param ell Degree of the spherical harmonic.
 * @param m Order of the spherical harmonic.
 * @returns Weighted power spectrum shot noise.
 *
 * @overload
 */
std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_powspec(
  ParticleCatalogue& particles, LineOfSight* los,
  double alpha, int ell, int m
);


// ***********************************************************************
// Full statistics
// ***********************************************************************

// STYLE: Standard naming convention is not always followed for
// intermediary quantities in the functions below.

/**
 * @brief Compute power spectrum from paired survey-type catalogues.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param catalogue_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param kbinning Wavenumber binning.
 * @param norm_factor Normalisation factor.
 * @returns Power spectrum measurements.
 */
trv::PowspecMeasurements compute_powspec(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
);

/**
 * @brief Compute two-point correlation function from paired
 *        survey-type catalogues.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param catalogue_rand (Random-source) particle catalogue.
 * @param los_data (Data-source) particle lines of sight.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbinning Separation binning.
 * @param norm_factor Normalisation factor.
 * @returns Two-point correlation function measurements.
 */
trv::TwoPCFMeasurements compute_corrfunc(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
);

/**
 * @brief Compute power spectrum in a periodic box in the global
 *        plane-parallel approximation.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param params Parameter set.
 * @param kbinning Wavenumber binning.
 * @param norm_factor Normalisation factor.
 * @returns Power spectrum measurements.
 */
trv::PowspecMeasurements compute_powspec_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning kbinning,
  double norm_factor
);

/**
 * @brief Compute two-point correlation function in a periodic box
 *        in the global plane-parallel approximation.
 *
 * @param catalogue_data (Data-source) particle catalogue.
 * @param params Parameter set.
 * @param rbinning Separation binning.
 * @param norm_factor Normalisation factor.
 * @returns Two-point correlation function measurements.
 */
trv::TwoPCFMeasurements compute_corrfunc_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
);

/**
 * @brief Compute two-point correlation function window from a random
 *        catalogue and optionally save the results.
 *
 * @param catalogue_rand (Random-source) particle catalogue.
 * @param los_rand (Random-source) particle lines of sight.
 * @param params Parameter set.
 * @param rbinning Separation binning.
 * @param alpha Alpha contrast.
 * @param norm_factor Normalisation factor.
 * @returns Two-point correlation function window measurements.
 */
trv::TwoPCFWindowMeasurements compute_corrfunc_window(
  trv::ParticleCatalogue& catalogue_rand, trv::LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning rbinning,
  double alpha, double norm_factor
);

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_TWOPT_HPP_INCLUDED_
