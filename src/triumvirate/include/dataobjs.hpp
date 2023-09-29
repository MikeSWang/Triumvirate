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
 * @file dataobjs.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Clustering measurement data objects.
 *
 * Clustering measurement data objects provided include:
 * - binning schemes;
 * - line of sight; and
 * - clustering statistics.
 */

#ifndef TRIUMVIRATE_INCLUDE_DATAOBJS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_DATAOBJS_HPP_INCLUDED_

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "monitor.hpp"
#include "parameters.hpp"

namespace trv {

// ***********************************************************************
// Binning schemes
// ***********************************************************************

/**
 * @brief Isotropic coordinate binning.
 *
 * This sets up isotropic wavenumber or separation bins
 * in configuration or Fourier space.
 *
 */
class Binning {
 public:
  std::string space;                ///< coordinate space
  std::string scheme;               ///< binning scheme
  double bin_min;                   ///< lowest bin edge
  double bin_max;                   ///< highest bin edge
  int num_bins;                     ///< number of bins
  std::vector<double> bin_edges;    ///< bin edges
  std::vector<double> bin_centres;  ///< bin centres
  std::vector<double> bin_widths;   ///< bin widths

  /**
   * @brief Construct binnng from bin specification.
   *
   * @param space Coordinate space, one of {"fourier", "config"}.
   * @param scheme Binning scheme, one of
   *               {"lin", "log", "linpad", "logpad"}.
   */
  Binning(std::string space, std::string scheme);

  /**
   * @brief Construct binning from a parameter set.
   *
   * @param params Paramater set.
   *
   * @overload
   */
  Binning(trv::ParameterSet& params);

  /**
   * @brief Set bins.
   *
   * @param coord_min Minimum coordinate in bin range.
   * @param coord_max Maximum coordinate in bin range.
   * @param nbin Number of bins.
   * @throws trv::sys::InvalidParameterError When @p coord_min is
   *                                         negative.
   * @throws trv::sys::InvalidParameterError When @p nbin is non-positive.
   *
   * @note If @p scheme is "lin" or "log", the bin edges are set linearly
   *       or log-linearly (i.e. exponentially) in the bin range.  If
   *       @p scheme is "linpad", "logpad", 5 linear bins are set from
   *       zero (leftmost edge) with width 1.e-3 ("fourier") or 10.
   *       ("config") with the remaining bin range divided linearly
   *       or log-linearly.
   */
  void set_bins(double coord_min, double coord_max, int nbin);

  /**
   * @brief Set bins.
   *
   * The bin properties are inferred from @ref trv::ParameterSet
   * when initialised with @ref trv::Binning::Binning(trv::ParameterSet&).
   *
   * @overload
   */
  void set_bins();

  /**
   * @brief Construct binning from a mesh grid.
   *
   * The bin width is given by the grid resolution in configuration
   * space or the fundamental wavenumber in Fourier space.  The bin
   * minimum is zero and the bin maximum is half the boxsize in
   * configuration space or the Nyquist wavenumber in Fourier space.
   *
   * @param boxsize_max (Maximum) box size.
   * @param ngrid_min (Minimum) grid number.
   *
   * @attention @ref trv::Binning.scheme is reset to "lin".
   *
   * @overload
   */
  void set_bins(double boxsize_max, int ngrid_min);

  /**
   * @brief Set bins from custom bin edges.
   *
   * @param bin_edges Custom bin edges (in ascending order).
   *
   * @overload
   */
  void set_bins(std::vector<double> bin_edges);

 private:
  // CAVEAT: Discretionary choices.
  int nbin_pad = 5;                 ///< number of padded bins
  double dbin_pad_fourier = 1.e-3;  ///< bin padding in Fourier space
  double dbin_pad_config = 10.;     ///< bin padding in configuration space

  /**
   * @brief Compute bin edges, centres and widths.
   *
   * @throws trv::sys::InvalidParameterError When @p scheme is not one of
   *                                         the allowed options.
   */
  void compute_binning();
};


// ***********************************************************************
// Mesh grids
// ***********************************************************************

/**
 * @brief Binned vectors.
 *
 */
struct BinnedVectors {
  int count = 0;                    ///< number of vectors
  int num_bins = 0;                 ///< number of bins
  std::vector<int> indices;         ///< bin indices
  std::vector<double> lower_edges;  ///< lower bin edges
  std::vector<double> upper_edges;  ///< upper bin edges
  std::vector<double> vecx;         ///< x-components
  std::vector<double> vecy;         ///< y-components
  std::vector<double> vecz;         ///< z-components
};


// ***********************************************************************
// Line of sight
// ***********************************************************************

/**
 * @brief Line-of-sight vector.
 *
 */
struct LineOfSight {
  double pos[3];  ///< 3-d position vector
};

// ***********************************************************************
// Clustering statistics
// ***********************************************************************

// -----------------------------------------------------------------------
// Two-point statistics
// -----------------------------------------------------------------------

/**
 * @brief Power spectrum measurements.
 *
 */
struct PowspecMeasurements {
  int dim = 0;               ///< dimension of data vector
  std::vector<double> kbin;  ///< central wavenumber in bins
  std::vector<double> keff;  ///< effective wavenumber in bins
  std::vector<int> nmodes;   ///< number of wavevectors in bins
  /// power spectrum raw measurements (with normalisation and shot noise)
  std::vector< std::complex<double> > pk_raw;
  /// power spectrum shot noise
  std::vector< std::complex<double> > pk_shot;
};

/**
 * @brief Two-point correlation function measurements.
 *
 */
struct TwoPCFMeasurements {
  int dim = 0;               ///< dimension of data vector
  std::vector<double> rbin;  ///< central separation in bins
  std::vector<double> reff;  ///< effective separation in bins
  std::vector<int> npairs;   ///< number of separation vectors in bins
  /// two-point correlation function measurements (with normalisation)
  std::vector< std::complex<double> > xi;
};

/**
 * @brief Two-point correlation function window measurements.
 *
 */
struct TwoPCFWindowMeasurements {
  int dim = 0;               ///< dimension of data vector
  std::vector<double> rbin;  ///< central separation in bins
  std::vector<double> reff;  ///< effective separation in bins
  std::vector<int> npairs;   ///< number of separation vectors in bins
  /// two-point correlation function window measurements
  /// (with normalisation)
  std::vector< std::complex<double> > xi;
};


// -----------------------------------------------------------------------
// Three-point statistics
// -----------------------------------------------------------------------

/**
 * @brief Bispectrum measurements.
 *
 */
struct BispecMeasurements {
  int dim = 0;                 ///< dimension of data vector
  std::vector<double> k1_bin;  ///< first central wavenumber in bins
  std::vector<double> k2_bin;  ///< second central wavenumber in bins
  std::vector<double> k1_eff;  ///< first effective wavenumber in bins
  std::vector<double> k2_eff;  ///< second effective wavenumber in bins
  std::vector<int> nmodes_1;   ///< number of first wavevectors in bins
  std::vector<int> nmodes_2;   ///< number of second wavevectors in bins
  /// bispectrum raw measurements (with normalisation and shot noise)
  std::vector< std::complex<double> > bk_raw;
  /// bispectrum shot noise
  std::vector< std::complex<double> > bk_shot;
};

/**
 * @brief Three-point correlation function measurements.
 *
 */
struct ThreePCFMeasurements {
  int dim = 0;                 ///< dimension of data vector
  std::vector<double> r1_bin;  ///< first central separation in bins
  std::vector<double> r2_bin;  ///< second central separation in bins
  std::vector<double> r1_eff;  ///< first effective separation in bins
  std::vector<double> r2_eff;  ///< second effective separation in bins
  /// number of first separation vectors in bins
  std::vector<int> npairs_1;
  /// number of second separation vectors in bins
  std::vector<int> npairs_2;
  /// three-point correlation function raw measurements
  /// (with normalisation and shot noise)
  std::vector< std::complex<double> > zeta_raw;
  /// three-point correlation function shot noise
  std::vector< std::complex<double> > zeta_shot;
};

/**
 * @brief Three-point correlation function window measurements.
 *
 */
struct ThreePCFWindowMeasurements {
  int dim = 0;                 ///< dimension of data vector
  std::vector<double> r1_bin;  ///< first central separation in bins
  std::vector<double> r2_bin;  ///< second central separation in bins
  std::vector<double> r1_eff;  ///< first effective separation in bins
  std::vector<double> r2_eff;  ///< second effective separation in bins
  /// number of first separation vectors in bins
  std::vector<int> npairs_1;
  /// number of second separation vectors in bins
  std::vector<int> npairs_2;
  /// three-point correlation function window raw measurements
  /// (with normalisation and shot noise)
  std::vector< std::complex<double> > zeta_raw;
  /// three-point correlation function window shot noise
  std::vector< std::complex<double> > zeta_shot;
};

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_DATAOBJS_HPP_INCLUDED_
