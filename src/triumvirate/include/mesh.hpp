// Triumvirate: Three-Point Clustering Measurements in LSS
//
// Copyright (C) 2025 Mike S Wang & Naonori S Sugiyama [GPL-3.0-or-later]
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
 * @file mesh.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang)
 * @brief Mesh grid (with one-point statistics) and
 *        pseudo two-point statistics.
 *
 * This module performs operations, transforms and calculations on a mesh
 * grid of discretely sampled fields, including the assignment of a
 * catalogue of particles with weights or the loading of external field
 * data onto the said mesh grid.  It provides the methods needed to
 * compute various constituent terms (one-point and pseudo two-point
 * statistics) in the estimators of two- and three-point statistics.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_MESH_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_MESH_HPP_INCLUDED_

#if defined(TRV_USE_HIP)
#include <hipfft/hipfftXt.h>
#elif defined(TRV_USE_CUDA)  // !TRV_USE_HIP && TRV_USE_CUDA
#include <cufftXt.h>
#endif                       // TRV_USE_HIP
#include <fftw3.h>

#include <chrono>
#include <cmath>
#include <complex>
#include <vector>

#include "arrayops.hpp"
#include "monitor.hpp"
#include "parameters.hpp"
#include "maths.hpp"
#include "dataobjs.hpp"
#include "io.hpp"
#include "particles.hpp"

/// @cond DOXYGEN_DOC_MACROS
#ifdef __GNUC__
#define PURE __attribute__((pure))
#else
#define PURE
#endif
/// @endcond

namespace trvm = trv::maths;

namespace trv {

// ***********************************************************************
// Mesh grid
// ***********************************************************************

/**
 * @brief Discretely sampled field on a mesh grid from particle catalogues.
 *
 */
class MeshGrid {
 public:
  trv::ParameterSet params;                ///< parameter set
  std::string name;                        ///< mesh grid name
  std::vector<double> source_field;        ///< source field on mesh grid
  std::vector<fftw_complex> field;         ///< field on mesh grid
  std::vector< std::vector<double> > los;  ///< lines of sight
  double dr[3];                            ///< grid cell sizes
  double dk[3];                            ///< fundamental wavenumbers
  double vol;                              ///< mesh volume
  double vol_cell;                         ///< mesh grid cell volume

  // ---------------------------------------------------------------------
  // Life cycle
  // ---------------------------------------------------------------------

  /**
   * @brief Construct the mesh grid.
   *
   * @param params Parameter set.
   * @param plan_ini Flag for FFTW plan initialisation
   *                 (default is `true`).
   * @param name Grid name (default is "mesh-grid").
   */
  explicit MeshGrid(
    trv::ParameterSet& params,
    bool plan_ini = true,
    const std::string& name = "mesh-grid"
  );

  /**
   * @brief Destruct the mesh grid.
   */
  ~MeshGrid();

  /**
   * @brief (Re-)initialise the complex field on mesh grid.
   *
   * This is an explicit method to reset values of
   * @ref trv::MeshGrid->field to zeros.
   */
  void reset_density_field();

  // ---------------------------------------------------------------------
  // Operators & reserved methods
  // ---------------------------------------------------------------------

  /**
   * @brief Return mesh grid cell field value.
   *
   * @param gid Grid index.
   * @returns Field value.
   */
  PURE const fftw_complex& operator[](long long gid);

  // ---------------------------------------------------------------------
  // Mesh assignment
  // ---------------------------------------------------------------------

  /**
   * @brief Assign a weighted field to a mesh by interpolation scheme.
   *
   * @param particles Particle catalogue.
   * @param weights Particle Weights (default is null).
   */
  void assign_field(
    ParticleCatalogue& particles,
    fftw_complex* weights = nullptr
  );

  /**
   * @brief Load external field data mesh grid.
   *
   * @param field Source field.
   * @param indexing Mesh grid indexing convention, either @c "ijk"
   *                 (default) for C-style row-major ordering, or @c "xyz"
   *                 for Fortran-style column-major ordering.
   */
  void load_field(
    std::vector<double>& field,
    std::string indexing = "ijk"
  );

  // ---------------------------------------------------------------------
  // Field computations
  // ---------------------------------------------------------------------

  /**
   * @brief Weight the field by the reduced spherical harmonic.
   *
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void weight_by_ylm(int ell, int m);

  /**
   * @brief Weight the quadratic field by the reduced spherical harmonic.
   *
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void square_weight_by_ylm(int ell, int m);

  // ---------------------------------------------------------------------
  // Field transforms
  // ---------------------------------------------------------------------

  /**
   * @brief Fourier transform the field.
   *
   */
  void fourier_transform();

  /**
   * @brief Inverse Fourier transform the (FFT-transformed) field.
   */
  void inv_fourier_transform();

  // ---------------------------------------------------------------------
  // Field operations
  // ---------------------------------------------------------------------

  /**
   * @brief Apply mesh assignment window compensation in Fourier space.
   */
  void apply_assignment_compensation(std::string assignment);

  // ---------------------------------------------------------------------
  // One-point statistics
  // ---------------------------------------------------------------------

  /**
   * @brief Inverse Fourier transform a field @f$ f @f$ weighted by the
   *        reduced spherical harmonics restricted to a wavenumber band.
   *
   * This method computes the quantity
   * @f[
   *   F_{LM}(\vec{x}; k) = \frac{(2\pi)^3}{4\pi k^2}
   *     \int \frac{\mathrm{d}^3\,\vec{k}'}{(2\pi)^3}
   *       \mathrm{e}^{\mathrm{i} \vec{k}' \cdot \vec{x}}
   *       \delta^{(\mathrm{D})}(k' - k)
   *       y_{LM}(\hat{\vec{k}}) f(\vec{k}) \,.
   * @f]
   *
   * @see Eq. (42) in Sugiyama et al. (2019)
   *      [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param[in] field_fourier A Fourier-space field.
   * @param[in] ylm Reduced spherical harmonic on a mesh.
   * @param[in] k_band Band wavenumber.
   * @param[in] dk_band Band wavenumber width.
   * @param[out] k_eff Effective band wavenumber.
   * @param[out] nmodes Number of wavevector modes in band.
   */
  void inv_fourier_transform_ylm_wgtd_field_band_limited(
    MeshGrid& field_fourier, std::vector< std::complex<double> >& ylm,
    double k_band, double dk_band,
    double& k_eff, int& nmodes
  );

  /**
   * @brief Inverse Fourier transform a field @f$ f @f$ weighted by the
   *        spherical Bessel function and reduced spherical harmonics.
   *
   * This method computes the quantity
   * @f[
   *   F_{LM}(\vec{x}; r) = \mathrm{i}^L
   *     \int \frac{\mathrm{d}^3\,\vec{k}}{(2\pi)^3}
   *       \mathrm{e}^{\mathrm{i} \vec{k} \cdot \vec{x}}
   *       j_L(kr) y_{LM}(\hat{\vec{k}}) f(\vec{k}) \,.
   * @f]
   *
   * @see Eq. (49) in Sugiyama et al. (2019)
   *      [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param field_fourier A Fourier-space field.
   * @param ylm Reduced spherical harmonic on a mesh.
   * @param sjl Spherical Bessel function interpolator.
   * @param r Separation in configuration space.
   */
  void inv_fourier_transform_sjl_ylm_wgtd_field(
    MeshGrid& field_fourier,
    std::vector< std::complex<double> >& ylm,
    trvm::SphericalBesselCalculator& sjl,
    double r
  );

  // ---------------------------------------------------------------------
  // Misc
  // ---------------------------------------------------------------------

  /**
   * Calculate the normalisation factor @f$ 1/I_N @f$ for <i>N</i>-point
   * (for two-point statistics), where
   * @f[
   *   I_N = \int\mathrm{d}^3\,\vec{x} w(\vec{x})^N \bar{n}(\vec{x})^N \,.
   * @f]
   *
   * @param order Order @f$ N @f$ of the <i>N</i>-point statistics.
   * @returns norm_factor Normalisation factor.
   */
  double calc_powlaw_norm(int order);

 private:

};

}  // namespace trv

#endif
