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
 * @file field.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Mesh field (with one-point statistics) and
 *        pseudo two-point statistics.
 *
 * This module performs the assignment of catalogue particles with
 * weights to a mesh grid to construct discretely sampled fields, and
 * the Fourier transform and its inverse of the said fields.  It also
 * provides the corrections needed due to sampling effects, interlaced
 * 'shadow' mesh grid for reducing aliasing effects, and methods to
 * compute various constituent terms (one-point and pseudo two-point
 * statistics) in the estimators of two- and three-point statistics.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_FIELD_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_FIELD_HPP_INCLUDED_

#include <fftw3.h>

#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <vector>

#include "arrayops.hpp"
#include "monitor.hpp"
#include "parameters.hpp"
#include "maths.hpp"
#include "dataobjs.hpp"
#include "io.hpp"
#include "particles.hpp"

namespace trvm = trv::maths;

namespace trv {

// ***********************************************************************
// Mesh field
// ***********************************************************************

/**
 * @brief Discretely sampled field on a mesh grid from particle catalogues.
 *
 */
class MeshField {
 public:
  trv::ParameterSet params;  ///< parameter set
  std::string name;          ///< field name
  fftw_complex* field;       ///< complex field on mesh
  double dr[3];              ///< grid size in each dimension
  double dk[3];              ///< fundamental wavenumber in each dimension
  double vol;                ///< mesh volume
  double vol_cell;           ///< mesh grid cell volume

  // ---------------------------------------------------------------------
  // Life cycle
  // ---------------------------------------------------------------------

  /**
   * @brief Construct the mesh field.
   *
   * @param params Parameter set.
   * @param plan_ini Flag for FFTW plan initialisation
   *                 (default is `true`).
   * @param name Field name (default is "mesh-field").
   */
  MeshField(
    trv::ParameterSet& params,
    bool plan_ini = true,
    const std::string name = "mesh-field"
  );

  /**
   * @brief Construct the mesh field with external FFTW plans.
   *
   * @param params Parameter set.
   * @param transform External FFTW plan for Fourier transform.
   * @param inv_transform External FFTW plan for inverse
   *                      Fourier transform.
   * @param name Field name (default is "mesh-field").
   *
   * @overload
   */
  MeshField(
    trv::ParameterSet& params,
    fftw_plan& transform, fftw_plan& inv_transform,
    const std::string name = "mesh-field"
  );

  /**
   * @brief Destruct the mesh field.
   */
  ~MeshField();

  /**
   * @brief (Re-)initialise the complex field (and its shadow) on mesh.
   *
   * This is an explicit method to reset values of
   * @ref trv::MeshField.field (and its interlaced counterpart) to zeros.
   */
  void reset_density_field();

  // ---------------------------------------------------------------------
  // Operators & reserved methods
  // ---------------------------------------------------------------------

  /**
   * @brief Return mesh field grid cell value.
   *
   * @param gid Grid index.
   * @returns Field value.
   */
  const fftw_complex& operator[](long long gid);

  // ---------------------------------------------------------------------
  // Mesh assignment
  // ---------------------------------------------------------------------

  /**
   * @brief Assign a weighted field to a mesh by interpolation scheme.
   *
   * @param particles Particle catalogue.
   * @param weights Weight field.
   */
  void assign_weighted_field_to_mesh(
    ParticleCatalogue& particles, fftw_complex* weights
  );

  // ---------------------------------------------------------------------
  // Field computations
  // ---------------------------------------------------------------------

  /**
   * @brief Compute the unweighted field.
   *
   * This is is the number density field
   * @f[
   *   n(\vec{x}) = \sum_i \delta^{(\mathrm{D})}(\vec{x} - \vec{x}_i) \,.
   * @f]
   *
   * @param particles Particle catalogue.
   */
  void compute_unweighted_field(ParticleCatalogue& particles);

  /**
   * @brief Compute the unweighted field fluctuations.
   *
   * This is typically used for the number density fluctuations in
   * a simulation snapshot in a periodic box, with a global field value
   * calculated from the particle number and box volume subtracted
   * to compute the fluctuations, i.e.
   * @f[
   *   \delta{n}(\vec{x}) =
   *     \sum_i \delta^{(\mathrm{D})}(\vec{x} - \vec{x}_i) - N / V \,,
   * @f]
   * where @f$ N @f$ is the total number of particles and @f$ V @f$ the
   * box volume.
   *
   * @param particles Particle catalogue.
   */
  void compute_unweighted_field_fluctuations_insitu(
    ParticleCatalogue& particles
  );

  /**
   * @brief Compute the weighted field (fluctuations) further weighted by
   *        the reduced spherical harmonics.
   *
   * For a field (or its fluctuations) @f$ f @f$, this is
   * @f[
   *   f_{LM}(\vec{x}) = {\sum_i}{\vphantom{\sum}}'
   *     y_{LM}^*(\hat{\vec{x}}) w(\vec{x})
   *     \delta^{(\mathrm{D})}(\vec{x} - \vec{x}_i) \,,
   * @f]
   * where if a pair of catalogues are provided,
   * @f[
   *   {\sum_i}{\vphantom{\sum}}' =
   *     \sum_{i \in \mathrm{data}} - \alpha \sum_{i \in \mathrm{rand}}
   * @f]
   * with @f$ f \equiv \delta{n} @f$, and otherwise
   * @f[
   *   {\sum_i}{\vphantom{\sum}}' = \sum_{i \in \mathrm{data\ or\ rand}}
   * @f]
   * with @f$ f \equiv n @f$.
   *
   * @see Eq. (34) in Sugiyama et al. (2019)
   *      [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param particles_data (Data-source) particle catalogue.
   * @param particles_rand (Random-source) particle catalogue.
   * @param los_data (Data-source) particle lines of sight.
   * @param los_rand (Random-source) particle lines of sight.
   * @param alpha Alpha contrast.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void compute_ylm_wgtd_field(
    ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
    LineOfSight* los_data, LineOfSight* los_rand,
    double alpha, int ell, int m
  );

  /**
   * @brief Compute the weighted field further weighted by the
   *        reduced spherical harmonics.
   *
   * @param particles Particle catalogue.
   * @param los Particle lines of sight.
   * @param alpha Alpha contrast.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   *
   * @overload
   */
  void compute_ylm_wgtd_field(
    ParticleCatalogue& particles, LineOfSight* los,
    double alpha, int ell, int m
  );

  /**
   * @brief Compute the quadratic weighted field (fluctuations) further
   *        weighted by the reduced spherical harmonics.
   *
   * For a quadratic field (of fluctuations) @f$ f @f$, this is
   * @f[
   *   f_{LM}(\vec{x}) = {\sum_i}{\vphantom{\sum}}'
   *     y_{LM}^*(\hat{\vec{x}}) w(\vec{x})^2
   *     \delta^{(\mathrm{D})}(\vec{x} - \vec{x}_i) \,,
   * @f]
   * where if a pair of catalogues are provided,
   * @f[
   *   {\sum_i}{\vphantom{\sum}}' =
   *     \sum_{i \in \mathrm{data}} + \alpha^2 \sum_{i \in \mathrm{rand}} \,,
   * @f]
   * and otherwise
   * @f[
   *   {\sum_i}{\vphantom{\sum}}' = \sum_{i \in \mathrm{data\ or\ rand}} \,.
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
   */
  void compute_ylm_wgtd_quad_field(
    ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
    LineOfSight* los_data, LineOfSight* los_rand,
    double alpha, int ell, int m
  );

  /**
   * @brief Compute the quadratic weighted field (fluctuations) further
   *        weighted by the reduced spherical harmonics.
   *
   * @param particles Particle catalogue.
   * @param los Particle lines of sight.
   * @param alpha Alpha contrast.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   *
   * @overload
   */
  void compute_ylm_wgtd_quad_field(
    ParticleCatalogue& particles, LineOfSight* los,
    double alpha, int ell, int m
  );

  // ---------------------------------------------------------------------
  // Field transforms
  // ---------------------------------------------------------------------

  /**
   * @brief Fourier transform the field.
   *
   * If @c trv::MeshField.params.interlace is set to "true", interlacing
   * is performed where a phase factor is multiplied into the 'shadow'
   * complex field before the average of the complex field and its shadow
   * is taken.
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
   * @brief Apply wide-angle correction kernel in configuration space.
   *
   * This multiplies the field by a power law @f$ r^{- i - j} @f$ at order
   * @f$ (i, j) @f$ where @f$ r @f$ is the grid cell radial distance.
   */
  void apply_wide_angle_pow_law_kernel();

  /**
   * @brief Apply compensation needed after Fourier transform for
   *        assignment schemes.
   */
  void apply_assignment_compensation();

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
    MeshField& field_fourier, std::vector< std::complex<double> >& ylm,
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
    MeshField& field_fourier,
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
   * @param particles (Typically random-source) particle catalogue.
   * @param order Order @f$ N @f$ of the <i>N</i>-point statistics.
   * @returns norm_factor Normalisation factor.
   */
  double calc_grid_based_powlaw_norm(ParticleCatalogue& particles, int order);

 private:
  /// assignment window on mesh
  std::vector<double> window;
  /// window assignment order (default -1 if unassigned)
  int window_assign_order = -1;

  /// half-grid shifted complex field on mesh
  fftw_complex* field_s = nullptr;

  /// FFTW plan for Fourier transform of the field
  fftw_plan transform;
  /// FFTW plan for Fourier transform of the shadow field
  fftw_plan transform_s;
  /// FFTW plan for inverse Fourier transform of the field
  fftw_plan inv_transform;

  bool plan_ini = false;  ///< FFTW plan initialisation flag
  bool plan_ext = false;  ///< FFTW plan externality flag

  friend class FieldStats;

  // ---------------------------------------------------------------------
  // Mesh grid properties
  // ---------------------------------------------------------------------

  /**
   * @brief Return the grid cell index.
   *
   * @param i, j, k Grid index in each dimension.
   * @returns Grid cell index.
   */
  long long ret_grid_index(int i, int j, int k);

  /**
   * @brief Shift the grid indices on a discrete Fourier mesh grid.
   *
   * @param i, j, k Grid index in each dimension.
   */
  void shift_grid_indices_fourier(int& i, int& j, int& k);

  /**
   * @brief Get the grid cell position vector.
   *
   * @param[in] i, j, k Grid index in each dimension.
   * @param[out] rvec Position vector to the grid cell centre.
   */
  void get_grid_pos_vector(int i, int j, int k, double rvec[3]);

  /**
   * @brief Get the grid cell wavevector.
   *
   * @param[in] i, j, k Grid index in each dimension.
   * @param[out] kvec Position vector to the grid cell centre.
   */
  void get_grid_wavevector(int i, int j, int k, double kvec[3]);

  // ---------------------------------------------------------------------
  // Mesh assignment
  // ---------------------------------------------------------------------

  /**
   * @brief Assign weighted field to a mesh by the nearest-grid-point
   *        (NGP) scheme.
   *
   * @param particles Particle catalogue.
   * @param weight Particle weights.
   */
  void assign_weighted_field_to_mesh_ngp(
    ParticleCatalogue& particles, fftw_complex* weight
  );

  /**
   * @brief Assign weighted field to a mesh by the cloud-in-cell
   *        (CIC) scheme.
   *
   * @param particles Particle catalogue.
   * @param weight Particle weights.
   */
  void assign_weighted_field_to_mesh_cic(
    ParticleCatalogue& particles, fftw_complex* weight
  );

  /**
   * @brief Assign weighted field to a mesh by the triangular-shaped-cloud
   *        (TSC) scheme.
   *
   * @param particles Particle catalogue.
   * @param weight Particle weights.
   */
  void assign_weighted_field_to_mesh_tsc(
    ParticleCatalogue& particles, fftw_complex* weight
  );

  /**
   * @brief Assign weighted field to a mesh by the piecewise cubib spline
   *        (PCS) scheme.
   *
   * @param particles Particle catalogue.
   * @param weight Particle weights.
   */
  void assign_weighted_field_to_mesh_pcs(
    ParticleCatalogue& particles, fftw_complex* weight
  );

  /**
   * @brief Calculate the interpolation window at each mesh grid
   *        in Fourier space for different assignment schemes.
   *
   * @param i, j, k Grid cell indices.
   * @param order Order of the assignment scheme
   *              (default is 0 as placeholder).
   * @returns Window value.
   */
  double calc_assignment_window_in_fourier(int i, int j, int k, int order = 0);

  /**
   * @brief Compute the interpolation window at each mesh grid
   *        in Fourier space for different assignment schemes and store
   *        the values in @ref trv::MeshField.window.
   *
   * @param order Order of the assignment scheme
   *              (default is 0 as placeholder).
   */
  void compute_assignment_window_in_fourier(int order = 0);
};


// ***********************************************************************
// Field statistics
// ***********************************************************************

/**
 * @brief Field (pseudo-two-point) statistics.
 *
 * This provides the computation of both binned and unbinned (pseudo)
 * two-point statistiscs that aid the computation of full two- and
 * three-point statistics.
 *
 */
class FieldStats {
 public:
  std::vector<int> nmodes;  ///< number of wavevector modes in bins
  std::vector<int> npairs;  ///< number of separation pairs in bins
  std::vector<double> k;    ///< average wavenumber in bins
  std::vector<double> r;    ///< average separation in bins
  /// shot-noise power in bins
  std::vector< std::complex<double> > sn;
  /// pseudo power spectrum in bins
  std::vector< std::complex<double> > pk;
  /// pseudo two-point correlation function in bins
  std::vector< std::complex<double> > xi;

  // ---------------------------------------------------------------------
  // Life cycle
  // ---------------------------------------------------------------------

  /**
   * @brief Construct pseudo two-point statistics.
   *
   * @param params Parameter set.
   * @param plan_ini Flag for FFTW plan initialisation
   *                 (default is `true`).
   */
  FieldStats(trv::ParameterSet& params, bool plan_ini = true);

  /**
   * @brief Destruct two-point statistics.
   */
  ~FieldStats();

  /**
   * @brief Reset two-point statistics.
   */
  void reset_stats();

  // ---------------------------------------------------------------------
  // Binned statistics
  // ---------------------------------------------------------------------

  /**
   * @brief Compute binned two-point statistics in Fourier space.
   *
   * For a pair of Fourier-space fields @f$ f_{a,b} @f$, this computes
   * @f[
   *   \int \frac{\mathrm{d}^2\,\hat{\vec{k}}}{4\pi} y_{LM}(\hat{\vec{k}})
   *     W(\vec{k})^{-2} [
   *       f_a(\vec{k}) f_b(\vec{k}) - P_\mathrm{shot} C_1(\vec{k})
   *     ] \,,
   * @f]
   * where @f$ W(\vec{k}) @f$ is the mesh assignment window in Fourier
   * space, @f$ P_\mathrm{shot} @f$ is the shot noise amplitude, and
   * @f$ C_1 @f$ is the mode-dependent aliasing function.
   *
   * @see Eq. (20) in Jing (2004)
   *      [<a href="https://arxiv.org/abs/astro-ph/0409240">astro-ph/0409240</a>].
   *
   * @param field_a First field.
   * @param field_b Second field.
   * @param shotnoise_amp Shot-noise amplitude.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   * @param kbinning Wavenumber binning.
   * @throws trv::sys::InvalidDataError When @p field_a and @p field_b
   *                                    have incompatible physical
   *                                    properties.
   *
   * @note @p field_a and @p field_b are Fourier-space fields and their
   *       entries are arranged in the FFTW convention (i.e. shifted).
   */
  void compute_ylm_wgtd_2pt_stats_in_fourier(
    MeshField& field_a, MeshField& field_b, std::complex<double> shotnoise_amp,
    int ell, int m, trv::Binning& kbinning
  );

  /**
   * @brief Compute binned two-point statistics in configuration space.
   *
   * This involves an inverse Fourier transform added to the algorithm
   * in @ref trv::FieldStats::compute_ylm_wgtd_2pt_stats_in_fourier
   * before binning.
   *
   * @param field_a First field.
   * @param field_b Second field.
   * @param shotnoise_amp Shot-noise amplitude.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   * @param rbinning Separation binning.
   *
   * @note @p field_a and @p field_b are Fourier-space fields and their
   *       entries are arranged in the FFTW convention (i.e. shifted).
   */
  void compute_ylm_wgtd_2pt_stats_in_config(
    MeshField& field_a, MeshField& field_b, std::complex<double> shotnoise_amp,
    int ell, int m, trv::Binning& rbinning
  );

  /**
   * @brief Compute binned uncoupled three-point correlation function
   *        shot noise.
   *
   * This computes the quantity
   * @f[
   *   \frac{(-1)^{\ell_a + \ell_b}}{N_\mathrm{mode}(r) V_\mathrm{cell}}
   *   \int \frac{\mathrm{d}^2\,\hat{\vec{r}}}{4\pi}
   *     y^*_{\ell_a m_a}(\hat{\vec{r}}) y^*_{\ell_b m_b}(\hat{\vec{r}})
   *   \int \frac{\mathrm{d}^3\,\vec{k}}{(2\pi)^3}
   *     \mathrm{e}^{\mathrm{i} \vec{k} \cdot \vec{r}}
   *     W(\vec{k})^{-2} [f_a(\vec{k}) f_b(\vec{k}) - S C_1(\vec{k})] \,,
   * @f]
   * where @f$ S @f$ is the shot noise amplitude (see
   * @ref trv::FieldStats::compute_ylm_wgtd_2pt_stats_in_fourier
   * for other notations).
   *
   * This method involves additional reduced spherical-harmonic weights
   * after the inverse Fourier transform as in
   * @ref trv::FieldStats::compute_ylm_wgtd_2pt_stats_in_config
   * before binning and additional normalisation.
   *
   * @see Eq. (51) in Sugiyama et al. (2019)
   *      [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param field_a First field.
   * @param field_b Second field.
   * @param ylm_a Reduced spherical harmonics over the first field mesh.
   * @param ylm_b Reduced spherical harmonics over the second field mesh.
   * @param shotnoise_amp Shot-noise amplitude.
   * @param rbinning Separation binning.
   */
  void compute_uncoupled_shotnoise_for_3pcf(
    MeshField& field_a, MeshField& field_b,
    std::vector< std::complex<double> >& ylm_a,
    std::vector< std::complex<double> >& ylm_b,
    std::complex<double> shotnoise_amp,
    trv::Binning& rbinning
  );

  /**
   * @brief Compute unbinned uncoupled bispectrum shot noise on the FFT
   *        mesh grid.
   *
   * This computes the quantity
   * @f[
   *   \int \mathrm{d}^3\,\vec{r}
   *     j_{\ell_a}(k_a r) j_{\ell_b}(k_b r)
   *     y^*_{\ell_a m_a}(\hat{\vec{x}}) y^*_{\ell_b m_b}(\hat{\vec{x}})
   *   \int \frac{\mathrm{d}^3\,\vec{k}}{(2\pi)^3}
   *     \mathrm{e}^{\mathrm{i} \vec{k} \cdot \vec{r}}
   *     W(\vec{k})^{-2} [f_a(\vec{k}) f_b(\vec{k}) - S C_1(\vec{k})]
   * @f]
   * (see
   * @ref trv::FieldStats::compute_ylm_wgtd_2pt_stats_in_fourier
   * for other notations).
   *
   * This method involves additional reduced spherical-harmonic weights
   * after the inverse Fourier transform as in
   * @ref trv::FieldStats::compute_ylm_wgtd_2pt_stats_in_config
   * before binning and additional normalisation.
   *
   * @see Eq. (45) in Sugiyama et al. (2019)
   *      [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param field_a First field.
   * @param field_b Second field.
   * @param ylm_a Reduced spherical harmonics over the first
   *              field mesh.
   * @param ylm_b Reduced spherical harmonics over the second
   *              field mesh.
   * @param sj_a First spherical Bessel function.
   * @param sj_b Second spherical Bessel function.
   * @param shotnoise_amp Shot-noise amplitude.
   * @param k_a, k_b Wavenumbers at which the shot noise is evaluated.
   * @returns Unbinned uncoupled bispectrum shot noise.
   */
  std::complex<double> compute_uncoupled_shotnoise_for_bispec_per_bin(
    MeshField& field_a, MeshField& field_b,
    std::vector< std::complex<double> >& ylm_a,
    std::vector< std::complex<double> >& ylm_b,
    trvm::SphericalBesselCalculator& sj_a,
    trvm::SphericalBesselCalculator& sj_b,
    std::complex<double> shotnoise_amp,
    double k_a, double k_b
  );

  /**
   * @brief Record binned vectors given a binning scheme.
   *
   * @param binning Binning.
   * @param save_file Saved filename if non-empty.
   */
  trv::BinnedVectors record_binned_vectors(
    trv::Binning& binning, const std::string& save_file
  );

 private:
  trv::ParameterSet params;  ///< parameter set
  double dr[3];              ///< grid size in each dimension
  double dk[3];              ///< fundamental wavenumber in each dimension
  double vol;                ///< mesh volume
  double vol_cell;           ///< mesh grid cell volume

  /// FFTW buffer array for pseudo-two-point statistics
  fftw_complex* twopt_3d = nullptr;
  /// FFTW plan for inverse Fourier transform
  fftw_plan inv_transform;
  /// FFTW plan initialisation flag
  bool plan_ini = false;

  /// shot-noise aliasing scale-dependence function
  std::vector<double> alias_sn;
  /// shot-noise aliasing function initialisation flag
  bool alias_ini = false;

  // ---------------------------------------------------------------------
  // Utilities
  // ---------------------------------------------------------------------

  /**
   * @brief Check if two mesh fields have compatible properties.
   *
   * @param field_a First field.
   * @param field_b Second field.
   * @returns { @c true , @c false }
   */
  bool if_fields_compatible(MeshField& field_a, MeshField& field_b);

  /**
   * @brief Resize binned two-point statistics.
   *
   * @param num_bins Number of bins.
   */
  void resize_stats(int num_bins);

  /**
   * @brief Return the grid cell index.
   *
   * @param i, j, k Grid index in each dimension.
   * @returns Grid cell index.
   *
   * @see trv::MeshField::ret_grid_index
   */
  long long ret_grid_index(int i, int j, int k);

  /**
   * @brief Shift the grid indices on a discrete Fourier mesh grid.
   *
   * @param i, j, k Grid index in each dimension.
   *
   * @see trv::MeshField::shift_grid_indices_fourier
   */
  void shift_grid_indices_fourier(int& i, int& j, int& k);

  // ---------------------------------------------------------------------
  // Sampling corrections
  // ---------------------------------------------------------------------

  /**
   * @brief Return the shot-noise aliasing scale-dependence function
   *        @f$ C_1(\vec{k}) @f$ at each mesh grid.
   *
   * @see Eqs. (45) and (46) in Sugiyama et al. (2019)
   *      [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>]
   *      and Jing (2004)
   *      [<a href="https://arxiv.org/abs/astro-ph/0409240">astro-ph/0409240</a>].
   *
   * @returns Aliasing function.
   */
  std::function<double(int, int, int)> ret_calc_shotnoise_aliasing();

  /**
   * @brief Get the square-sine arguments for the shot-noise aliasing
   *        function.
   *
   * @param i, j, k Grid index in each dimension.
   * @param cx2, cy2, cz2 Square-sine arguments.
   */
  void get_shotnoise_aliasing_sin2(
    int i, int j, int k, double& cx2, double& cy2, double& cz2
  );

  /**
   * Calculate the shot-noise aliasing function for the
   * nearest-grid-point (NGP) assignment scheme.
   *
   * @param i, j, k Grid indices.
   * @returns Function value.
   */
  double calc_shotnoise_aliasing_ngp(int i, int j, int k);

  /**
   * Calculate the shot-noise aliasing function for the
   * cloud-in-cell (CIC) assignment scheme.
   *
   * @param i, j, k Grid indices.
   * @returns Function value.
   */
  double calc_shotnoise_aliasing_cic(int i, int j, int k);

  /**
   * Calculate the shot-noise aliasing function for the
   * triangular-shaped-cloud (TSC) assignment scheme.
   *
   * @param i, j, k Grid indices.
   * @returns Function value.
   */
  double calc_shotnoise_aliasing_tsc(int i, int j, int k);

  /**
   * Calculate the shot-noise aliasing function for the
   * piecewise-cubic-spline (PCS) assignment scheme.
   *
   * @param i, j, k Grid indices.
   * @returns Function value.
   */
  double calc_shotnoise_aliasing_pcs(int i, int j, int k);

  /**
   * Compute the shot-noise aliasing function.
   *
   */
  void compute_shotnoise_aliasing();
};

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_FIELD_HPP_INCLUDED_
