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
 * @file particles.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Particle containers with I/O methods and operations.
 *
 * This module defines a particle catalogue object with I/O methods,
 * summary information and its computations, and methods to offset
 * particle coordinates (in particular in a mesh grid box).
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "monitor.hpp"

namespace trv {

/**
 * @brief Particle data container.
 *
 */
struct ParticleData {
  double pos[3];  ///< particle position vector
  double nz;      ///< redshift-dependent expected number density
  double ws;      ///< particle sample weight
  double wc;      ///< particle clustering weight
  double w;       ///< particle overall weight
};

/**
 * @brief Particle catalogue.
 *
 * The catalogue object contains particle data and summary information,
 * as well as methods for computing its attributes.
 *
 */
class ParticleCatalogue {
 public:
  std::string source;  ///< catalogue source

  ParticleData* pdata;  ///< particle data

  int ntotal;      ///< total number of particles
  double wtotal;   ///< total overall weight of particles
  double wstotal;  ///< total sample weight of particles

  double pos_min[3];   ///< minimum values of particle coordinates
  double pos_max[3];   ///< maximum values of particle coordinates
  double pos_span[3];  ///< span of particle coordinates

  // ---------------------------------------------------------------------
  // Life cycle
  // ---------------------------------------------------------------------

  /**
   * @brief Construct the particle catalogue with initial values.
   *
   * @param verbose Logger verbosity (default is -1 for no change).
   *
   * @note @p verbose is optional and used to control the logging level
   *       when an instance is created as a standalone without
   *       @c ParameterSet .
   */
  ParticleCatalogue(int verbose = -1);

  /**
   * @brief Destruct the particle catalogue.
   */
  ~ParticleCatalogue();

  /**
   * @brief Initialise particle data container.
   *
   * @attention This method does not set the values of
   *            @ref trv::ParticleCatalogue.pdata,
   *            @ref trv::ParticleCatalogue.wtotal,
   *            @ref trv::ParticleCatalogue.wstotal,
   *            @ref trv::ParticleCatalogue.pos_min or
   *            @ref trv::ParticleCatalogue.pos_max, or
   *            @ref trv::ParticleCatalogue.pos_span.
   *
   * @param num Number of data units (i.e. particles).
   */
  void initialise_particles(const int num);

  /**
   * @brief Finalise particle data container.
   *
   */
  void finalise_particles();

  /**
   * @brief Reset particle data container.
   *
   * This is an explicit method to free the resources occupied by
   * @ref trv::ParticleCatalogue.pdata and may be called outside the
   * class destructor.
   */
  void reset_particles();

  // ---------------------------------------------------------------------
  // Operators & reserved methods
  // ---------------------------------------------------------------------

  /**
   * @brief Return individual particle information.
   *
   * @param pid Particle index.
   * @returns Individual particle data.
   */
  ParticleData& operator[](const int pid);

  // ---------------------------------------------------------------------
  // Data I/O
  // ---------------------------------------------------------------------

  /**
   * @brief Read in a catalogue file.
   *
   * @param catalogue_filepath Catalogue file path.
   * @param catalogue_columns Catalogue data column names
   *                          (comma-separated without space).
   * @param volume Catalogue volume (default is 0.) used for computing
   *               the default 'nz' value when the field is missing.
   * @returns Exit status.
   */
  int load_catalogue_file(
    const std::string& catalogue_filepath,
    const std::string& catalogue_columns,
    double volume = 0.
  );

  /**
   * @brief Read in particle data.
   *
   * @param x, y, z, nz, ws, wc Particle data by column.
   * @returns Exit status.
   */
  int load_particle_data(
    std::vector<double> x, std::vector<double> y, std::vector<double> z,
    std::vector<double> nz, std::vector<double> ws, std::vector<double> wc
  );

  // ---------------------------------------------------------------------
  // Catalogue properties
  // ---------------------------------------------------------------------

  /**
   * @brief Calculate total overall weight of particles.
   *
   * @attention This method resets @ref trv::ParticleCatalogue::wtotal
   *            and @ref trv::ParticleCatalogue::wstotal.
   */
  void calc_total_weights();

  /**
   * @brief Calculate the extents of particle positions.
   *
   * @param init Initialisation flag (default is `true`).
   *
   * @attention This method resets @ref trv::ParticleCatalogue::pos_min,
   *            @ref trv::ParticleCatalogue::pos_max and
   *            @ref trv::ParticleCatalogue::pos_span.
   */
  void calc_pos_extents(bool init = true);

  // ---------------------------------------------------------------------
  // Catalogue operations
  // ---------------------------------------------------------------------

  /**
   * @brief Offset particle positions by a given vector.
   *
   * The position specified by the input vector is the new origin.
   *
   * @param dpos (Subtractive) offset position vector.
   */
  void offset_coords(const double dpos[3]);

  /**
   * @brief Offset particle positions for periodic boundary conditions.
   *
   * @param boxsize Periodic box size in each dimension.
   */
  void offset_coords_for_periodicity(const double boxsize[3]);

  /**
   * @brief Centre a catalogue in a box.
   *
   * @param catalogue Particle catalogue.
   * @param boxsize Box size in each dimension.
   */
  static void centre_in_box(
    ParticleCatalogue& catalogue,
    const double boxsize[3]
  );

  /**
   * @brief Centre a pair of catalogues in a box.
   *
   * The secondary catalogue's centre is used as the reference point
   * to also offset the primary catalogue's particle positions.
   *
   * @param catalogue Primary particle catalogue.
   * @param catalogue_ref Secondary reference particle catalogue.
   * @param boxsize Box size in each dimension.
   *
   * @overload
   */
  static void centre_in_box(
    ParticleCatalogue& catalogue, ParticleCatalogue& catalogue_ref,
    const double boxsize[3]
  );

  /**
   * @brief Pad a catalogue in a box.
   *
   * The amount of padding is a fraction of the box size.
   *
   * @param catalogue Particle catalogue.
   * @param boxsize Box size in each dimension.
   * @param boxsize_pad Box size padding factor in each dimension.
   */
  static void pad_in_box(
    ParticleCatalogue& catalogue,
    const double boxsize[3], const double boxsize_pad[3]
  );

  /**
   * @brief Pad a pair of catalogues in a box.
   *
   * The amount of padding is a fraction of the box size.
   *
   * The secondary catalogue's extents are used as the reference points
   * to also offset the primary catalogue's particle positions.
   *
   * @param catalogue Primary particle catalogue.
   * @param catalogue_ref Secondary reference particle catalogue.
   * @param boxsize Box size in each dimension.
   * @param boxsize_pad Box size padding factor in each dimension.
   *
   * @overload
   */
  static void pad_in_box(
    ParticleCatalogue& catalogue, ParticleCatalogue& catalogue_ref,
    const double boxsize[3], const double boxsize_pad[3]
  );

  /**
   * @brief Pad a catalogue in a box.
   *
   * The amount of padding is a multiple of the mesh grid size.
   *
   * @param catalogue Particle catalogue.
   * @param boxsize Box size in each dimension.
   * @param ngrid Grid number in each dimension.
   * @param ngrid_pad Grid number factor for padding.
   */
  static void pad_grids(
    ParticleCatalogue& catalogue,
    const double boxsize[3], const int ngrid[3], const double ngrid_pad[3]
  );

  /**
   * @brief Pad a catalogue in a box.
   *
   * The amount of padding is a multiple of the mesh grid size.
   *
   * The secondary catalogue's extents are used as the reference points
   * to also offset the primary catalogue's particle positions.
   *
   * @param catalogue Primary particle catalogue.
   * @param catalogue_ref Secondary reference particle catalogue.
   * @param boxsize Box size in each dimension.
   * @param ngrid Grid number in each dimension.
   * @param ngrid_pad Grid number factor for padding.
   *
   * @overload
   */
  static void pad_grids(
    ParticleCatalogue& catalogue, ParticleCatalogue& catalogue_ref,
    const double boxsize[3], const int ngrid[3], const double ngrid_pad[3]
  );
};

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_
