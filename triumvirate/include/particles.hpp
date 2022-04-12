/**
 * @file particles.hpp
 * @brief Particle containers with I/O methods and operations, as well as
 *        line-of-sight vectors.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "common.hpp"

/// TODO: Remove in public releases.
/// Set placeholder 'nz' default value.
#ifdef DBGNZ
  double NZ_DEFAULT = 5.e-4;
#else
  double NZ_DEFAULT = 0.;
#endif

/**
 * Line-of-sight vector.
 *
 */
struct LineOfSight {
  double pos[3];  ///< 3-d position vector
};

/**
 * Particle catalogue containing particle data and summary information.
 *
 */
class ParticleCatalogue {
 public:
  struct ParticleData {
    double pos[3];  ///< particle position vector
    double nz;  ///< expected redshift-dependent number density in-situ
    double ws;  ///< particle systematic weight
    double wc;  ///< particle clustering weight
    double w;  ///< particle overall weight
  }* pdata;  ///< particle data

  int ntotal;  ///< total number of particles

  double wtotal;  ///< total systematic weight of particles

  double pos_min[3];  ///< minimum values of particle positions
  double pos_max[3];  ///< maximum values of particle positions

  /**
   * Initialise the particle container with default values.
   */
  ParticleCatalogue () {
    /// Set default values.
    this->pdata = NULL;
    this->ntotal = 0;
    this->wtotal = 0.;
    this->pos_min[0] = 0.; this->pos_min[1] = 0.; this->pos_min[2] = 0.;
    this->pos_max[0] = 0.; this->pos_max[1] = 0.; this->pos_max[2] = 0.;
  }

  /**
   * Destruct the particle container.
   */
  ~ParticleCatalogue() {
    this->finalise_particles();
  }

  /**
   * Return individual particle information.
   *
   * @param pid Particle index.
   * @returns Individual particle data.
   */
  ParticleData& operator[](int pid) {
    return this->pdata[pid];
  }

  /**
   * Initialise particle data container.
   *
   * @param num Number of particles.
   */
  void _initialise_particles(const int num) {
    /// Check the total number of particles.
    if (num <= 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[WARN] (+%s) Number of particles is negative.\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
      );
      return;
    }
    this->ntotal = num;

    /// Renew particle data.
    delete[] this->pdata; this->pdata = NULL;
    this->pdata = new ParticleData[this->ntotal];

    /// Determine memory usage.
    gbytesMem += double(this->ntotal)
      * sizeof(struct ParticleData) / BYTES_PER_GBYTES;
  }

  /**
   * Finalise particle data and information.
   */
  void finalise_particles() {
    /// Free particle usage.
    if (this->pdata != NULL) {
      delete[] this->pdata; this->pdata = NULL;
      gbytesMem -= double(this->ntotal)
        * sizeof(struct ParticleData) / BYTES_PER_GBYTES;
    }
  }

  /**
   * Read in particle data from a file.
   *
   * @param particles_file Particle data file path.
   * @param names Field names, comma-separated without space, in the file.
   * @returns Exit status.
   */
  int read_particle_data_from_file(
    std::string& particles_file,
    const std::string& names
  ) {
    this->source = particles_file;

    /// Initialise the counter for the number of lines/particles.
    int num_lines = 0;

    /// Extract field names and ordering.
    std::vector<std::string> fields;
    std::istringstream iss(names);
    std::string name;
    while (std::getline(iss, name, ',')) {
      fields.push_back(name);
    }

    /// CAVEAT: Hard-wired ordered field names.
    std::vector<std::string> names_ordered = {"x", "y", "z", "nz", "ws", "wc"};
    std::vector<int> name_indices(names_ordered.size(), -1);
    for (int idx_name = 0; idx_name < names_ordered.size(); idx_name++) {
      ptrdiff_t col_idx = std::distance(
        fields.begin(),
        find(fields.begin(), fields.end(), names_ordered[idx_name])
      );
      if (0 <= col_idx && col_idx < fields.size()) {
        name_indices[idx_name] = col_idx;
      }
    }

    if (name_indices[3] == -1) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[WARN] (+%s) Catalogue 'nz' field is unavailable, "
          "which may raise errors in some computations (source=%s).\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          this->source.c_str()
        );
      }
    }

    /// Check and size up data from the file.
    std::ifstream fin;

    fin.open(particles_file.c_str(), std::ios::in);

    if (fin.fail()) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Failed to open file '%s'.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          this->source.c_str()
        );
      }
      fin.close();
      return -1;
    }

    std::string str_line;
    while (getline(fin, str_line)) {
      /// Terminate at the end of file.
      if (!fin) {break;}

      /// Skip empty lines or comment lines.
      if (str_line[0] == '#' || str_line.empty()) {continue;}

      /// Count the line as a particle.
      num_lines++;
    }

    fin.close();

    /// Fill in particle data.
    this->_initialise_particles(num_lines);

    fin.open(particles_file.c_str(), std::ios::in);

    num_lines = 0;  // reset

    double nz, ws, wc;
    double entry;
    while (getline(fin, str_line)) {
      /// Terminate at the end of file.
      if (!fin) {break;}

      /// Skip empty lines or comment lines.
      if (str_line[0] == '#' || str_line.empty()) {continue;}

      /// Extract row entries.
      std::stringstream ss(
        str_line, std::ios_base::out | std::ios_base::in | std::ios_base::binary
      );

      std::vector<double> row;
      while (ss >> entry) {
        row.push_back(entry);
      }

      /// Add the current line as a particle.
      this->pdata[num_lines].pos[0] = row[name_indices[0]];  // x
      this->pdata[num_lines].pos[1] = row[name_indices[1]];  // y
      this->pdata[num_lines].pos[2] = row[name_indices[2]];  // z

      if (name_indices[3] != -1) {
        nz = row[name_indices[3]];
      } else {
        nz = NZ_DEFAULT;  // default value
      }

      if (name_indices[4] != -1) {
        ws = row[name_indices[4]];
      } else {
        ws = 1.;  // default value
      }

      if (name_indices[5] != -1) {
        wc = row[name_indices[5]];
      } else {
        wc = 1.;  // default value
      }

      this->pdata[num_lines].nz = nz;
      this->pdata[num_lines].ws = ws;
      this->pdata[num_lines].wc = wc;
      this->pdata[num_lines].w = ws * wc;

      num_lines++;
    }

    fin.close();

    /// Calculate weight sum.
    this->_calc_weighted_total();

    /// Calculate extreme particle positions.
    this->_calc_pos_min_and_max(true);

    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[INFO] (+%s) Catalogue loaded: %d particles with "
        "total systematic weights %.3f (source=%s).\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
        this->ntotal, this->wtotal, this->source.c_str()
      );
    }

    return 0;
  }

  /**
   * Read in particle data from external data.
   *
   * @param x Particle x-positions.
   * @param y Particle y-positions.
   * @param z Particle z-positions.
   * @param nz Particle redshift-dependent mean number density in-situ.
   * @param ws Particle systematic weights.
   * @param wc Particle clustering weights.
   * @returns Exit status.
   */
  int read_particle_data(
    std::vector<double> x, std::vector<double> y, std::vector<double> z,
    std::vector<double> nz, std::vector<double> ws, std::vector<double> wc
  ) {
    this->source = "runtime";

    /// Check array sizes.
    if (!(
      x.size() == y.size() && y.size() == z.size() && z.size() == nz.size()
      && nz.size() == ws.size() && ws.size() == wc.size()
    )) {
      return -1;
    }

    /// Fill in particle data.
    int ntotal = nz.size();

    this->_initialise_particles(ntotal);

    for (int pid = 0; pid < ntotal; pid++) {
      this->pdata[pid].pos[0] = x[pid];
      this->pdata[pid].pos[1] = y[pid];
      this->pdata[pid].pos[2] = z[pid];
      this->pdata[pid].nz = nz[pid];
      this->pdata[pid].ws = ws[pid];
      this->pdata[pid].wc = wc[pid];
      this->pdata[pid].w = ws[pid] * wc[pid];
    }

    /// Calculate weight sum.
    this->_calc_weighted_total();

    /// Calculate extreme particle positions.
    this->_calc_pos_min_and_max(true);

    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[INFO] (+%s) Catalogue constructed: %d particles with "
        "total systematic weights %.3f (source=%s).\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
        this->ntotal, this->wtotal, this->source.c_str()
      );
    }

    return 0;
  }

  /**
   * Calculate total systematic weights of particles.
   */
  void _calc_weighted_total() {
    if (this->pdata == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    double wtotal = 0.;
    for (int pid = 0; pid < this->ntotal; pid++) {
      wtotal += this->pdata[pid].ws;
    }

    this->wtotal = wtotal;
  }

  /**
   * Calculate extreme particle positions.
   *
   * @param verbose Print out particle coordinate extents if `true`.
   */
  void _calc_pos_min_and_max(bool verbose=false) {
    if (this->pdata == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    /// Initialise minimum and maximum values with the first
    /// data entry/row.
    double min[3], max[3];
    for (int iaxis = 0; iaxis < 3; iaxis++) {
      min[iaxis] = this->pdata[0].pos[iaxis];
      max[iaxis] = this->pdata[0].pos[iaxis];
    }

    /// Update minimum and maximum values line-by-line.
    for (int id = 0; id < this->ntotal; id++) {
      for (int iaxis = 0; iaxis < 3; iaxis++) {
        if (min[iaxis] > this->pdata[id].pos[iaxis]) {
          min[iaxis] = this->pdata[id].pos[iaxis];
        }
        if (max[iaxis] < this->pdata[id].pos[iaxis]) {
          max[iaxis] = this->pdata[id].pos[iaxis];
        }
      }
    }

    for (int iaxis = 0; iaxis < 3; iaxis++) {
      this->pos_min[iaxis] = min[iaxis];
      this->pos_max[iaxis] = max[iaxis];
    }

    if (currTask == 0 && verbose) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[INFO] (+%s) Extents of particle coordinates: "
        "{'x': (%.3f, %.3f), 'y': (%.3f, %.3f), 'z': (%.3f, %.3f)} "
        "(source=%s).\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
        this->pos_min[0], this->pos_max[0],
        this->pos_min[1], this->pos_max[1],
        this->pos_min[2], this->pos_max[2],
        this->source.c_str()
      );
    }
  }

  /**
   * Offset particle positions by a given vector (as the new origin).
   *
   * @param dpos (Subtractive) offset position vector.
   */
  void offset_coords(const double dpos[3]) {
    if (this->pdata == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    /// Make adjustments required.
    for (int pid = 0; pid < this->ntotal; pid++) {
      for (int iaxis = 0; iaxis < 3; iaxis++) {
        this->pdata[pid].pos[iaxis] -= dpos[iaxis];
      }
    }

    /// Recalculate extreme particle positions.
    this->_calc_pos_min_and_max(true);
  }

  /**
   * Offset particle positions by centring the catalogue inside the
   * specified box.
   *
   * @param boxsize Boxsize in each dimension.
   */
  void offset_coords_for_centring(const double boxsize[3]) {
    if (this->pdata == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    /// Calculate adjustments needed.
    this->_calc_pos_min_and_max();

    double xmid = (this->pos_min[0] + this->pos_max[0]) / 2.;
    double ymid = (this->pos_min[1] + this->pos_max[1]) / 2.;
    double zmid = (this->pos_min[2] + this->pos_max[2]) / 2.;

    double dvec[3] = {
      boxsize[0]/2. - xmid, boxsize[1]/2. - ymid, boxsize[2]/2. - zmid
    };

    /// Centre the catalogue in box.
    for (int pid = 0; pid < this->ntotal; pid++) {
      for (int iaxis = 0; iaxis < 3; iaxis++) {
        this->pdata[pid].pos[iaxis] += dvec[iaxis];
      }
    }

    /// Recalculate extreme particle positions.
    this->_calc_pos_min_and_max(true);
  }

  /**
   * Offset particle positions for periodic boundary conditions.
   *
   * @param boxsize Periodic boxsize in each dimension.
   */
  void offset_coords_for_periodicity(const double boxsize[3]) {
    /// Adjust the box.
    for (int pid = 0; pid < this->ntotal; pid++) {
      for (int iaxis = 0; iaxis < 3; iaxis++) {
        if (
          this->pdata[pid].pos[iaxis] >= boxsize[iaxis]
        ) {
          this->pdata[pid].pos[iaxis] -= boxsize[iaxis];
        } else if (
          this->pdata[pid].pos[iaxis] < 0.
        ) {
          this->pdata[pid].pos[iaxis] += boxsize[iaxis];
        }
      }
    }

    /// Recalculate extreme particle positions.
    this->_calc_pos_min_and_max(true);
  }

  /**
   * Centre a pair of catalogues in a box, with the secondary catalogue
   * as the reference.
   *
   * @param particles_data (Data-source) particle catalogue.
   * @param particles_rand (Random-source) particle catalogue.
   * @param boxsize Box size in each dimension.
   */
  static void centre_pair_in_box(
    ParticleCatalogue& particles_data,
    ParticleCatalogue& particles_rand,
    const double boxsize[3]
  ) {
    /// Calculate adjustments needed using the random-source catalogue.
    particles_rand._calc_pos_min_and_max();

    double xmid = (particles_rand.pos_min[0] + particles_rand.pos_max[0]) / 2.;
    double ymid = (particles_rand.pos_min[1] + particles_rand.pos_max[1]) / 2.;
    double zmid = (particles_rand.pos_min[2] + particles_rand.pos_max[2]) / 2.;

    double dvec[3] = {
      xmid - boxsize[0]/2., ymid - boxsize[1]/2., zmid - boxsize[2]/2.
    };

    /// Centre the catalogue in box.
    particles_data.offset_coords(dvec);
    particles_rand.offset_coords(dvec);
  }

  /**
   * Align a pair of catalogues in a box for FFT by grid shift.
   *
   * @param particles_data (Data-source) particle catalogue.
   * @param particles_rand (Random-source) particle catalogue.
   * @param boxsize Box size in each dimension.
   * @param ngrid Grid number in each dimension.
   * @param ngrid_pad Grid number factor for padding.
   */
  static void pad_pair_in_box(
    ParticleCatalogue& particles_data,
    ParticleCatalogue& particles_rand,
    const double boxsize[3],
    const int ngrid[3],
    const double ngrid_pad[3]
  ) {
    /// Calculate adjustments needed using the random-source catalogue.
    particles_rand._calc_pos_min_and_max();

    double dvec[3] = {
      particles_rand.pos_min[0],
      particles_rand.pos_min[1],
      particles_rand.pos_min[2],
    };

    dvec[0] -= ngrid_pad[0] * boxsize[0] / double(ngrid[0]);
    dvec[1] -= ngrid_pad[1] * boxsize[1] / double(ngrid[1]);
    dvec[2] -= ngrid_pad[2] * boxsize[2] / double(ngrid[2]);

    /// Shift mesh grid and recalculate extreme particle positions.
    particles_data.offset_coords(dvec);
    particles_rand.offset_coords(dvec);
  }

  /**
   * Calculate particle-based power spectrum normalisation.
   *
   * @param alpha Alpha ratio.
   * @returns norm Power spectrum normalisation constant.
   */
  double _calc_powspec_normalisation(double alpha=1.) {
    if (this->pdata == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    double vol_norm = 0.;
    for (int pid = 0; pid < this->ntotal; pid++) {
      vol_norm += this->pdata[pid].nz
        * this->pdata[pid].ws * this->pdata[pid].wc * this->pdata[pid].wc;
    }

    if (vol_norm == 0.) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Particle 'nz' values appear to be all zeros. "
          "Check the input catalogue contains valid 'nz' field.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    double norm = 1. / (alpha * vol_norm);

    return norm;
  }

  /**
   * Calculate particle-based power spectrum shot noise.
   *
   * @param alpha Alpha ratio.
   * @returns shotnoise Power spectrum shot noise.
   */
  double _calc_powspec_shotnoise(double alpha=1.) {
    if (this->pdata == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    double shotnoise = 0.;
    for (int pid = 0; pid < this->ntotal; pid++) {
      shotnoise += this->pdata[pid].ws * this->pdata[pid].ws
        * this->pdata[pid].wc * this->pdata[pid].wc;
    }

    shotnoise *= alpha * alpha;

    return shotnoise;
  }

 private:
  std::string source;
};

#endif  // TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_
