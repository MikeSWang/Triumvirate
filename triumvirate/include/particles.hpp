/**
 * @file particles.hpp
 * @brief Particle containers with I/O methods and operations, as well as
 *        line-of-sight vectors.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "common.hpp"
#include "tools.hpp"

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
    double w;  ///< particle weight
  }* particles;  ///< particle data

  int ntotal;  ///< total number of particles

  double wtotal;  ///< total weight of particles

  double pos_min[3];  ///< minimum values of particle positions
  double pos_max[3];  ///< maximum values of particle positions

  /**
   * Initialise the particle container with default values.
   */
  ParticleCatalogue () {
    /// Set default values.
    this->particles = NULL;
    this->ntotal = 0;
    this->wtotal = 0.;
    this->pos_min[0] = 0.; this->pos_min[1] = 0.; this->pos_min[2] = 0.;
    this->pos_max[0] = 0.; this->pos_max[1] = 0.; ; this->pos_max[2] = 0.;
  }

  /**
   * Destruct the particle container.
   */
  ~ParticleCatalogue() {
    finalise_particles();
  }

  /**
   * Return individual particle information.
   *
   * @param pid Particle index.
   * @returns Individual particle data.
   */
  ParticleData& operator[](int pid) {
    return this->particles[pid];
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
        "[Warning] (+%s) Number of particles is negative.\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
      );
      return;
    }
    this->ntotal = num;

    /// Renew particle data.
    delete[] this->particles; this->particles = NULL;
    this->particles = new ParticleData[this->ntotal];

    /// Determine memory usage.
    bytesMem += double(this->ntotal) * sizeof(struct ParticleData)
      / 1024. / 1024. / 1024.;
  }

  /**
   * Finalise particle data and information.
   */
  void finalise_particles() {
    /// Free particle usage.
    if (this->particles != NULL) {
      delete[] this->particles; this->particles = NULL;
      bytesMem -= double(this->ntotal) * sizeof(struct ParticleData)
        / 1024. / 1024. / 1024.;
    }
  }

  /**
   * Read in particle data from a file.
   *
   * @param particles_file Particle data file path.
   * @returns Exit status.
   */
  int read_particle_data_from_file(std::string& particles_file) {
    /// Initialise the counter for the number of lines/particles.
    int num_lines = 0;

    /// Check and size up data from the file.
    std::ifstream fin;

    fin.open(particles_file.c_str(), std::ios::in);

    if (fin.fail()) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[Error] (+%s) Cannot open file '%s'.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          particles_file.c_str()
        );
      }
      fin.close();
      return -1;
    }

    std::string str_line;
    double x, y, z, w;
    while (getline(fin, str_line)) {
      /// Check if the line conforms to the expected format.
      if (sscanf(str_line.c_str(), "%lf %lf %lf %lf", &x, &y, &z, &w) != 4) {
        continue;
      }

      /// Count the line as a particle.
      num_lines++;
    }

    fin.close();

    /// Fill in particle data.
    this->_initialise_particles(num_lines);

    fin.open(particles_file.c_str(), std::ios::in);

    num_lines = 0;  // reset
    while (getline(fin, str_line)) {
      /// Check if the line conforms to the expected format.
      if (sscanf(str_line.c_str(), "%lf %lf %lf %lf", &x, &y, &z, &w) != 4) {
        continue;
      }

      /// Add the current line as a particle.
      this->particles[num_lines].pos[0] = x;
      this->particles[num_lines].pos[1] = y;
      this->particles[num_lines].pos[2] = z;
      this->particles[num_lines].w = w;

      num_lines++;
    }

    fin.close();

    /// Calculate weight sum.
    this->_calc_weighted_total();

    /// Calculate extreme particle positions.
    this->_calc_pos_min_and_max();

    return 0;
  }

  /**
   * Read in particle data from external data.
   *
   * @param x Particle x-positions.
   * @param y Particle y-positions.
   * @param z Particle z-positions.
   * @param w Particle weights.
   * @returns Exit status.
   */
  int read_particle_data(
    std::vector<double> x, std::vector<double> y, std::vector<double> z,
    std::vector<double> w
  ) {
    /// Check array sizes.
    if (
      !(x.size() == y.size() && y.size() == z.size() && z.size() == w.size())
    ) {
      return -1;
    }

    /// Fill in particle data.
    int ntotal = w.size();

    this->_initialise_particles(ntotal);

    for (int pid = 0; pid < ntotal; pid++) {
      this->particles[pid].pos[0] = x[pid];
      this->particles[pid].pos[1] = y[pid];
      this->particles[pid].pos[2] = z[pid];
      this->particles[pid].w = w[pid];
    }

    /// Calculate weight sum.
    this->_calc_weighted_total();

    /// Calculate extreme particle positions.
    this->_calc_pos_min_and_max();

    return 0;
  }

  /**
   * Calculate the weighted number of particles, i.e. the total weights.
   */
  void _calc_weighted_total() {
    if (this->particles == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[Error] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    double wtotal = 0.;
    for (int pid = 0; pid < this->ntotal; pid++) {
      wtotal += this->particles[pid].w;
    }

    this->wtotal = wtotal;
  }

  /**
   * Calculate extreme particle positions.
   */
  void _calc_pos_min_and_max() {
    if (this->particles == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[Error] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    /// Initialise minimum and maximum values with the first
    /// data entry/row.
    double min[3], max[3];
    for (int axis = 0; axis < 3; axis++) {
      min[axis] = this->particles[0].pos[axis];
      max[axis] = this->particles[0].pos[axis];
    }

    /// Update minimum and maximum values line-by-line.
    for (int id = 0; id < this->ntotal; id++) {
      for (int axis = 0; axis < 3; axis++) {
        if (min[axis] > this->particles[id].pos[axis]) {
          min[axis] = this->particles[id].pos[axis];
        }
        if (max[axis] < this->particles[id].pos[axis]) {
          max[axis] = this->particles[id].pos[axis];
        }
      }
    }

    for (int axis = 0; axis < 3; axis++) {
      this->pos_min[axis] = min[axis];
      this->pos_max[axis] = max[axis];
    }
  }

  /**
   * Offset particle positions by a given vector (as the new origin).
   *
   * @param dpos (Subtractive) offset position vector.
   */
  void offset_coords(const double dpos[3]) {
    if (this->particles == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[Error] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    for (int pid = 0; pid < this->ntotal; pid++) {
      for (int axis = 0; axis < 3; axis++) {
        this->particles[pid].pos[axis] -= dpos[axis];
      }
    }
  }

  /**
   * Offset particle positions by centring the catalogue inside the
   * specified box.
   *
   * @param boxsize Boxsize in each dimension.
   */
  void offset_coords_for_centring(const double boxsize[3]) {
    if (this->particles == NULL) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[Error] (+%s) Particle data are uninitialised.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }

    /// Calculate adjustments needed.
    double xmid = this->pos_min[0]
      + (this->pos_max[0] - this->pos_min[0]) / 2.;
    double ymid = this->pos_min[1]
      + (this->pos_max[1] - this->pos_min[1]) / 2.;
    double zmid = this->pos_min[2]
      + (this->pos_max[2] - this->pos_min[2]) / 2.;

    double dvec[3] = {
      boxsize[0]/2. - xmid,
      boxsize[1]/2. - ymid,
      boxsize[2]/2. - zmid,
    };

    /// Centre the catalogue in box.
    for (int pid = 0; pid < this->ntotal; pid++) {
      for (int axis = 0; axis < 3; axis++) {
        this->particles[pid].pos[axis] += dvec[axis];
      }
    }

    /// Recalculate extreme particle positions.
    this->_calc_pos_min_and_max();
  }

  /**
   * Offset particle positions for periodic boundary conditions.
   *
   * @param boxsize Periodic boxsize in each dimension.
   */
  void offset_coords_for_periodicity(const double boxsize[3]) {
    /// Re-adjust the grid.
    for (int pid = 0; pid < this->ntotal; pid++) {
      for (int axis = 0; axis < 3; axis++) {
        if (
          this->particles[pid].pos[axis] >= boxsize[axis]
        ) {
          this->particles[pid].pos[axis] -= boxsize[axis];
        } else if (
          this->particles[pid].pos[axis] < 0.
        ) {
          this->particles[pid].pos[axis] += boxsize[axis];
        }
      }
    }

    /// Recalculate extreme particle positions.
    this->_calc_pos_min_and_max();
  }

  /**
   * Align a pair of catalogues by offsetting particle positions
   * for FFT (though mesh grid shift).
   *
   * @param particles_data (Data-source) particle catalogue.
   * @param particles_rand (Random-source) particle catalogue.
   * @param boxsize Box size in each dimension.
   * @param ngrid Grid number in each dimension.
   * @param factor Offset grid shift factor (default is 3.).
   * @returns Exit status.
   */
  static void align_catalogues_for_fft(
    ParticleCatalogue& particles_data,
    ParticleCatalogue& particles_rand,
    const double boxsize[3],
    const int ngrid[3],
    double shift_factor=3.  // CAVEAT: discretionary choice
  ) {
    /// Calculate adjustments needed.
    particles_data._calc_pos_min_and_max();
    particles_rand._calc_pos_min_and_max();

    double dpos[3] = {
      particles_rand.pos_min[0],
      particles_rand.pos_min[1],
      particles_rand.pos_min[2],
    };

    dpos[0] -= shift_factor * boxsize[0] / double(ngrid[0]);
    dpos[1] -= shift_factor * boxsize[1] / double(ngrid[1]);
    dpos[2] -= shift_factor * boxsize[2] / double(ngrid[2]);

    /// Shift mesh grid and recalculate extreme particle positions.
    particles_data.offset_coords(dpos);
    particles_rand.offset_coords(dpos);

    particles_data._calc_pos_min_and_max();
    particles_rand._calc_pos_min_and_max();
  }
};

#endif  // TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_
