#ifndef TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_PARTICLES_HPP_INCLUDED_

#include <cstdio>
#include <string>
#include <vector>

#include "parameters.hpp"

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

  int nparticles;  ///< total number of particles

  double pos_min[3];  ///< minimum values of particle positions
  double pos_max[3];  ///< maximum values of particle positions

  /**
   * Initialise the particle container with default values.
   */
  ParticleCatalogue () {
    /// Set default values.
    this->particles = NULL;
    this->nparticles = 0;
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
  void initialise_particles(const int num) {
    /// Check the total number of particles.
    if (num <= 0) {
      printf("[Warning] :: Number of particles is negative.\n");
      return;
    }
    this->nparticles = num;

    /// Renew particle data.
    delete[] this->particles; this->particles = NULL;
    this->particles = new ParticleData[this->nparticles];

    /// Determine memory usage.
    bytesMem += double(this->nparticles) * sizeof(struct ParticleData)
      / 1024. / 1024. / 1024.;

    // /// Initialise particle data default values.
    // for (int id = 0; id < this->nparticles; id++) {
    //   this->particles[id].pos[0] = 0.;
    //   this->particles[id].pos[1] = 0.;
    //   this->particles[id].pos[2] = 0.;
    //   this->particles[id].w = 0.;
    // }
  }

  /**
   * Finalise particle data and information.
   */
  void finalise_particles() {
    /// Free particle usage.
    if (this->particles != NULL) {
      delete[] this->particles; this->particles = NULL;
      bytesMem -= double(this->nparticles) * sizeof(struct ParticleData)
        / 1024. / 1024. / 1024.;
    }
  }

  /**
   * Read in particle data from a file.
   *
   * @param particles_file Particle data file path.
   * @returns Exit status.
   */
  int read_particles_catalogue(std::string& particles_file) {
    /// Initialise the counter for the number of lines/particles.
    int num_lines = 0;

    /// Check and size up data from the file.
    std::ifstream fin;

    fin.open(particles_file.c_str(), std::ios::in);

    if (fin.fail()) {
      printf("[Error] :: Cannot open file '%s'.\n", particles_file.c_str());
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
    this->initialise_particles(num_lines);

    fin.open(particles_file.c_str(), std::ios::in);

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
    }

    fin.close();

    /// Calculate extreme data values.
    this->calc_min_and_max();

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
  int read_particles_catalogue(
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
    int nparticles = w.size();

    this->initialise_particles(nparticles);

    for (int pid = 0; pid < nparticles; pid++) {
      this->particles[pid].pos[0] = x[pid];
      this->particles[pid].pos[1] = y[pid];
      this->particles[pid].pos[2] = z[pid];
      this->particles[pid].w = w[pid];
    }

    /// Calculate extreme data values.
    this->calc_min_and_max();

    return 0;
  }

  /**
   * Calculate extreme data values.
   *
   * @returns Exit status.
   */
  int calc_min_and_max() {
    if (this->particles == NULL) {
      return -1;
    }

    /// Initialise minimum and maximum values with the first
    /// data entry/row.
    double min[3], max[3];

    min[0] = this->particles[0].pos[0]; max[0] = this->particles[0].pos[0];
    min[1] = this->particles[0].pos[1]; max[1] = this->particles[0].pos[1];
    min[2] = this->particles[0].pos[2]; max[2] = this->particles[0].pos[2];

    /// Find/update minimum and maximum values line-by-line.
    for (int id = 0; id < this->nparticles; id++) {
      if (min[0] > this->particles[id].pos[0]) {
        min[0] = this->particles[id].pos[0];
      }
      if (min[1] > this->particles[id].pos[1]) {
        min[1] = this->particles[id].pos[1];
      }
      if (min[2] > this->particles[id].pos[2]) {
        min[2] = this->particles[id].pos[2];
      }

      if (max[0] < this->particles[id].pos[0]) {
        max[0] = this->particles[id].pos[0];
      }
      if (max[1] < this->particles[id].pos[1]) {
        max[1] = this->particles[id].pos[1];
      }
      if (max[2] < this->particles[id].pos[2]) {
        max[2] = this->particles[id].pos[2];
      }
    }

    this->pos_min[0] = min[0]; this->pos_max[0] = max[0];
    this->pos_min[1] = min[1]; this->pos_max[1] = max[1];
    this->pos_min[2] = min[2]; this->pos_max[2] = max[2];

    return 0;
  }

  /**
   * Calculate the alpha ratio (of weighted number counts or
   * number densities).
   *
   * @param particles_data Data-source particle catalogue.
   * @param particles_rand Random-source particle catalogue.
   * @returns alpha Alpha ratio.
   */
  static double calc_alpha_ratio(
      ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand
    ) {
    double num_wgt_data = 0.;
    for (int id = 0; id < particles_data.nparticles; id++) {
      num_wgt_data += particles_data[id].w;
    }

    double num_wgt_rand = 0.;
    for (int id = 0; id < particles_rand.nparticles; id++) {
      num_wgt_rand += particles_rand[id].w;
    }

    double alpha = num_wgt_data / num_wgt_rand;

    return alpha;
  }

  /**
   * Calculate power spectrum normalisation.
   *
   * @param particles (Data-source) particle catalogue.
   * @param vol_norm Volume normalisation constant.
   * @returns norm Normalisation factor.
   */
  static double calc_norm_for_powspec(
      ParticleCatalogue& particles, double vol_norm
  ) {
    double num_wgt = 0.;
    for (int id = 0; id < particles.nparticles; id++) {
      num_wgt += particles[id].w;
    }

    // This is equivalent to 1/I_2 where I_2 = ∫d^3x \bar{n}(x)^2.
    double norm = vol_norm / num_wgt / num_wgt;

    return norm;
  }

  /**
   * Calculate bispectrum normalisation.
   *
   * @param particles (Data-source) particle catalogue.
   * @param vol_norm Volume normalisation constant.
   * @returns norm Normalisation factor.
   */
  static double calc_norm_for_bispec(
      ParticleCatalogue& particles, double vol_norm
  ) {
    double num_wgt = 0.;
    for (int id = 0; id < particles.nparticles; id++) {
      num_wgt += particles[id].w;
    }

    double norm = vol_norm / num_wgt / num_wgt;

    // This is NOT equivalent to 1/I_3 where I_3 = ∫d^3x \bar{n}(x)^3.
    norm *= vol_norm / num_wgt;

    return norm;
  }

  /**
   * Offset particle positions by a given vector to subtract.
   *
   * @param dpos (Subtractive) position offset vector .
   * @returns Exit status.
   */
  int offset_particles(const double* dpos) {
    if (this->particles == NULL) {
      if (currTask == 0) {
        printf("[Error] :: Particle data uninitialised.\n");
      }
      return -1;
    }

    for (int pid = 0; pid < this->nparticles; pid++) {
      this->particles[pid].pos[0] -= dpos[0];
      this->particles[pid].pos[1] -= dpos[1];
      this->particles[pid].pos[2] -= dpos[2];
    }

    return 0;
  }

  /**
   * Offset particle positions for FFT (by grid adjustment).
   *
   * @param particles_data (Data-source) particle catalogue.
   * @param particles_rand (Random-soruce) particle catalogue.
   * @param params Parameter set.
   * @param factor Offset grid adjustment factor (default is 3.).
   * @returns Exit status.
   */
  static int offset_particles_for_fft(
      ParticleCatalogue& particles_data,
      ParticleCatalogue& particles_rand,
      Parameters& params,
      double factor=3.  // CAVEAT: discretionary choice
    ) {
    /// Calculate adjustments needed.
    particles_data.calc_min_and_max();
    particles_rand.calc_min_and_max();

    double dpos[3] = {
      particles_rand.pos_min[0],
      particles_rand.pos_min[1],
      particles_rand.pos_min[2]
    };

    /// Re-adjust the grid.
    dpos[0] -= factor * params.boxsize[0] / double(params.nmesh[0]);
    dpos[1] -= factor * params.boxsize[1] / double(params.nmesh[1]);
    dpos[2] -= factor * params.boxsize[2] / double(params.nmesh[2]);

    particles_data.offset_particles(dpos);
    particles_rand.offset_particles(dpos);

    /// Recalculate extreme data values.
    particles_data.calc_min_and_max();
    particles_rand.calc_min_and_max();

    return 0;
  }

  /**
   * Offset particle positions for window function calculations (by
   * centring them inside the specified box).
   *
   * @param params Parameter set.
   * @returns Exit status.
   */
  int offset_particles_for_window(Parameters& params) {
    /// Calculate grid centre.
    double xmid = this->pos_min[0]
      + (this->pos_max[0] - this->pos_min[0]) / 2.;
    double ymid = this->pos_min[1]
      + (this->pos_max[1] - this->pos_min[1]) / 2.;
    double zmid = this->pos_min[2]
      + (this->pos_max[2] - this->pos_min[2]) / 2.;

    /// Calculate adjustments needed.
    double dx[3] = {0., 0., 0.};
    dx[0] = params.boxsize[0]/2. - xmid;
    dx[1] = params.boxsize[1]/2. - ymid;
    dx[2] = params.boxsize[2]/2. - zmid;

    /// Centre the grid.
    for (int pid = 0; pid < this->nparticles; pid++) {
      for (int axis = 0; axis < 3; axis++) {
        this->particles[pid].pos[axis] += dx[axis];
      }
    }

    /// Recalculate extreme data values.
    this->calc_min_and_max();

    return 0;
  }

  /**
   * Offset particle positions for periodic boundary conditions.
   *
   * @param params Parameter set.
   * @returns Exit status.
   */
  int offset_particles_for_periodicity(Parameters& params) {
    /// Re-adjust the grid.
    for (int pid = 0; pid < this->nparticles; pid++) {
      for (int axis = 0; axis < 3; axis++) {
        if (this->particles[pid].pos[axis] >= params.boxsize[axis]) {
          this->particles[pid].pos[axis] -= params.boxsize[axis];
        } else if (this->particles[pid].pos[axis] < 0.) {
          this->particles[pid].pos[axis] += params.boxsize[axis];
        }
      }
    }

    /// Recalculate extreme data values.
    this->calc_min_and_max();

    return 0;
  }
};

#endif
