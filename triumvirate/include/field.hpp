#ifndef TRIUMVIRATE_INCLUDE_FIELD_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_FIELD_HPP_INCLUDED_

#include <complex>
#include <cstring>

#include <fftw3.h>

#include "common.hpp"
#include "parameters.hpp"

/**
 * Density field object instantiated from a particle container.
 *
 */
template <class ParticleContainer>
class DensityField {
 public:
  fftw_complex* field;  ///> gridded complex field

  /**
   * Construct density field.
   *
   * @param params Parameter set.
   */
  DensityField(ParameterSet& params) {
    /// Attach external parameters.
    this->params = params;

    /// Initialise complex field and increase allocated memory.
    this->field = fftw_alloc_complex(this->params.nmesh);
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] = 0.;
      this->field[gid][1] = 0.;
    }

    bytesMem += double(this->params.nmesh) * sizeof(fftw_complex)
      / 1024. / 1024. / 1024.;
  }

  /**
   * Destruct density field.
   */
  ~DensityField() {
    finalise_density_field();
  }

  /**
   * Return individual grid value of the field.
   *
   * @param gid Grid ID/index.
   * @returns Field value.
   */
  const fftw_complex& operator[](int gid) {
    return this->field[gid];
  }

  /**
   * Finalise density field data.
   */
  void finalise_density_field() {
    /// Free memory usage.
    if (this->field != NULL) {
      fftw_free(this->field); this->field = NULL;
      bytesMem -= double(this->params.nmesh) * sizeof(fftw_complex)
        / 1024. / 1024. / 1024.;
    }
  }

  /**
   * Assign weighted density field to a grid by interpolation scheme.
   *
   * @param particles Particle container.
   * @param weights Weight field.
   * @returns Exit status.
   */
  void assign_weighted_field_to_grid(
    ParticleContainer& particles,
    fftw_complex* weights
  ) {
    if (
      this->params.assignment == "NGP"
    ) {
      this->assign_weighted_field_to_grid_ngp(particles, weights);
    } else if (
      this->params.assignment == "CIC"
    ) {
      this->assign_weighted_field_to_grid_cic(particles, weights);
    } else if (
      this->params.assignment == "TSC"
    ) {
      this->assign_weighted_field_to_grid_tsc(particles, weights);
    } else {
      if (currTask == 0) {
        printf(
          "[Error] :: Unsupported grid assignment scheme: '%s'.",
          this->params.assignment.c_str()
        );
      }
      exit(1);
    }
  }

  /**
   * Calculate weighted density field fluctuation(s) weighted by the
   * reduced spherical harmonics from data and random sources,
   *
   * f@[
   *   {\delta n_{LM}}(\vec{x}) = \left(
   *     \sum_{i \in \mathrm{data}} - \alpha \sum_{i \in \mathrm{rand}}
   *   \right) \delta_{\mathrm{D}}(\vec{x} - \vec{x}_i)
   *   y_{LM}^*(\hat{\vec{x}}_i) w(\vec{x}_i) \,.
   * f@]
   *
   * See eq. (34) in Sugiyama et al. (2018)
   * [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param particles_data (Data-source) particle container.
   * @param particles_rand (Random-source) particle container.
   * @param los_data (Data-source) particle lines of sight.
   * @param los_rand (Random-source) particle lines of sight.
   * @param alpha Alpha ratio.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void calc_ylm_wgtd_fluctuation(
    ParticleContainer& particles_data, ParticleContainer& particles_rand,
    LineOfSight* los_data, LineOfSight* los_rand,
    double alpha,
    int ell, int m
  ) {
    /// Initialise the random density field and the kernel weight field.
    DensityField<ParticleCatalogue> density_rand(this->params);
    fftw_complex* weight_kern = NULL;

    /// Calculate the transformed data-source field.
    weight_kern = fftw_alloc_complex(particles_data.ntotal);
    for (int pid = 0; pid < particles_data.ntotal; pid++) {
      double los_[3] = {
        los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
          ell, m, los_
        );

      weight_kern[pid][0] = ylm.real() * particles_data[pid].w;
      weight_kern[pid][1] = ylm.imag() * particles_data[pid].w;
    }

    this->assign_weighted_field_to_grid(particles_data, weight_kern);

    fftw_free(weight_kern); weight_kern = NULL;

    /// Calculate the random-source transformed weighted field.
    weight_kern = fftw_alloc_complex(particles_rand.ntotal);
    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      double los_[3] = {
        los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
          ell, m, los_
        );

      weight_kern[pid][0] = ylm.real() * particles_rand[pid].w;
      weight_kern[pid][1] = ylm.imag() * particles_rand[pid].w;
    }

    density_rand.assign_weighted_field_to_grid(particles_rand, weight_kern);

    fftw_free(weight_kern); weight_kern = NULL;

    /// Subtract to compute fluctuations, i.e. δn_LM.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] -= alpha * density_rand.field[gid][0];
      this->field[gid][1] -= alpha * density_rand.field[gid][1];
    }
  }

  /**
   * Calculate the weighted density field weighted by the
   * reduced spherical harmonics from a (random) source,
   *
   * f@[
   *   \bar{n}_{LM}(\vec{x}) = \alpha \sum_{i \in \mathrm{rand}}
   *     \delta_{\mathrm{D}}(\vec{x} - \vec{x}_i)
   *     y_{LM}^*(\hat{\vec{x}}_i) w(\vec{x}_i) \,.
   * f@]
   *
   * This is analogous to DensityField::calc_ylm_wgtd_fluctuation
   * but for the (mean-)density field represented by a (random) particle
   * container.
   *
   * @param particles_rand (Random-source) particle container.
   * @param los_rand (Random-source) particle lines of sight.
   * @param alpha Alpha ratio.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void calc_ylm_wgtd_density(
    ParticleContainer& particles_rand,
    LineOfSight* los_rand,
    double alpha,
    int ell, int m
  ) {
    /// Initialise the kernel weight field.
    fftw_complex* weight_kern = NULL;

    /// Calculate the transformed weighted field.
    weight_kern = fftw_alloc_complex(particles_rand.ntotal);
    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      double los_[3] = {
        los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
          ell, m, los_
        );

      weight_kern[pid][0] = ylm.real() * particles_rand[pid].w;
      weight_kern[pid][1] = ylm.imag() * particles_rand[pid].w;
    }

    this->assign_weighted_field_to_grid(particles_rand, weight_kern);

    fftw_free(weight_kern); weight_kern = NULL;

    /// Apply mean-density matching normalisation (i.e. alpha ratio)
    /// to compute \bar{n}_LM.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] *= alpha;
      this->field[gid][1] *= alpha;
    }
  }

  /**
   * Calculate two-point self-contribution component to shot noise,
   * weighted by the reduced spherical harmonics, from data and random
   * sources,
   *
   * f@[
   *   N_{LM}(\vec{x}) = \left(
   *     \sum_{i \in \mathrm{data}} + \alpha^2 \sum_{i \in \mathrm{rand}}
   *   \right) \delta_{\mathrm{D}}(\vec{x} - \vec{x}_i)
   *   y_{LM}(\vec{x}) w(\vec{x}_i)^2 \,.
   * f@]
   *
   * See eq. (46) in Sugiyama et al. (2018)
   * [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param particles_data (Data-source) particle container.
   * @param particles_rand (Random-source) particle container.
   * @param los_data (Data-source) particle lines of sight.
   * @param los_rand (Random-source) particle lines of sight.
   * @param alpha Alpha ratio.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void calc_ylm_wgtd_2pt_self_contrib_for_shotnoise(
    ParticleContainer& particles_data, ParticleContainer& particles_rand,
    LineOfSight* los_data, LineOfSight* los_rand,
    double alpha,
    int ell, int m
  ) {
    /// Initialise the random density field and the kernel weight field.
    DensityField<ParticleCatalogue> density_rand(this->params);
    fftw_complex* weight_kern = NULL;

    /// Calculate the transformed data-source field.
    weight_kern = fftw_alloc_complex(particles_data.ntotal);
    for (int pid = 0; pid < particles_data.ntotal; pid++) {
      double los_[3] = {
        los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
          ell, m, los_
        );

      ylm = std::conj(ylm);  // NOTE: cojugation is essential

      weight_kern[pid][0] = ylm.real() * pow(particles_data[pid].w, 2);
      weight_kern[pid][1] = ylm.imag() * pow(particles_data[pid].w, 2);
    }

    this->assign_weighted_field_to_grid(particles_data, weight_kern);

    fftw_free(weight_kern); weight_kern = NULL;

    /// Calculate the random-source transformed weighted field.
    weight_kern = fftw_alloc_complex(particles_rand.ntotal);
    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      double los_[3] = {
        los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
          ell, m, los_
        );

      ylm = std::conj(ylm);  // NOTE: cojugation is essential

      weight_kern[pid][0] = ylm.real() * pow(particles_rand[pid].w, 2);
      weight_kern[pid][1] = ylm.imag() * pow(particles_rand[pid].w, 2);
    }

    density_rand.assign_weighted_field_to_grid(particles_rand, weight_kern);

    fftw_free(weight_kern); weight_kern = NULL;

    /// Add to compute total shot noise contribution, i.e. N_LM.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] += alpha * alpha * density_rand.field[gid][0];
      this->field[gid][1] += alpha * alpha * density_rand.field[gid][1];
    }
  }

  /**
   * Calculate two-point self-contribution component to shot noise,
   * weighted by the reduced spherical harmonics, from a (random) source.
   *
   * This is analogous to
   * DensityField::calc_ylm_wgtd_2pt_self_contrib_for_shotnoise(),
   * but used for calculating shot noise in window functions.
   *
   * This is analogous to
   * DensityField::calc_ylm_wgtd_2pt_self_contrib_for_shotnoise
   * but for the (mean-)density field represented by a (random) particle
   * container only.
   *
   * @param particles_rand (Random-source) particle container.
   * @param los_rand (Random-source) particle lines of sight.
   * @param alpha Alpha ratio.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void calc_ylm_wgtd_2pt_self_contrib_for_shotnoise(
    ParticleContainer& particles_rand,
    LineOfSight* los_rand,
    double alpha,
    int ell, int m
  ) {
    /// Initialise the kernel weight field.
    fftw_complex* weight_kern = NULL;

    /// Calculate the transformed weighted field.
    weight_kern = fftw_alloc_complex(particles_rand.ntotal);
    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      double los_[3] = {
        los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
          ell, m, los_
        );

      ylm = std::conj(ylm);  // NOTE: cojugation is essential

      weight_kern[pid][0] = ylm.real() * pow(particles_rand[pid].w, 2);
      weight_kern[pid][1] = ylm.imag() * pow(particles_rand[pid].w, 2);
    }

    this->assign_weighted_field_to_grid(particles_rand, weight_kern);

    fftw_free(weight_kern); weight_kern = NULL;

    /// Apply mean-density matching normalisation (i.e. alpha ratio)
    /// to compute N_LM.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] *= alpha * alpha;
      this->field[gid][1] *= alpha * alpha;
    }
  }

  /**
   * Calculate unweighted density field fluctuation(s).
   *
   * This is typically used in the case of reconstruction.
   *
   * @param particles_data (Data-source) particle container.
   * @param particles_rand (Random-source) particle container.
   * @param alpha Alpha ratio.
   */
  void calc_unweighted_fluctuation(
    ParticleContainer& particles_data, ParticleContainer& particles_rand,
    double alpha
  ) {
    /// Initialise the random field and the weight field.
    DensityField<ParticleCatalogue> density_rand(this->params);
    fftw_complex* weight = NULL;

    /// Calculate the data-source field.
    weight = fftw_alloc_complex(particles_data.ntotal);
    for (int pid = 0; pid < particles_data.ntotal; pid++) {
      weight[pid][0] = 1.;
      weight[pid][1] = 0.;
    }

    this->assign_weighted_field_to_grid(particles_data, weight);

    fftw_free(weight); weight = NULL;

    /// Calculate the random-source field.
    weight = fftw_alloc_complex(particles_rand.ntotal);
    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      weight[pid][0] = 1.;
      weight[pid][1] = 0.;
    }

    density_rand.assign_weighted_field_to_grid(particles_rand, weight);

    fftw_free(weight); weight = NULL;

    /// Subtract to compute fluctuations, i.e. δn.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] -= alpha * density_rand.field[gid][0];
      this->field[gid][1] -= alpha * density_rand.field[gid][1];
    }
  }

  /**
   * Calculate unweighted density field fluctuation(s) within a source.
   *
   * This is typically used for simulations in a box.
   *
   * @param particles_data (Data-source) particle container.
   * @param volume Box volume.
   */
  void calc_unweighted_fluctuation_insitu(
    ParticleContainer& particles_data,
    double volume
  ) {
    /// Initialise the weight field.
    fftw_complex* weight = NULL;

    /// Calculate the unit-weighted field.
    weight = fftw_alloc_complex(particles_data.ntotal);
    for (int pid = 0; pid < particles_data.ntotal; pid++) {
      weight[pid][0] = 1.;
      weight[pid][1] = 0.;
    }

    this->assign_weighted_field_to_grid(particles_data, weight);

    fftw_free(weight); weight = NULL;

    /// Subtract the global mean density to compute fluctuations, i.e. δn.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] -= double(particles_data.ntotal) / params.volume;
      this->field[gid][1] -= 0.;
    }
  }

  /**
   * Calculate the unweighted density field.
   *
   * This is typically used for simulations in a periodic box and
   * for bispectrum calculations.
   *
   * @param particles_data (Data-source) particle container.
   */
  void calc_unweighted_density(ParticleContainer& particles_data) {
    /// Initialise the weight field.
    fftw_complex* weight = NULL;

    /// Calculate the unit-weighted field, i.e. n.
    weight = fftw_alloc_complex(particles_data.ntotal);
    for (int pid = 0; pid < particles_data.ntotal; pid++) {
      weight[pid][0] = 1.;
      weight[pid][1] = 0.;
    }

    this->assign_weighted_field_to_grid(particles_data, weight);

    fftw_free(weight); weight = NULL;
  }

  /**
   * Calculate Fourier transform of the field.
   */
  void fourier_transform() {
    /// Apply FFT volume normalisation, where ∫d^3x corresponds to
    /// dV Σ_i, dV =: `vol_cell`.
    double vol_cell = this->params.volume / double(this->params.nmesh);

    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] *= vol_cell;
      this->field[gid][1] *= vol_cell;
    }

    /// Perform FFT.
    fftw_plan transform = fftw_plan_dft_3d(
      this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
      this->field, this->field,
      FFTW_FORWARD, FFTW_ESTIMATE
    );

    fftw_execute(transform);
    fftw_destroy_plan(transform);
  }

  /**
   * Calculate inverse Fourier transform of the (FFT-transformed) field.
   */
  void inv_fourier_transform() {
    /// Apply inverse FFT volume normalisation, where ∫d^3k / (2\pi)^3
    /// corresponds to (1/V) Σ_i, V =: `vol`.
    double vol = this->params.volume;

    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] /= vol;
      this->field[gid][1] /= vol;
    }

    /// Perform inverse FFT.
    fftw_plan inv_transform = fftw_plan_dft_3d(
      this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
      this->field, this->field,
      FFTW_BACKWARD, FFTW_ESTIMATE
    );

    fftw_execute(inv_transform);
    fftw_destroy_plan(inv_transform);
  }

  /**
   * Calculate inverse Fourier transform of the (FFT-transformed)
   * (density fluctuation) field weighted by the reduce spherical
   * harmonics in a wavenumber bin,
   *
   * This is used to compute
   *
   * f@[
   *   F_{LM}(\vec{x}; k) = \frac{(2\pi)^3}{4\pi k^2}
   *     \int \frac{\mathrm{d}^3\,\vec{k}'}{(2\pi)^3}
   *       \mathrm{e}^{\mathrm{i} \vec{k}' \cdot \vec{x}}
   *       \delta_\mathrm{D}(k' - k) y_{LM}(\hat{\vec{k}})
   *       \delta n_{LM}(\vec{k}) \,,
   * f@]
   *
   * where f@$ \delta n_{LM} f@$ evaluated using FFT is corrected for
   * the interpolation window.
   *
   * See eq. (42) in Sugiyama et al. (2018)
   * [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param density (FFT-transformed) density (fluctuation) field.
   * @param ylm Reduced spherical harmonic on a grid.
   * @param k_in Wavenumber bin wavenumber.
   * @param dk_in Wavenumber bin width.
   */
  void inv_fourier_transform_for_ylm_wgtd_field_in_wavenum_bin(
    DensityField& density, std::complex<double>* ylm,
    double k_in, double dk_in
  ) {
    /// Reset field with zero values.
    for (int i = 0; i < this->params.nmesh; i++) {
      this->field[i][0] = 0.;
      this->field[i][1] = 0.;
    }

    /// Set up mode sampling within the wavenumber bin.
    double k_lower = (k_in > dk_in/2) ? (k_in - dk_in/2) : 0.;
    double k_upper = k_in + dk_in/2;

    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    double kv[3];
    int nmode = 0;
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          /// Calculate the wave-vector/-number representing the grid.
          kv[0] = (i < this->params.ngrid[0]/2) ?
            i * dk[0] : (i - this->params.ngrid[0]) * dk[0];
          kv[1] = (j < this->params.ngrid[1]/2) ?
            j * dk[1] : (j - this->params.ngrid[1]) * dk[1];
          kv[2] = (k < this->params.ngrid[2]/2) ?
            k * dk[2] : (k - this->params.ngrid[2]) * dk[2];

          double k_ = sqrt(kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2]);

          /// Determine grid contribution to the bin.
          if (k_ > k_lower && k_ <= k_upper) {
            std::complex<double> delta_n(
              density[idx_grid][0], density[idx_grid][1]
            );

            /// Apply assignment compensation.
            double win = this->calc_assignment_window_in_fourier(kv);
            delta_n /= win;

            /// Weight the field.
            this->field[idx_grid][0] = (ylm[idx_grid] * delta_n).real();
            this->field[idx_grid][1] = (ylm[idx_grid] * delta_n).imag();

            nmode++;
          } else {
            this->field[idx_grid][0] = 0.;
            this->field[idx_grid][1] = 0.;
          }
        }
      }
    }

    fftw_plan inv_transform = fftw_plan_dft_3d(
      this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
      this->field, this->field,
      FFTW_BACKWARD, FFTW_ESTIMATE
    );

    fftw_execute(inv_transform);
    fftw_destroy_plan(inv_transform);

    /// Apply the 4π-factor by equivalently mode averaging.
    for (int i = 0; i < this->params.nmesh; i++) {
      this->field[i][0] /= double(nmode);
      this->field[i][1] /= double(nmode);
    }
  }

  /**
   * Calculate the inverse Fourier transform of a transformed field
   * for three-point correlation functions.
   *
   * @param density FFT-transformed density field.
   * @param ylm Reduced spherical harmonic on a grid.
   * @param sjl Spherical Bessel function interpolator.
   * @param r_in Separation magnitude.
   * @returns Exit status.
   */
  void inv_fourier_transform_for_sjl_ylm_wgtd_field(
    DensityField& density,
    std::complex<double>* ylm,
    SphericalBesselCalculator& sjl,
    double r_in
  ) {
    /// Reset field with zero values.
    for (int i = 0; i < this->params.nmesh; i++) {
      this->field[i][0] = 0.;
      this->field[i][1] = 0.;
    }

    /// Set inverse FFT volume normalisation, where ∫d^3k / (2\pi)^3
    /// corresponds to (1/V) Σ_i, V =: `vol`.
    double vol = this->params.volume;

    /// Set up grid modes.
    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    double kv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          /// Calculate the wave-vector/-number representing the grid.
          kv[0] = (i < this->params.ngrid[0]/2) ?
            i * dk[0] : (i - this->params.ngrid[0]) * dk[0];
          kv[1] = (j < this->params.ngrid[1]/2) ?
            j * dk[1] : (j - this->params.ngrid[1]) * dk[1];
          kv[2] = (k < this->params.ngrid[2]/2) ?
            k * dk[2] : (k - this->params.ngrid[2]) * dk[2];

          double k_ = sqrt(kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2]);

          /// Apply assignment compensation.
          std::complex<double> den(
            density[idx_grid][0], density[idx_grid][1]
          );

          double win = this->calc_assignment_window_in_fourier(kv);
          den /= win;

          /// Weight the field.  Apply inverse FFT volume normalisation.
          this->field[idx_grid][0] =
            sjl.eval(k_ * r_in) * (ylm[idx_grid] * den).real() / vol;
          this->field[idx_grid][1] =
            sjl.eval(k_ * r_in) * (ylm[idx_grid] * den).imag() / vol;
        }
      }
    }

    fftw_plan inv_transform = fftw_plan_dft_3d(
      this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
      this->field, this->field,
      FFTW_BACKWARD, FFTW_ESTIMATE
    );

    fftw_execute(inv_transform);
    fftw_destroy_plan(inv_transform);
  }

  /**
   * Apply compensation needed in Fourier transform for
   * assignment schemes.
   */
  void apply_assignment_compensation() {
    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    double kv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          kv[0] = (i < this->params.ngrid[0]/2) ?
            i * dk[0] : (i - this->params.ngrid[0]) * dk[0];
          kv[1] = (j < this->params.ngrid[1]/2) ?
            j * dk[1] : (j - this->params.ngrid[1]) * dk[1];
          kv[2] = (k < this->params.ngrid[2]/2) ?
            k * dk[2] : (k - this->params.ngrid[2]) * dk[2];

          double win = this->calc_assignment_window_in_fourier(kv);

          this->field[idx_grid][0] /= win;
          this->field[idx_grid][1] /= win;
        }
      }
    }
  }

  /**
   * Apply separation power-law weight for wide-angle corrections.
   *
   * This weights the field by f@$ r^{- i - j} f@$, where f@$i, jf@$
   * are wide-angle correction orders.
   */
  void apply_power_law_weight_for_wide_angle() {
    const double eps = 1.e-10;

    double dr[3];
    dr[0] = this->params.boxsize[0] / double(this->params.ngrid[0]);
    dr[1] = this->params.boxsize[1] / double(this->params.ngrid[1]);
    dr[2] = this->params.boxsize[2] / double(this->params.ngrid[2]);

    double rv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;
          rv[0] = (i < this->params.ngrid[0]/2) ?
            i * dr[0] : (i - this->params.ngrid[0]) * dr[0];
          rv[1] = (j < this->params.ngrid[1]/2) ?
            j * dr[1] : (j - this->params.ngrid[1]) * dr[1];
          rv[2] = (k < this->params.ngrid[2]/2) ?
            k * dr[2] : (k - this->params.ngrid[2]) * dr[2];

          double r_ = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);

          if (r_ < eps) {
            this->field[idx_grid][0] *= 0.;
            this->field[idx_grid][1] *= 0.;
          } else {
            this->field[idx_grid][0] *=
              pow(r_, - this->params.i_wa - this->params.j_wa);
            this->field[idx_grid][1] *=
              pow(r_, - this->params.i_wa - this->params.j_wa);
          }
        }
      }
    }
  }

  /**
   * Calculate (effective) volume normalisation for two-point correlators,
   *
   * f@[
   *   I_2 = \int \mathrm{d}^3\,\vec{x} \bar{n}(\vec{x})^2 \,.   *
   * f@]
   *
   * @param particles_rand (Random-source) particle container.
   * @returns vol_norm Survey volume normalisation.
   */
  double calc_volume_normalisation(ParticleContainer& particles_rand) {
    /// Initialise the weight field.
    fftw_complex* weight = NULL;

    weight = fftw_alloc_complex(particles_rand.ntotal);
    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      weight[pid][0] = particles_rand[pid].w;
      weight[pid][1] = 0.;
    }

    /// Calculate the weighted field.
    this->assign_weighted_field_to_grid(particles_rand, weight);

    fftw_free(weight); weight = NULL;

    /// Compute volume normalisation integral, where ∫d^3x corresponds
    /// to dV Σ_i, dV =: `vol_cell`.
    double vol_cell = this->params.volume / double(this->params.nmesh);

    double vol_factor = 0.;
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      vol_factor += vol_cell * this->field[gid][0] * this->field[gid][0];
    }

    /// Eliminate dependence on total particle number.
    /// CAVEAT: `double` needed to prevent int overflow.
    double vol_norm =
      double(particles_rand.ntotal) * double(particles_rand.ntotal)
      / vol_factor;

    return vol_norm;
  }

 private:
  ParameterSet params;

  /**
   * Assign weighted density field to a mesh grid by the
   * nearest-grid-point (NGP) scheme.
   *
   * @param particles Particle container.
   * @param weight Particle weight.
   */
  void assign_weighted_field_to_grid_ngp(
    ParticleContainer& particles,
    fftw_complex* weight
  ) {
    /// Set interpolation order, i.e. number of mesh grid, per dimension,
    /// to which a single particle is assigned
    int order = 1;

    /// Reset field with zero values.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] = 0.;
      this->field[gid][1] = 0.;
    }

    /// Here the density field is given by Σ_i w_i δ_D(x - x_i),
    /// where δ_D corresponds to δ_K / dV, dV =: `vol_cell`.
    double vol_cell = this->params.volume / double(this->params.nmesh);
    double inv_vol_cell = 1. / vol_cell;

    /// Perform assignment.
    for (int pid = 0; pid < particles.ntotal; pid++) {
      int ijk[order][3];  // coordinates of covered mesh grids
      double win[order][3];  // interpolation window
      for (int axis = 0; axis < 3; axis++) {
        double loc_grid = this->params.ngrid[axis]
          * particles[pid].pos[axis] / this->params.boxsize[axis];

        /// Set only 0th element as `order == 1`.
        ijk[0][axis] = int(loc_grid + 0.5);
        win[0][axis] = 1.;
      }

      for (int iloc = 0; iloc < order; iloc++) {
        for (int jloc = 0; jloc < order; jloc++) {
          for (int kloc = 0; kloc < order; kloc++) {
            long long idx_grid = (
              ijk[iloc][0] * this->params.ngrid[1] + ijk[jloc][1]
            ) * this->params.ngrid[2] + ijk[kloc][2];

            if (idx_grid >= 0 && idx_grid < this->params.nmesh) {
              this->field[idx_grid][0] += inv_vol_cell
                * weight[pid][0]
                * win[iloc][0] * win[jloc][1] * win[kloc][2];
              this->field[idx_grid][1] += inv_vol_cell
                * weight[pid][1]
                * win[iloc][0] * win[jloc][1] * win[kloc][2];
            }
          }
        }
      }
    }
  }

  /**
   * Assign weighted density field to a mesh grid by the cloud-in-cell
   * (CIC) scheme.
   *
   * @param particles Particle container.
   * @param weight Particle weight.
   */
  void assign_weighted_field_to_grid_cic(
    ParticleContainer& particles,
    fftw_complex* weight
  ) {
    /// Set interpolation order, i.e. number of mesh grid, per dimension,
    /// to which a single particle is assigned
    int order = 2;

    /// Reset field with zero values.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] = 0.;
      this->field[gid][1] = 0.;
    }

    /// Here the density field is given by Σ_i w_i δ_D(x - x_i),
    /// where δ_D corresponds to δ_K / dV, dV =: `vol_cell`.
    double vol_cell = this->params.volume / double(this->params.nmesh);
    double inv_vol_cell = 1. / vol_cell;

    /// Perform assignment.
    for (int pid = 0; pid < particles.ntotal; pid++) {
      int ijk[order][3];  // coordinates of covered mesh grids
      double win[order][3];  // interpolation window
      for (int axis = 0; axis < 3; axis++) {
        double loc_grid = this->params.ngrid[axis]
          * particles[pid].pos[axis] / this->params.boxsize[axis];
        double s = loc_grid - double(int(loc_grid));

        /// Set up to 1st element as `order == 2`.
        ijk[0][axis] = int(loc_grid);
        ijk[1][axis] = int(loc_grid) + 1;

        win[0][axis] = 1. - s;
        win[1][axis] = s;
      }

      for (int iloc = 0; iloc < order; iloc++) {
        for (int jloc = 0; jloc < order; jloc++) {
          for (int kloc = 0; kloc < order; kloc++) {
            long long idx_grid = (
              ijk[iloc][0] * this->params.ngrid[1] + ijk[jloc][1]
            ) * this->params.ngrid[2] + ijk[kloc][2];

            if (idx_grid >= 0 && idx_grid < this->params.nmesh) {
              this->field[idx_grid][0] += inv_vol_cell
                * weight[pid][0]
                * win[iloc][0] * win[jloc][1] * win[kloc][2];
              this->field[idx_grid][1] += inv_vol_cell
                * weight[pid][1]
                * win[iloc][0] * win[jloc][1] * win[kloc][2];
            }
          }
        }
      }
    }
  }

  /**
   * Assign weighted density field to a grid by the
   * triangular-shaped-cloud (TSC) scheme.
   *
   * @param particles Particle container.
   * @param weight Particle weight.
   */
  void assign_weighted_field_to_grid_tsc(
    ParticleContainer& particles,
    fftw_complex* weight
  ) {
    /// Set interpolation order, i.e. number of mesh grid, per dimension,
    /// to which a single particle is assigned
    int order = 3;

    /// Reset field with zero values.
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field[gid][0] = 0.;
      this->field[gid][1] = 0.;
    }

    /// Here the over-density is given by Σ_i w_i δ_D(x - x_i),
    /// where δ_D corresponds to δ_K / dV, dV =: `vol_cell`.
    double vol_cell = this->params.volume / double(this->params.nmesh);
    double inv_vol_cell = 1. / vol_cell;

    /// Perform assignment.
    for (int pid = 0; pid < particles.ntotal; pid++) {
      int ijk[order][3];  // coordinates of covered mesh grids
      double win[order][3];  // interpolation window
      for (int axis = 0; axis < 3; axis++) {
        double loc_grid = this->params.ngrid[axis]
          * particles[pid].pos[axis] / this->params.boxsize[axis];
        double s = loc_grid - double(int(loc_grid + 0.5));

        /// Set up to 2nd element as `order == 3`.
        ijk[0][axis] = int(loc_grid + 0.5) - 1;
        ijk[1][axis] = int(loc_grid + 0.5);
        ijk[2][axis] = int(loc_grid + 0.5) + 1;

        win[0][axis] = 0.5 * (0.5 - s) * (0.5 - s);
        win[1][axis] = 0.75 - s * s;
        win[2][axis] = 0.5 * (0.5 + s) * (0.5 + s);
      }

      for (int iloc = 0; iloc < order; iloc++) {
        for (int jloc = 0; jloc < order; jloc++) {
          for (int kloc = 0; kloc < order; kloc++) {
            long long idx_grid = (
              ijk[iloc][0] * this->params.ngrid[1] + ijk[jloc][1]
            ) * this->params.ngrid[2] + ijk[kloc][2];

            if (idx_grid >= 0 && idx_grid < this->params.nmesh) {
              this->field[idx_grid][0] += inv_vol_cell
                * weight[pid][0]
                * win[iloc][0] * win[jloc][1] * win[kloc][2];
              this->field[idx_grid][1] += inv_vol_cell
                * weight[pid][1]
                * win[iloc][0] * win[jloc][1] * win[kloc][2];
            }
          }
        }
      }
    }
  }

  /**
   * Calculate the interpolarion window in Fourier space for
   * assignment schemes.
   *
   * @param kvec Wavevector.
   * @returns Window value in Fourier space.
   */
  double calc_assignment_window_in_fourier(double* kvec) {
    const double tol = 1.e-5;

    int order;
    if (
      this->params.assignment == "NGP"
    ) {
      order = 1;
    } else if (
      this->params.assignment == "CIC"
    ) {
      order = 2;
    } else if (
      this->params.assignment == "TSC"
    ) {
      order = 3;
    }

    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    int i = int(kvec[0] / dk[0] + tol);
    int j = int(kvec[1] / dk[1] + tol);
    int k = int(kvec[2] / dk[2] + tol);

    double k_x = M_PI * i / double(this->params.ngrid[0]);
    double k_y = M_PI * j / double(this->params.ngrid[1]);
    double k_z = M_PI * k / double(this->params.ngrid[2]);

    /// Note sin(u) / u -> 1 as u -> 0.
    double wk_x = (i != 0) ? sin(k_x) / k_x : 1.;
    double wk_y = (j != 0) ? sin(k_y) / k_y : 1.;
    double wk_z = (k != 0) ? sin(k_z) / k_z : 1.;

    double wk = wk_x * wk_y * wk_z;

    return pow(wk, order);
  }
};

/**
 * Binned two-point statistics from a particle container.
 *
 */
template <class ParticleContainer>
class TwoPointStatistics {
 public:
  double* k;  ///< central wavenumber in bins/shells
  std::complex<double>* sn;  ///< shot-noise statistics in Fourier space
  std::complex<double>* pk;  ///< power spectrum statistics
  std::complex<double>* xi;  ///< two-point correlation function statistics

  /**
   * Construct two-point statistics.
   *
   * @param params Parameter set.
   */
  TwoPointStatistics(ParameterSet& params){
    this->params = params;

    /// Set up binned power spectrum and mode counter.
    this->k = new double[params.num_kbin];
    this->sn = new std::complex<double>[params.num_kbin];
    this->pk = new std::complex<double>[params.num_kbin];
    this->nmode = new int[params.num_kbin];

    for (int ibin = 0; ibin < params.num_kbin; ibin++) {
      this->k[ibin] = 0.;
      this->pk[ibin] = 0.;
      this->sn[ibin] = 0.;
      this->nmode[ibin] = 0;
    }

    /// Set up binned two-point correlation function and pair counter.
    this->xi = new std::complex<double>[params.num_rbin];
    this->npair = new int[params.num_rbin];

    for (int ibin = 0; ibin < params.num_rbin; ibin++) {
      this->xi[ibin] = 0.;
      this->npair[ibin] = 0;
    }
  }

  /**
   * Destruct two-point statistics.
   */
  ~TwoPointStatistics(){
    finalise_2pt_stats();
  }

  /**
   * Finalise two-point statistics.
   */
  void finalise_2pt_stats() {
    delete[] this->k; this->k = NULL;
    delete[] this->sn; this->sn = NULL;
    delete[] this->pk; this->pk = NULL;
    delete[] this->nmode; this->nmode = NULL;
    delete[] this->xi; this->xi = NULL;
    delete[] this->npair; this->npair = NULL;
  }

  /**
   * Calculate two-point statistics in Fourier space.
   *
   * @param field_a First density field.
   * @param field_b Second density field.
   * @param shotnoise Shot noise amplitude.
   * @param kbin Wavenumber bins.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void calc_2pt_stats_in_fourier(
    DensityField<ParticleContainer> & field_a,
    DensityField<ParticleContainer> & field_b,
    std::complex<double> shotnoise,
    double* kbin,
    int ell, int m
  ) {
    /// Set up fine-sampling of wavenumbers.
    /// CAVEAT: discretionary choices.
    const double dk_sample = 1.e-4;
    const int n_sample = 1e5;

    double* ks_sample = new double[n_sample];
    int* nmode_sample = new int[n_sample];
    std::complex<double>* sn_sample = new std::complex<double>[n_sample];
    std::complex<double>* pk_sample = new std::complex<double>[n_sample];
    for (int i = 0; i < n_sample; i++) {
      ks_sample[i] = 0.;
      sn_sample[i] = 0.;
      pk_sample[i] = 0.;
      nmode_sample[i] = 0;
    }

    /// Set up binned statistics.
    for (int ibin = 0; ibin < this->params.num_kbin; ibin++) {
      this->k[ibin] = 0.;
      this->pk[ibin] = 0.;
      this->sn[ibin] = 0.;
      this->nmode[ibin] = 0;
    }

    /// Perform fine sampling.
    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    double kv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          kv[0] = (i < this->params.ngrid[0]/2) ?
            i * dk[0] : (i - this->params.ngrid[0]) * dk[0];
          kv[1] = (j < this->params.ngrid[1]/2) ?
            j * dk[1] : (j - this->params.ngrid[1]) * dk[1];
          kv[2] = (k < this->params.ngrid[2]/2) ?
            k * dk[2] :  (k - this->params.ngrid[2]) * dk[2];

          double k_ = sqrt(kv[0] * kv[0] + kv[1] * kv[1] + kv[2] * kv[2]);

          int idx_k = int(k_ / dk_sample + 0.5);
          if (idx_k < n_sample) {
            std::complex<double> delta_a(
              field_a[idx_grid][0], field_a[idx_grid][1]
            );
            std::complex<double> delta_b(
              field_b[idx_grid][0], field_b[idx_grid][1]
            );

            std::complex<double> mode_power = delta_a * conj(delta_b);
            std::complex<double> mode_sn =
              shotnoise * calc_shotnoise_scale_dependence(kv);

            /// Apply assignment compensation.
            double win = this->calc_assignment_window_in_fourier(kv);

            mode_power /= pow(win, 2);
            mode_sn /= pow(win, 2);

            /// Weight by reduced spherical harmonics.
            std::complex<double> ylm =
              SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
                ell, m, kv
              );

            mode_power *= ylm;
            mode_sn *= ylm;

            /// Add contribution.
            ks_sample[idx_k] += k_;
            pk_sample[idx_k] += mode_power;
            sn_sample[idx_k] += mode_sn;
            nmode_sample[idx_k]++;
          }
        }
      }
    }

    /// Perform binning.
    double dkbin[this->params.num_rbin - 1];
    for (int ibin = 0; ibin < this->params.num_rbin - 1; ibin++) {
      dkbin[ibin] = kbin[ibin + 1] - kbin[ibin];
    }

    for (int ibin = 0; ibin < this->params.num_kbin; ibin++) {
      double k_lower = 0.;
      if (ibin == 0) {
        k_lower = (kbin[ibin] > dkbin[ibin]/2) ?
          (kbin[ibin] - dkbin[ibin]/2) : 0.;
      } else {
        k_lower = kbin[ibin] - dkbin[ibin - 1]/2;
      }

      double k_upper = 0.;
      if (ibin == this->params.num_rbin - 1) {
        k_upper = kbin[ibin] + dkbin[ibin - 1]/2;
      } else {
        k_upper = kbin[ibin] + dkbin[ibin]/2;
      }

      for (int i = 0; i < n_sample; i++) {
        double k_sample = i * dk_sample;
        if (k_lower < k_sample && k_sample <= k_upper) {
          this->k[ibin] += ks_sample[i];
          this->pk[ibin] += pk_sample[i];
          this->sn[ibin] += sn_sample[i];
          this->nmode[ibin] += nmode_sample[i];
        }
      }
    }

    for (int ibin = 0; ibin < this->params.num_kbin; ibin++) {
      if (this->nmode[ibin] != 0) {
        this->k[ibin] /= double(this->nmode[ibin]);
        this->pk[ibin] /= double(this->nmode[ibin]);
        this->sn[ibin] /= double(this->nmode[ibin]);
      } else {
        this->k[ibin] = (kbin[ibin] + kbin[ibin + 1]) / 2.;
        this->pk[ibin] = 0.;
        this->sn[ibin] = 0.;
      }
    }

    delete[] ks_sample;
    delete[] pk_sample;
    delete[] sn_sample;
    delete[] nmode_sample;
  }

  /**
   * Calculate two-point statistics in configuration space.
   *
   * @param field_a First density field.
   * @param field_b Second density field.
   * @param shotnoise Shot noise amplitude.
   * @param rbin Separation bins.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void calc_2pt_stats_in_config(
    DensityField<ParticleContainer>& field_a,
    DensityField<ParticleContainer>& field_b,
    std::complex<double> shotnoise,
    double* rbin,
    int ell, int m
  ) {
    /// Set up grid sampling (before inverse Fourier transform).
    /// FIXME: Change variable name.
    fftw_complex* twopt3d_sample = fftw_alloc_complex(this->params.nmesh);
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      twopt3d_sample[gid][0] = 0.;
      twopt3d_sample[gid][1] = 0.;
    }

    /// Set inverse FFT volume normalisation, where ∫d^3k / (2\pi)^3
    /// corresponds to (1/V) Σ_i, V =: `vol`.
    double vol = this->params.volume;

    /// Compute gridded statistics.
    /// FIXME: "See also `calc_powspec` above."
    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    double kv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          kv[0] = (i < this->params.ngrid[0]/2) ?
            i * dk[0] : (i - this->params.ngrid[0]) * dk[0];
          kv[1] = (j < this->params.ngrid[1]/2) ?
            j * dk[1] : (j - this->params.ngrid[1]) * dk[1];
          kv[2] = (k < this->params.ngrid[2]/2) ?
            k * dk[2] : (k - this->params.ngrid[2]) * dk[2];

          std::complex<double> delta_a(
            field_a[idx_grid][0], field_a[idx_grid][1]
          );
          std::complex<double> delta_b(
            field_b[idx_grid][0], field_b[idx_grid][1]
          );

          std::complex<double> mode_power = delta_a * conj(delta_b);

          /// Subtract shot noise component.
          mode_power -= shotnoise * calc_shotnoise_scale_dependence(kv);

          /// Apply assignment compensation.
          double win = this->calc_assignment_window_in_fourier(kv);

          mode_power /= pow(win, 2);

          twopt3d_sample[idx_grid][0] = mode_power.real() / vol;
          twopt3d_sample[idx_grid][1] = mode_power.imag() / vol;
        }
      }
    }

    /// Inverse Fourier transform.
    fftw_plan inv_transform = fftw_plan_dft_3d(
      this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
      twopt3d_sample, twopt3d_sample,
      FFTW_BACKWARD, FFTW_ESTIMATE
    );

    fftw_execute(inv_transform);
    fftw_destroy_plan(inv_transform);

    /// Set up fine-sampling of separations.
    /// CAVEAT: discretionary choices.
    const double dr_sample = 0.5;
    const int n_sample = 10000;

    std::complex<double>* xi_sample = new std::complex<double>[n_sample];
    int* npair_sample = new int[n_sample];
    for (int i = 0; i < n_sample; i++) {
      xi_sample[i] = 0.;
      npair_sample[i] = 0;
    }

    /// Set up binned two-point statistics.
    for (int ibin = 0; ibin < this->params.num_rbin; ibin++) {
      this->xi[ibin] = 0.;
      this->npair[ibin] = 0;
    }

    /// Perform fine sampling.
    double dr[3];
    dr[0] = this->params.boxsize[0] / double(this->params.ngrid[0]);
    dr[1] = this->params.boxsize[1] / double(this->params.ngrid[1]);
    dr[2] = this->params.boxsize[2] / double(this->params.ngrid[2]);

    double rv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          rv[0] = (i < this->params.ngrid[0]/2) ?
            i * dr[0] : (i - this->params.ngrid[0]) * dr[0];
          rv[1] = (j < this->params.ngrid[1]/2) ?
            j * dr[1] : (j - this->params.ngrid[1]) * dr[1];
          rv[2] = (k < this->params.ngrid[2]/2) ?
            k * dr[2] : (k - this->params.ngrid[2]) * dr[2];

          double r_ = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);

          int idx_r = int(r_ / dr_sample + 0.5);
          if (idx_r < n_sample) {
            std::complex<double> pair_corr(
              twopt3d_sample[idx_grid][0], twopt3d_sample[idx_grid][1]
            );

            /// Weight by reduced spherical harmonics.
            std::complex<double> ylm =
              SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
                ell, m, rv
              );

            pair_corr *= ylm;

            /// Add contribution.
            xi_sample[idx_r] += pair_corr;
            npair_sample[idx_r]++;
          }
        }
      }
    }

    /// Perform binning.
    double drbin[this->params.num_rbin - 1];
    for (int ibin = 0; ibin < this->params.num_rbin - 1; ibin++) {
      drbin[ibin] = rbin[ibin + 1] - rbin[ibin];
    }

    for (int ibin = 0; ibin < this->params.num_rbin; ibin++) {
      double r_lower = 0.;
      if (ibin == 0) {
        r_lower = (rbin[ibin] > drbin[ibin]/2) ?
          (rbin[ibin] - drbin[ibin]/2) : 0.;
      } else {
        r_lower = rbin[ibin] - drbin[ibin - 1]/2;
      }

      double r_upper = 0.;
      if (ibin == this->params.num_rbin - 1) {
        r_upper = rbin[ibin] + drbin[ibin - 1]/2;
      } else {
        r_upper = rbin[ibin] + drbin[ibin]/2;
      }

      for (int i = 0; i < n_sample; i++) {
        double r_sample = i * dr_sample;
        if (r_lower < r_sample && r_sample <= r_upper) {
          this->xi[ibin] += xi_sample[i];
          this->npair[ibin] += npair_sample[i];
        }
      }
    }

    for (int ibin = 0; ibin < this->params.num_rbin; ibin++) {
      if (this->npair[ibin] != 0) {
        this->xi[ibin] /= double(this->npair[ibin]);
      } else {
        this->xi[ibin] = 0.;
      }
    }

    delete[] twopt3d_sample;
    delete[] xi_sample;
    delete[] npair_sample;
  }

  /**
   * Calculate two-point statistics for calculating shot-noise components
   * in the three-point correlation function.
   *
   * See eq. (51) in Sugiyama et al. (2018)
   * [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>].
   *
   * @param field_a First field.
   * @param field_b Second field.
   * @param ylm_a Reduced spherical harmonics over the first field grid.
   * @param ylm_b Reduced spherical harmonics over the second field grid.
   * @param shotnoise Shot noise amplitude.
   * @param rbin Separation bins.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   */
  void calc_2pt_stats_for_3pcf(
    DensityField<ParticleContainer>& field_a,
    DensityField<ParticleContainer>& field_b,
    std::complex<double>* ylm_a, std::complex<double>* ylm_b,
    std::complex<double> shotnoise,
    double* rbin,
    int ell, int m
  ) {
    /// Set up grid sampling (before inverse Fourier transform).
    /// FIXME: Change variable name.
    fftw_complex* twopt3d_sample = fftw_alloc_complex(this->params.nmesh);
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      twopt3d_sample[gid][0] = 0.;
      twopt3d_sample[gid][1] = 0.;
    }

    /// Set inverse FFT volume normalisation, where ∫d^3k / (2\pi)^3
    /// corresponds to (1/V) Σ_i, V =: `vol`.
    double vol = this->params.volume;

    /// Compute gridded statistics.
    /// FIXME: "See also `calc_powspec` above."
    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    double kv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          kv[0] = (i < this->params.ngrid[0]/2) ?
            i * dk[0] : (i - this->params.ngrid[0]) * dk[0];
          kv[1] = (j < this->params.ngrid[1]/2) ?
            j * dk[1] : (j - this->params.ngrid[1]) * dk[1];
          kv[2] = (k < this->params.ngrid[2]/2) ?
            k * dk[2] : (k - this->params.ngrid[2]) * dk[2];

          std::complex<double> delta_a(
            field_a[idx_grid][0], field_a[idx_grid][1]
          );
          std::complex<double> delta_b(
            field_b[idx_grid][0], field_b[idx_grid][1]
          );

          std::complex<double> mode_power = delta_a * conj(delta_b);

          /// Subtract shot noise component.
          mode_power -= shotnoise * calc_shotnoise_scale_dependence(kv);

          /// Apply assignment compensation.
          double win = this->calc_assignment_window_in_fourier(kv);

          mode_power /= pow(win, 2);

          twopt3d_sample[idx_grid][0] = mode_power.real() / vol;
          twopt3d_sample[idx_grid][1] = mode_power.imag() / vol;
        }
      }
    }

    /// Inverse Fourier transform.
    fftw_plan inv_transform = fftw_plan_dft_3d(
      this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
      twopt3d_sample, twopt3d_sample,
      FFTW_BACKWARD, FFTW_ESTIMATE
    );

    fftw_execute(inv_transform);
    fftw_destroy_plan(inv_transform);

    /// Set up fine-sampling of separations.
    /// CAVEAT: discretionary choices.
    const double dr_sample = 0.5;
    const int n_sample = 10000;

    std::complex<double>* xi_sample = new std::complex<double>[n_sample];
    int* npair_sample = new int[n_sample];
    for (int i = 0; i < n_sample; i++) {
      xi_sample[i] = 0.;
      npair_sample[i] = 0;
    }

    for (int ibin = 0; ibin < this->params.num_rbin; ibin++) {
      this->xi[ibin] = 0.;
      this->npair[ibin] = 0;
    }

    /// Perform fine sampling.
    double dr[3];
    dr[0] = this->params.boxsize[0] / double(this->params.ngrid[0]);
    dr[1] = this->params.boxsize[1] / double(this->params.ngrid[1]);
    dr[2] = this->params.boxsize[2] / double(this->params.ngrid[2]);

    double rv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          rv[0] = (i < this->params.ngrid[0]/2) ?
            i * dr[0] : (i - this->params.ngrid[0]) * dr[0];
          rv[1] = (j < this->params.ngrid[1]/2) ?
            j * dr[1] : (j - this->params.ngrid[1]) * dr[1];
          rv[2] = (k < this->params.ngrid[2]/2) ?
            k * dr[2] : (k - this->params.ngrid[2]) * dr[2];

          double r_ = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);

          int idx_r = int(r_ / dr_sample + 0.5);
          if (idx_r < n_sample) {
            std::complex<double> pair_corr(
              twopt3d_sample[idx_grid][0], twopt3d_sample[idx_grid][1]
            );

            /// Weight by reduced spherical harmonics.
            pair_corr *= ylm_a[idx_grid] * ylm_b[idx_grid];

            /// Add contribution.
            xi_sample[idx_r] += pair_corr;
            npair_sample[idx_r]++;
          }
        }
      }
    }

    /// Perform binning.
    double drbin[this->params.num_rbin - 1];
    for (int ibin = 0; ibin < this->params.num_rbin - 1; ibin++) {
      drbin[ibin] = rbin[ibin + 1] - rbin[ibin];
    }

    for (int ibin = 0; ibin < this->params.num_rbin; ibin++) {
      double r_lower = 0.;
      if (ibin == 0) {
        r_lower = (rbin[ibin] > drbin[ibin]/2) ?
          (rbin[ibin] - drbin[ibin]/2) : 0.;
      } else {
        r_lower = rbin[ibin] - drbin[ibin - 1]/2;
      }

      double r_upper = 0.;
      if (ibin == this->params.num_rbin - 1) {
        r_upper = rbin[ibin] + drbin[ibin - 1]/2;
      } else {
        r_upper = rbin[ibin] + drbin[ibin]/2;
      }

      for (int i = 0; i < n_sample; i++) {
        double r_sample = i * dr_sample;
        if (r_lower < r_sample && r_sample <= r_upper) {
          this->xi[ibin] += xi_sample[i];
          this->npair[ibin] += npair_sample[i];
        }
      }
    }

    /// Apply normalisation constants.
    double vol_cell = this->params.volume / double(this->params.nmesh);

    for (int ibin = 0; ibin < this->params.num_rbin; ibin++) {
      if (this->npair[ibin] != 0) {
        this->xi[ibin] *=
          pow(-1., this->params.ell1 + this->params.ell2)
          / vol_cell
          / double(this->npair[ibin]) / double(this->npair[ibin]);
      } else {
        this->xi[ibin] = 0.;
      }
    }

    delete[] twopt3d_sample;
    delete[] xi_sample;
    delete[] npair_sample;
  }

  /**
   * Calculate shot noise for power spectrum.
   *
   * @param particles_data (Data-source) particle container.
   * @param particles_rand (Random-source) particle container.
   * @param los_data (Data-source) particle lines of sight.
   * @param los_rand (Random-source) particle lines of sight.
   * @param alpha Alpha ratio.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   * @returns Shot noise for power spectrum.
   */
  std::complex<double> calc_shotnoise_for_powspec(
    ParticleContainer& particles_data,
    ParticleContainer& particles_rand,
    LineOfSight* los_data,
    LineOfSight* los_rand,
    double alpha,
    int ell, int m
  ) {
    /// Initialise shot noise contributions.
    std::complex<double> sum_data = 0.;
    std::complex<double> sum_rand = 0.;

    /// Perform direct summation with spherical harmonic weighting.
    /// See Jing (2005) [arXiv: astro-ph/0409240].
    for (int pid = 0; pid < particles_data.ntotal; pid++) {
      double los[3] = {
        los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(ell, m, los);

      sum_data += ylm * pow(particles_data[pid].w, 2);
    }

    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      double los[3] = {
        los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(ell, m, los);

      sum_rand += ylm * pow(particles_rand[pid].w, 2);
    }

    return sum_data + pow(alpha, 2) * sum_rand;
  }

  /**
   * Calculate shot noise for power spectrum in a periodic box
   * for reconstruction.
   *
   * @param particles_data (Data-source) particle container.
   * @param particles_rand (Random-source) particle container.
   * @param alpha Alpha ratio.
   * @returns Shot noise for power spectrum.
   */
  std::complex<double> calc_shotnoise_for_powspec_in_box_for_recon(
    ParticleContainer& particles_data,
    ParticleContainer& particles_rand,
    double alpha
  ) {
    std::complex<double> sum_data = double(particles_data.ntotal);
    std::complex<double> sum_rand = double(particles_rand.ntotal);

    return sum_data + pow(alpha, 2) * sum_rand;
  }

  /**
   * Calculate shot noise for correlation function window.
   *
   * @param particles_rand (Random-source) particle container.
   * @param los_rand (Random-source) particle lines of sight.
   * @param alpha Alpha ratio.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   * @returns Shot noise for correlation function window.
   */
  std::complex<double> calc_shotnoise_for_corrfunc_window(
    ParticleContainer& particles_rand,
    LineOfSight* los_rand,
    double alpha,
    int ell, int m
  ) {
    std::complex<double> sum_rand = 0.;
    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      double los[3] = {
        los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(ell, m, los);

      sum_rand += ylm * pow(particles_rand[pid].w, 2);
    }

    return pow(alpha, 2) * sum_rand;
  }

  /**
   * Calculate shot noise for bispectrum from pure self-contributions.
   *
   * @param particles_data (Data-source) particle container.
   * @param particles_rand (Random-source) particle container.
   * @param los_data (Data-source) particle lines of sight.
   * @param los_rand (Random-source) particle lines of sight.
   * @param alpha Alpha ratio.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   * @returns Shot noise for bispectrum.
   */
  std::complex<double> calc_shotnoise_for_bispec_from_self(
    ParticleContainer& particles_data, ParticleContainer& particles_rand,
    LineOfSight* los_data, LineOfSight* los_rand,
    double alpha,
    int ell, int m
  ) {
    /// Initialise shot noise contributions.
    std::complex<double> sum_data = 0.;
    std::complex<double> sum_rand = 0.;

    /// Perform direct summation with spherical harmonic weighting.
    for (int pid = 0; pid < particles_data.ntotal; pid++) {
      double los[3] = {
        los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(ell, m, los);

      sum_data += ylm * pow(particles_data[pid].w, 3);
    }

    for (int pid = 0; pid < particles_rand.ntotal; pid++) {
      double los[3] = {
        los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
      };

      std::complex<double> ylm =
        SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(ell, m, los);

      sum_rand += ylm * pow(particles_rand[pid].w, 3);
    }

    /// Calculate \bar{S}_LM in eq. (46) in arXiv:1803.02132.
    return sum_data - pow(alpha, 3) * sum_rand;
  }

  /**
   * Calculate shot noise for bispectrum for each mesh grid.
   *
   * @param field_a First density field.
   * @param field_b Second density field.
   * @param shotnoise Shot noise.
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   * @param three_pt Three-point function placeholder.
   * @returns Exit status.
   */
  int calc_shotnoise_for_bispec_on_grid(
    DensityField<ParticleContainer>& field_a,
    DensityField<ParticleContainer>& field_b,
    std::complex<double> shotnoise,
    int ell, int m,
    fftw_complex* three_pt
  ) {
    /// Set inverse FFT volume normalisation, where ∫d^3k / (2\pi)^3
    /// corresponds to (1/V) Σ_i, V =: `vol`.
    double vol = this->params.volume;

    /// Compute gridded power spectrum.  See also `calc_powspec` above.
    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    double kv[3];
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid =
            (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;

          kv[0] = (i < this->params.ngrid[0]/2) ?
            i * dk[0] : (i - this->params.ngrid[0]) * dk[0];
          kv[1] = (j < this->params.ngrid[1]/2) ?
            j * dk[1] : (j - this->params.ngrid[1]) * dk[1];
          kv[2] = (k < this->params.ngrid[2]/2) ?
            k * dk[2] : (k - this->params.ngrid[2]) * dk[2];

          std::complex<double> delta_a(
            field_a[idx_grid][0], field_a[idx_grid][1]
          );
          std::complex<double> delta_b(
            field_b[idx_grid][0], field_b[idx_grid][1]
          );
          std::complex<double> mode_power = delta_a * conj(delta_b);

          mode_power -= shotnoise * calc_shotnoise_scale_dependence(kv);

          double win = this->calc_assignment_window_in_fourier(kv);
          mode_power /= pow(win, 2);

          three_pt[idx_grid][0] = mode_power.real() / vol;
          three_pt[idx_grid][1] = mode_power.imag() / vol;
        }
      }
    }

    /// Inverse Fourier transform.  Calculate
    /// S_{\ell_1 \ell_2 L; i = j != k} in eq. (45) in arXiv:1803.02132.
    fftw_plan fft_plan_backward = fftw_plan_dft_3d(
      this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
      three_pt, three_pt,
      FFTW_BACKWARD, FFTW_ESTIMATE
    );
    fftw_execute(fft_plan_backward);
    fftw_destroy_plan(fft_plan_backward);

    return 0;
  }

 private:
  ParameterSet params;
  int* nmode;  ///< number of wavevector modes
  int* npair;  ///< number of separation pairs

  /**
   * Calculate the interpolarion window in Fourier space for
   * assignment schemes.
   *
   * @param kvec Wavevector.
   * @returns Window value in Fourier space.
   */
  double calc_assignment_window_in_fourier(double* kvec) {
    const double tol = 1.e-5;

    int order;
    if (
      this->params.assignment == "NGP"
    ) {
      order = 1;
    } else if (
      this->params.assignment == "CIC"
    ) {
      order = 2;
    } else if (
      this->params.assignment == "TSC"
    ) {
      order = 3;
    }

    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    int i = int(kvec[0] / dk[0] + tol);
    int j = int(kvec[1] / dk[1] + tol);
    int k = int(kvec[2] / dk[2] + tol);

    double k_x = M_PI * i / double(this->params.ngrid[0]);
    double k_y = M_PI * j / double(this->params.ngrid[1]);
    double k_z = M_PI * k / double(this->params.ngrid[2]);

    /// Note sin(u) / u -> 1 as u -> 0.
    double wk_x = (i != 0) ? sin(k_x) / k_x : 1.;
    double wk_y = (j != 0) ? sin(k_y) / k_y : 1.;
    double wk_z = (k != 0) ? sin(k_z) / k_z : 1.;

    double wk = wk_x * wk_y * wk_z;

    return pow(wk, order);
  }

  /**
   * Calculate the shot noise scale-dependent function f@$ C(\vec{k}) f@$.
   *
   * See eqs. (45) and (46) in Sugiyama et al. (2018)
   * [<a href="https://arxiv.org/abs/1803.02132">1803.02132</a>]
   * and Jing (2004)
   * [<a href="https://arxiv.org/abs/astro-ph/0409240">astro-ph/0409240</a>].
   *
   * @param kvec Wavevector.
   * @returns Value of the shot noise function.
   */
  double calc_shotnoise_scale_dependence(double* kvec) {
    if (
      this->params.assignment == "NGP"
    ) {
      return this->calc_shotnoise_scale_dependence_ngp(kvec);
    } else if (
      this->params.assignment == "CIC"
    ) {
      return this->calc_shotnoise_scale_dependence_cic(kvec);
    } else if (
      this->params.assignment == "TSC"
    ) {
      return this->calc_shotnoise_scale_dependence_tsc(kvec);
    }
  }

  /**
   * Calculate the shot noise scale-dependent function for the
   * nearest-grid-point assignment scheme.
   *
   * @param kvec Wavevector.
   * @returns Function value.
   */
  double calc_shotnoise_scale_dependence_ngp(const double* kvec) {
    return 1.;
  }

  /**
   * Calculate the shot noise scale-dependent function for the
   * cloud-in-cell assignment scheme.
   *
   * @param kvec Wavevector.
   * @returns Function value.
   */
  double calc_shotnoise_scale_dependence_cic(const double* kvec) {
    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    int i = int(kvec[0] / dk[0] + tol);
    int j = int(kvec[1] / dk[1] + tol);
    int k = int(kvec[2] / dk[2] + tol);

    double k_x = M_PI * i / double(this->params.ngrid[0]);
    double k_y = M_PI * j / double(this->params.ngrid[1]);
    double k_z = M_PI * k / double(this->params.ngrid[2]);

    double cx = (i != 0) ? sin(k_x): 0.;
    double cy = (j != 0) ? sin(k_y): 0.;
    double cz = (k != 0) ? sin(k_z): 0.;

    double val =
      (1. - 2./3. * cx * cx)
      * (1. - 2./3. * cy * cy)
      * (1. - 2./3. * cz * cz);

    return val;
  }

  /**
   * Calculate the shot noise scale-dependent function for the
   * triangular-shaped-cloud assignment scheme.
   *
   * @param kvec Wavevector.
   * @returns Function value.
   */
  double calc_shotnoise_scale_dependence_tsc(const double* kvec) {
    const double tol = 1.e-5;

    double dk[3];
    dk[0] = 2.*M_PI / this->params.boxsize[0];
    dk[1] = 2.*M_PI / this->params.boxsize[1];
    dk[2] = 2.*M_PI / this->params.boxsize[2];

    int i = int(kvec[0] / dk[0] + tol);
    int j = int(kvec[1] / dk[1] + tol);
    int k = int(kvec[2] / dk[2] + tol);

    double k_x = M_PI * i / double(this->params.ngrid[0]);
    double k_y = M_PI * j / double(this->params.ngrid[1]);
    double k_z = M_PI * k / double(this->params.ngrid[2]);

    double cx = (i != 0) ? sin(k_x): 0.;
    double cy = (j != 0) ? sin(k_y): 0.;
    double cz = (k != 0) ? sin(k_z): 0.;

    double val =
      (1. - cx * cx + 2./15. * pow(cx, 4))
      * (1. - cy * cy + 2./15. * pow(cy, 4))
      * (1. - cz * cz + 2./15. * pow(cz, 4));

    return val;
  }
};

#endif
