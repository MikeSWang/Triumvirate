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
 * @file field.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "field.hpp"

namespace trvs = trv::sys;
namespace trvm = trv::maths;

namespace trv {

// ***********************************************************************
// Mesh field
// ***********************************************************************

// -----------------------------------------------------------------------
// Life cycle
// -----------------------------------------------------------------------

MeshField::MeshField(trv::ParameterSet& params) {
  // Attach the full parameter set to @ref trv::MeshField.
  this->params = params;

  trvs::logger.reset_level(params.verbose);

  // Initialise the field (and its shadow field if interlacing is used)
  // and increase allocated memory.
  this->field = fftw_alloc_complex(this->params.nmesh);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(this->params.nmesh);
  trv::sys::update_maxmem();

  if (this->params.interlace == "true") {
    this->field_s = fftw_alloc_complex(this->params.nmesh);

    trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(this->params.nmesh);
    trv::sys::update_maxmem();
  }

  this->initialise_density_field();  // likely redundant but safe

  // Calculate grid sizes in configuration space.
  this->dr[0] = this->params.boxsize[0] / this->params.ngrid[0];
  this->dr[1] = this->params.boxsize[1] / this->params.ngrid[1];
  this->dr[2] = this->params.boxsize[2] / this->params.ngrid[2];

  // Calculate fundamental wavenumbers in Fourier space.
  this->dk[0] = 2.*M_PI / this->params.boxsize[0];
  this->dk[1] = 2.*M_PI / this->params.boxsize[1];
  this->dk[2] = 2.*M_PI / this->params.boxsize[2];

  // Calculate mesh volume and mesh grid cell volume.
  this->vol = this->params.volume;
  this->vol_cell = this->vol / double(this->params.nmesh);
}

MeshField::~MeshField() {this->finalise_density_field();}

void MeshField::initialise_density_field() {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] = 0.;
    this->field[gid][1] = 0.;
  }
  if (this->params.interlace == "true") {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field_s[gid][0] = 0.;
      this->field_s[gid][1] = 0.;
    }
  }
}

void MeshField::finalise_density_field() {
  // Free memory usage.
  if (this->field != nullptr) {
    fftw_free(this->field); this->field = nullptr;
    trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(this->params.nmesh);
  }
  if (this->field_s != nullptr) {
    fftw_free(this->field_s); this->field_s = nullptr;
    trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(this->params.nmesh);
  }
}


// -----------------------------------------------------------------------
// Operators & reserved methods
// -----------------------------------------------------------------------

const fftw_complex& MeshField::operator[](int gid) {return this->field[gid];}


// -----------------------------------------------------------------------
// Mesh grid properties
// -----------------------------------------------------------------------

long long MeshField::get_grid_index(int i, int j, int k) {
  long long idx_grid =
    (i * this->params.ngrid[1] + j) * this->params.ngrid[2] + k;
  return idx_grid;
}

void MeshField::get_grid_pos_vector(int i, int j, int k, double rvec[3]) {
  rvec[0] = (i < this->params.ngrid[0]/2) ?
    i * this->dr[0] : (i - this->params.ngrid[0]) * this->dr[0];
  rvec[1] = (j < this->params.ngrid[1]/2) ?
    j * this->dr[1] : (j - this->params.ngrid[1]) * this->dr[1];
  rvec[2] = (k < this->params.ngrid[2]/2) ?
    k * this->dr[2] : (k - this->params.ngrid[2]) * this->dr[2];
}

void MeshField::get_grid_wavevector(int i, int j, int k, double kvec[3]) {
  kvec[0] = (i < this->params.ngrid[0]/2) ?
    i * this->dk[0] : (i - this->params.ngrid[0]) * this->dk[0];
  kvec[1] = (j < this->params.ngrid[1]/2) ?
    j * this->dk[1] : (j - this->params.ngrid[1]) * this->dk[1];
  kvec[2] = (k < this->params.ngrid[2]/2) ?
    k * this->dk[2] : (k - this->params.ngrid[2]) * this->dk[2];
}


// -----------------------------------------------------------------------
// Mesh assignment
// -----------------------------------------------------------------------

void MeshField::assign_weighted_field_to_mesh(
  ParticleCatalogue& particles, fftw_complex* weights
) {
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    double extent = particles.pos_max[iaxis] - particles.pos_min[iaxis];
    if (params.boxsize[iaxis] < extent) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Box size in dimension %d is smaller than catalogue extents: "
          "%.3f < %.3f.",
          iaxis, params.boxsize[iaxis], extent
        );
      }
    }
  }

  if (this->params.assignment == "ngp") {
    this->assign_weighted_field_to_mesh_ngp(particles, weights);
  } else
  if (this->params.assignment == "cic") {
    this->assign_weighted_field_to_mesh_cic(particles, weights);
  } else
  if (this->params.assignment == "tsc") {
    this->assign_weighted_field_to_mesh_tsc(particles, weights);
  } else
  if (this->params.assignment == "pcs") {
    this->assign_weighted_field_to_mesh_pcs(particles, weights);
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Unsupported mesh assignment scheme: '%s'.",
        this->params.assignment.c_str()
      );
      throw trvs::InvalidParameterError(
        "Unsupported mesh assignment scheme: '%s'.\n",
        this->params.assignment.c_str()
      );
    };
  }
}

void MeshField::assign_weighted_field_to_mesh_ngp(
  ParticleCatalogue& particles, fftw_complex* weight
) {
  // Set interpolation order, i.e. number of grids, per dimension,
  // to which a single particle is assigned.
  const int order = 1;

  // Here the field is given by Σᵢ wᵢ δᴰ(x - xᵢ), where δᴰ ↔ δᴷ / dV,
  // dV =: `vol_cell`.
  const double inv_vol_cell = 1 / this->vol_cell;

  // Reset field values to zero.
  this->initialise_density_field();

  // Assign particles to grid cells.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    int ijk[order][3];     // grid index coordinates of covered grid cells
    double win[order][3];  // sampling window
    long long gid = 0;     // flattened grid cell index

    for (int iaxis = 0; iaxis < 3; iaxis++) {
      // Carefully set covered sampling window grid indices.
      double loc_grid = this->params.ngrid[iaxis] *
        particles[pid].pos[iaxis] / this->params.boxsize[iaxis];

      int idx_grid = int(loc_grid);
      if (loc_grid - idx_grid >= 0.5) {
        idx_grid = (idx_grid == this->params.ngrid[iaxis] - 1)
          ? 0 : idx_grid + 1;
      }

      ijk[0][iaxis] = idx_grid;

      // Set sampling window value (only 0th element as ``order == 1``).
      win[0][iaxis] = 1.;
    }

    for (int iloc = 0; iloc < order; iloc++) {
      for (int jloc = 0; jloc < order; jloc++) {
        for (int kloc = 0; kloc < order; kloc++) {
          gid = this->get_grid_index(
            ijk[iloc][0], ijk[jloc][1], ijk[kloc][2]
          );
          if (0 <= gid && gid < this->params.nmesh) {
OMP_ATOMIC
            this->field[gid][0] += inv_vol_cell
              * weight[pid][0] * win[iloc][0] * win[jloc][1] * win[kloc][2];
OMP_ATOMIC
            this->field[gid][1] += inv_vol_cell
              * weight[pid][1] * win[iloc][0] * win[jloc][1] * win[kloc][2];
          }
        }
      }
    }
  }

  // Perform interlacing if needed.
  if (this->params.interlace == "true") {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int pid = 0; pid < particles.ntotal; pid++) {
      int ijk[order][3];
      double win[order][3];
      long long gid = 0;

      for (int iaxis = 0; iaxis < 3; iaxis++) {
        // Apply a half-grid shift and impose the periodic boundary condition.
        double loc_grid = this->params.ngrid[iaxis]
          * particles[pid].pos[iaxis] / this->params.boxsize[iaxis] + 0.5;

        if (loc_grid > this->params.ngrid[iaxis]) {
          loc_grid -= this->params.ngrid[iaxis];
        }

        int idx_grid = int(loc_grid);
        if (loc_grid - idx_grid >= 0.5) {
          idx_grid = (idx_grid == this->params.ngrid[iaxis] - 1)
            ? 0 : idx_grid + 1;
        }

        ijk[0][iaxis] = idx_grid;

        win[0][iaxis] = 1.;
      }

      for (int iloc = 0; iloc < order; iloc++) {
        for (int jloc = 0; jloc < order; jloc++) {
          for (int kloc = 0; kloc < order; kloc++) {
            gid = this->get_grid_index(
              ijk[iloc][0], ijk[jloc][1], ijk[kloc][2]
            );
            if (0 <= gid && gid < this->params.nmesh) {
OMP_ATOMIC
              this->field_s[gid][0] += inv_vol_cell
                * weight[pid][0] * win[iloc][0] * win[jloc][1] * win[kloc][2];
OMP_ATOMIC
              this->field_s[gid][1] += inv_vol_cell
                * weight[pid][1] * win[iloc][0] * win[jloc][1] * win[kloc][2];
            }
          }
        }
      }
    }
  }
}

void MeshField::assign_weighted_field_to_mesh_cic(
  ParticleCatalogue& particles, fftw_complex* weight
) {
  // Set interpolation order, i.e. number of grids, per dimension,
  // to which a single particle is assigned.
  const int order = 2;

  // Here the field is given by Σᵢ wᵢ δᴰ(x - xᵢ), where δᴰ ↔ δᴷ / dV,
  // dV =: `vol_cell`.
  const double inv_vol_cell = 1 / this->vol_cell;

  // Reset field values to zero.
  this->initialise_density_field();

  // Assign particles to grid cells.

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    int ijk[order][3];     // grid index coordinates of covered grid cells
    double win[order][3];  // sampling window
    long long gid = 0;     // flattened grid cell index

    for (int iaxis = 0; iaxis < 3; iaxis++) {
      // Carefully set covered sampling window grid indices.
      double loc_grid = this->params.ngrid[iaxis]
        * particles[pid].pos[iaxis] / this->params.boxsize[iaxis];

      int idx_grid = int(loc_grid);

      ijk[0][iaxis] = idx_grid;
      ijk[1][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
        ? 0 : idx_grid + 1;

      // Set sampling window value (up to the 1st element as `order == 2`).
      double s = loc_grid - idx_grid;  // particle-to-grid grid-index distance

      win[0][iaxis] = 1. - s;
      win[1][iaxis] = s;
    }

    for (int iloc = 0; iloc < order; iloc++) {
      for (int jloc = 0; jloc < order; jloc++) {
        for (int kloc = 0; kloc < order; kloc++) {
          gid = this->get_grid_index(
            ijk[iloc][0], ijk[jloc][1], ijk[kloc][2]
          );
          if (0 <= gid && gid < this->params.nmesh) {
OMP_ATOMIC
            this->field[gid][0] += inv_vol_cell
              * weight[pid][0] * win[iloc][0] * win[jloc][1] * win[kloc][2];
OMP_ATOMIC
            this->field[gid][1] += inv_vol_cell
              * weight[pid][1] * win[iloc][0] * win[jloc][1] * win[kloc][2];
          }
        }
      }
    }
  }

  // Perform interlacing if needed.
  if (this->params.interlace == "true") {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int pid = 0; pid < particles.ntotal; pid++) {
      int ijk[order][3];
      double win[order][3];
      long long gid = 0;

      for (int iaxis = 0; iaxis < 3; iaxis++) {
        // Apply a half-grid shift and impose the periodic boundary condition.
        double loc_grid = this->params.ngrid[iaxis]
          * particles[pid].pos[iaxis] / this->params.boxsize[iaxis] + 0.5;

        if (loc_grid > this->params.ngrid[iaxis]) {
          loc_grid -= this->params.ngrid[iaxis];
        }

        int idx_grid = int(loc_grid);

        ijk[0][iaxis] = idx_grid;
        ijk[1][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
          ? 0 : idx_grid + 1;

        double s = loc_grid - idx_grid;
        win[0][iaxis] = 1. - s;
        win[1][iaxis] = s;
      }

      for (int iloc = 0; iloc < order; iloc++) {
        for (int jloc = 0; jloc < order; jloc++) {
          for (int kloc = 0; kloc < order; kloc++) {
            gid = this->get_grid_index(
              ijk[iloc][0], ijk[jloc][1], ijk[kloc][2]
            );
            if (0 <= gid && gid < this->params.nmesh) {
OMP_ATOMIC
              this->field_s[gid][0] += inv_vol_cell
                * weight[pid][0] * win[iloc][0] * win[jloc][1] * win[kloc][2];
OMP_ATOMIC
              this->field_s[gid][1] += inv_vol_cell
                * weight[pid][1] * win[iloc][0] * win[jloc][1] * win[kloc][2];
            }
          }
        }
      }
    }
  }
}

void MeshField::assign_weighted_field_to_mesh_tsc(
  ParticleCatalogue& particles, fftw_complex* weight
) {
  // Set interpolation order, i.e. number of grids, per dimension,
  // to which a single particle is assigned.
  const int order = 3;

  // Here the field is given by Σᵢ wᵢ δᴰ(x - xᵢ), where δᴰ ↔ δᴷ / dV,
  // dV =: `vol_cell`.
  const double inv_vol_cell = 1 / this->vol_cell;

  // Reset field values to zero.
  this->initialise_density_field();

  // Perform assignment.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    int ijk[order][3];     // grid index coordinates of covered grid cells
    double win[order][3];  // sampling window
    long long gid = 0;     // flattened grid cell index

    for (int iaxis = 0; iaxis < 3; iaxis++) {
      // Carefully set covered sampling window grid indices.
      double loc_grid = this->params.ngrid[iaxis]
        * particles[pid].pos[iaxis] / this->params.boxsize[iaxis];

      int idx_grid = int(loc_grid);

      if (loc_grid - idx_grid < 0.5) {
        ijk[0][iaxis] = (idx_grid == 0)
          ? this->params.ngrid[iaxis] - 1 : idx_grid - 1;
        ijk[1][iaxis] = idx_grid;
        ijk[2][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
          ? 0 : idx_grid + 1;
      } else {
        ijk[0][iaxis] = idx_grid;
        ijk[1][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
          ? 0 : idx_grid + 1;
        ijk[2][iaxis] = (ijk[1][iaxis] == this->params.ngrid[iaxis] - 1)
          ? 0 : ijk[1][iaxis] + 1;
      }

      // Set sampling window value (up to the 2nd element as `order == 3`).
      double s = loc_grid - idx_grid;

      if (s < 0.5) {
        win[0][iaxis] = 1./2 * (1./2 - s) * (1./2 - s);
        win[1][iaxis] = 3./4 - s * s;
        win[2][iaxis] = 1./2 * (1./2 + s) * (1./2 + s);
      } else {
        s = 1 - s;
        win[0][iaxis] = 1./2 * (1./2 + s) * (1./2 + s);
        win[1][iaxis] = 3./4 - s * s;
        win[2][iaxis] = 1./2 * (1./2 - s) * (1./2 - s);
      }
    }

    for (int iloc = 0; iloc < order; iloc++) {
      for (int jloc = 0; jloc < order; jloc++) {
        for (int kloc = 0; kloc < order; kloc++) {
          gid = this->get_grid_index(
            ijk[iloc][0], ijk[jloc][1], ijk[kloc][2]
          );
          if (0 <= gid && gid < this->params.nmesh) {
OMP_ATOMIC
            this->field[gid][0] += inv_vol_cell
              * weight[pid][0] * win[iloc][0] * win[jloc][1] * win[kloc][2];
OMP_ATOMIC
            this->field[gid][1] += inv_vol_cell
              * weight[pid][1] * win[iloc][0] * win[jloc][1] * win[kloc][2];
          }
        }
      }
    }
  }

  // Perform interlacing if needed.
  if (this->params.interlace == "true") {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int pid = 0; pid < particles.ntotal; pid++) {
      int ijk[order][3];
      double win[order][3];
      long long gid = 0;

      for (int iaxis = 0; iaxis < 3; iaxis++) {
        // Apply a half-grid shift and impose the periodic boundary condition.
        double loc_grid = this->params.ngrid[iaxis]
          * particles[pid].pos[iaxis] / this->params.boxsize[iaxis] + 0.5;

        if (loc_grid > this->params.ngrid[iaxis]) {
          loc_grid -= this->params.ngrid[iaxis];
        }

        int idx_grid = int(loc_grid);

        if (loc_grid - idx_grid < 0.5) {
          ijk[0][iaxis] = (idx_grid == 0)
            ? this->params.ngrid[iaxis] - 1 : idx_grid - 1;
          ijk[1][iaxis] = idx_grid;
          ijk[2][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
            ? 0 : idx_grid + 1;
        } else {
          ijk[0][iaxis] = idx_grid;
          ijk[1][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
            ? 0 : ijk[0][iaxis] + 1;
          ijk[2][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
            ? 0 : ijk[1][iaxis] + 1;
        }

        double s = loc_grid - idx_grid;

        if (s < 0.5) {
          win[0][iaxis] = 1./2 * (1./2 - s) * (1./2 - s);
          win[1][iaxis] = 3./4 - s * s;
          win[2][iaxis] = 1./2 * (1./2 + s) * (1./2 + s);
        } else {
          s = 1 - s;
          win[0][iaxis] = 1./2 * (1./2 + s) * (1./2 + s);
          win[1][iaxis] = 3./4 - s * s;
          win[2][iaxis] = 1./2 * (1./2 - s) * (1./2 - s);
        }
      }

      for (int iloc = 0; iloc < order; iloc++) {
        for (int jloc = 0; jloc < order; jloc++) {
          for (int kloc = 0; kloc < order; kloc++) {
            gid = this->get_grid_index(
              ijk[iloc][0], ijk[jloc][1], ijk[kloc][2]
            );
            if (0 <= gid && gid < this->params.nmesh) {
OMP_ATOMIC
              this->field_s[gid][0] += inv_vol_cell
                * weight[pid][0] * win[iloc][0] * win[jloc][1] * win[kloc][2];
OMP_ATOMIC
              this->field_s[gid][1] += inv_vol_cell
                * weight[pid][1] * win[iloc][0] * win[jloc][1] * win[kloc][2];
            }
          }
        }
      }
    }
  }
}

void MeshField::assign_weighted_field_to_mesh_pcs(
  ParticleCatalogue& particles, fftw_complex* weight
) {
  // Set interpolation order, i.e. number of grids, per dimension,
  // to which a single particle is assigned.
  const int order = 4;

  // Here the field is given by Σᵢ wᵢ δᴰ(x - xᵢ),
  // where δᴰ corresponds to δᴷ / dV, dV =: `vol_cell`.
  const double inv_vol_cell = 1 / this->vol_cell;

  // Reset field values to zero.
  this->initialise_density_field();

  // Perform assignment.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    int ijk[order][3];     // grid index coordinates of covered grid cells
    double win[order][3];  // sampling window
    long long gid = 0;     // flattened grid cell index

    for (int iaxis = 0; iaxis < 3; iaxis++) {
      // Carefully set covered sampling window grid indices.
      double loc_grid = this->params.ngrid[iaxis]
        * particles[pid].pos[iaxis] / this->params.boxsize[iaxis];

      int idx_grid = int(loc_grid);

      ijk[0][iaxis] = (idx_grid == 0)
        ? this->params.ngrid[iaxis] - 1 : idx_grid - 1;
      ijk[1][iaxis] = idx_grid;
      ijk[2][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
        ? 0 : idx_grid + 1;
      ijk[3][iaxis] = (ijk[2][iaxis] == this->params.ngrid[iaxis] - 1)
        ? 0 : ijk[2][iaxis] + 1;

      // Set sampling window value (up to the 2nd element as `order == 3`).
      double s = loc_grid - idx_grid;

      win[0][iaxis] = 1./6 * (1. - s) * (1. - s) * (1. - s);
      win[1][iaxis] = 1./6 * (4. - 6. * s * s + 3. * s * s * s);
      win[2][iaxis] = 1./6 * (
        4. - 6. * (1. - s) * (1. - s) + 3. * (1. - s) * (1. - s) * (1. - s)
      );
      win[3][iaxis] = 1./6 * s * s * s;
    }

    for (int iloc = 0; iloc < order; iloc++) {
      for (int jloc = 0; jloc < order; jloc++) {
        for (int kloc = 0; kloc < order; kloc++) {
          gid = this->get_grid_index(
            ijk[iloc][0], ijk[jloc][1], ijk[kloc][2]
          );
          if (0 <= gid && gid < this->params.nmesh) {
OMP_ATOMIC
            this->field[gid][0] += inv_vol_cell
              * weight[pid][0] * win[iloc][0] * win[jloc][1] * win[kloc][2];
OMP_ATOMIC
            this->field[gid][1] += inv_vol_cell
              * weight[pid][1] * win[iloc][0] * win[jloc][1] * win[kloc][2];
          }
        }
      }
    }
  }

  // Perform interlacing if needed.
  if (this->params.interlace == "true") {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int pid = 0; pid < particles.ntotal; pid++) {
      int ijk[order][3];
      double win[order][3];
      long long gid = 0;

      for (int iaxis = 0; iaxis < 3; iaxis++) {
        // Apply a half-grid shift and impose the periodic boundary condition.
        double loc_grid = this->params.ngrid[iaxis]
          * particles[pid].pos[iaxis] / this->params.boxsize[iaxis] + 0.5;

        if (loc_grid > this->params.ngrid[iaxis]) {
          loc_grid -= this->params.ngrid[iaxis];
        }

        int idx_grid = int(loc_grid);

        ijk[0][iaxis] = (idx_grid == 0)
          ? this->params.ngrid[iaxis] - 1 : idx_grid - 1;
        ijk[1][iaxis] = idx_grid;
        ijk[2][iaxis] = (idx_grid == this->params.ngrid[iaxis] - 1)
          ? 0 : idx_grid + 1;
        ijk[3][iaxis] = (ijk[2][iaxis] == this->params.ngrid[iaxis] - 1)
          ? 0 : ijk[2][iaxis] + 1;
        double s = loc_grid - idx_grid;

        win[0][iaxis] = 1./6 * (1. - s) * (1. - s) * (1. - s);
        win[1][iaxis] = 1./6 * (4. - 6. * s * s + 3. * s * s * s);
        win[2][iaxis] = 1./6 * (
          4. - 6. * (1. - s) * (1. - s) + 3. * (1. - s) * (1. - s) * (1. - s)
        );
        win[3][iaxis] = 1./6 * s * s * s;
      }

      for (int iloc = 0; iloc < order; iloc++) {
        for (int jloc = 0; jloc < order; jloc++) {
          for (int kloc = 0; kloc < order; kloc++) {
            gid = this->get_grid_index(
              ijk[iloc][0], ijk[jloc][1], ijk[kloc][2]
            );
            if (0 <= gid && gid < this->params.nmesh) {
OMP_ATOMIC
              this->field_s[gid][0] += inv_vol_cell
                * weight[pid][0] * win[iloc][0] * win[jloc][1] * win[kloc][2];
OMP_ATOMIC
              this->field_s[gid][1] += inv_vol_cell
                * weight[pid][1] * win[iloc][0] * win[jloc][1] * win[kloc][2];
            }
          }
        }
      }
    }
  }
}

double MeshField::calc_assignment_window_in_fourier(int i, int j, int k) {
  int order = 0;
  if (this->params.assignment == "ngp") {
    order = 1;
  } else
  if (this->params.assignment == "cic") {
    order = 2;
  } else
  if (this->params.assignment == "tsc") {
    order = 3;
  } else
  if (this->params.assignment == "pcs") {
    order = 4;
  } else {
    order = 0;
  }

  i = (i < this->params.ngrid[0]/2) ? i : i - this->params.ngrid[0];
  j = (j < this->params.ngrid[1]/2) ? j : j - this->params.ngrid[1];
  k = (k < this->params.ngrid[2]/2) ? k : k - this->params.ngrid[2];

  double u_x = M_PI * i / double(this->params.ngrid[0]);
  double u_y = M_PI * j / double(this->params.ngrid[1]);
  double u_z = M_PI * k / double(this->params.ngrid[2]);

  // Note sin(u) / u -> 1 as u -> 0.
  double wk_x = (i != 0) ? std::sin(u_x) / u_x : 1.;
  double wk_y = (j != 0) ? std::sin(u_y) / u_y : 1.;
  double wk_z = (k != 0) ? std::sin(u_z) / u_z : 1.;

  double wk = wk_x * wk_y * wk_z;

  return std::pow(wk, order);
}


// -----------------------------------------------------------------------
// Field computations
// -----------------------------------------------------------------------

void MeshField::compute_unweighted_field(ParticleCatalogue& particles) {
  fftw_complex* unit_weight = nullptr;

  unit_weight = fftw_alloc_complex(particles.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles.ntotal);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    unit_weight[pid][0] = 1.;
    unit_weight[pid][1] = 0.;
  }

  this->assign_weighted_field_to_mesh(particles, unit_weight);

  fftw_free(unit_weight); unit_weight = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles.ntotal);
}

void MeshField::compute_unweighted_field_fluctuations_insitu(
  ParticleCatalogue& particles
) {
  this->compute_unweighted_field(particles);

  // Subtract the global mean density to compute fluctuations, i.e. δn.
  double nbar = double(particles.ntotal) / this->vol;

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] -= nbar;
    // this->field[gid][1] -= 0.; (unused)
  }
}

void MeshField::compute_ylm_wgtd_field(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  double alpha, int ell, int m
) {
  fftw_complex* weight_kern = nullptr;

  // Compute the weighted data-source field.
  weight_kern = fftw_alloc_complex(particles_data.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles_data.ntotal);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles_data.ntotal; pid++) {
    double los_[3] = {
      los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    weight_kern[pid][0] = ylm.real() * particles_data[pid].w;
    weight_kern[pid][1] = ylm.imag() * particles_data[pid].w;
  }

  this->assign_weighted_field_to_mesh(particles_data, weight_kern);

  fftw_free(weight_kern); weight_kern = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles_data.ntotal);

  // Compute the weighted random-source field.
  weight_kern = fftw_alloc_complex(particles_rand.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles_rand.ntotal);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles_rand.ntotal; pid++) {
    double los_[3] = {
      los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    weight_kern[pid][0] = ylm.real() * particles_rand[pid].w;
    weight_kern[pid][1] = ylm.imag() * particles_rand[pid].w;
  }

  MeshField field_rand(this->params);
  field_rand.assign_weighted_field_to_mesh(particles_rand, weight_kern);

  fftw_free(weight_kern); weight_kern = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles_rand.ntotal);

  // Subtract to compute fluctuations, i.e. δn_LM.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] -= alpha * field_rand[gid][0];
    this->field[gid][1] -= alpha * field_rand[gid][1];
  }

  if (this->params.interlace == "true") {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field_s[gid][0] -= alpha * field_rand.field_s[gid][0];
      this->field_s[gid][1] -= alpha * field_rand.field_s[gid][1];
    }
  }
}

void MeshField::compute_ylm_wgtd_field(
  ParticleCatalogue& particles, LineOfSight* los,
  double alpha, int ell, int m
) {
  fftw_complex* weight_kern = nullptr;

  // Compute the weighted field.
  weight_kern = fftw_alloc_complex(particles.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles.ntotal);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    double los_[3] = {los[pid].pos[0], los[pid].pos[1], los[pid].pos[2]};

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    weight_kern[pid][0] = ylm.real() * particles[pid].w;
    weight_kern[pid][1] = ylm.imag() * particles[pid].w;
  }

  this->assign_weighted_field_to_mesh(particles, weight_kern);

  fftw_free(weight_kern); weight_kern = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles.ntotal);

  // Apply the normalising alpha contrast.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] *= alpha;
    this->field[gid][1] *= alpha;
  }
}

void MeshField::compute_ylm_wgtd_quad_field(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  double alpha,
  int ell, int m
) {
  fftw_complex* weight_kern = nullptr;

  // Compute the quadratic weighted data-source field.
  weight_kern = fftw_alloc_complex(particles_data.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles_data.ntotal);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles_data.ntotal; pid++) {
    double los_[3] = {
      los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    ylm = std::conj(ylm);  // additional conjugation

    weight_kern[pid][0] = ylm.real() * std::pow(particles_data[pid].w, 2);
    weight_kern[pid][1] = ylm.imag() * std::pow(particles_data[pid].w, 2);
  }

  this->assign_weighted_field_to_mesh(particles_data, weight_kern);

  fftw_free(weight_kern); weight_kern = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles_data.ntotal);

  // Compute the quadratic weighted random-source field.
  weight_kern = fftw_alloc_complex(particles_rand.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles_rand.ntotal);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles_rand.ntotal; pid++) {
    double los_[3] = {
      los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    ylm = std::conj(ylm);  // additional conjugation

    weight_kern[pid][0] = ylm.real() * std::pow(particles_rand[pid].w, 2);
    weight_kern[pid][1] = ylm.imag() * std::pow(particles_rand[pid].w, 2);
  }

  MeshField field_rand(this->params);
  field_rand.assign_weighted_field_to_mesh(particles_rand, weight_kern);

  fftw_free(weight_kern); weight_kern = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles_rand.ntotal);

  // Add to compute quadratic fluctuations, i.e. N_LM.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] += std::pow(alpha, 2) * field_rand[gid][0];
    this->field[gid][1] += std::pow(alpha, 2) * field_rand[gid][1];
  }

  if (this->params.interlace == "true") {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field_s[gid][0] += std::pow(alpha, 2) * field_rand.field_s[gid][0];
      this->field_s[gid][1] += std::pow(alpha, 2) * field_rand.field_s[gid][1];
    }
  }
}

void MeshField::compute_ylm_wgtd_quad_field(
  ParticleCatalogue& particles, LineOfSight* los,
  double alpha, int ell, int m
) {
  fftw_complex* weight_kern = nullptr;

  // Compute the quadratic weighted field.
  weight_kern = fftw_alloc_complex(particles.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles.ntotal);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    double los_[3] = {los[pid].pos[0], los[pid].pos[1], los[pid].pos[2]};

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    ylm = std::conj(ylm);  // conjugation is essential

    weight_kern[pid][0] = ylm.real() * std::pow(particles[pid].w, 2);
    weight_kern[pid][1] = ylm.imag() * std::pow(particles[pid].w, 2);
  }

  this->assign_weighted_field_to_mesh(particles, weight_kern);

  fftw_free(weight_kern); weight_kern = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles.ntotal);

  // Apply mean-density matching normalisation (i.e. alpha contrast)
  // to compute N_LM.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] *= std::pow(alpha, 2);
    this->field[gid][1] *= std::pow(alpha, 2);
  }
}


// -----------------------------------------------------------------------
// Field transforms
// -----------------------------------------------------------------------

void MeshField::fourier_transform() {
  // Apply FFT volume normalisation, where ∫d³x ↔ dV Σᵢ, dV =: `vol_cell`.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] *= this->vol_cell;
    this->field[gid][1] *= this->vol_cell;
  }

  // Perform FFT.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
  fftw_plan transform = fftw_plan_dft_3d(
    this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
    this->field, this->field,
    FFTW_FORWARD, FFTW_ESTIMATE
  );

  fftw_execute(transform);
  fftw_destroy_plan(transform);

  // Interlace with the shadow field.
  if (this->params.interlace == "true") {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int gid = 0; gid < this->params.nmesh; gid++) {
      this->field_s[gid][0] *= this->vol_cell;
      this->field_s[gid][1] *= this->vol_cell;
    }

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
    fftw_plan transform_s = fftw_plan_dft_3d(
      this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
      this->field_s, this->field_s,
      FFTW_FORWARD, FFTW_ESTIMATE
    );

    fftw_execute(transform_s);
    fftw_destroy_plan(transform_s);

#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
    for (int i = 0; i < this->params.ngrid[0]; i++) {
      for (int j = 0; j < this->params.ngrid[1]; j++) {
        for (int k = 0; k < this->params.ngrid[2]; k++) {
          long long idx_grid = this->get_grid_index(i, j, k);

          // Calculate the index vector representing the grid cell.
          double m[3];
          m[0] = (i < this->params.ngrid[0]/2)
            ? double(i) / this->params.ngrid[0]
            : double(i) / this->params.ngrid[0] - 1;
          m[1] = (j < this->params.ngrid[1]/2)
            ? double(j) / this->params.ngrid[1]
            : double(j) / this->params.ngrid[1] - 1;
          m[2] = (k < this->params.ngrid[2]/2)
            ? double(k) / this->params.ngrid[2]
            : double(k) / this->params.ngrid[2] - 1;

          // Multiply by the phase factor from the half-grid shift and
          // add the shadow mesh field contribution.  Note the positive
          // sign of `arg`.
          double arg = M_PI * (m[0] + m[1] + m[2]);

          this->field[idx_grid][0] +=
            std::cos(arg) * this->field_s[idx_grid][0]
            - std::sin(arg) * this->field_s[idx_grid][1]
          ;
          this->field[idx_grid][1] +=
            std::sin(arg) * this->field_s[idx_grid][0]
            + std::cos(arg) * this->field_s[idx_grid][1]
          ;

          this->field[idx_grid][0] /= 2.;
          this->field[idx_grid][1] /= 2.;
        }
      }
    }
  }
}

void MeshField::inv_fourier_transform() {
  // Apply inverse FFT volume normalisation, where ∫d³k/(2π)³ ↔ (1/V) Σᵢ,
  // V =: `vol`.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] /= this->vol;
    this->field[gid][1] /= this->vol;
  }

  // Perform inverse FFT.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
  fftw_plan inv_transform = fftw_plan_dft_3d(
    this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
    this->field, this->field,
    FFTW_BACKWARD, FFTW_ESTIMATE
  );

  fftw_execute(inv_transform);
  fftw_destroy_plan(inv_transform);
}


// -----------------------------------------------------------------------
// Field operations
// -----------------------------------------------------------------------

void MeshField::apply_wide_angle_pow_law_kernel() {
  // CAVEAT: Discretionary choice such that eps_r / r = O(1.e-9).
  const double eps_r = 1.e-6;

#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = this->get_grid_index(i, j, k);

        double rv[3];
        this->get_grid_pos_vector(i, j, k, rv);

        double r_ = trvm::get_vec3d_magnitude(rv);

        if (r_ < eps_r) {
          // this->field[idx_grid][0] *= 0.; (unused)
          // this->field[idx_grid][1] *= 0.; (unused)
        } else {
          this->field[idx_grid][0] *=
            std::pow(r_, - this->params.i_wa - this->params.j_wa);
          this->field[idx_grid][1] *=
            std::pow(r_, - this->params.i_wa - this->params.j_wa);
        }
      }
    }
  }
}

void MeshField::apply_assignment_compensation() {
#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = this->get_grid_index(i, j, k);

        double win = this->calc_assignment_window_in_fourier(i, j, k);

        this->field[idx_grid][0] /= win;
        this->field[idx_grid][1] /= win;
      }
    }
  }
}


// -----------------------------------------------------------------------
// One-point statistics
// -----------------------------------------------------------------------

void MeshField::inv_fourier_transform_ylm_wgtd_field_band_limited(
  MeshField& field_fourier, std::vector< std::complex<double> >& ylm,
  double k_lower, double k_upper,
  double& k_eff, int& nmodes
) {
  // Reset field values to zero.
  this->initialise_density_field();

  // Reset effective wavenumber and wavevector modes.
  k_eff = 0.;
  nmodes = 0;

  // Perform wavevector mode binning in the band.
#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3) reduction(+:k_eff, nmodes)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = this->get_grid_index(i, j, k);

        double kv[3];
        this->get_grid_wavevector(i, j, k, kv);

        double k_ = trvm::get_vec3d_magnitude(kv);

        // Determine the grid cell contribution to the band.
        if (k_lower < k_ && k_ <= k_upper) {
          std::complex<double> fk(
            field_fourier[idx_grid][0], field_fourier[idx_grid][1]
          );

          // Apply assignment compensation.
          double win = this->calc_assignment_window_in_fourier(i, j, k);
          fk /= win;

          // Weight the field.
          this->field[idx_grid][0] = (ylm[idx_grid] * fk).real();
          this->field[idx_grid][1] = (ylm[idx_grid] * fk).imag();

          k_eff += k_;
          nmodes++;
        } else {
          // this->field[idx_grid][0] = 0.; (unused)
          // this->field[idx_grid][1] = 0.; (unused)
        }
      }
    }
  }

  // Perform inverse FFT.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
  fftw_plan inv_transform = fftw_plan_dft_3d(
    this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
    this->field, this->field,
    FFTW_BACKWARD, FFTW_ESTIMATE
  );

  fftw_execute(inv_transform);
  fftw_destroy_plan(inv_transform);

  // Average over wavevector modes in the band.
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    this->field[gid][0] /= double(nmodes);
    this->field[gid][1] /= double(nmodes);
  }

  k_eff /= double(nmodes);
}

void MeshField::inv_fourier_transform_sjl_ylm_wgtd_field(
    MeshField& field_fourier,
    std::vector< std::complex<double> >& ylm,
    trvm::SphericalBesselCalculator& sjl,
    double r
) {
  // Reset field values to zero.
  this->initialise_density_field();

  // Compute the field weighted by the spherical Bessel function and
  // reduced spherical harmonics.
#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = this->get_grid_index(i, j, k);

        double kv[3];
        this->get_grid_wavevector(i, j, k, kv);

        double k_ = trvm::get_vec3d_magnitude(kv);

        // Apply assignment compensation.
        std::complex<double> fk(
          field_fourier[idx_grid][0], field_fourier[idx_grid][1]
        );

        double win = this->calc_assignment_window_in_fourier(i, j, k);
        fk /= win;

        // Weight the field including the volume normalisation,
        // where ∫d³k/(2π)³ ↔ (1/V) Σᵢ, V =: `vol`.
        this->field[idx_grid][0] =
          sjl.eval(k_ * r) * (ylm[idx_grid] * fk).real() / this->vol;
        this->field[idx_grid][1] =
          sjl.eval(k_ * r) * (ylm[idx_grid] * fk).imag() / this->vol;
      }
    }
  }

  // Perform inverse FFT.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
  fftw_plan inv_transform = fftw_plan_dft_3d(
    this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
    this->field, this->field,
    FFTW_BACKWARD, FFTW_ESTIMATE
  );

  fftw_execute(inv_transform);
  fftw_destroy_plan(inv_transform);
}


// -----------------------------------------------------------------------
// Misc
// -----------------------------------------------------------------------

double MeshField::calc_grid_based_powlaw_norm(
  ParticleCatalogue& particles, int order
) {
  // Initialise the weight field.
  fftw_complex* weight = nullptr;

  weight = fftw_alloc_complex(particles.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles.ntotal);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    weight[pid][0] = particles[pid].w;
    weight[pid][1] = 0.;
  }

  // Compute the weighted field.
  this->assign_weighted_field_to_mesh(particles, weight);

  fftw_free(weight); weight = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles.ntotal);

  // Compute normalisation volume integral, where ∫d³x ↔ dV Σᵢ,
  // dV =: `vol_cell`.
  double vol_int = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:vol_int)
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    vol_int += std::pow(this->field[gid][0], order);
  }

  vol_int *= this->vol_cell;

  double norm_factor = 1. / vol_int;

  return norm_factor;
}


// ***********************************************************************
// Field statistics
// ***********************************************************************

// -----------------------------------------------------------------------
// Life cycle
// -----------------------------------------------------------------------

FieldStats::FieldStats(trv::ParameterSet& params){
  this->params = params;

  this->reset_stats();

  // Calculate grid sizes in configuration space.
  this->dr[0] = this->params.boxsize[0] / this->params.ngrid[0];
  this->dr[1] = this->params.boxsize[1] / this->params.ngrid[1];
  this->dr[2] = this->params.boxsize[2] / this->params.ngrid[2];

  // Calculate fundamental wavenumbers in Fourier space.
  this->dk[0] = 2.*M_PI / this->params.boxsize[0];
  this->dk[1] = 2.*M_PI / this->params.boxsize[1];
  this->dk[2] = 2.*M_PI / this->params.boxsize[2];

  // Calculate mesh volume and mesh grid cell volume.
  this->vol = this->params.volume;
  this->vol_cell = this->vol / double(this->params.nmesh);
}

// An empty destructor is redundant but left here for future implementations.
FieldStats::~FieldStats() {}

void FieldStats::reset_stats() {
  std::fill(this->nmodes.begin(), this->nmodes.end(), 0);
  std::fill(this->npairs.begin(), this->npairs.end(), 0);
  std::fill(this->k.begin(), this->k.end(), 0.);
  std::fill(this->r.begin(), this->r.end(), 0.);
  std::fill(this->sn.begin(), this->sn.end(), 0.);
  std::fill(this->pk.begin(), this->pk.end(), 0.);
  std::fill(this->xi.begin(), this->xi.end(), 0.);
}

void FieldStats::resize_stats(int num_bins){
  this->nmodes.resize(num_bins);
  this->npairs.resize(num_bins);
  this->k.resize(num_bins);
  this->r.resize(num_bins);
  this->sn.resize(num_bins);
  this->pk.resize(num_bins);
  this->xi.resize(num_bins);
}


// -----------------------------------------------------------------------
// Utilities
// -----------------------------------------------------------------------

bool FieldStats::if_fields_compatible(
  MeshField& field_a, MeshField& field_b
) {
  bool flag_compatible = true;

  // Check if physical dimensions match.
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    if (
      this->params.boxsize[iaxis] != field_a.params.boxsize[iaxis]
      || this->params.boxsize[iaxis] != field_b.params.boxsize[iaxis]
      || this->params.ngrid[iaxis] != field_a.params.ngrid[iaxis]
      || this->params.ngrid[iaxis] != field_b.params.ngrid[iaxis]
    ) {
      flag_compatible = false;
    }
  }

  // Check if derived physical dimensions match.  This is usually
  // redundant but parameters may have been altered.
  if (
    this->params.nmesh != field_a.params.nmesh
    || this->params.nmesh != field_b.params.nmesh
    || this->params.volume != field_b.params.volume
    || this->params.volume != field_b.params.volume
  ) {
    flag_compatible = false;
  }

  return flag_compatible;
}


// -----------------------------------------------------------------------
// Binned statistics
// -----------------------------------------------------------------------

void FieldStats::compute_ylm_wgtd_2pt_stats_in_fourier(
  MeshField& field_a, MeshField& field_b, std::complex<double> shotnoise_amp,
  int ell, int m, trv::Binning& kbinning
) {
  this->resize_stats(kbinning.num_bins);

  // Check mesh fields compatibility and reuse methods of the first mesh field.
  if (!this->if_fields_compatible(field_a, field_b)) {
    trvs::logger.error(
      "Input mesh fields have incompatible physical properties."
    );
    throw trvs::InvalidDataError(
      "Input mesh fields have incompatible physical properties.\n"
    );
  }

  auto ret_grid_index = [&field_a](int i, int j, int k) {
    return field_a.get_grid_index(i, j, k);
  };

  auto ret_grid_wavevector = [&field_a](int i, int j, int k, double kvec[3]) {
    field_a.get_grid_wavevector(i, j, k, kvec);
  };

  // Perform fine binning.
  // NOTE: Dynamically allocate owing to size.
  // CAVEAT: Discretionary choices such that 0.0 < k < 10.0.
  const int n_sample = 1e6;
  const double dk_sample = 1.e-5;
  if (kbinning.bin_max > n_sample * dk_sample) {
    trvs::logger.warn(
      "Input binning range exceeds sampled range. "
      "Statistics in bins beyond sampled range are uncomputed."
    );
  }

  int* nmodes_sample = new int[n_sample];
  double* k_sample = new double[n_sample];
  double* pk_sample_real = new double[n_sample];
  double* pk_sample_imag = new double[n_sample];
  double* sn_sample_real = new double[n_sample];
  double* sn_sample_imag = new double[n_sample];
  std::complex<double>* pk_sample = new std::complex<double>[n_sample];
  std::complex<double>* sn_sample = new std::complex<double>[n_sample];
  for (int i = 0; i < n_sample; i++) {
    nmodes_sample[i] = 0;
    k_sample[i] = 0.;
    pk_sample_real[i] = 0.;
    pk_sample_imag[i] = 0.;
    sn_sample_real[i] = 0.;
    sn_sample_imag[i] = 0.;
  }

  this->reset_stats();

#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = ret_grid_index(i, j, k);

        double kv[3];
        ret_grid_wavevector(i, j, k, kv);

        double k_ = trvm::get_vec3d_magnitude(kv);

        int idx_k = int(k_ / dk_sample);
        if (0 <= idx_k && idx_k < n_sample) {
          std::complex<double> fa(field_a[idx_grid][0], field_a[idx_grid][1]);
          std::complex<double> fb(field_b[idx_grid][0], field_b[idx_grid][1]);

          std::complex<double> pk_mode = fa * std::conj(fb);
          std::complex<double> sn_mode =
            shotnoise_amp * this->calc_shotnoise_aliasing(i, j, k);

          // Apply grid corrections.
          double win_pk = 1., win_sn = 1.;
          if (this->params.interlace == "true") {
            win_pk = field_a.calc_assignment_window_in_fourier(i, j, k) *
              field_b.calc_assignment_window_in_fourier(i, j, k);
            win_sn = win_pk;
          } else
          if (this->params.interlace == "false") {
#ifndef DBG_FLAG_NOAC
            win_sn = calc_shotnoise_aliasing(i, j, k);
            win_pk = win_sn;
#else   // !DBG_FLAG_NOAC
            win_pk = field_a.calc_assignment_window_in_fourier(i, j, k) *
              field_b.calc_assignment_window_in_fourier(i, j, k);
            win_sn = calc_shotnoise_aliasing(i, j, k);
#endif  // !DBG_FLAG_NOAC
          }

          pk_mode /= win_pk;
          sn_mode /= win_sn;

          // Weight by reduced spherical harmonics.
          std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
            calc_reduced_spherical_harmonic(ell, m, kv);

          pk_mode *= ylm;
          sn_mode *= ylm;

          double pk_mode_real = pk_mode.real();
          double pk_mode_imag = pk_mode.imag();
          double sn_mode_real = sn_mode.real();
          double sn_mode_imag = sn_mode.imag();

          // Add contribution.
OMP_ATOMIC
          nmodes_sample[idx_k]++;
OMP_ATOMIC
          k_sample[idx_k] += k_;
OMP_ATOMIC
          pk_sample_real[idx_k] += pk_mode_real;
OMP_ATOMIC
          pk_sample_imag[idx_k] += pk_mode_imag;
OMP_ATOMIC
          sn_sample_real[idx_k] += sn_mode_real;
OMP_ATOMIC
          sn_sample_imag[idx_k] += sn_mode_imag;
        }
      }
    }
  }

  for (int i = 0; i < n_sample; i++) {
    pk_sample[i] = pk_sample_real[i] + trvm::M_I * pk_sample_imag[i];
    sn_sample[i] = sn_sample_real[i] + trvm::M_I * sn_sample_imag[i];
  }

  // Perform binning.
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    double k_lower = kbinning.bin_edges[ibin];
    double k_upper = kbinning.bin_edges[ibin + 1];
    for (int i = 0; i < n_sample; i++) {
      double k_ = i * dk_sample;
      if (k_lower < k_ && k_ <= k_upper) {
        this->nmodes[ibin] += nmodes_sample[i];
        this->k[ibin] += k_sample[i];
        this->pk[ibin] += pk_sample[i];
        this->sn[ibin] += sn_sample[i];
      }
    }

    if (this->nmodes[ibin] != 0) {
      this->k[ibin] /= double(this->nmodes[ibin]);
      this->pk[ibin] /= double(this->nmodes[ibin]);
      this->sn[ibin] /= double(this->nmodes[ibin]);
    } else {
      this->k[ibin] = kbinning.bin_centres[ibin];
      this->pk[ibin] = 0.;
      this->sn[ibin] = 0.;
    }
  }

  delete[] nmodes_sample;
  delete[] k_sample;
  delete[] pk_sample_real;
  delete[] pk_sample_imag;
  delete[] sn_sample_real;
  delete[] sn_sample_imag;
  delete[] pk_sample;
  delete[] sn_sample;
}

void FieldStats::compute_ylm_wgtd_2pt_stats_in_config(
  MeshField& field_a, MeshField& field_b, std::complex<double> shotnoise_amp,
  int ell, int m, trv::Binning& rbinning
) {
  this->resize_stats(rbinning.num_bins);

  // Check mesh fields compatibility and reuse properties and methods of
  // the first mesh field.
  if (!this->if_fields_compatible(field_a, field_b)) {
    trvs::logger.error(
      "Input mesh fields have incompatible physical properties."
    );
    throw trvs::InvalidDataError(
      "Input mesh fields have incompatible physical properties.\n"
    );
  }

  auto ret_grid_index = [&field_a](int i, int j, int k) {
    return field_a.get_grid_index(i, j, k);
  };

  auto ret_grid_pos_vector = [&field_a](int i, int j, int k, double rvec[3]) {
    field_a.get_grid_pos_vector(i, j, k, rvec);
  };

  // auto ret_grid_wavevector = [&field_a](
  //   int i, int j, int k, double kvec[3]
  // ) {
  //   field_a.get_grid_wavevector(i, j, k, kvec);
  // };

  // Set up 3-d two-point statistics mesh grids (before inverse
  // Fourier transform).
  fftw_complex* twopt_3d = fftw_alloc_complex(this->params.nmesh);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(this->params.nmesh);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    twopt_3d[gid][0] = 0.;
    twopt_3d[gid][1] = 0.;
  }  // likely redundant but safe

  // Compute shot noise--subtracted mode powers on mesh grids.
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = ret_grid_index(i, j, k);

        std::complex<double> fa(field_a[idx_grid][0], field_a[idx_grid][1]);
        std::complex<double> fb(field_b[idx_grid][0], field_b[idx_grid][1]);

        std::complex<double> pk_mode = fa * std::conj(fb);
        std::complex<double> sn_mode =
          shotnoise_amp * this->calc_shotnoise_aliasing(i, j, k);

        // Apply grid corrections.
        double win_pk = 1., win_sn = 1.;
        if (this->params.interlace == "true") {
          win_pk = field_a.calc_assignment_window_in_fourier(i, j, k) *
            field_b.calc_assignment_window_in_fourier(i, j, k);
          win_sn = win_pk;
        } else
        if (this->params.interlace == "false") {
#ifndef DBG_FLAG_NOAC
          win_sn = calc_shotnoise_aliasing(i, j, k);
          win_pk = win_sn;
#else   // !DBG_FLAG_NOAC
          win_pk = field_a.calc_assignment_window_in_fourier(i, j, k) *
            field_b.calc_assignment_window_in_fourier(i, j, k);
          win_sn = calc_shotnoise_aliasing(i, j, k);
#endif  // !DBG_FLAG_NOAC
        }

        pk_mode /= win_pk;
        sn_mode /= win_sn;

        pk_mode -= sn_mode;

        twopt_3d[idx_grid][0] = pk_mode.real() / this->vol;
        twopt_3d[idx_grid][1] = pk_mode.imag() / this->vol;
      }
    }
  }

  // Inverse Fourier transform.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
  fftw_plan inv_transform = fftw_plan_dft_3d(
    this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
    twopt_3d, twopt_3d,
    FFTW_BACKWARD, FFTW_ESTIMATE
  );

  fftw_execute(inv_transform);
  fftw_destroy_plan(inv_transform);

  // Perform fine binning.
  // NOTE: Dynamically allocate owing to size.
  // CAVEAT: Discretionary choices such that 0 < r < 100k.
  const int n_sample = 1e6;
  const double dr_sample = 1.e-1;
  if (rbinning.bin_max > n_sample * dr_sample) {
    trvs::logger.warn(
      "Input binning range exceeds sampled range. "
      "Statistics in bins beyond sampled range are uncomputed."
    );
  }

  int* npairs_sample = new int[n_sample];
  double* r_sample = new double[n_sample];
  double* xi_sample_real = new double[n_sample];
  double* xi_sample_imag = new double[n_sample];
  std::complex<double>* xi_sample = new std::complex<double>[n_sample];
  for (int i = 0; i < n_sample; i++) {
    npairs_sample[i] = 0;
    r_sample[i] = 0.;
    xi_sample_real[i] = 0.;
    xi_sample_imag[i] = 0.;
  }

  this->reset_stats();

#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = ret_grid_index(i, j, k);

        double rv[3];
        ret_grid_pos_vector(i, j, k, rv);

        double r_ = trvm::get_vec3d_magnitude(rv);

        int idx_r = int(r_ / dr_sample);
        if (0 <= idx_r && idx_r < n_sample) {
          std::complex<double> xi_pair(
            twopt_3d[idx_grid][0], twopt_3d[idx_grid][1]
          );

          // Weight by reduced spherical harmonics.
          std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
            calc_reduced_spherical_harmonic(ell, m, rv);

          xi_pair *= ylm;

          double xi_pair_real = xi_pair.real();
          double xi_pair_imag = xi_pair.imag();

          // Add contribution.
OMP_ATOMIC
          npairs_sample[idx_r]++;
OMP_ATOMIC
          r_sample[idx_r] += r_;
OMP_ATOMIC
          xi_sample_real[idx_r] += xi_pair_real;
OMP_ATOMIC
          xi_sample_imag[idx_r] += xi_pair_imag;
        }
      }
    }
  }

  for (int i = 0; i < n_sample; i++) {
    xi_sample[i] = xi_sample_real[i] + trvm::M_I * xi_sample_imag[i];
  }

  // Perform binning.
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    double r_lower = rbinning.bin_edges[ibin];
    double r_upper = rbinning.bin_edges[ibin + 1];
    for (int i = 0; i < n_sample; i++) {
      double r_ = i * dr_sample;
      if (r_lower <= r_ && r_ < r_upper) {
        this->npairs[ibin] += npairs_sample[i];
        this->r[ibin] += r_sample[i];
        this->xi[ibin] += xi_sample[i];
      }
    }

    if (this->npairs[ibin] != 0) {
      this->r[ibin] /= double(this->npairs[ibin]);
      this->xi[ibin] /= double(this->npairs[ibin]);
      this->npairs[ibin] /= 2;  // reality condition
    } else {
      this->r[ibin] = rbinning.bin_centres[ibin];
      this->xi[ibin] = 0.;
    }
  }

  fftw_free(twopt_3d); twopt_3d = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(this->params.nmesh);

  delete[] npairs_sample;
  delete[] r_sample;
  delete[] xi_sample_real;
  delete[] xi_sample_imag;
  delete[] xi_sample;
}

void FieldStats::compute_uncoupled_shotnoise_for_3pcf(
  MeshField& field_a, MeshField& field_b,
  std::vector< std::complex<double> >& ylm_a,
  std::vector< std::complex<double> >& ylm_b,
  std::complex<double> shotnoise_amp,
  trv::Binning& rbinning
) {
  this->resize_stats(rbinning.num_bins);

  // Check mesh fields compatibility and reuse properties and methods of
  // the first mesh field.
  if (!this->if_fields_compatible(field_a, field_b)) {
    trvs::logger.error(
      "Input mesh fields have incompatible physical properties."
    );
    throw trvs::InvalidDataError(
      "Input mesh fields have incompatible physical properties.\n"
    );
  }

  auto ret_grid_index = [&field_a](int i, int j, int k) {
    return field_a.get_grid_index(i, j, k);
  };

  auto ret_grid_pos_vector = [&field_a](int i, int j, int k, double rvec[3]) {
    field_a.get_grid_pos_vector(i, j, k, rvec);
  };

  // auto ret_grid_wavevector = [&field_a](
  //   int i, int j, int k, double kvec[3]
  // ) {
  //   field_a.get_grid_wavevector(i, j, k, kvec);
  // };

  // Set up 3-d two-point statistics mesh grids (before inverse
  // Fourier transform).
  fftw_complex* twopt_3d = fftw_alloc_complex(this->params.nmesh);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(this->params.nmesh);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    twopt_3d[gid][0] = 0.;
    twopt_3d[gid][1] = 0.;
  }  // likely redundant but safe

  // Compute meshed statistics.
#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = ret_grid_index(i, j, k);

        std::complex<double> fa(field_a[idx_grid][0], field_a[idx_grid][1]);
        std::complex<double> fb(field_b[idx_grid][0], field_b[idx_grid][1]);

        std::complex<double> pk_mode = fa * std::conj(fb);
        std::complex<double> sn_mode =
          shotnoise_amp * this->calc_shotnoise_aliasing(i, j, k);

        // Apply grid corrections.
        double win_pk = 1., win_sn = 1.;
        if (this->params.interlace == "true") {
          win_pk = field_a.calc_assignment_window_in_fourier(i, j, k) *
            field_b.calc_assignment_window_in_fourier(i, j, k);
          win_sn = win_pk;
        } else
        if (this->params.interlace == "false") {
#ifndef DBG_FLAG_NOAC
          win_sn = calc_shotnoise_aliasing(i, j, k);
          win_pk = win_sn;
#else   // !DBG_FLAG_NOAC
          win_pk = field_a.calc_assignment_window_in_fourier(i, j, k) *
            field_b.calc_assignment_window_in_fourier(i, j, k);
          win_sn = calc_shotnoise_aliasing(i, j, k);
#endif  // !DBG_FLAG_NOAC
        }

        pk_mode /= win_pk;
        sn_mode /= win_sn;

        pk_mode -= sn_mode;

        twopt_3d[idx_grid][0] = pk_mode.real() / this->vol;
        twopt_3d[idx_grid][1] = pk_mode.imag() / this->vol;
      }
    }
  }

  // Inverse Fourier transform.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
  fftw_plan inv_transform = fftw_plan_dft_3d(
    this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
    twopt_3d, twopt_3d,
    FFTW_BACKWARD, FFTW_ESTIMATE
  );

  fftw_execute(inv_transform);
  fftw_destroy_plan(inv_transform);

  // Perform fine binning.
  // NOTE: Dynamically allocate owing to size.
  // CAVEAT: Discretionary choices such that 0 < r < 100k.
  const int n_sample = 1e5;
  const double dr_sample = 1.;

  int* npairs_sample = new int[n_sample];
  double* r_sample = new double[n_sample];
  double* xi_sample_real = new double[n_sample];
  double* xi_sample_imag = new double[n_sample];
  std::complex<double>* xi_sample = new std::complex<double>[n_sample];
  for (int i = 0; i < n_sample; i++) {
    npairs_sample[i] = 0;
    r_sample[i] = 0.;
    xi_sample_real[i] = 0.;
    xi_sample_imag[i] = 0.;
  }

  this->reset_stats();

#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = ret_grid_index(i, j, k);

        double rv[3];
        ret_grid_pos_vector(i, j, k, rv);

        double r_ = trvm::get_vec3d_magnitude(rv);

        int idx_r = int(r_ / dr_sample);
        if (0 <= idx_r && idx_r < n_sample) {
          std::complex<double> xi_pair(
            twopt_3d[idx_grid][0], twopt_3d[idx_grid][1]
          );

          // Weight by reduced spherical harmonics.
          xi_pair *= ylm_a[idx_grid] * ylm_b[idx_grid];

          double xi_pair_real = xi_pair.real();
          double xi_pair_imag = xi_pair.imag();

          // Add contribution.
OMP_ATOMIC
          npairs_sample[idx_r]++;
OMP_ATOMIC
          r_sample[idx_r] += r_;
OMP_ATOMIC
          xi_sample_real[idx_r] += xi_pair_real;
OMP_ATOMIC
          xi_sample_imag[idx_r] += xi_pair_imag;
        }
      }
    }
  }

  for (int i = 0; i < n_sample; i++) {
    xi_sample[i] = xi_sample_real[i] + trvm::M_I * xi_sample_imag[i];
  }

  // Perform binning.
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    double r_lower = rbinning.bin_edges[ibin];
    double r_upper = rbinning.bin_edges[ibin + 1];
    for (int i = 0; i < n_sample; i++) {
      double r_ = i * dr_sample;
      if (r_lower <= r_ && r_ < r_upper) {
        this->npairs[ibin] += npairs_sample[i];
        this->r[ibin] += r_sample[i];
        this->xi[ibin] += xi_sample[i];
      }
    }

    if (this->npairs[ibin] != 0) {
      this->r[ibin] /= double(this->npairs[ibin]);
      this->xi[ibin] /= double(this->npairs[ibin]);
      this->npairs[ibin] /= 2;  // reality condition
    } else {
      this->r[ibin] = rbinning.bin_centres[ibin];
      this->xi[ibin] = 0.;
    }
  }

  // Apply normalisation factors.
  double norm_factors = 1 / this->vol_cell
    * std::pow(-1, this->params.ell1 + this->params.ell2);

  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    if (this->npairs[ibin] != 0) {
      this->xi[ibin] *= norm_factors / this->npairs[ibin];
    }
  }

  fftw_free(twopt_3d); twopt_3d = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(this->params.nmesh);

  delete[] npairs_sample;
  delete[] r_sample;
  delete[] xi_sample_real;
  delete[] xi_sample_imag;
  delete[] xi_sample;
}

std::complex<double> \
FieldStats::compute_uncoupled_shotnoise_for_bispec_per_bin(
  MeshField& field_a, MeshField& field_b,
  std::vector< std::complex<double> >& ylm_a,
  std::vector< std::complex<double> >& ylm_b,
  trvm::SphericalBesselCalculator& sj_a, trvm::SphericalBesselCalculator& sj_b,
  std::complex<double> shotnoise_amp,
  double k_a, double k_b
) {
  // Check mesh fields compatibility and reuse properties and methods of
  // the first mesh field.
  if (!this->if_fields_compatible(field_a, field_b)) {
    trvs::logger.error(
      "Input mesh fields have incompatible physical properties."
    );
    throw trvs::InvalidDataError(
      "Input mesh fields have incompatible physical properties.\n"
    );
  }

  auto ret_grid_index = [&field_a](int i, int j, int k) {
    return field_a.get_grid_index(i, j, k);
  };

  auto ret_grid_pos_vector = [&field_a](int i, int j, int k, double rvec[3]) {
    field_a.get_grid_pos_vector(i, j, k, rvec);
  };

  // auto ret_grid_wavevector = [&field_a](
  //   int i, int j, int k, double kvec[3]
  // ) {
  //   field_a.get_grid_wavevector(i, j, k, kvec);
  // };

  // Set up 3-d two-point statistics mesh grids (before inverse
  // Fourier transform).
  fftw_complex* twopt_3d = fftw_alloc_complex(this->params.nmesh);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(this->params.nmesh);
  trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int gid = 0; gid < this->params.nmesh; gid++) {
    twopt_3d[gid][0] = 0.;
    twopt_3d[gid][1] = 0.;
  }  // likely redundant but safe

  // Compute meshed statistics.
#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = ret_grid_index(i, j, k);

        std::complex<double> fa(field_a[idx_grid][0], field_a[idx_grid][1]);
        std::complex<double> fb(field_b[idx_grid][0], field_b[idx_grid][1]);

        std::complex<double> pk_mode = fa * std::conj(fb);
        std::complex<double> sn_mode =
          shotnoise_amp * this->calc_shotnoise_aliasing(i, j, k);

        // Apply grid corrections.
        double win_pk = 1., win_sn = 1.;
        if (this->params.interlace == "true") {
          win_pk = field_a.calc_assignment_window_in_fourier(i, j, k) *
            field_b.calc_assignment_window_in_fourier(i, j, k);
          win_sn = win_pk;
        } else
        if (this->params.interlace == "false") {
#ifndef DBG_FLAG_NOAC
          win_sn = calc_shotnoise_aliasing(i, j, k);
          win_pk = win_sn;
#else   // !DBG_FLAG_NOAC
          win_pk = field_a.calc_assignment_window_in_fourier(i, j, k) *
            field_b.calc_assignment_window_in_fourier(i, j, k);
          win_sn = calc_shotnoise_aliasing(i, j, k);
#endif  // !DBG_FLAG_NOAC
        }

        pk_mode /= win_pk;
        sn_mode /= win_sn;

        pk_mode -= sn_mode;

        twopt_3d[idx_grid][0] = pk_mode.real() / this->vol;
        twopt_3d[idx_grid][1] = pk_mode.imag() / this->vol;
      }
    }
  }

  // Inverse Fourier transform.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
  fftw_plan inv_transform = fftw_plan_dft_3d(
    this->params.ngrid[0], this->params.ngrid[1], this->params.ngrid[2],
    twopt_3d, twopt_3d,
    FFTW_BACKWARD, FFTW_ESTIMATE
  );

  fftw_execute(inv_transform);
  fftw_destroy_plan(inv_transform);

  // Weight by spherical Bessel functions and harmonics before summing
  // over the configuration-space grids.
  double S_ij_k_real = 0., S_ij_k_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3) reduction(+:S_ij_k_real, S_ij_k_imag)
#endif  // TRV_USE_OMP
  for (int i = 0; i < this->params.ngrid[0]; i++) {
    for (int j = 0; j < this->params.ngrid[1]; j++) {
      for (int k = 0; k < this->params.ngrid[2]; k++) {
        long long idx_grid = ret_grid_index(i, j, k);

        double rv[3];
        ret_grid_pos_vector(i, j, k, rv);

        double r_ = trvm::get_vec3d_magnitude(rv);

        double ja = sj_a.eval(k_a * r_);
        double jb = sj_b.eval(k_b * r_);

        std::complex<double> S_ij_k_3d(
          twopt_3d[idx_grid][0], twopt_3d[idx_grid][1]
        );

        S_ij_k_3d *= ja * jb * ylm_a[idx_grid] * ylm_b[idx_grid];

        double S_ij_k_3d_real = S_ij_k_3d.real();
        double S_ij_k_3d_imag = S_ij_k_3d.imag();

        S_ij_k_real += S_ij_k_3d_real;
        S_ij_k_imag += S_ij_k_3d_imag;
      }
    }
  }

  std::complex<double> S_ij_k(S_ij_k_real, S_ij_k_imag);

  S_ij_k *= this->vol_cell;

  fftw_free(twopt_3d); twopt_3d = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(this->params.nmesh);

  return S_ij_k;
}


// -----------------------------------------------------------------------
// Sampling corrections
// -----------------------------------------------------------------------

double FieldStats::calc_shotnoise_aliasing(int i, int j, int k) {
  i = (i < this->params.ngrid[0]/2) ? i : i - this->params.ngrid[0];
  j = (j < this->params.ngrid[1]/2) ? j : j - this->params.ngrid[1];
  k = (k < this->params.ngrid[2]/2) ? k : k - this->params.ngrid[2];

  if (this->params.assignment == "ngp") {
    return this->calc_shotnoise_aliasing_ngp(i, j, k);
  }
  if (this->params.assignment == "cic") {
    return this->calc_shotnoise_aliasing_cic(i, j, k);
  }
  if (this->params.assignment == "tsc") {
    return this->calc_shotnoise_aliasing_tsc(i, j, k);
  }
  if (this->params.assignment == "pcs") {
    return this->calc_shotnoise_aliasing_pcs(i, j, k);
  }

  return 1.;  // default
}

double FieldStats::calc_shotnoise_aliasing_ngp(int i, int j, int k) {
  return 1.;
}

double FieldStats::calc_shotnoise_aliasing_cic(int i, int j, int k) {
  double u_x = M_PI * i / double(this->params.ngrid[0]);
  double u_y = M_PI * j / double(this->params.ngrid[1]);
  double u_z = M_PI * k / double(this->params.ngrid[2]);

  double cx2 = (i != 0) ? std::sin(u_x) * std::sin(u_x) : 0.;
  double cy2 = (j != 0) ? std::sin(u_y) * std::sin(u_y) : 0.;
  double cz2 = (k != 0) ? std::sin(u_z) * std::sin(u_z) : 0.;

  return (1. - 2./3. * cx2) * (1. - 2./3. * cy2) * (1. - 2./3. * cz2);
}

double FieldStats::calc_shotnoise_aliasing_tsc(int i, int j, int k) {
  double u_x = M_PI * i / double(this->params.ngrid[0]);
  double u_y = M_PI * j / double(this->params.ngrid[1]);
  double u_z = M_PI * k / double(this->params.ngrid[2]);

  double cx2 = (i != 0) ? std::sin(u_x) * std::sin(u_x) : 0.;
  double cy2 = (j != 0) ? std::sin(u_y) * std::sin(u_y) : 0.;
  double cz2 = (k != 0) ? std::sin(u_z) * std::sin(u_z) : 0.;

  return (1. - cx2 + 2./15. * cx2 * cx2)
    * (1. - cy2 + 2./15. * cy2 * cy2)
    * (1. - cz2 + 2./15. * cz2 * cz2);
}

double FieldStats::calc_shotnoise_aliasing_pcs(int i, int j, int k) {
  double u_x = M_PI * i / double(this->params.ngrid[0]);
  double u_y = M_PI * j / double(this->params.ngrid[1]);
  double u_z = M_PI * k / double(this->params.ngrid[2]);

  double cx2 = (i != 0) ? std::sin(u_x) * std::sin(u_x) : 0.;
  double cy2 = (j != 0) ? std::sin(u_y) * std::sin(u_y) : 0.;
  double cz2 = (k != 0) ? std::sin(u_z) * std::sin(u_z) : 0.;

  return (1. - 4./3. * cx2 + 2./5. * cx2 * cx2 - 4./315. * cx2 * cx2 * cx2)
    * (1. - 4./3. * cy2 + 2./5. * cy2 * cy2 - 4./315. * cy2 * cy2 * cy2)
    * (1. - 4./3. * cz2 + 2./5. * cz2 * cz2 - 4./315. * cz2 * cz2 * cz2);
}

}  // namespace trv
