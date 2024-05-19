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
 * @file particles.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "particles.hpp"

namespace trvs = trv::sys;

namespace trv {

// ***********************************************************************
// Life cycle
// ***********************************************************************

ParticleCatalogue::ParticleCatalogue(int verbose) {
  if (verbose >= 0) {
    trvs::logger.reset_level(verbose);
  }

  // Set default values (likely redundant but safe).
  this->pdata = nullptr;
  this->ntotal = 0;
  this->wtotal = 0.;
  this->wstotal = 0.;
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    this->pos_min[iaxis] = 0.;
    this->pos_max[iaxis] = 0.;
    this->pos_span[iaxis] = 0.;
  }
}

ParticleCatalogue::~ParticleCatalogue() {this->finalise_particles();}

void ParticleCatalogue::initialise_particles(const int num) {
  // Check the total number of particles.
  if (num <= 0) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Number of particles is non-positive.");
    }
    throw trvs::InvalidParameterError(
      "Number of particles is non-positive.\n"
    );
  }

  this->ntotal = num;

  // Renew particle data.
  this->reset_particles();

  this->pdata = new ParticleData[this->ntotal];
  trvs::gbytesMem += trvs::size_in_gb<struct ParticleData>(this->ntotal);
  trvs::update_maxmem();
}

void ParticleCatalogue::finalise_particles() {
  this->reset_particles();
}

void ParticleCatalogue::reset_particles() {
  // Free particle data.
  if (this->pdata != nullptr) {
    delete[] this->pdata; this->pdata = nullptr;
    trvs::gbytesMem -= trvs::size_in_gb<struct ParticleData>(this->ntotal);
  }
}


// ***********************************************************************
// Operators & reserved methods
// ***********************************************************************

ParticleData& ParticleCatalogue::operator[](const int pid) {
  return this->pdata[pid];
}


// ***********************************************************************
// Data I/O
// ***********************************************************************

int ParticleCatalogue::load_catalogue_file(
  const std::string& catalogue_filepath,
  const std::string& catalogue_columns,
  double volume
) {
  if (!(this->source.empty())) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Catalogue already loaded from another source: %s.", this->source.c_str()
      );
    }
    throw trvs::InvalidDataError(
      "Catalogue already loaded from another source: %s.\n",
      this->source.c_str()
    );
  }
  this->source = "extfile:" + catalogue_filepath;

  // ---------------------------------------------------------------------
  // Columns & fields
  // ---------------------------------------------------------------------

  // CAVEAT: Hard-coded ordered column names.
  const std::vector<std::string> names_ordered = {
    "x", "y", "z", "nz", "ws", "wc"
  };

  std::istringstream iss(catalogue_columns);
  std::vector<std::string> colnames;
  std::string name;
  while (std::getline(iss, name, ',')) {
    colnames.push_back(name);
  }

  // CAVEAT: Default -1 index as a flag for unfound column names.
  std::vector<int> name_indices(names_ordered.size(), -1);
  for (int iname = 0; iname < int(names_ordered.size()); iname++) {
    std::ptrdiff_t col_idx = std::distance(
      colnames.begin(),
      std::find(colnames.begin(), colnames.end(), names_ordered[iname])
    );
    if (0 <= col_idx && col_idx < int(colnames.size())) {
      name_indices[iname] = col_idx;
    }
  }

  // Check for the 'nz' column.
  if (name_indices[3] == -1) {
    if (trvs::currTask == 0) {
      trvs::logger.warn(
        "Catalogue 'nz' field is unavailable and "
        "will be set to the mean density in the bounding box (source=%s).",
        this->source.c_str()
      );
    }
  }

  // ---------------------------------------------------------------------
  // Data reading
  // ---------------------------------------------------------------------

  std::ifstream fin;

  fin.open(catalogue_filepath.c_str(), std::ios::in);

  if (fin.fail()) {
    fin.close();
    if (trvs::currTask == 0) {
      trvs::logger.error("Failed to open file: %s", this->source.c_str());
    }
    throw trvs::IOError("Failed to open file: %s\n", this->source.c_str());
  }

  // Initialise particle data.
  int num_lines = 0;
  std::string line_str;
  while (std::getline(fin, line_str)) {
    // Terminate at the end of file.
    if (!fin) {break;}

    // Skip empty lines or comment lines.
    if (line_str.empty() || line_str[0] == '#') {continue;}

    // Count the line as valid otherwise.
    num_lines++;
  }

  fin.close();

  this->initialise_particles(num_lines);

  // Set particle data.
  double nz_box_default = 0.;
  if (volume > 0.) {
    nz_box_default = this->ntotal / volume;
  }

  fin.open(catalogue_filepath.c_str(), std::ios::in);

  int idx_line = 0;   // current line number
  double nz, ws, wc;  // placeholder variables
  double entry;       // data entry (per column per row)
  while (std::getline(fin, line_str)) {  // std::string line_str;
    // Terminate at the end of file.
    if (!fin) {break;}

    // Skip empty lines or comment lines.
    if (line_str.empty() || line_str[0] == '#') {continue;}

    // Extract row entries.
    std::vector<double> row;

    std::stringstream ss(
      line_str, std::ios_base::out | std::ios_base::in | std::ios_base::binary
    );
    while (ss >> entry) {row.push_back(entry);}

    // Add the current line as a particle.
    this->pdata[idx_line].pos[0] = row[name_indices[0]];  // x
    this->pdata[idx_line].pos[1] = row[name_indices[1]];  // y
    this->pdata[idx_line].pos[2] = row[name_indices[2]];  // z

    if (name_indices[3] != -1) {
      nz = row[name_indices[3]];
    } else {
      nz = nz_box_default;  // default value
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

    this->pdata[idx_line].nz = nz;
    this->pdata[idx_line].ws = ws;
    this->pdata[idx_line].wc = wc;
    this->pdata[idx_line].w = ws * wc;

    idx_line++;
  }

  fin.close();

  // ---------------------------------------------------------------------
  // Catalogue properties
  // ---------------------------------------------------------------------

  // Calculate total weights.
  this->calc_total_weights();

  // Calculate particle extents.
  this->calc_pos_extents();

  return 0;
}

int ParticleCatalogue::load_particle_data(
  std::vector<double> x, std::vector<double> y, std::vector<double> z,
  std::vector<double> nz, std::vector<double> ws, std::vector<double> wc
) {
  this->source = "extdata";

  // Check data array sizes.
  int ntotal = x.size();
  if (!(
    ntotal == int(y.size())
    && ntotal == int(z.size())
    && ntotal == int(nz.size())
    && ntotal == int(ws.size())
    && ntotal == int(wc.size())
  )) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Inconsistent particle data dimensions (source=%s).",
        this->source.c_str()
      );
    }
    throw trvs::InvalidDataError(
      "Inconsistent particle data dimensions (source=%s).\n",
      this->source.c_str()
    );
  }

  // Fill in particle data.
  this->initialise_particles(ntotal);

#ifdef TRV_USE_OMP
#pragma omp parallel for simd
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < ntotal; pid++) {
    this->pdata[pid].pos[0] = x[pid];
    this->pdata[pid].pos[1] = y[pid];
    this->pdata[pid].pos[2] = z[pid];
    this->pdata[pid].nz = nz[pid];
    this->pdata[pid].ws = ws[pid];
    this->pdata[pid].wc = wc[pid];
    this->pdata[pid].w = ws[pid] * wc[pid];
  }

  // Calculate sample weight sum.
  this->calc_total_weights();

  // Calculate the extents of particles.
  this->calc_pos_extents();

  return 0;
}


// ***********************************************************************
// Catalogue properties
// ***********************************************************************

void ParticleCatalogue::calc_total_weights() {
  if (this->pdata == nullptr) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Particle data are uninitialised.");
    }
    throw trvs::InvalidDataError("Particle data are uninitialised.\n");
  }

  double wtotal = 0., wstotal = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for simd reduction(+:wtotal, wstotal)
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < this->ntotal; pid++) {
    wtotal += this->pdata[pid].w;
    wstotal += this->pdata[pid].ws;
  }

  this->wtotal = wtotal;
  this->wstotal = wstotal;

  if (trvs::currTask == 0) {
    trvs::logger.info(
      "Catalogue loaded: "
      "ntotal = %d, wtotal = %.3f, wstotal = %.3f (source=%s).",
      this->ntotal, this->wtotal, this->wstotal, this->source.c_str()
    );
  }
}

void ParticleCatalogue::calc_pos_extents(bool init) {
  if (this->pdata == nullptr) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Particle data are uninitialised.");
    }
    throw trvs::InvalidDataError("Particle data are uninitialised.\n");
  }

  // Initialise minimum and maximum values with the 0th particle's.
  double pos_min[3], pos_max[3];
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    pos_min[iaxis] = this->pdata[0].pos[iaxis];
    pos_max[iaxis] = this->pdata[0].pos[iaxis];
  }

  // Update minimum and maximum values partice by particle.
#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(min:pos_min) reduction(max:pos_max)
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < this->ntotal; pid++) {
    for (int iaxis = 0; iaxis < 3; iaxis++) {
      pos_min[iaxis] = (pos_min[iaxis] < this->pdata[pid].pos[iaxis]) ?
        pos_min[iaxis] : this->pdata[pid].pos[iaxis];
      pos_max[iaxis] = (pos_max[iaxis] > this->pdata[pid].pos[iaxis]) ?
        pos_max[iaxis] : this->pdata[pid].pos[iaxis];
    }
  }

  for (int iaxis = 0; iaxis < 3; iaxis++) {
    this->pos_min[iaxis] = pos_min[iaxis];
    this->pos_max[iaxis] = pos_max[iaxis];
    this->pos_span[iaxis] = pos_max[iaxis] - pos_min[iaxis];
  }

  if (trvs::currTask == 0 && init) {
    trvs::logger.info(
      "Extents of particle coordinates: "
      "{'x': (%.3f, %.3f | %.3f),"
      " 'y': (%.3f, %.3f | %.3f),"
      " 'z': (%.3f, %.3f | %.3f)} "
      "(source=%s).",
      this->pos_min[0], this->pos_max[0], this->pos_span[0],
      this->pos_min[1], this->pos_max[1], this->pos_span[1],
      this->pos_min[2], this->pos_max[2], this->pos_span[2],
      this->source.c_str()
    );
  }
}


// ***********************************************************************
// Catalogue operations
// ***********************************************************************

void ParticleCatalogue::offset_coords(const double dpos[3]) {
  if (this->pdata == nullptr) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Particle data are uninitialised.");
    }
    throw trvs::InvalidDataError("Particle data are uninitialised.\n");
  }

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < this->ntotal; pid++) {
    for (int iaxis = 0; iaxis < 3; iaxis++) {
      this->pdata[pid].pos[iaxis] -= dpos[iaxis];
    }
  }

  this->calc_pos_extents();
}

void ParticleCatalogue::\
offset_coords_for_periodicity(const double boxsize[3]) {
#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < this->ntotal; pid++) {
    for (int iaxis = 0; iaxis < 3; iaxis++) {
      if (this->pdata[pid].pos[iaxis] >= boxsize[iaxis]) {
        this->pdata[pid].pos[iaxis] = std::fmod(
          this->pdata[pid].pos[iaxis], boxsize[iaxis]
        );
      } else
      if (this->pdata[pid].pos[iaxis] < 0.) {
        this->pdata[pid].pos[iaxis] = std::fmod(
          this->pdata[pid].pos[iaxis], boxsize[iaxis]
        ) + boxsize[iaxis];
      }
    }
  }

  this->calc_pos_extents();
}

void ParticleCatalogue::centre_in_box(
  ParticleCatalogue& catalogue,
  const double boxsize[3]
) {
  catalogue.calc_pos_extents(false);  // likely redundant but safe
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    if (catalogue.pos_span[iaxis] > boxsize[iaxis]) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Catalogue extent exceeds the box size along axis %d: "
          "span = %.3f, boxsize = %.3f (source=%s). "
          "Some partcles may lie outside the box after centring.",
          iaxis,
          catalogue.pos_span[iaxis], boxsize[iaxis],
          catalogue.source.c_str()
        );
      }
    }
  }

  double xmid = (catalogue.pos_min[0] + catalogue.pos_max[0]) / 2.;
  double ymid = (catalogue.pos_min[1] + catalogue.pos_max[1]) / 2.;
  double zmid = (catalogue.pos_min[2] + catalogue.pos_max[2]) / 2.;

  double dvec[3] = {
    xmid - boxsize[0]/2., ymid - boxsize[1]/2., zmid - boxsize[2]/2.
  };

  catalogue.offset_coords(dvec);
}

void ParticleCatalogue::centre_in_box(
  ParticleCatalogue& catalogue, ParticleCatalogue& catalogue_ref,
  const double boxsize[3]
) {
  catalogue.calc_pos_extents(false);  // likely redundant but safe
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    if (catalogue.pos_span[iaxis] > boxsize[iaxis]) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Catalogue extent exceeds the box size along axis %d: "
          "span = %.3f, boxsize = %.3f (source=%s). "
          "Some partcles may lie outside the box after centring.",
          iaxis,
          catalogue.pos_span[iaxis], boxsize[iaxis],
          catalogue.source.c_str()
        );
      }
    }
  }

  catalogue_ref.calc_pos_extents(false);  // likely redundant but safe
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    if (catalogue_ref.pos_span[iaxis] > boxsize[iaxis]) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Catalogue extent exceeds the box size along axis %d: "
          "span = %.3f, boxsize = %.3f (source=%s). "
          "Some partcles may lie outside the box after centring.",
          iaxis,
          catalogue_ref.pos_span[iaxis], boxsize[iaxis],
          catalogue_ref.source.c_str()
        );
      }
    }
  }

  double xmid = (catalogue_ref.pos_min[0] + catalogue_ref.pos_max[0]) / 2.;
  double ymid = (catalogue_ref.pos_min[1] + catalogue_ref.pos_max[1]) / 2.;
  double zmid = (catalogue_ref.pos_min[2] + catalogue_ref.pos_max[2]) / 2.;

  double dvec[3] = {
    xmid - boxsize[0]/2., ymid - boxsize[1]/2., zmid - boxsize[2]/2.
  };

  catalogue_ref.offset_coords(dvec);
  catalogue.offset_coords(dvec);
}

void ParticleCatalogue::pad_in_box(
  ParticleCatalogue& catalogue,
  const double boxsize[3], const double boxsize_pad[3]
) {
  catalogue.calc_pos_extents(false);  // likely redundant but safe
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    if (catalogue.pos_span[iaxis] > boxsize[iaxis]) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Catalogue extent exceeds the box size along axis %d: "
          "span = %.3f, boxsize = %.3f (source=%s). "
          "Some partcles may lie outside the box after padding.",
          iaxis,
          catalogue.pos_span[iaxis], boxsize[iaxis],
          catalogue.source.c_str()
        );
      }
    }
  }

  double dvec[3] = {
    catalogue.pos_min[0] - boxsize_pad[0] * boxsize[0],
    catalogue.pos_min[1] - boxsize_pad[1] * boxsize[1],
    catalogue.pos_min[2] - boxsize_pad[2] * boxsize[2]
  };

  catalogue.offset_coords(dvec);
}

void ParticleCatalogue::pad_in_box(
    ParticleCatalogue& catalogue, ParticleCatalogue& catalogue_ref,
    const double boxsize[3], const double boxsize_pad[3]
) {
  catalogue.calc_pos_extents(false);  // likely redundant but safe
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    if (catalogue.pos_span[iaxis] > boxsize[iaxis]) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Catalogue extent exceeds the box size along axis %d: "
          "span = %.3f, boxsize = %.3f (source=%s). "
          "Some partcles may lie outside the box after padding.",
          iaxis,
          catalogue.pos_span[iaxis], boxsize[iaxis],
          catalogue.source.c_str()
        );
      }
    }
  }

  catalogue_ref.calc_pos_extents(false);  // likely redundant but safe
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    if (catalogue_ref.pos_span[iaxis] > boxsize[iaxis]) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Catalogue extent exceeds the box size along axis %d: "
          "span = %.3f, boxsize = %.3f (source=%s). "
          "Some partcles may lie outside the box after padding.",
          iaxis,
          catalogue_ref.pos_span[iaxis], boxsize[iaxis],
          catalogue_ref.source.c_str()
        );
      }
    }
  }

  double dvec[3] = {
    catalogue_ref.pos_min[0] - boxsize_pad[0] * boxsize[0],
    catalogue_ref.pos_min[1] - boxsize_pad[1] * boxsize[1],
    catalogue_ref.pos_min[2] - boxsize_pad[2] * boxsize[2]
  };

  catalogue_ref.offset_coords(dvec);
  catalogue.offset_coords(dvec);
}

void ParticleCatalogue::pad_grids(
  ParticleCatalogue& catalogue,
  const double boxsize[3], const int ngrid[3], const double ngrid_pad[3]
) {
  catalogue.calc_pos_extents(false);  // likely redundant but safe

  double dvec[3] = {
    catalogue.pos_min[0], catalogue.pos_min[1], catalogue.pos_min[2]
  };

  dvec[0] -= ngrid_pad[0] * boxsize[0] / double(ngrid[0]);
  dvec[1] -= ngrid_pad[1] * boxsize[1] / double(ngrid[1]);
  dvec[2] -= ngrid_pad[2] * boxsize[2] / double(ngrid[2]);

  catalogue.offset_coords(dvec);
}

void ParticleCatalogue::pad_grids(
    ParticleCatalogue& catalogue, ParticleCatalogue& catalogue_ref,
    const double boxsize[3], const int ngrid[3], const double ngrid_pad[3]
) {
  catalogue_ref.calc_pos_extents(false);  // likely redundant but safe

  double dvec[3] = {
    catalogue_ref.pos_min[0] - ngrid_pad[0] * boxsize[0] / double(ngrid[0]),
    catalogue_ref.pos_min[1] - ngrid_pad[1] * boxsize[1] / double(ngrid[1]),
    catalogue_ref.pos_min[2] - ngrid_pad[2] * boxsize[2] / double(ngrid[2])
  };

  catalogue_ref.offset_coords(dvec);
  catalogue.offset_coords(dvec);
}

}  // namespace trv
