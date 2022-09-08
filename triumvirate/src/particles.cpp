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
 * @author Mike S Wang (https://github.com/MikeSWang)
 *
 */

#include "particles.hpp"

namespace trv {

/// **********************************************************************
/// Life cycle
/// **********************************************************************

ParticleCatalogue::ParticleCatalogue () {
  /// Set default values.
  /// NOTE: These assignments are not needed but kept for completeness.
  this->pdata = nullptr;
  this->ntotal = 0;
  this->wtotal = 0.;
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    this->pos_min[iaxis] = 0.;
    this->pos_max[iaxis] = 0.;
  }
}

ParticleCatalogue::~ParticleCatalogue() {this->finalise_particles();}

void ParticleCatalogue::initialise_particles(const int num) {
  /// Check the total number of particles.
  if (num <= 0) {
    throw trv::sys::InvalidParameter(
      "[%s WARN] Number of particles is non-positive.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  this->ntotal = num;

  /// Renew particle data.
  delete[] this->pdata; this->pdata = nullptr;

  this->pdata = new ParticleData[this->ntotal];

  trv::sys::gbytesMem += double(this->ntotal)
    * sizeof(struct ParticleData) / BYTES_PER_GBYTES;
}

void ParticleCatalogue::finalise_particles() {
  /// Free particle data.
  if (this->pdata != nullptr) {
    delete[] this->pdata; this->pdata = nullptr;
    trv::sys::gbytesMem -= double(this->ntotal)
      * sizeof(struct ParticleData) / BYTES_PER_GBYTES;
  }
}


/// **********************************************************************
/// Operators & reserved methods
/// **********************************************************************

ParticleCatalogue::ParticleData& ParticleCatalogue::operator[](const int pid) {
  return this->pdata[pid];
}


/// **********************************************************************
/// Data I/O
/// **********************************************************************

int ParticleCatalogue::load_catalogue_file(
  const std::string& catalogue_filepath,
  const std::string& catalogue_columns,
  double volume
) {
  this->source = "extfile:" + catalogue_filepath;

  /// --------------------------------------------------------------------
  /// Columns & fields
  /// --------------------------------------------------------------------

  /// CAVEAT: Hard-coded ordered column names.
  const std::vector<std::string> names_ordered = {
    "x", "y", "z", "nz", "ws", "wc"
  };

  std::istringstream iss(catalogue_columns);
  std::vector<std::string> colnames;
  std::string name;
  while (std::getline(iss, name, ',')) {
    colnames.push_back(name);
  }

  /// CAVEAT: Default -1 index as a flag for unfound column names.
  std::vector<int> name_indices(names_ordered.size(), -1);
  for (int iname = 0; iname < names_ordered.size(); iname++) {
    std::ptrdiff_t col_idx = std::distance(
      colnames.begin(),
      std::find(colnames.begin(), colnames.end(), names_ordered[iname])
    );
    if (0 <= col_idx && col_idx < colnames.size()) {
      name_indices[iname] = col_idx;
    }
  }

  /// Check for the 'nz' column.
  if (name_indices[3] == -1) {
    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s WARN] Catalogue 'nz' field is unavailable, "
        "which may raise errors in some computations (source=%s).\n",
        trv::sys::show_timestamp().c_str(),
        this->source.c_str()
      );
    }
  }

  /// --------------------------------------------------------------------
  /// Data reading
  /// --------------------------------------------------------------------

  std::ifstream fin;

  fin.open(catalogue_filepath.c_str(), std::ios::in);

  if (fin.fail()) {
    fin.close();
    if (trv::sys::currTask == 0) {
      throw trv::sys::IOError(
        "[%s ERRO] Failed to open file '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->source.c_str()
      );
    }
  }

  /// Initialise particle data.
  int num_lines = 0;
  std::string line_str;
  while (std::getline(fin, line_str)) {
    /// Terminate at the end of file.
    if (!fin) {break;}

    /// Skip empty lines or comment lines.
    if (line_str.empty() || line_str[0] == '#') {continue;}

    /// Count the line as valid otherwise.
    num_lines++;
  }

  fin.close();

  this->initialise_particles(num_lines);

  /// Set particle data.
  double nz_box_default;
  if (volume > 0.) {
    nz_box_default = this->ntotal / volume;
  }

  fin.open(catalogue_filepath.c_str(), std::ios::in);

  int idx_line = 0;   // current line number
  double nz, ws, wc;  // placeholder variables
  double entry;       // data entry (per column per row)
  while (std::getline(fin, line_str)) {  // std::string line_str;
    /// Terminate at the end of file.
    if (!fin) {break;}

    /// Skip empty lines or comment lines.
    if (line_str.empty() || line_str[0] == '#') {continue;}

    /// Extract row entries.
    std::vector<double> row;

    std::stringstream ss(
      line_str, std::ios_base::out | std::ios_base::in | std::ios_base::binary
    );
    while (ss >> entry) {row.push_back(entry);}

    /// Add the current line as a particle.
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

  /// --------------------------------------------------------------------
  /// Catalogue properties
  /// --------------------------------------------------------------------

  /// Calculate systematic weight sum.
  this->calc_wtotal();

  /// Calculate the extents of particles.
  this->calc_pos_min_and_max();

  return 0;
}

int ParticleCatalogue::load_particle_data(
  std::vector<double> x, std::vector<double> y, std::vector<double> z,
  std::vector<double> nz, std::vector<double> ws, std::vector<double> wc
) {
  this->source = "extdata";

  /// Check data array sizes.
  int ntotal = x.size();
  if (!(
    ntotal == y.size()
    && ntotal == z.size()
    && ntotal == nz.size()
    && ntotal == ws.size()
    && ntotal == wc.size()
  )) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Inconsistent particle data dimensions (source=%s).\n",
        trv::sys::show_timestamp().c_str(),
        this->source.c_str()
      );
    }
  }

  /// Fill in particle data.
  this->initialise_particles(ntotal);

  for (int pid = 0; pid < ntotal; pid++) {
    this->pdata[pid].pos[0] = x[pid];
    this->pdata[pid].pos[1] = y[pid];
    this->pdata[pid].pos[2] = z[pid];
    this->pdata[pid].nz = nz[pid];
    this->pdata[pid].ws = ws[pid];
    this->pdata[pid].wc = wc[pid];
    this->pdata[pid].w = ws[pid] * wc[pid];
  }

  /// Calculate systematic weight sum.
  this->calc_wtotal();

  /// Calculate the extents of particles.
  this->calc_pos_min_and_max();

  return 0;
}


/// **********************************************************************
/// Catalogue properties
/// **********************************************************************

void ParticleCatalogue::calc_wtotal() {
  if (this->pdata == nullptr) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle data are uninitialised.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  double wtotal = 0.;
  for (int pid = 0; pid < this->ntotal; pid++) {
    wtotal += this->pdata[pid].ws;
  }

  this->wtotal = wtotal;

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s INFO] Catalogue loaded: %d particles with "
      "total systematic weights %.3f (source=%s).\n",
      trv::sys::show_timestamp().c_str(),
      this->ntotal, this->wtotal, this->source.c_str()
    );
  }
}

void ParticleCatalogue::calc_pos_min_and_max() {
  if (this->pdata == nullptr) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle data are uninitialised.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  /// Initialise minimum and maximum values with the 0th particle's.
  double pos_min[3], pos_max[3];
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    pos_min[iaxis] = this->pdata[0].pos[iaxis];
    pos_max[iaxis] = this->pdata[0].pos[iaxis];
  }

  /// Update minimum and maximum values partice by particle.
  for (int pid = 0; pid < this->ntotal; pid++) {
    for (int iaxis = 0; iaxis < 3; iaxis++) {
      if (pos_min[iaxis] > this->pdata[pid].pos[iaxis]) {
        pos_min[iaxis] = this->pdata[pid].pos[iaxis];
      }
      if (pos_max[iaxis] < this->pdata[pid].pos[iaxis]) {
        pos_max[iaxis] = this->pdata[pid].pos[iaxis];
      }
    }
  }

  for (int iaxis = 0; iaxis < 3; iaxis++) {
    this->pos_min[iaxis] = pos_min[iaxis];
    this->pos_max[iaxis] = pos_max[iaxis];
  }

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s INFO] Extents of particle coordinates: "
      "{'x': (%.3f, %.3f), 'y': (%.3f, %.3f), 'z': (%.3f, %.3f)} "
      "(source=%s).\n",
      trv::sys::show_timestamp().c_str(),
      this->pos_min[0], this->pos_max[0],
      this->pos_min[1], this->pos_max[1],
      this->pos_min[2], this->pos_max[2],
      this->source.c_str()
    );
  }
}


/// **********************************************************************
/// Catalogue operations
/// **********************************************************************

void ParticleCatalogue::offset_coords(const double dpos[3]) {
  if (this->pdata == nullptr) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidData(
        "[%s ERRO] Particle data are uninitialised.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  for (int pid = 0; pid < this->ntotal; pid++) {
    for (int iaxis = 0; iaxis < 3; iaxis++) {
      this->pdata[pid].pos[iaxis] -= dpos[iaxis];
    }
  }

  this->calc_pos_min_and_max();
}

void ParticleCatalogue::offset_coords_for_periodicity(const double boxsize[3]) {
  for (int pid = 0; pid < this->ntotal; pid++) {
    for (int iaxis = 0; iaxis < 3; iaxis++) {
      if (this->pdata[pid].pos[iaxis] >= boxsize[iaxis]) {
        this->pdata[pid].pos[iaxis] -= boxsize[iaxis];
      } else
      if (this->pdata[pid].pos[iaxis] < 0.) {
        this->pdata[pid].pos[iaxis] += boxsize[iaxis];
      }
    }
  }

  this->calc_pos_min_and_max();
}

void ParticleCatalogue::centre_in_box(
  ParticleCatalogue& catalogue,
  const double boxsize[3]
) {
  catalogue.calc_pos_min_and_max();  // likely redundant, but safe

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
  catalogue_ref.calc_pos_min_and_max();  // likely redundant, but safe

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
  catalogue.calc_pos_min_and_max();  // likely redundant, but safe

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
  catalogue_ref.calc_pos_min_and_max();  // likely redundant, but safe

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
  catalogue.calc_pos_min_and_max();  // likely redundant, but safe

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
  catalogue_ref.calc_pos_min_and_max();  // likely redundant, but safe

  double dvec[3] = {
    catalogue_ref.pos_min[0] - ngrid_pad[0] * boxsize[0] / double(ngrid[0]),
    catalogue_ref.pos_min[1] - ngrid_pad[1] * boxsize[1] / double(ngrid[1]),
    catalogue_ref.pos_min[2] - ngrid_pad[2] * boxsize[2] / double(ngrid[2])
  };

  catalogue_ref.offset_coords(dvec);
  catalogue.offset_coords(dvec);
}

}  // namespace trv
