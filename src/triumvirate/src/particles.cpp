// Triumvirate: Three-Point Clustering Measurements in LSS
//
// Copyright (C) 2023 Mike S Wang & Naonori S Sugiyama [GPL-3.0-or-later]
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
 *          Naonori S Sugiyama (https://github.com/naonori)
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
}

ParticleCatalogue::~ParticleCatalogue() {this->finalise_particles();}

void ParticleCatalogue::initialise_particles(const int num) {
  // Check the total number of particles.
  if (num <= 0) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Number of particles is non-positive.");
    }
    throw trvs::InvalidParameterError(
      "Number of particles is non-positive."
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

PURE ParticleData& ParticleCatalogue::operator[](const int pid) {
  return this->pdata[pid];
}


// ***********************************************************************
// Data I/O
// ***********************************************************************

int ParticleCatalogue::load_catalogue_file(
  const std::string& catalogue_filepath,
  const std::string& catalogue_columns,
  const std::string& catalogue_dataset,
  double volume
) {
  if (!(this->source.empty())) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Catalogue already loaded from another source: %s.", this->source.c_str()
      );
    }
    throw trvs::InvalidDataError(
      "Catalogue already loaded from another source: %s.",
      this->source.c_str()
    );
  }
  this->source = "extfile:" + catalogue_filepath;

  // Check for HDF5 file extension.
  bool is_hdf5 = (
    trvs::has_extension(catalogue_filepath, ".h5")
    || trvs::has_extension(catalogue_filepath, ".hdf5")
  );
#ifndef TRV_USE_H5
  if (is_hdf5) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "HDF5 file format is not supported in this build: %s",
        catalogue_filepath.c_str()
      );
    }
    throw trvs::InvalidDataError(
      "HDF5 file format is not supported in this build: %s",
      catalogue_filepath.c_str()
    );
  }
#endif  // !TRV_USE_H5

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
  for (std::size_t iname = 0; iname < names_ordered.size(); iname++) {
    std::size_t col_idx = static_cast<std::size_t>(std::distance(
      colnames.begin(),
      std::find(colnames.begin(), colnames.end(), names_ordered[iname])
    ));
    if (col_idx < colnames.size()) {
      name_indices[iname] = col_idx;
    }
  }

  // ---------------------------------------------------------------------
  // Data reading
  // ---------------------------------------------------------------------

  std::vector<std::string> catalogue_subfilepaths = trvs::split_string(
    catalogue_filepath, trvs::fn_delimiter
  );

  if (is_hdf5) {
#ifdef TRV_USE_H5
    try {
      int nentry = 0;
      std::vector<std::vector<double>> ctlg_data;
      for (const auto& ctlg_subfilepath : catalogue_subfilepaths) {
        HighFive::File ctlg_subfile(
          ctlg_subfilepath, HighFive::File::ReadOnly
        );

        // Get dataset (assumed to be the first if unspecified)
        // from the file.
        std::vector<std::string> obj_names = ctlg_subfile.listObjectNames();
        if (obj_names.empty()) {
          if (trvs::currTask == 0) {
            trvs::logger.error(
              "No objects found in HDF5 file: %s",
              catalogue_filepath.c_str()
            );
          }
          throw trvs::InvalidDataError(
            "No objects found in HDF5 file: %s",
            catalogue_filepath.c_str()
          );
        }

        std::string ctlg_dset_name = catalogue_dataset;
        if (ctlg_dset_name.empty()) {
          for (const auto& obj_name : obj_names) {
            HighFive::ObjectType obj_type =
              ctlg_subfile.getObjectType(obj_name);
            if (obj_type == HighFive::ObjectType::Dataset) {
              ctlg_dset_name = obj_name;
              trvs::logger.info(
                "Catalogue dataset name inferred from HDF5 file: %s",
                ctlg_dset_name.c_str()
              );
              break;
            }
          }
          if (ctlg_dset_name.empty()) {
            if (trvs::currTask == 0) {
              trvs::logger.error(
                "No datasets found in HDF5 file: %s",
                catalogue_filepath.c_str()
              );
            }
            throw trvs::InvalidDataError(
              "No datasets found in HDF5 file: %s",
              catalogue_filepath.c_str()
            );
          }
        }

        HighFive::DataSet ctlg_dset = ctlg_subfile.getDataSet(ctlg_dset_name);

        // Get dataset data type, which must be float or double and the same
        // for all columns.
        HighFive::DataType ctlg_dtype = ctlg_dset.getDataType();
        if (
          ctlg_dtype != HighFive::create_datatype<float>()
          && ctlg_dtype != HighFive::create_datatype<double>()
        ) {
          if (trvs::currTask == 0) {
            trvs::logger.error(
              "Unsupported or mixed data type in HDF5 file: %s",
              catalogue_filepath.c_str()
            );
          }
          throw trvs::InvalidDataError(
            "Unsupported or mixed data type in HDF5 file: %s",
            catalogue_filepath.c_str()
          );
        }

        // Get dataset column names and indices.
        std::vector<std::string> ctlg_colnames;
        if (colnames.empty()) {
          std::vector<std::string> dset_attrs = ctlg_dset.listAttributeNames();
          if (dset_attrs.empty()) {
            if (trvs::currTask == 0) {
              trvs::logger.error(
                "Catalogue column names are not specified, and "
                "no attributes found in the dataset: %s",
                ctlg_dset_name.c_str()
              );
            }
            throw trvs::InvalidDataError(
              "Catalogue column names are not specified, and "
              "no attributes found in the dataset: %s",
              ctlg_dset_name.c_str()
            );
          }

          ctlg_dset.getAttribute(dset_attrs[0]).read(ctlg_colnames);
          if (trvs::currTask == 0) {
            trvs::logger.info(
              "Catalogue column names inferred from HDF5 file: %s",
              catalogue_filepath.c_str()
            );
          }
          for (std::size_t iname = 0; iname < names_ordered.size(); iname++) {
            std::size_t col_idx = static_cast<std::size_t>(std::distance(
              ctlg_colnames.begin(),
              std::find(
                ctlg_colnames.begin(), ctlg_colnames.end(), names_ordered[iname]
              )
            ));
            if (col_idx < ctlg_colnames.size()) {
              name_indices[iname] = col_idx;
            } else {
              name_indices[iname] = -1;
            }
          }
        } else
        if (colnames.size() == 1 and colnames[0].find("attr::") == 0) {
          std::string attr_name = colnames[0].substr(6);
          ctlg_dset.getAttribute(attr_name).read(ctlg_colnames);
          if (trvs::currTask == 0) {
            trvs::logger.info(
              "Catalogue column names inferred from dataset attribute: %s",
              attr_name.c_str()
            );
          }
          for (std::size_t iname = 0; iname < names_ordered.size(); iname++) {
            std::size_t col_idx = static_cast<std::size_t>(std::distance(
              ctlg_colnames.begin(),
              std::find(
                ctlg_colnames.begin(), ctlg_colnames.end(), names_ordered[iname]
              )
            ));
            if (col_idx < ctlg_colnames.size()) {
              name_indices[iname] = col_idx;
            } else {
              name_indices[iname] = -1;
            }
          }
        } else {
          ctlg_colnames = colnames;
        }

        // Get dataset dimensions.
        std::vector<std::size_t> ctlg_dims = ctlg_dset.getDimensions();
        int num_rows = static_cast<int>(ctlg_dims[0]);

        // Get dataset data.
        std::vector<std::vector<double>> ctlg_subdata(num_rows);
        ctlg_dset.read(ctlg_subdata);

        ctlg_data.insert(
          ctlg_data.end(), ctlg_subdata.begin(), ctlg_subdata.end()
        );
        nentry += num_rows;
      }

      // Initialise particle data.
      this->initialise_particles(nentry);

      // Set particle data.
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

      double nz_box_default = 0.;
      if (volume > 0.) {
        nz_box_default = this->ntotal / volume;
      }

      int idx_row = 0;    // current row number
      double nz, ws, wc;  // placeholder variables
      for (const auto& row : ctlg_data) {
        this->pdata[idx_row].pos[0] = row[name_indices[0]];  // x
        this->pdata[idx_row].pos[1] = row[name_indices[1]];  // y
        this->pdata[idx_row].pos[2] = row[name_indices[2]];  // z

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

        this->pdata[idx_row].nz = nz;
        this->pdata[idx_row].ws = ws;
        this->pdata[idx_row].wc = wc;
        this->pdata[idx_row].w = ws * wc;

        idx_row++;
      }
    } catch (const HighFive::Exception& e) {
      if (trvs::currTask == 0) {
        trvs::logger.error("[HighFive] %s", e.what());
      }
      throw trvs::IOError("[HighFive] %s", e.what());
    }
#endif  // TRV_USE_H5
  } else {
    int nentry = 0;
    for (const auto& ctlg_subfilepath : catalogue_subfilepaths) {
      std::ifstream fin;

      fin.open(ctlg_subfilepath.c_str(), std::ios::in);

      if (fin.fail()) {
        fin.close();
        if (trvs::currTask == 0) {
          trvs::logger.error(
            "Failed to open file: %s", ctlg_subfilepath.c_str()
          );
        }
        throw trvs::IOError(
          "Failed to open file: %s", ctlg_subfilepath.c_str()
        );
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

      nentry += num_lines;
    }

    this->initialise_particles(nentry);

    // Set particle data.
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

    double nz_box_default = 0.;
    if (volume > 0.) {
      nz_box_default = this->ntotal / volume;
    }

    int idx_entry = 0;     // current entry number
    std::string entry_str;  // current line string
    for (const auto& ctlg_subfilepath : catalogue_subfilepaths) {
      std::ifstream fin;
      fin.open(ctlg_subfilepath.c_str(), std::ios::in);

      double nz, ws, wc;  // placeholder variables
      double entry;       // data entry (per column per row)
      while (std::getline(fin, entry_str)) {
        // Terminate at the end of file.
        if (!fin) {break;}

        // Skip empty lines or comment lines.
        if (entry_str.empty() || entry_str[0] == '#') {continue;}

        // Extract row entries.
        std::vector<double> row;

        std::stringstream ss(
          entry_str,
          std::ios_base::out | std::ios_base::in | std::ios_base::binary
        );
        while (ss >> entry) {row.push_back(entry);}

        // Add the current line as a particle.
        this->pdata[idx_entry].pos[0] = row[name_indices[0]];  // x
        this->pdata[idx_entry].pos[1] = row[name_indices[1]];  // y
        this->pdata[idx_entry].pos[2] = row[name_indices[2]];  // z

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

        this->pdata[idx_entry].nz = nz;
        this->pdata[idx_entry].ws = ws;
        this->pdata[idx_entry].wc = wc;
        this->pdata[idx_entry].w = ws * wc;

        idx_entry++;
      }

      fin.close();
    }
  }

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
  std::size_t ntotal_t = x.size();
  int ntotal = static_cast<int>(ntotal_t);
  if (!(
    ntotal_t == y.size()
    && ntotal_t == z.size()
    && ntotal_t == nz.size()
    && ntotal_t == ws.size()
    && ntotal_t == wc.size()
  )) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Inconsistent particle data dimensions (source=%s).",
        this->source.c_str()
      );
    }
    throw trvs::InvalidDataError(
      "Inconsistent particle data dimensions (source=%s).",
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
    throw trvs::InvalidDataError("Particle data are uninitialised.");
  }

  double wtotal = 0., wstotal = 0.;

#ifdef TRV_USE_OMP
#if defined(__GNUC__) && !defined(__clang__)
#pragma omp parallel for simd reduction(+:wtotal, wstotal)
#else  // !__GNUC__ || __clang__
#pragma omp parallel for reduction(+:wtotal, wstotal)
#endif  // __GNUC__ && !__clang__
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
    throw trvs::InvalidDataError("Particle data are uninitialised.");
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
    throw trvs::InvalidDataError("Particle data are uninitialised.");
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
          "Some particles may lie outside the box after centring.",
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
          "Some particles may lie outside the box after centring.",
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
          "Some particles may lie outside the box after centring.",
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
          "Some particles may lie outside the box after padding.",
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
          "Some particles may lie outside the box after padding.",
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
          "Some particles may lie outside the box after padding.",
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
