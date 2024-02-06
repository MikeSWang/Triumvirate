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
 * @file dataobjs.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "dataobjs.hpp"

namespace trvs = trv::sys;

namespace trv {

// ***********************************************************************
// Binning schemes
// ***********************************************************************

Binning::Binning(std::string space, std::string scheme) {
  this->space = space;
  this->scheme = scheme;
}

Binning::Binning(trv::ParameterSet& params) : Binning::Binning(
  params.space, params.binning
) {
  this->bin_min = params.bin_min;
  this->bin_max = params.bin_max;
  this->num_bins = params.num_bins;

  // Change default padding to grid scale.
  this->dbin_pad_config = (1 + 5.e-3)
    * *std::max_element(params.boxsize, params.boxsize + 3)
    / *std::min_element(params.ngrid, params.ngrid + 3);
  this->dbin_pad_fourier = (1 + 5.e-3)
    * (2 * M_PI)
    / *std::max_element(params.boxsize, params.boxsize + 3);
}

void Binning::set_bins(double coord_min, double coord_max, int nbin) {
  if (coord_min < 0.) {
    throw trvs::InvalidParameterError(
      "Bin range must be non-negative."
    );
  }
  if (nbin <= 0) {
    throw trvs::InvalidParameterError(
      "Number of bins must be positive."
    );
  }

  this->bin_min = coord_min;
  this->bin_max = coord_max;
  this->num_bins = nbin;

  this->set_bins();
}

void Binning::set_bins() {
  // Reset vector attributes.
  this->bin_edges.clear();
  this->bin_centres.clear();
  this->bin_widths.clear();

  this->compute_binning();
}

void Binning::set_bins(double boxsize_max, int ngrid_min) {
  this->scheme = "lin";

  this->bin_min = 0.;

  double bin_width = 0.;
  if (this->space == "config") {
    bin_width = boxsize_max / double(ngrid_min);
    this->bin_max = boxsize_max / 2.;
  } else
  if (this->space == "fourier") {
    bin_width = 2*M_PI / boxsize_max;
    this->bin_max = M_PI * double(ngrid_min) / boxsize_max;
  }
  this->bin_max += bin_width / 2.;

  this->num_bins = ngrid_min / 2;

  // Reset vector attributes.
  this->bin_edges.clear();
  this->bin_centres.clear();
  this->bin_widths.clear();

  this->compute_binning();
}

void Binning::set_bins(std::vector<double> bin_edges) {
  this->bin_min = bin_edges.front();
  this->bin_max = bin_edges.back();
  this->num_bins = bin_edges.size() - 1;

  // Reset attributes.
  this->bin_edges.clear();
  this->bin_centres.clear();
  this->bin_widths.clear();
  this->scheme = "custom";

  this->bin_edges = bin_edges;

  for (int ibin = 0; ibin < this->num_bins; ibin++) {
    double centre = (this->bin_edges[ibin] + this->bin_edges[ibin + 1]) / 2.;
    double width = this->bin_edges[ibin + 1] - this->bin_edges[ibin];

    this->bin_centres.push_back(centre);
    this->bin_widths.push_back(width);
  }
}

void Binning::compute_binning() {
  // Set up padding parameters.
  double dbin_pad;
  if (this->space == "fourier") {
    dbin_pad = this->dbin_pad_fourier;
  } else
  if (this->space == "config") {
    dbin_pad = this->dbin_pad_config;
  }

  // Implement binning scheme.

  // ---------------------------------------------------------------------
  // Linear binning
  // ---------------------------------------------------------------------
  if (scheme == "lin") {
    double dbin = (this->bin_max - this->bin_min) / double(this->num_bins);

    for (int ibin = 0; ibin < this->num_bins; ibin++) {
      double edge_left = this->bin_min + dbin * ibin;
      double centre = edge_left + dbin / 2.;

      this->bin_edges.push_back(edge_left);
      this->bin_centres.push_back(centre);
      this->bin_widths.push_back(dbin);
    }
    this->bin_edges.push_back(this->bin_max);
  } else
  // ---------------------------------------------------------------------
  // Logarithmic binning
  // ---------------------------------------------------------------------
  if (scheme == "log")  {
    if (this->bin_min == 0.) {
      throw trvs::InvalidParameterError(
        "Cannot use logarithmic binning when the lowest edge is zero."
      );
    }

    double dlnbin = (std::log(this->bin_max) - std::log(this->bin_min))
      / double(this->num_bins);

    for (int ibin = 0; ibin < this->num_bins; ibin++) {
      double edge_left = this->bin_min * std::exp(dlnbin * ibin);
      double edge_right = this->bin_min * std::exp(dlnbin * (ibin + 1));
      double centre = (edge_left + edge_right) / 2.;

      this->bin_edges.push_back(edge_left);
      this->bin_centres.push_back(centre);
      this->bin_widths.push_back(edge_right - edge_left);
    }
    this->bin_edges.push_back(this->bin_max);
  } else
  // ---------------------------------------------------------------------
  // Padded linear binning
  // ---------------------------------------------------------------------
  if (scheme == "linpad") {
    for (int ibin = 0; ibin < this->nbin_pad; ibin++) {
      double edge_left = dbin_pad * ibin;
      double centre = edge_left + dbin_pad / 2.;

      this->bin_edges.push_back(edge_left);
      this->bin_centres.push_back(centre);
      this->bin_widths.push_back(dbin_pad);
    }

    double bin_min = dbin_pad * this->nbin_pad;

    double dbin = (this->bin_max - bin_min)
      / double(this->num_bins - this->nbin_pad);

    for (int ibin = this->nbin_pad; ibin < this->num_bins; ibin++) {
      double edge_left = bin_min + dbin * (ibin - this->nbin_pad);
      double centre = edge_left + dbin / 2.;

      this->bin_edges.push_back(edge_left);
      this->bin_centres.push_back(centre);
      this->bin_widths.push_back(dbin);
    }
    this->bin_edges.push_back(this->bin_max);
  } else
  // ---------------------------------------------------------------------
  // Padded logarithmic binning
  // ---------------------------------------------------------------------
  if (scheme == "logpad") {
    for (int ibin = 0; ibin < this->nbin_pad; ibin++) {
      double edge_left = dbin_pad * ibin;
      double centre = edge_left + dbin_pad / 2.;

      this->bin_edges.push_back(edge_left);
      this->bin_centres.push_back(centre);
      this->bin_widths.push_back(dbin_pad);
    }

    double bin_min = dbin_pad * this->nbin_pad;

    double dlnbin = (std::log(this->bin_max) - std::log(bin_min))
      / double(this->num_bins - this->nbin_pad);

    for (int ibin = this->nbin_pad; ibin < this->num_bins; ibin++) {
      double edge_left =
        bin_min * std::exp(dlnbin * (ibin - this->nbin_pad));
      double edge_right =
        bin_min * std::exp(dlnbin * (ibin - this->nbin_pad + 1));
      double centre = (edge_left + edge_right) / 2.;

      this->bin_edges.push_back(edge_left);
      this->bin_centres.push_back(centre);
      this->bin_widths.push_back(edge_right - edge_left);
    }
    this->bin_edges.push_back(this->bin_max);
  } else {
    throw trvs::InvalidParameterError(
      "Unrecognised/unsupported binning `scheme`: %s.", scheme.c_str()
    );
  }
}


// ***********************************************************************
// Mesh grids
// ***********************************************************************


// ***********************************************************************
// Line of sight
// ***********************************************************************


// ***********************************************************************
// Clustering statistics
// ***********************************************************************

}  // namespace trv
