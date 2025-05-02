// Triumvirate: Three-Point Clustering Measurements in LSS
//
// Copyright (C) 2025 Mike S Wang & Naonori S Sugiyama [GPL-3.0-or-later]
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
 * @file mesh.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang)
 * @brief Mesh grid (with one-point statistics) and
 *        pseudo two-point statistics.
 *
 * This module performs operations, transforms and calculations on a mesh
 * grid of discretely sampled fields, including the assignment of a
 * catalogue of particles with weights or the loading of external field
 * data onto the said mesh grid.  It provides the methods needed to
 * compute various constituent terms (one-point and pseudo two-point
 * statistics) in the estimators of two- and three-point statistics.
 *
 */

#include "mesh.hpp"

namespace trva = trv::array;
namespace trvs = trv::sys;
namespace trvm = trv::maths;

namespace trv {

}  // namespace trv
