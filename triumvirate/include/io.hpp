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
 * @file io.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang)
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief I/O support including custom exceptions and utility functions.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_IO_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_IO_HPP_INCLUDED_

#include <cstdarg>
#include <stdexcept>
#include <string>

namespace trv {
namespace sys {

/**
 * @brief Exception raised when an input/output operation fails.
 *
 */
class IOError: public std::runtime_error {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @ref trv::sys::IOError exception.
   *
   * @param fmt_string Error message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  IOError(const char* fmt_string, ...);

  /**
   * @brief Exception string representation.
   *
   * @returns String representation of the exception.
   */
  virtual const char* what() const noexcept;
};

/**
 * @brief Check if a file path is set.
 *
 * @param pathstr File path string.
 * @returns { @c true , @c false }
 */
bool if_filepath_is_set(std::string pathstr);

}  // namespace trv::sys
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_IO_HPP_INCLUDED_
