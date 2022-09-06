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
 * @file monitor.hpp
 * @author Mike S Wang (https://github.com/MikeSWang)
 * @brief Provide tracking of program time and memory usage
 *        and exceptions.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

/// FIXME: To be removed.
extern const double BYTES_PER_GBYTES;  ///< bytes per gibibyte

// const double BYTES_PER_GBYTES = 1073741824.;  ///< 1024^3 bytes
//                                               ///< per gibibyte

namespace trv {
namespace mon {

/// **********************************************************************
/// Program tracking
/// **********************************************************************

/// RFE: Sort out MPI implementation.

extern int currTask;  ///< current task
extern int numTasks;  ///< number of tasks (in a batch)

extern double gbytesMem;  ///< memory usage in gibibytes

// /**
//  * @brief Return size in gibibytes.
//  *
//  * @tparam T A @c typename.
//  * @returns Size in gibibytes.
//  */
// template <typename T>
// double size_in_gb() {return sizeof(T) / BYTES_PER_GBYTES;}

/**
 * @brief Return the current date-time string in
 *        "YYYY-MM-DD HH:MM:SS" format.
 *
 * @returns Timestamp string.
 */
std::string show_current_datetime();

/**
 * @brief Return the elapsed-time string in "HH:MM:SS" format.
 *
 * @param duration_in_seconds Duration in seconds.
 * @returns Elapsed-time string.
 */
std::string show_elapsed_time(double duration_in_seconds);

/**
 * @brief Return the timestamp string including the elapsed time.
 *
 * @returns Timestamp string.
 */
std::string show_timestamp();

/// **********************************************************************
/// Program I/O
/// **********************************************************************

/**
 * @brief Check if a file path is set.
 *
 * @param pathstr File path string.
 * @returns {true, false}
 */
bool if_filepath_is_set(std::string pathstr);

/// **********************************************************************
/// Program exceptions
/// **********************************************************************

/**
 * @brief Exception raised when an input/output operation fails.
 *
 */
class IOError: public std::runtime_error {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @c IOError exception.
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
 * @brief Exception raised when a function or method is unimplemented.
 *
 */
class UnimplementedError: public std::logic_error {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @c UnimplementedError exception.
   *
   * @param fmt_string Error message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  UnimplementedError(const char* fmt_string, ...);

  /**
   * @brief Exception string representation.
   *
   * @returns String representation of the exception.
   */
  virtual const char* what() const noexcept;
};

/**
 * @brief Exception raised when parameters are invalid.
 *
 */
class InvalidParameter: public std::invalid_argument {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @c InvalidParameter exception.
   *
   * @param fmt_string Error message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  InvalidParameter(const char* fmt_string, ...);

  /**
   * @brief Exception string representation.
   *
   * @returns String representation of the exception.
   */
  virtual const char* what() const noexcept;
};

/**
 * @brief Exception raised when the data to be operated on are invalid.
 *
 */
class InvalidData: public std::runtime_error {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @c InvalidData exception.
   *
   * @param fmt_string Error message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  InvalidData(const char* fmt_string, ...);

  /**
   * @brief Exception string representation.
   *
   * @returns String representation of the exception.
   */
  virtual const char* what() const noexcept;
};

}  // namespace trv::mon
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_
