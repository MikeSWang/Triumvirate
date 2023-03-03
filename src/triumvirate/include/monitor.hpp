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
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Provide tracking of program resources and exceptions.
 *
 * This module defines the global variables and functions to track
 * time and memory usage of the program.  It also defines custom
 * exceptions used in the program.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>
#include <string>

/// @cond DOXYGEN_DOC_MACROS
// Declares OMP macros.
#ifdef TRV_USE_OMP
#include <omp.h>
#define OMP_ATOMIC _Pragma("omp atomic")
#define OMP_CRITICAL _Pragma("omp critical")
#else  // !TRV_USE_OMP
#define OMP_ATOMIC
#define OMP_CRITICAL
#endif  // TRV_USE_OMP
/// @endcond

// Enter debugging mode.
#ifdef DBG_MODE
#include <iostream>
#endif  // DBG_MODE

namespace trv {
namespace sys {

// ***********************************************************************
// Program tracking
// ***********************************************************************

// RFE: Implement MPI.
extern int currTask;  ///< current task

extern double gbytesMem;     ///< current memory usage in gibibytes
extern double gbytesMaxMem;  ///< maximum memory usage in gibibytes

/**
 * @brief Return size in gibibytes.
 *
 * @tparam T A @c typename.
 * @param num Number of type-double elements.
 * @returns Size in gibibytes.
 */
template <typename T>
double size_in_gb(int num) {
  const double BYTES_PER_GBYTES = 1073741824.;  // 1024³ bytes per gibibyte
  return double(num) * sizeof(T) / BYTES_PER_GBYTES;
}

/**
 * @brief Update the maximum memory usage estimate.
 *
 */
void update_maxmem();

/**
 * @brief Return the current date-time string in
 *        'YYYY-MM-DD HH:MM:SS' format.
 *
 * @returns Timestamp string.
 */
std::string show_current_datetime();

/**
 * @brief Return the elapsed-time string in 'HH:MM:SS' format.
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

/**
 * @brief Logging levels.
 *
 */
enum LogLevel {
  // CAVEAT: Discretionary choices.
  NSET = 0,   ///<  0: unset
  DBUG = 10,  ///< 10: debugging
  STAT = 20,  ///< 20: status
  INFO = 30,  ///< 30: info
  WARN = 40,  ///< 40: warning
  ERRO = 50   ///< 50: error/critical
};

/**
 * @brief Logger with logging level differentiation.
 *
 */
class Logger {
 public:
  int level_limit;  ///< logger threshold level

  /**
   * @brief Construct the logger with the specified threshold level.
   *
   * @param level Threshold level (from enumerated options) (default
   *              is `NSET`).
   */
  Logger(LogLevel level);

  /**
   * @brief Construct the logger with the specified threshold level.
   *
   * @param level Threshold level (as a non-negative integer) (default
   *              is 0).
   *
   * @overload
   */
  Logger(int level);

  /**
   * @brief Reset the logger threshold level.
   *
   * @param level Threshold level (from enumerated options).
   */
  void reset_level(LogLevel level);

  /**
   * @brief Reset the logger threshold level.
   *
   * @param level Threshold level (as a non-negative integer).
   *
   * @overload
   */
  void reset_level(int level);

  /**
   * @brief Log a message at the specified level.
   *
   * If the specified level is below the threshold, no messages will be
   * logged.
   *
   * @param level_entry Message level (from enumerated options).
   * @param fmt_string Log message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  void log(LogLevel level_entry, const char* fmt_string, ...);

  /**
   * @brief Log a message at the specified level.
   *
   * If the specified level is below the threshold, no messages will be
   * logged.
   *
   * @param level_entry Message level (as a non-negative integer).
   * @param fmt_string Log message format string.
   * @param ... An arbitrary number of substitution arguments.
   *
   * @overload
   */
  void log(int level_entry, const char* fmt_string, ...);

  /**
   * @brief Emit a debugging-level message.
   *
   * If the threshold level is higher, no messages will be emitted.
   *
   * @param fmt_string Log message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  void debug(const char* fmt_string, ...);

  /**
   * @brief Emit a status-level message.
   *
   * If the threshold level is higher, no messages will be emitted.
   *
   * @param fmt_string Log message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  void stat(const char* fmt_string, ...);

  /**
   * @brief Emit a information-level message.
   *
   * If the threshold level is higher, no messages will be emitted.
   *
   * @param fmt_string Log message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  void info(const char* fmt_string, ...);

  /**
   * @brief Emit a warning-level message.
   *
   * If the threshold level is higher, no messages will be emitted.
   *
   * @param fmt_string Log message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  void warn(const char* fmt_string, ...);

  /**
   * @brief Emit a warning-level message.
   *
   * If the threshold level is higher, no messages will be emitted.
   *
   * @param fmt_string Log message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  void error(const char* fmt_string, ...);

 private:
  void emit(std::string log_type, const char* fmt_string, std::va_list args);
};

extern Logger logger;  ///< default logger (at `NSET` logging level)


// ***********************************************************************
// Program exceptions
// ***********************************************************************

/**
 * @brief Exception raised when a function or method is unimplemented.
 *
 */
class UnimplementedError: public std::logic_error {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @ref trv::sys::UnimplementedError exception.
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
 * @brief Exception raised when parameters are invalid.
 *
 */
class InvalidParameterError: public std::invalid_argument {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @ref trv::sys::InvalidParameterError exception.
   *
   * @param fmt_string Error message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  InvalidParameterError(const char* fmt_string, ...);

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
class InvalidDataError: public std::runtime_error {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @ref trv::sys::InvalidDataError exception.
   *
   * @param fmt_string Error message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  InvalidDataError(const char* fmt_string, ...);

  /**
   * @brief Exception string representation.
   *
   * @returns String representation of the exception.
   */
  virtual const char* what() const noexcept;
};


// ***********************************************************************
// Program notices
// ***********************************************************************

/**
 * @brief Display program notice including logo in @c stdout.
 *
 */
void display_prog_notice();

}  // namespace trv::sys
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_
