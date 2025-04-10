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
 * @file monitor.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori S Sugiyama (https://github.com/naonori)
 * @brief Provide tracking of program resources and exceptions.
 *
 * This module defines the global variables and functions to track
 * time and memory usage of the program.  It also defines custom
 * exceptions used in the program.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_

#include <gsl/gsl_version.h>

#if defined(TRV_USE_HIP)
#include <hipfft/hipfftXt.h>
#elif defined(TRV_USE_CUDA)  // !TRV_USE_HIP && TRV_USE_CUDA
#include <cufftXt.h>
#endif                       // TRV_USE_HIP
#include <fftw3.h>

#ifdef TRV_USE_OMP
#include <omp.h>
#endif  // TRV_USE_OMP

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef TRV_USE_H5
#include "highfive/H5File.hpp"
#endif  // TRV_USE_H5

// Enter debugging mode.
#ifdef DBG_MODE
#include <iostream>
#endif  // DBG_MODE

/// @cond DOXYGEN_DOC_MACROS
// Declares OMP macros.
#ifdef TRV_USE_OMP
#define OMP_ATOMIC _Pragma("omp atomic")
#define OMP_CRITICAL _Pragma("omp critical")
// NOTE: OpenMP is unavailabled in separate GPU compilation pass.
#ifndef __HIP_DEVICE_COMPILE__
#define _OMP_VERSION std::to_string(_OPENMP)
#else   // !__HIP_DEVICE_COMPILE__
#define _OMP_VERSION std::string("unknown")
#endif  // __HIP_DEVICE_COMPILE__
#define _OMP_NTHREADS omp_get_max_threads()
#else   // !TRV_USE_OMP
#define OMP_ATOMIC
#define OMP_CRITICAL
#define _OMP_VERSION std::string("unknown")
#define _OMP_NTHREADS 1
#endif  // TRV_USE_OMP

#ifndef GSL_VERSION
#define GSL_VERSION "unknown"
#endif  // !GSL_VERSION

// Record (fallback) version number.
#ifndef __TRV_VERSION__
#define __TRV_VERSION__ "0.6+"
#endif  // !__TRV_VERSION__

// Indicate external interface call.
#ifdef TRV_EXTCALL
#define SHOW_CPPSTATE " C++"
#else   // !TRV_EXTCALL
#define SHOW_CPPSTATE ""
#endif  // TRV_EXTCALL

// Record (fallback) build timezone offset.
#ifndef __TZOFFSET__
#define __TZOFFSET__ "+0000"
#endif  // !__TZOFFSET__

#if defined(TRV_USE_HIP)
#ifndef HIP_EXEC
/**
 * @brief Check if a HIP function call succeeded.
 *
 * @param hip_call A HIP function call.
 */
#define HIP_EXEC(hip_call) {                                             \
  auto ret_status = static_cast<hipError_t>(hip_call);                   \
  if (ret_status != hipSuccess) {                                        \
    std::fprintf(                                                        \
      stderr,                                                            \
      "HIP error: %s. "                                                  \
      "Function <%s> returned error code %d "                            \
      "in file \"%s\", line %d.\n",                                      \
      hipGetErrorString(ret_status),                                     \
      #hip_call, ret_status,                                             \
      __FILE__, __LINE__                                                 \
    );                                                                   \
  }                                                                      \
}
#endif  // !HIP_EXEC

#ifndef HIPFFT_EXEC
/**
 * @brief Check if a hipFFT function call succeeded.
 *
 * @param hipfft_call A hipFFT function call.
 */
#define HIPFFT_EXEC(hipfft_call) {                                       \
  auto ret_status = static_cast<hipfftResult>(hipfft_call);              \
  if (ret_status != HIPFFT_SUCCESS) {                                    \
    std::fprintf(                                                        \
      stderr,                                                            \
      "hipFFT error: function <%s> returned error code %d "              \
      "in file \"%s\", line %d.\n",                                      \
      #hipfft_call, ret_status,                                          \
      __FILE__, __LINE__                                                 \
    );                                                                   \
  }                                                                      \
}
#endif  // !HIPFFT_EXEC

#elif defined(TRV_USE_CUDA)
#ifndef CUDA_EXEC
/**
 * @brief Check if a CUDA function call succeeded.
 *
 * @param cuda_call A CUDA function call.
 */
#define CUDA_EXEC(cuda_call) {                                           \
  auto ret_status = static_cast<cudaError_t>(cuda_call);                 \
  if (ret_status != cudaSuccess) {                                       \
    std::fprintf(                                                        \
      stderr,                                                            \
      "CUDA error: %s. "                                                 \
      "Function <%s> returned error code %d "                            \
      "in file \"%s\", line %d.\n",                                      \
      cudaGetErrorString(ret_status),                                    \
      #cuda_call, ret_status,                                            \
      __FILE__, __LINE__                                                 \
    );                                                                   \
  }                                                                      \
}
#endif  // !CUDA_EXEC

#ifndef CUFFT_EXEC
/**
 * @brief Check if a cuFFT function call succeeded.
 *
 * @param cufft_call A cuFFT function call.
 */
#define CUFFT_EXEC(cufft_call) {                                         \
  auto ret_status = static_cast<cufftResult>(cufft_call);                \
  if (ret_status != CUFFT_SUCCESS) {                                     \
    std::fprintf(                                                        \
      stderr,                                                            \
      "cuFFT error: function <%s> returned error code %d "               \
      "in file \"%s\", line %d.\n",                                      \
      #cufft_call, ret_status,                                           \
      __FILE__, __LINE__                                                 \
    );                                                                   \
  }                                                                      \
}
#endif // !CUFFT_EXEC

#endif  // TRV_USE_HIP
/// @endcond

namespace trv {
namespace sys {

// ***********************************************************************
// Program helpers
// ***********************************************************************

/**
 * @brief Check if a file has a given extension.
 *
 * @param fname File name.
 * @param fext File extension (with the dot).
 * @returns {`true`, `false`}
 */
bool has_extension(const std::string& fname, const std::string& fext);

/**
 * @brief Join a vector of strings with a delimiter.
 *
 * @param strings A vector of strings.
 * @param delimiter Delimiter string.
 * @return std::string Joined string.
 */
std::string join_strings(
  const std::vector<std::string>& strings,
  const std::string& delimiter
);


// ***********************************************************************
// Program tracking
// ***********************************************************************

// STYLE: Standard naming convention is not followed below.

// RFE: Implement MPI.
extern int currTask;  ///< current task

extern double gbytesMem;     ///< current memory usage in gibibytes
extern double gbytesMaxMem;  ///< maximum memory usage in gibibytes

extern double gbytesMemGPU;     ///< current (GPU) memory usage in gibibytes
extern double gbytesMaxMemGPU;  ///< maximum (GPU) memory usage in gibibytes

extern int count_rgrid;       ///< number of 3-d real grids
extern int count_cgrid;       ///< number of 3-d complex grids
extern float count_grid;      ///< number of grids
extern int max_count_rgrid;   ///< maximum number of 3-d real grids
extern int max_count_cgrid;   ///< maximum number of 3-d complex grids
extern float max_count_grid;  ///< maximum number of grids

extern int count_fft;   ///< number of FFTs
extern int count_ifft;  ///< number of IFFTs

/// wisdom import status for forward transform
extern bool fftw_wisdom_f_imported;
/// wisdom import status for backward transform
extern bool fftw_wisdom_b_imported;

/**
 * @brief Return size in gibibytes.
 *
 * @tparam T A @c typename.
 * @param num Number of type-double elements.
 * @returns Size in gibibytes.
 */
template <typename T>
double size_in_gb(long long num) {
  const double BYTES_PER_GBYTES = 1073741824.;  // 1024³ bytes per gibibyte
  return double(num) * sizeof(T) / BYTES_PER_GBYTES;
}

/**
 * @brief Return size in gibibytes.
 *
 * @tparam T A @c typename.
 * @param num Number of type-double elements.
 * @returns Size in gibibytes.
 *
 * @overload
 */
template <typename T>
double size_in_gb(int num) {
  const double BYTES_PER_GBYTES = 1073741824.;  // 1024³ bytes per gibibyte
  return double(num) * sizeof(T) / BYTES_PER_GBYTES;
}

#if defined(TRV_USE_HIP)
/**
 * @brief Return the work area size in gibibytes for a 3D FFT plan.
 *
 * @param plan 3D FFT plan handle.
 * @param ffttype FFT type.
 * @param nx, ny, nz Transform dimensions.
 * @param gpus GPU device indices.
 * @returns double Work area size in gibibytes.
 */
double worksize_in_gb(
  hipfftHandle plan,
  hipfftType ffttype,
  int nx, int ny, int nz,
  std::vector<int> gpus
);
#elif defined(TRV_USE_CUDA)  // !TRV_USE_HIP && TRV_USE_CUDA
/**
 * @brief Return the work area size in gibibytes for a 3D FFT plan.
 *
 * @param plan 3D FFT plan handle.
 * @param ffttype FFT type.
 * @param nx, ny, nz Transform dimensions.
 * @param worksizes FFT plan work sizes.
 * @param gpus GPU device indices.
 * @returns double Work area size in gibibytes.
 */
double worksize_in_gb(
  cufftHandle plan,
  cufftType ffttype,
  int nx, int ny, int nz,
  std::size_t* worksizes,
  std::vector<int> gpus
);
#endif                       // TRV_USE_HIP

/**
 * @brief Update the maximum memory usage estimate.
 *
 * @param gpu If @c true, update the GPU memory usage estimate.
 */
void update_maxmem(bool gpu = false);

/**
 * @brief Update the maximum 3-d grid counts.
 *
 */
void update_maxcntgrid();

/**
 * @brief Return the current datetime string.
 *
 * @param utc If @c true, return UTC time in ISO 8601--like format, else
 *            return local time in 'YYYY-MM-DD HH:MM:SS' format.
 * @returns Timestamp string.
 */
std::string show_current_datetime(bool utc = false);

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
 * @brief Get the number of GPUs available.
 *
 * @param sys If @c true, return the number of GPU devices in the system,
 *            else return the number of GPU devices available for use.
 * @returns int Number of GPUs.
 *
 * @note This function checks both the number of GPU devices present and
 *       the maximum number available for use as limited by the
 *       environmental variable @c TRV_GPU_MAXNUM.
 */
int get_gpu_count(bool sys = false);

/**
 * @brief Get the indices of GPUs available for use.
 *
 * @returns GPU indices.
 */
std::vector<int> get_gpu_ids();

/**
 * @brief Check if GPUs are available in the system.
 *
 * @returns {`true`, `false`}
 */
bool is_gpu_available();

/**
 * @brief Check if GPU mode is enabled.
 *
 * @returns {`true`, `false`}
 *
 * @note This function checks if at least one GPU is available for use.
 *       It also checks for any override by the environmental
 *       variable @c TRV_GPU_MODE: when it is set to @c false, @c no,
 *       @c off or @c 0, the function will return @c false.
 */
bool is_gpu_enabled();

/**
 * @brief Check if in single-GPU mode.
 *
 * @returns {`true`, `false`}
 */
bool is_gpu_single();

// Determine CUDA stream support for multi-GPU cuFFT (in cuFFT 10.4.0+).
#if defined(TRV_USE_CUDA) && CUFFT_VERSION >= 10400
#define _CUDA_STREAM
#endif  // TRV_USE_CUDA && CUFFT_VERSION >= 10400

/**
 * @brief Terminate the program with exit status `EXIT_FAILURE`.
 *
 */
void exit_fatal(const std::string& msg);

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

/**
 * @brief Check if the program @c stdout is colourable.
 *
 * @returns {`true`, `false`}
 * @note `true` if the environmental variable @c TRV_INTERACTIVE is set
 *       to @c "true" or @c "yes" or @c "1" or @c "on" and @c TERM
 *       contains @c "color", or `false` otherwise.
 */
bool is_colourable();

extern Logger logger;  ///< default logger (at `NSET` logging level)

/**
 * @brief Progress bar for tracking tasks.
 *
 */
class ProgressBar {
 public:
  std::string name;     ///< progress bar name
  int task_count = 1;   ///< total task count
  int task_idx = 0;     ///< task index
  float progress = 0.;  ///< progress value in [0, 1] interval

  /**
   * @brief Construct a progress bar.
   *
   * @param task_count Total task count (at least 1).
   * @param name Progress bar name (default is an empty string).
   */
  explicit ProgressBar(int task_count, std::string name = "");

  /**
   * @brief Set progress bar width.
   *
   * @param bar_width Progress bar width (at least 1).
   */
  void set_bar_width(int bar_width);

  /**
   * @brief Set progress bar width.
   *
   * @param nodes Progress values (strictly increasing in the interval
   *              [0, 1]) at which the bar is updated.
   */
  void set_nodes(std::vector<float> nodes);

  /**
   * @brief Set current (or initial) task index.
   *
   * @param task_idx Current task index within the total task count.
   *
   * @note When the tas index is zero, it means no task has completed yet.
   */
  void set_task_idx(int task_idx);

  /**
   * @brief Set current (or initial) progress value.
   *
   * @param progress Current progress value in the interval [0, 1].
   */
  void set_progress(float progress);

  /**
   * @brief Update the progress bar.
   *
   * @param task_idx_now Latest completed task index.
   *
   * @note When the tas index is zero, it means no task has completed yet.
   */
  void update(int task_idx_now);

  /**
   * @brief Update the progress bar.
   *
   * @param progress_now Latest progress value.
   *
   * @overload
   */
  void update(float progress_now);

 private:
  int bar_width = 50;        ///< progress bar width

  std::vector<float> nodes;  ///< progress nodes
  int next_node_idx = 0;     ///< next progress node index

  /**
   * @brief Set default progress nodes at percentage points.
   *
   */
  void set_default_pcpt_nodes();
};

/**
 * @brief Set a node list possibly from a string.
 *
 * If the string corresponds to a number between 0 and 100, that is then
 * the percentage-point interval at which the nodes are set as by
 * ProgressBar::set_nodes().
 *
 * @param interval_str Interval string.
 * @returns Node list.
 * @throws trv::sys::InvalidParameterError If the interval string
 *                                         is invalid.
 */
std::vector<float> set_nodes_by_str(std::string interval_str);


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
 * @brief Return the build datetime string in ISO 8601--like format.
 *
 * @returns Build datetime string.
 */
std::string get_build_datetime();

/**
 * @brief Display help message in @c stdout.
 *
 */
void display_help();

/**
 * @brief Display program logo in @c stdout.
 *
 */
void display_prog_logo();

/**
 * @brief Display program licence in @c stdout.
 *
 * @param brief Display brief notice only (default is @c false).
 */
void display_prog_licence(bool brief = false);

/**
 * @brief Display program information in @c stdout.
 *
 * @param runtime Display runtime information (default is @c false).
 */
void display_prog_info(bool runtime = false);

/**
 * @brief Display program log bars in @c stdout.
 *
 * @param endpoint Progress bar endpoint, either 0 (start) or 1 (finish).
 */
void display_prog_logbars(int endpoint);


// ***********************************************************************
// Program utilities
// ***********************************************************************

/**
 * @brief Expand environment variables in a path string.
 *
 * @param[in,out] path_str Path string.
 */
void expand_envar_in_path(std::string& path_str);

}  // namespace trv::sys
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_
