/**
 * @file monitor.hpp
 * @brief Program global variables and classes for monitoring.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <string>

/// Include potentially unused libraries.
#include <iostream>

const double BYTES_PER_GBYTES = 1073741824.;  ///< 1024^3 bytes per gibibyte

namespace trv {
namespace runtime {

/// //////////////////////////////////////////////////////////////////////
/// Program tracking
/// //////////////////////////////////////////////////////////////////////

/// RFE: Sort out MPI implementation.

int currTask = 0;  ///< current task
int numTasks = 1;  ///< number of tasks (in a batch)

double gbytesMem = 0.;  ///< memory usage in gibibytes

double clockStart;  ///< program start clock
double clockElapsed;  ///< program elapsed clock

/**
 * Return string showing the current date-time in "YYYY-MM-DD HH:MM:SS"
 * format.
 *
 * @returns timestamp Timestamp string.
 */
std::string show_current_datetime() {
  auto now = std::chrono::system_clock::now();
  auto timenow = std::chrono::system_clock::to_time_t(now);

  char buffer[64];
  std::strftime(
    buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&timenow)
  );

  std::string timestamp = std::string(buffer);

  return timestamp;
}

/// Provide program time tracking.

/**
 * Return string showing elapsed time in "HH:MM:SS" format.
 *
 * @param clock_duration Clock time duration.
 * @returns elapsed_time Elapsed-time string.
 */
std::string show_elapsed_time(double clock_duration) {
  int time = int(clock_duration / CLOCKS_PER_SEC);

  std::string h = std::to_string(time / 3600);
  std::string m = std::to_string((time % 3600) / 60);
  std::string s = std::to_string(time % 60);

  std::string hh;
  if (h.length() < 2) {
    hh = std::string(2 - h.length(), '0') + h;
  } else {
    hh = h;
  }

  std::string mm = std::string(2 - m.length(), '0') + m;
  std::string ss = std::string(2 - s.length(), '0') + s;

  std::string elapsed_time = hh + ":" + mm + ":" + ss;

  return elapsed_time;
}

/**
 * Return logging string showing the timestamp including the elapsed time.
 *
 * @param elapsed_time Elapsed time.
 * @returns timestamp Timestamp string.
 */
std::string show_timestamp() {
  double elapsed_time = double(clock() - clockStart);

  char timestamp_[128];
  std::sprintf(
    timestamp_, "%s (+%s)",
    show_current_datetime().c_str(), show_elapsed_time(elapsed_time).c_str()
  );

  std::string timestamp(timestamp_);

  return timestamp;
}

/// //////////////////////////////////////////////////////////////////////
/// Program I/O
/// //////////////////////////////////////////////////////////////////////

/**
 * Check if a path is set (after trimming).
 *
 * This is used for checking whether a path is set.
 *
 * @param pathstr
 * @returns true, false
 */
bool if_path_is_set(std::string pathstr){
  /// Check if the string is empty.
  if (pathstr.empty()) {return false;}

  /// Check if the string contains non-whitespace characters.  If so,
  /// the path is set, otherwise not.
  for (int ichar = 0; ichar < pathstr.length(); ichar++){
    if (!std::isspace(pathstr[ichar])) {return true;}
  }
  return false;
}

/// //////////////////////////////////////////////////////////////////////
/// Program exceptions
/// //////////////////////////////////////////////////////////////////////

/**
 * Exception raised when an input/output operation fails.
 *
 */
class IOError: public std::runtime_error {
 public:
  std::string err_mesg;

  IOError(const char* fmt_string, ...): std::runtime_error(
    "I/O error."  // default error message; essential
  ) {
    std::va_list args;
    char err_mesg_buf[4096];

    va_start(args, fmt_string);
    std::vsprintf(err_mesg_buf, fmt_string, args);
    va_end(args);

    this->err_mesg = std::string(err_mesg_buf);
  }

  virtual const char* what() const noexcept {
    return err_mesg.c_str();
  }
};

/**
 * Exception raised when parameters are invalid.
 *
 */
class InvalidParameter: public std::invalid_argument {
 public:
  std::string err_mesg;

  InvalidParameter(const char* fmt_string, ...): std::invalid_argument(
    "Invalid parameter error."  // default error message; essential
  ) {
    std::va_list args;
    char err_mesg_buf[4096];

    va_start(args, fmt_string);
    std::vsprintf(err_mesg_buf, fmt_string, args);
    va_end(args);

    this->err_mesg = std::string(err_mesg_buf);
  }

  virtual const char* what() const noexcept {
    return err_mesg.c_str();
  }
};

/**
 * Exception raised when the data to be operated on are invalid.
 *
 */
class InvalidData: public std::runtime_error {
 public:
  std::string err_mesg;

  InvalidData(const char* fmt_string, ...): std::runtime_error(
    "Invalid data error."  // default error message; essential
  ) {
    std::va_list args;
    char err_mesg_buf[4096];

    va_start(args, fmt_string);
    std::vsprintf(err_mesg_buf, fmt_string, args);
    va_end(args);

    this->err_mesg = std::string(err_mesg_buf);
  }

  virtual const char* what() const noexcept {
    return err_mesg.c_str();
  }
};

}  // trv::runtime::
}  // trv::

#endif  // TRIUMVIRATE_INCLUDE_MONITOR_HPP_INCLUDED_
