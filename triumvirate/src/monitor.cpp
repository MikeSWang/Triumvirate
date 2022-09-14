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
 * @file monitor.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang)
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "monitor.hpp"

/// FIXME: To be removed.
const double BYTES_PER_GBYTES = 1073741824.;  // 1024^3 bytes per gibibyte

namespace trv {
namespace sys {

/// **********************************************************************
/// Program tracking
/// **********************************************************************

int currTask = 0;

double gbytesMem = 0.;

auto clockStart = std::chrono::steady_clock::now();  ///< program
                                                     ///< starting time

std::string show_current_datetime() {
  /// Get current time.
  auto now = std::chrono::system_clock::now();
  auto timenow = std::chrono::system_clock::to_time_t(now);

  /// Format timestamp.
  char buffer[64];
  std::strftime(
    buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&timenow)
  );

  /// Print timestamp to string.
  std::string timestamp = std::string(buffer);

  return timestamp;
}

std::string show_elapsed_time(double elapsed_time_in_seconds) {
  /// Round time to an integer number of seconds.
  int time = int(elapsed_time_in_seconds);

  /// Calculate the hour, minute and second and convert to strings.
  std::string h = std::to_string(time / 3600);
  std::string m = std::to_string((time % 3600) / 60);
  std::string s = std::to_string(time % 60);

  /// Format the strings.
  std::string hh;
  if (h.length() < 2) {
    hh = std::string(2 - h.length(), '0') + h;
  } else {
    hh = h;
  }

  std::string mm = std::string(2 - m.length(), '0') + m;
  std::string ss = std::string(2 - s.length(), '0') + s;

  /// Form the elapsed time string.
  std::string elapsed_time = hh + ":" + mm + ":" + ss;

  return elapsed_time;
}

std::string show_timestamp() {
  /// Calculate the elapsed time in seconds.
  double elapsed_time = double(
    std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - clockStart
    ).count()
  );

  /// Format the timestamp.
  char timestamp_[128];
  std::sprintf(
    timestamp_, "%s (+%s)",
    show_current_datetime().c_str(), show_elapsed_time(elapsed_time).c_str()
  );

  std::string timestamp(timestamp_);

  return timestamp;
}


/// **********************************************************************
/// Program exceptions
/// **********************************************************************

/// STYLE: Column limit exceeded here.
UnimplementedError::UnimplementedError(const char* fmt_string, ...): std::logic_error(
    "Unimplemented error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsprintf(err_mesg_buf, fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* UnimplementedError::what() const noexcept {
  return this->err_mesg.c_str();
}

/// STYLE: Column limit exceeded here.
InvalidParameter::InvalidParameter(const char* fmt_string, ...): std::invalid_argument(
  "Invalid parameter error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsprintf(err_mesg_buf, fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* InvalidParameter::what() const noexcept {
  return this->err_mesg.c_str();
}

InvalidData::InvalidData(const char* fmt_string, ...): std::runtime_error(
  "Invalid data error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsprintf(err_mesg_buf, fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* InvalidData::what() const noexcept {return this->err_mesg.c_str();}

}  // namespace trv::sys
}  // namespace trv
