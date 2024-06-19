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
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "monitor.hpp"

/// @cond DOXYGEN_DOC_MACROS
#ifdef TRV_EXTCALL
#define SHOW_CPPSTATE "C++"
#else  // !TRV_EXTCALL
#define SHOW_CPPSTATE "\b"
#endif  // TRV_EXTCALL
/// @endcond

namespace trv {
namespace sys {

// ***********************************************************************
// Program tracking
// ***********************************************************************

int currTask = 0;

double gbytesMem = 0.;
double gbytesMaxMem = 0.;

int count_rgrid = 0;
int count_cgrid = 0;
float count_grid = 0.;
int max_count_rgrid = 0;
int max_count_cgrid = 0;
float max_count_grid = 0.;

int count_fft = 0;
int count_ifft = 0;

bool fftw_wisdom_f_imported = false;
bool fftw_wisdom_b_imported = false;

auto clockStart = std::chrono::steady_clock::now();  ///< program starting time

/// @cond DOXYGEN_DOC_MISC
Logger logger(LogLevel::NSET);  ///< default logger at `NSET` logging level
/// @endcond

void update_maxmem() {
  trv::sys::gbytesMaxMem = (trv::sys::gbytesMem > trv::sys::gbytesMaxMem) ?
    trv::sys::gbytesMem : trv::sys::gbytesMaxMem;
}

void update_maxcntgrid() {
  trv::sys::max_count_rgrid =
    (trv::sys::count_rgrid > trv::sys::max_count_rgrid) ?
    trv::sys::count_rgrid : trv::sys::max_count_rgrid;
  trv::sys::max_count_cgrid =
    (trv::sys::count_cgrid > trv::sys::max_count_cgrid) ?
    trv::sys::count_cgrid : trv::sys::max_count_cgrid;
  trv::sys::max_count_grid =
    (trv::sys::count_grid > trv::sys::max_count_grid) ?
    trv::sys::count_grid : trv::sys::max_count_grid;
}

std::string show_current_datetime() {
  // Get current time.
  auto now = std::chrono::system_clock::now();
  auto timenow = std::chrono::system_clock::to_time_t(now);

  // Format timestamp.
  char buffer[64];
  std::strftime(
    buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&timenow)
  );

  // Print timestamp to string.
  std::string timestamp = std::string(buffer);

  return timestamp;
}

std::string show_elapsed_time(double elapsed_time_in_seconds) {
  // Round time to an integer number of seconds.
  int time = int(elapsed_time_in_seconds);

  // Calculate the hour, minute and second and convert to strings.
  std::string h = std::to_string(time / 3600);
  std::string m = std::to_string((time % 3600) / 60);
  std::string s = std::to_string(time % 60);

  // Format the strings.
  std::string hh;
  if (h.length() < 2) {
    hh = std::string(2 - h.length(), '0') + h;
  } else {
    hh = h;
  }

  std::string mm = std::string(2 - m.length(), '0') + m;
  std::string ss = std::string(2 - s.length(), '0') + s;

  // Form the elapsed time string.
  std::string elapsed_time = hh + ":" + mm + ":" + ss;

  return elapsed_time;
}

std::string show_timestamp() {
  // Calculate the elapsed time in seconds.
  double elapsed_time = double(
    std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - clockStart
    ).count()
  );

  // Format the timestamp.
  char timestamp_[128];
  std::snprintf(
    timestamp_, sizeof(timestamp_), "%s (+%s)",
    show_current_datetime().c_str(), show_elapsed_time(elapsed_time).c_str()
  );

  std::string timestamp(timestamp_);

  return timestamp;
}

Logger::Logger(LogLevel level) {
  Logger::reset_level(level);
}

Logger::Logger(int level) {
  Logger::reset_level(level);
}

void Logger::reset_level(LogLevel level) {
  this->level_limit = level;
}

void Logger::reset_level(int level) {
  this->level_limit = level;
}

void Logger::emit(
  std::string log_type, const char* fmt_string, std::va_list args
) {
  char log_mesg_buf[4096];
  std::vsnprintf(log_mesg_buf, sizeof(log_mesg_buf), fmt_string, args);

  std::printf(
    "[%s %s %s] %s\n",
    trv::sys::show_timestamp().c_str(), log_type.c_str(), SHOW_CPPSTATE,
    log_mesg_buf
  );
}

void Logger::log(LogLevel entry_level, const char* fmt_string, ...) {
  if (entry_level >= this->level_limit) {
    std::string log_type;
    switch (entry_level) {
      case LogLevel::NSET:
        log_type = "\b";
      case LogLevel::DBUG:
        log_type = "DBUG";
      case LogLevel::STAT:
        log_type = "STAT";
      case LogLevel::INFO:
        log_type = "INFO";
      case LogLevel::WARN:
        log_type = "WARN";
      case LogLevel::ERRO:
        log_type = "ERRO";
    }

    std::va_list args;
    va_start(args, fmt_string);
    Logger::emit(log_type, fmt_string, args);
    va_end(args);
  }
}

void Logger::log(int level_entry, const char* fmt_string, ...) {
  if (level_entry >= this->level_limit) {
    // CAVEAT: See @ref trv::sys::LogLevel.
    level_entry /= 10;

    std::string log_type;
    switch (level_entry) {
      case LogLevel::DBUG:
        log_type = "DBUG";
      case LogLevel::STAT:
        log_type = "STAT";
      case LogLevel::INFO:
        log_type = "INFO";
      case LogLevel::WARN:
        log_type = "WARN";
      case LogLevel::ERRO:
        log_type = "ERRO";
    }

    std::va_list args;
    va_start(args, fmt_string);
    Logger::emit(log_type, fmt_string, args);
    va_end(args);
  }
}

void Logger::debug(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::DBUG) {
    std::va_list args;
    va_start(args, fmt_string);
    Logger::emit("DBUG", fmt_string, args);
    va_end(args);
  }
}

void Logger::stat(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::STAT) {
    std::va_list args;
    va_start(args, fmt_string);
    Logger::emit("STAT", fmt_string, args);
    va_end(args);
  }
}

void Logger::info(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::INFO) {
    std::va_list args;
    va_start(args, fmt_string);
    Logger::emit("INFO", fmt_string, args);
    va_end(args);
  }
}

void Logger::warn(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::WARN) {
    std::va_list args;
    va_start(args, fmt_string);
    Logger::emit("WARN", fmt_string, args);
    va_end(args);
  }
}

void Logger::error(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::ERRO) {
    std::va_list args;
    va_start(args, fmt_string);
    Logger::emit("ERRO", fmt_string, args);
    va_end(args);
  }
}

ProgressBar::ProgressBar(int task_count, std::string name) {
  this->name = name;

  if (task_count < 1) {
    throw InvalidParameterError(
      "Progress bar must count at least one task in total.\n"
    );
  }
  this->task_count = task_count;

  this->set_default_pcpt_nodes();
}

void ProgressBar::set_bar_width(int width) {
  if (width < 1) {
    throw InvalidParameterError("Progress bar width must be at least 1.\n");
  }
  this->bar_width = width;
}

void ProgressBar::set_nodes(std::vector<float> nodes) {
  if (nodes.size() < 1) {
    throw InvalidParameterError(
      "Progress bar nodes must have at least one element.\n"
    );
  }
  this->nodes = nodes;
}

void ProgressBar::set_task_idx(int task_idx) {
  if (task_idx < 0 or task_idx > this->task_count) {
    throw InvalidParameterError(
      "Progress bar task index must be non-negative and "
      "within the total task count.\n"
    );
  }
  this->task_idx = task_idx;

  float progress = float(this->task_idx) / float(this->task_count);
  this->set_progress(progress);
}

void ProgressBar::set_progress(float progress) {
  if (progress < 0. or progress > 1.) {
    throw InvalidParameterError(
      "Progress bar progress must be within the range [0, 1].\n"
    );
  }
  this->progress = progress;

  this->next_node_idx = 0;
  while (this->progress > this->nodes[this->next_node_idx]) {
    this->next_node_idx += 1;
  }
}

void ProgressBar::update(int task_idx_now) {
  this->task_idx = task_idx_now;
  float progress_now = float(this->task_idx) / float(this->task_count);
  this->update(progress_now);
}

void ProgressBar::update(float progress_now) {
  this->progress = progress_now;

  if (this->progress <= 1.) {
    if (this->progress >= this->nodes[this->next_node_idx]) {
      int pos = this->bar_width * this->progress;

      if (this->name != "") {
        std::cout << this->name << " [";
      } else {
        std::cout << "[";
      }

      for (int ichar = 0; ichar < this->bar_width; ichar++) {
          if (ichar < pos) {
            std::cout << "=";
          } else
          if (ichar == pos) {
            std::cout << ">";
          } else {
            std::cout << " ";
          }
      }
      std::cout << "] " << int(progress * 100.) << " %\r";
      std::cout.flush();

      this->next_node_idx += 1;
    }
  } else {
    throw InvalidDataError(
      "Progress bar has already completed: progress %f > 1.\n", this->progress
    );
  }
  if (this->progress == 1.) {std::cout << std::endl;}
}

void ProgressBar::set_default_pcpt_nodes() {
  this->nodes.resize(101, 0.);
  for (int pcpt = 0; pcpt <= 100; ++pcpt) {
    this->nodes.push_back(float(pcpt) / 100.);
  }
}


// ***********************************************************************
// Program exceptions
// ***********************************************************************

UnimplementedError::UnimplementedError(const char* fmt_string, ...):
  std::logic_error(
    "Unimplemented error."  // mandatory default error message
  ) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsnprintf(err_mesg_buf, sizeof(err_mesg_buf), fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* UnimplementedError::what() const noexcept {
  return this->err_mesg.c_str();
}

IOError::IOError(const char* fmt_string, ...) : std::runtime_error(
  "I/O error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsnprintf(err_mesg_buf, sizeof(err_mesg_buf), fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* IOError::what() const noexcept {return this->err_mesg.c_str();}

InvalidParameterError::InvalidParameterError(const char* fmt_string, ...):
std::invalid_argument(
  "Invalid parameter error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsnprintf(err_mesg_buf, sizeof(err_mesg_buf), fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* InvalidParameterError::what() const noexcept {
  return this->err_mesg.c_str();
}

InvalidDataError::InvalidDataError(const char* fmt_string, ...):
std::runtime_error(
  "Invalid data error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsnprintf(err_mesg_buf, sizeof(err_mesg_buf), fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* InvalidDataError::what() const noexcept {
  return this->err_mesg.c_str();
}


// ***********************************************************************
// Program notices
// ***********************************************************************

void display_prog_notice() {
  std::printf(
    "     //\\        ___  __                      __       ___  ___      \n"
    "    //  \\        |  |__) | |  | |\\/| \\  / | |__)  /\\   |  |__    \n"
    "   //    \\       |  |  \\ | \\__/ |  |  \\/  | |  \\ /~~\\  |  |___ \n"
    "  //      \\                                                         \n"
    " //________\\    (C) 2023 Mike S Wang & Naonori Sugiyama [GPLv3]     \n"
    "                                                                     \n"
    " •••  Three-Point Clustering Measurements in LSS  •••              \n\n"
  );
}

void display_prog_licence() {
  std::printf(
    "This program is free software: you can redistribute it and/or modify \n"
    "it under the terms of the GNU General Public License as published by \n"
    "the Free Software Foundation, either version 3 of the License, or    \n"
    "(at your option) any later version.                                  \n"
    "                                                                     \n"
    "This program is distributed in the hope that it will be useful, but  \n"
    "WITHOUT ANY WARRANTY; without even the implied warranty of           \n"
    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     \n"
    "General Public License for more details.                             \n"
    "                                                                     \n"
    "You should have received a copy of the GNU General Public License    \n"
    "along with this program. If not, see <https://www.gnu.org/licenses/>.\n"
    "                                                                     \n"
  );
}

}  // namespace trv::sys
}  // namespace trv
