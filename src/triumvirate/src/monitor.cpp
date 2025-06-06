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
 * @file monitor.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori S Sugiyama (https://github.com/naonori)
 *
 */

#include "monitor.hpp"

namespace trv {
namespace sys {

// ***********************************************************************
// Program helpers
// ***********************************************************************

bool has_extension(const std::string& fname, const std::string& fext) {
  if (fname.length() >= fext.length()) {
    return (fname.rfind(fext) == (fname.length() - fext.length()));
  }
  return false;
}

std::string join_strings(
  const std::vector<std::string>& strings,
  const std::string& delimiter
) {
  std::ostringstream oss;
  for (auto it = strings.begin(); it != strings.end(); ++it) {
      oss << *it;
      if (std::next(it) != strings.end()) {
          oss << delimiter;
      }
  }
  return oss.str();
}

std::vector<std::string> split_string(
  const std::string& str,
  const std::string& delimiter
) {
  if (str.empty()) {
    return {};
  }
  if (delimiter.empty()) {
    return {str};
  }

  std::string fullstr = {str};
  std::string substr;
  std::size_t dlpos = 0;

  std::vector<std::string> strings;
  while ((dlpos = fullstr.find(delimiter)) != std::string::npos) {
    substr = fullstr.substr(0, dlpos);
    fullstr.erase(0, dlpos + delimiter.length());
    strings.push_back(substr);
  }
  if (!fullstr.empty()) {
    strings.push_back(fullstr);
  }

  return strings;
}


// ***********************************************************************
// Program tracking
// ***********************************************************************

int currTask = 0;

double gbytesMem = 0.;
double gbytesMaxMem = 0.;

double gbytesMemGPU = 0.;
double gbytesMaxMemGPU = 0.;

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

#if defined(TRV_USE_HIP)
double worksize_in_gb(
  hipfftHandle plan,
  hipfftType ffttype,
  int nx, int ny, int nz,
  std::vector<int> gpus
) {
  std::vector<std::size_t> worksizes(gpus.size());
  HIPFFT_EXEC(hipfftGetSize3d(plan, nx, ny, nz, ffttype, worksizes.data()));

  std::size_t totsize = 0;
  for (std::size_t id = 0; id < gpus.size(); ++id) {
    totsize += worksizes[id];
  }

  const double BYTES_PER_GBYTES = 1073741824.;  // 1024³ bytes per gibibyte
  double worksize_gb = double(totsize) / BYTES_PER_GBYTES;

  return worksize_gb;
}
#elif defined(TRV_USE_CUDA)  // !TRV_USE_HIP && TRV_USE_CUDA
double worksize_in_gb(
  cufftHandle plan,
  cufftType ffttype,
  int nx, int ny, int nz,
  std::size_t* worksizes,
  std::vector<int> gpus
) {
  std::size_t totsize = 0;
  if (gpus.size() == 1) {
    CUFFT_EXEC(cufftGetSize3d(plan, nx, ny, nz, ffttype, &totsize));
  } else {
    CUFFT_EXEC(cufftGetSize3d(plan, nx, ny, nz, ffttype, worksizes));

    for (std::size_t id = 0; id < gpus.size(); ++id) {
      totsize += worksizes[id];
    }
  }

  const double BYTES_PER_GBYTES = 1073741824.;  // 1024³ bytes per gibibyte
  double worksize_gb = double(totsize) / BYTES_PER_GBYTES;

  return worksize_gb;
}
#endif                       // TRV_USE_HIP

void update_maxmem(bool gpu) {
  if (gpu) {
    trv::sys::gbytesMaxMemGPU =
      (trv::sys::gbytesMemGPU > trv::sys::gbytesMaxMemGPU) ?
      trv::sys::gbytesMemGPU : trv::sys::gbytesMaxMemGPU;
  } else {
    trv::sys::gbytesMaxMem =
      (trv::sys::gbytesMem > trv::sys::gbytesMaxMem) ?
      trv::sys::gbytesMem : trv::sys::gbytesMaxMem;
  }
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

std::string show_current_datetime(bool utc) {
  // Get current time.
  auto now = std::chrono::system_clock::now();
  auto timenow = std::chrono::system_clock::to_time_t(now);

  // Format timestamp.
  char buffer[21];
  if (utc) {
    std::strftime(
      buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&timenow)
    );
  } else {
    std::strftime(
      buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&timenow)
    );
  };

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

int get_gpu_count(bool sys) {
  int num_gpus = 0;
#if defined(TRV_USE_HIP)
  HIP_EXEC(hipGetDeviceCount(&num_gpus));
#elif defined(TRV_USE_CUDA)  // !TRV_USE_HIP && TRV_USE_CUDA
  CUDA_EXEC(cudaGetDeviceCount(&num_gpus));
#endif  // TRV_USE_HIP

  if (sys) {
    return num_gpus;
  }

  int max_gpus = num_gpus;

  char* env_gpu_maxnum = std::getenv("TRV_GPU_MAXNUM");
  if (env_gpu_maxnum != nullptr) {
    std::string gpu_maxnum(env_gpu_maxnum);
    std::regex re_int("^[0-9]+$");
    if (std::regex_match(gpu_maxnum, re_int)) {
      max_gpus = std::stoi(gpu_maxnum);
    }
  }

  return num_gpus > max_gpus ? max_gpus : num_gpus;
}

std::vector<int> get_gpu_ids() {
  std::vector<int> gpu_ids;
  for (int i = 0; i < get_gpu_count(); ++i) {
    gpu_ids.push_back(i);
  }
  return gpu_ids;
}

bool is_gpu_available() {
  bool sys = true;
  return (get_gpu_count(sys) > 0);
}

bool is_gpu_enabled() {
  if (!(get_gpu_count() > 0)) {
    return false;
  }

  // Check for environmental variable override.
  char* env_gpu_mode = std::getenv("TRV_GPU_MODE");
  if (env_gpu_mode != nullptr) {
    std::string gpu_mode(env_gpu_mode);
    // If ``TRV_GPU_MODE`` is set to "false", "no", "off" or "0",
    // then GPU mode is disabled.
    if (
      gpu_mode == "false"
      || gpu_mode == "no"
      || gpu_mode == "off"
      || gpu_mode == "0"
    ) {
      return false;
    }
  }

  // By default, GPU mode is enabled.
  return true;
}

bool is_gpu_single() {
  return (get_gpu_count() == 1);
}

void exit_fatal(const std::string& msg) {
  if (is_colourable()) {
    std::cout << "\n\033[1;37;41mFATAL\033[0m: " << msg << std::endl;
  } else {
    std::cout << "\nFATAL: " << msg << std::endl;
  }

  std::exit(EXIT_FAILURE);
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

  if (log_type.empty()) {
    std::printf(
      "[%s%s] %s\n",
      trv::sys::show_timestamp().c_str(), SHOW_CPPSTATE,
      log_mesg_buf
    );
  } else {
    std::printf(
      "[%s %s%s] %s\n",
      trv::sys::show_timestamp().c_str(), log_type.c_str(), SHOW_CPPSTATE,
      log_mesg_buf
    );
  }
}

void Logger::log(LogLevel entry_level, const char* fmt_string, ...) {
  Logger::log(static_cast<int>(entry_level), fmt_string);
}

void Logger::log(int level_entry, const char* fmt_string, ...) {
  if (level_entry >= this->level_limit) {
    // CAVEAT: See @ref trv::sys::LogLevel.
    level_entry /= 10;

    std::string log_type;

    if (is_colourable()) {
      // Use colours in interactive mode.
      switch (level_entry) {
        case LogLevel::NSET:
          log_type = "";
          break;
        case LogLevel::DBUG:
          log_type = "\033[0;36mDBUG\033[0m";  // cyan
          break;
        case LogLevel::STAT:
          log_type = "\033[0;34mSTAT\033[0m";  // blue
          break;
        case LogLevel::INFO:
          log_type = "\033[0;32mINFO\033[0m";  // green
          break;
        case LogLevel::WARN:
          log_type = "\033[0;33mWARN\033[0m";  // yellow
          break;
        case LogLevel::ERRO:
          log_type = "\033[0;31mERRO\033[0m";  // red
          break;
        default:
          throw InvalidParameterError("Unsupported log level.");
      }
    } else {
      // No colours otherwise.
      switch (level_entry) {
        case LogLevel::NSET:
          log_type = "";
          break;
        case LogLevel::DBUG:
          log_type = "DBUG";
          break;
        case LogLevel::STAT:
          log_type = "STAT";
          break;
        case LogLevel::INFO:
          log_type = "INFO";
          break;
        case LogLevel::WARN:
          log_type = "WARN";
          break;
        case LogLevel::ERRO:
          log_type = "ERRO";
          break;
        default:
          throw InvalidParameterError("Unsupported log level.");
      }
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
    if (is_colourable()) {
      Logger::emit("\033[0;36mDBUG\033[0m", fmt_string, args);
    } else {
      Logger::emit("DBUG", fmt_string, args);
    }
    va_end(args);
  }
}

void Logger::stat(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::STAT) {
    std::va_list args;
    va_start(args, fmt_string);
    if (is_colourable()) {
      Logger::emit("\033[0;34mSTAT\033[0m", fmt_string, args);
    } else {
      Logger::emit("STAT", fmt_string, args);
    }
    va_end(args);
  }
}

void Logger::info(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::INFO) {
    std::va_list args;
    va_start(args, fmt_string);
    if (is_colourable()) {
      Logger::emit("\033[0;32mINFO\033[0m", fmt_string, args);
    } else {
      Logger::emit("INFO", fmt_string, args);
    }
    va_end(args);
  }
}

void Logger::warn(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::WARN) {
    std::va_list args;
    va_start(args, fmt_string);
    if (is_colourable()) {
      Logger::emit("\033[0;33mWARN\033[0m", fmt_string, args);
    } else {
      Logger::emit("WARN", fmt_string, args);
    }
    va_end(args);
  }
}

void Logger::error(const char* fmt_string, ...) {
  if (this->level_limit <= LogLevel::ERRO) {
    std::va_list args;
    va_start(args, fmt_string);
    if (is_colourable()) {
      Logger::emit("\033[1;31mERRO\033[0m", fmt_string, args);
    } else {
      Logger::emit("ERRO", fmt_string, args);
    }
    va_end(args);
  }
}

bool is_colourable() {
  char* ev_term = std::getenv("TERM");
  char* ev_interactive = std::getenv("TRV_INTERACTIVE");
  if (ev_term == nullptr || ev_interactive == nullptr) {
    return false;
  }
  if (std::strstr(ev_term, "color") == nullptr) {
    return false;
  }
  std::string str_interactive = std::string(ev_interactive);
  return (
    str_interactive == "true" ||
    str_interactive == "yes" ||
    str_interactive == "1" ||
    str_interactive == "on"
  );
}

ProgressBar::ProgressBar(int task_count, std::string name) {
  this->name = name;

  if (task_count < 1) {
    throw InvalidParameterError(
      "Progress bar must count at least one task in total."
    );
  }
  this->task_count = task_count;

  this->set_default_pcpt_nodes();
}

void ProgressBar::set_bar_width(int width) {
  if (width < 1) {
    throw InvalidParameterError("Progress bar width must be at least 1.");
  }
  this->bar_width = width;
}

void ProgressBar::set_nodes(std::vector<float> nodes) {
  if (nodes.size() < 1) {
    throw InvalidParameterError(
      "Progress bar nodes must have at least one element."
    );
  }
  this->nodes = nodes;
}

void ProgressBar::set_task_idx(int task_idx) {
  if (task_idx < 0 or task_idx > this->task_count) {
    throw InvalidParameterError(
      "Progress bar task index must be non-negative and "
      "within the total task count."
    );
  }
  this->task_idx = task_idx;

  float progress = float(this->task_idx) / float(this->task_count);
  this->set_progress(progress);
}

void ProgressBar::set_progress(float progress) {
  if (progress < 0. or progress > 1.) {
    throw InvalidParameterError(
      "Progress bar progress must be within the range [0, 1]."
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

      // Print the timestamped status.
      std::cout << "[";
      std::cout << trv::sys::show_timestamp();
      if (is_colourable()) {
        std::cout << " \033[0;34mSTAT\033[0m";
      } else {
        std::cout << " STAT";
      }
      std::cout << "] ";

      // Print the progress bar.
      std::cout << "[";
      if (is_colourable()) {
        std::cout << "\033[0;34m";
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
      if (is_colourable()) {
        std::cout << "\033[0m";
      }
      std::cout << "] " << int(progress * 100.) << "%";

      // Print the progress bar name.
      if (this->name != "") {
        std::cout << " < " << this->name;
      }

      // Flush the progress-bar line.
      std::cout << "\r";
      std::cout.flush();

      this->next_node_idx += 1;
    }
  } else {
    throw InvalidDataError(
      "Progress bar has already completed: progress %f > 1.", this->progress
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

std::vector<float> set_nodes_by_str(std::string interval_str) {
  try {
    std::stof(interval_str);
  } catch (const std::invalid_argument& e) {
    throw InvalidParameterError(
      "Progress bar interval must be a float number."
    );
  }

  float interval = std::stof(interval_str);
  if (!(0. < interval && interval < 100.)) {
    throw InvalidParameterError(
      "Progress bar interval must be in (0, 100) interval."
    );
  }
  interval /= 100.;

  std::vector<float> nodes = {interval};
  while (nodes.back() < 1.) {
    if (nodes.back() + interval < 1.) {
      nodes.push_back(nodes.back() + interval);
    } else {
      nodes.push_back(1.);
    }
  }

  return nodes;
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

std::string get_build_datetime() {
  // Define helper functions.
  auto convert_month_from_name_to_num = [](const std::string& month) -> int {
    if (month == "Jan") return 1;
    if (month == "Feb") return 2;
    if (month == "Mar") return 3;
    if (month == "Apr") return 4;
    if (month == "May") return 5;
    if (month == "Jun") return 6;
    if (month == "Jul") return 7;
    if (month == "Aug") return 8;
    if (month == "Sep") return 9;
    if (month == "Oct") return 10;
    if (month == "Nov") return 11;
    if (month == "Dec") return 12;
    return 0;
  };

  // Parse macros.
  std::string month_str = std::string(__DATE__).substr(0, 3);
  int month = convert_month_from_name_to_num(month_str);

  int day = std::stoi(std::string(__DATE__).substr(4, 2));
  int year = std::stoi(std::string(__DATE__).substr(7, 4));

  int hour = std::stoi(std::string(__TIME__).substr(0, 2));
  int minute = std::stoi(std::string(__TIME__).substr(3, 2));
  int second = std::stoi(std::string(__TIME__).substr(6, 2));

  int hour_tzoffset = 0;
  int minute_tzoffset = 0;
  int status_tzoffset = 0;
  try {
    if (
      std::strlen(__TZOFFSET__) == 5
      && (__TZOFFSET__[0] == '+' || __TZOFFSET__[0] == '-')
    ) {
      int _tzos_sign = (__TZOFFSET__[0] == '-') ? -1 : 1;
      int _hr_tzoffset = std::stoi(std::string(__TZOFFSET__).substr(1, 2));
      int _min_tzoffset = std::stoi(std::string(__TZOFFSET__).substr(3, 2));
      if (
        0 <= hour_tzoffset && hour_tzoffset < 24
        && 0 <= minute_tzoffset && minute_tzoffset < 60
      ) {
        hour_tzoffset = _tzos_sign * _hr_tzoffset;
        minute_tzoffset = _tzos_sign * _min_tzoffset;
      } else {
        status_tzoffset = 1;
      }
    } else {
      status_tzoffset = 1;
    }
  } catch (const std::invalid_argument& e) {
    status_tzoffset = 1;
  } catch (const std::out_of_range& e) {
    status_tzoffset = 1;
  }

  if (status_tzoffset == 0) {
    hour -= hour_tzoffset;
    minute -= minute_tzoffset;
    if (minute < 0) {
      minute += 60;
      hour -= 1;
    } else
    if (minute >= 60) {
      minute -= 60;
      hour += 1;
    }
    if (hour < 0) {
      hour += 24;
      day -= 1;
    } else
    if (hour >= 24) {
      hour -= 24;
      day += 1;
    }
  }

  // Convert to timestamp.
  std::tm bdtm = {};
  bdtm.tm_year = year - 1900;
  bdtm.tm_mon = month - 1;
  bdtm.tm_mday = day;
  bdtm.tm_hour = hour;
  bdtm.tm_min = minute;
  bdtm.tm_sec = second;

  // Format string.
  char bdt[21];
  if (status_tzoffset == 0) {
    std::strftime(bdt, sizeof(bdt), "%Y-%m-%dT%H:%M:%SZ", &bdtm);
  } else {
    std::strftime(bdt, sizeof(bdt), "%Y-%m-%dT%H:%M:%S", &bdtm);
  }

  return std::string(bdt);
}

void display_help() {
  std::printf(
    "Triumvirate: Three-Point Clustering Measurements in LSS\n"
    "\n"
    "\033[1mUsage:\033[0m triumvirate [-h] [-V] <parameter-ini-file>\n"
    "\n"
    "\033[1mPositional arguments:\033[0m\n"
    "  <parameter-ini-file>  path to the parameter INI file\n"
    "\n"
    "\033[1mOptions:\033[0m\n"
    "  -h, --help     show help message and exit\n"
    "  -V, --version  show version and licensing information and exit\n"
  );
}

void display_prog_logo() {
  std::printf("\n");
  std::printf(
    "     //\\        ___  __                      __       ___  ___      \n"
    "    //  \\        |  |__) | |  | |\\/| \\  / | |__)  /\\   |  |__    \n"
    "   //    \\       |  |  \\ | \\__/ |  |  \\/  | |  \\ /~~\\  |  |___ \n"
    "  //      \\                                                         \n"
    " //________\\    "
  );
  if (is_colourable()) {
    std::printf(
      "• \033[1mThree-Point Clustering Measurements in LSS\033[0m •      \n"
    );
  } else {
    std::printf("• Three-Point Clustering Measurements in LSS •      \n");
  }
  std::printf("\n");
}

void display_prog_licence(bool brief) {
  std::printf("Copyright (C) 2023 Mike S Wang & Naonori S Sugiyama\n");
  std::printf("\n");
  if (brief) {
    if (is_colourable()) {
      std::printf("\033[1mLICENCE NOTICE\033[0m >\n\n");
    } else {
      std::printf("LICENCE NOTICE >\n\n");
    }
    std::printf(
      "This program comes with ABSOLUTELY NO WARRANTY. This is     \n"
      "free software, and you are welcome to redistribute it under \n"
      "certain conditions; run `triumvirate --version` for details.\n"
    );
  } else {
    if (is_colourable()) {
      std::printf("\033[1mLICENCE\033[0m >\n\n");
      std::printf("\033[2mGPL-3.0-or-later\033[0m\n\n");
    } else {
      std::printf("LICENCE >\n\n");
      std::printf("GPL-3.0-or-later\n\n");
    }
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
    );
  }
  std::printf("\n");
}

void display_prog_info(bool runtime) {
  if (runtime) {
    if (is_colourable()) {
      std::printf("\033[1mRUNTIME INFORMATION\033[0m >\n\n");
    } else {
      std::printf("RUNTIME INFORMATION >\n\n");
    }
  } else {
    if (is_colourable()) {
      std::printf("\033[1mPROGRAM INFORMATION\033[0m >\n\n");
    } else {
      std::printf("PROGRAM INFORMATION >\n\n");
    }
  }

  std::printf("Triumvirate version: %s\n", __TRV_VERSION__);

  std::printf("GSL version: %s\n", GSL_VERSION);

#if defined(TRV_USE_HIP)
  std::string hipfft_version =
    std::to_string(HIPFFT_MAJOR_VERSION) + "." +
    std::to_string(HIPFFT_MINOR_VERSION) + "." +
    std::to_string(HIPFFT_PATCH_LEVEL);
  std::printf("hipFFT version: %s\n", hipfft_version.c_str());
#elif defined(TRV_USE_CUDA)  // !TRV_USE_HIP && TRV_USE_CUDA
  std::string cufft_version =
    std::to_string(CUFFT_VER_MAJOR) + "." +
    std::to_string(CUFFT_VER_MINOR) + "." +
    std::to_string(CUFFT_VER_PATCH) + "." +
    std::to_string(CUFFT_VER_BUILD);
  std::printf("cuFFT version: %s\n", cufft_version.c_str());
#else   // !TRV_USE_HIP && !TRV_USE_CUDA
  std::printf("FFTW version: %s\n", fftw_version);
#endif  // TRV_USE_HIP

  std::printf("OpenMP version: %s\n", _OMP_VERSION.c_str());

#ifdef TRV_USE_H5
  unsigned hdf5_major, hdf5_minor, hdf5_release;
  H5get_libversion(&hdf5_major, &hdf5_minor, &hdf5_release);
  std::printf(
    "HDF5 version: %d.%d.%d\n", hdf5_major, hdf5_minor, hdf5_release
  );
#endif  // TRV_USE_H5

  std::printf("Build timestamp: %s\n", get_build_datetime().c_str());

  if (runtime) {
    std::printf("CPU thread count: %d\n", _OMP_NTHREADS);

    std::string gpu_mode = " (unavailable)";
    if (is_gpu_available()) {
      gpu_mode = is_gpu_enabled() ? "" : " (disabled)";
    }
    std::printf(
      "GPU device count: %d%s\n",
      trv::sys::get_gpu_count(), gpu_mode.c_str()
    );
  }

  std::printf("\n");
}

void display_prog_logbars(int endpoint) {
  if (endpoint == 0) {
    // std::printf("%s\n", std::string(80, '>').c_str());
    if (is_colourable()) {
      std::printf("\033[1mPROGRAM LOG\033[0m >\n\n");
    } else {
      std::printf("PROGRAM LOG >\n\n");
    }
  } else
  if (endpoint == 1) {
    // std::printf("%s\n", std::string(80, '<').c_str());
  } else {
    throw InvalidParameterError(
      "Invalid endpoint for log bars: %d.", endpoint
    );
  }
}

void expand_envar_in_path(std::string& path_str) {
  std::regex envar_pattern(R"(\$\{([^}]+)\})");
  std::smatch matches;

  auto search_start = path_str.cbegin();
  while (
    std::regex_search(search_start, path_str.cend(), matches, envar_pattern)
  ) {
    std::string var_name = matches[1].str();
    const char* var_value = std::getenv(var_name.c_str());

    if (var_value != nullptr) {
      // Calculate offset before replacing.
      auto start_pos = matches.prefix().second - path_str.cbegin();
      path_str.replace(
        matches.prefix().second, matches.suffix().first, var_value
      );

      // Update `search_start` using offset.
      search_start = path_str.cbegin() + start_pos + std::strlen(var_value);
    } else {
      // If the variable isn't found, skip past this match.
      search_start = matches.suffix().first;
    }
  }
}

}  // namespace trv::sys
}  // namespace trv
