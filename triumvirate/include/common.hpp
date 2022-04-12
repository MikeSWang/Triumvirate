/**
 * @file common.hpp
 * @brief Program global variables.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_COMMON_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_COMMON_HPP_INCLUDED_

#include <cstring>
#include <ctime>
#include <string>

/// Initialise program trackers.

int currTask = 0;  ///< current task
int numTasks = 1;  ///< number of tasks (in a batch)

const double BYTES_PER_GBYTES = 1024. * 1024. * 1024.;  ///< bytes per gibibyte
double gbytesMem = 0;  ///< memory usage in gibibytes

double clockStart;  ///< program start clock
double clockElapsed;  ///< program elapsed clock

/**
 * Return string showing elapsed time.
 *
 * @param clocktime Clock time.
 * @return elapsed_time Elapsed-time string.
 */
std::string calc_elapsed_time_in_hhmmss(double clocktime) {
  int time = int(clocktime / CLOCKS_PER_SEC);

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

#endif  // TRIUMVIRATE_INCLUDE_COMMON_HPP_INCLUDED_
