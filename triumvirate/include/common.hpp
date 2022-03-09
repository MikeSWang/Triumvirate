/**
 * @file common.hpp
 * @brief Program global variables.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_COMMON_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_COMMON_HPP_INCLUDED_

/// Initialise program trackers.

int currTask = 0;  ///< current task
int numTasks = 1;  ///< number of tasks (in a batch)

double bytesMem = 0;  ///< memory usage in bytes

double clockStart;  ///< program start clock
double clockElapsed;  ///< program elapsed clock

#endif  // TRIUMVIRATE_INCLUDE_COMMON_HPP_INCLUDED_
