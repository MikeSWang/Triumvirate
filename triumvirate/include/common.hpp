#ifndef TRIUMVIRATE_INCLUDE_COMMON_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_COMMON_HPP_INCLUDED_

#include <ctime>

/// Initialise program trackers.

/// FIXME: Remove out of header files.
int currTask = 0;  ///< current task
int numTasks = 1;  ///< number of tasks (in a batch)

double bytesMem = 0;  ///< memory usage in bytes

double timeStart;  ///< program start time
double durationInSec;  ///< program duration in seconds

#endif
