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
 * @file io.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang)
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "io.hpp"

namespace trv {
namespace sys {

IOError::IOError(const char* fmt_string, ...): std::runtime_error(
    "I/O error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsprintf(err_mesg_buf, fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* IOError::what() const noexcept {return this->err_mesg.c_str();}

bool if_filepath_is_set(std::string pathstr) {
  /// Check if the string is empty.
  if (pathstr.empty()) {return false;}

  /// Check if the path is a directory not file.
  std::string endchar = "/";
  int strcomp = pathstr.compare(
    pathstr.length() - endchar.length(), endchar.length(), endchar
  );
  if (strcomp == 0) {return false;}  // `pathstr` ends in "/"

  /// Check if the string contains non-whitespace characters.  If so,
  /// the path is set, otherwise not.
  for (int ichar = 0; ichar < pathstr.length(); ichar++){
    if (!std::isspace(pathstr[ichar])) {return true;}
  }

  return false;
}

}  // namespace trv::sys
}  // namespace trv
