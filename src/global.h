//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: global.h
// Copyright (C) 2009 Quantum Simulation Technologies, Inc.
//
// Author: Toru Shiozaki <shiozaki@qsimulate.com>
// Maintainer: QSimulate
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __SRC_GLOBAL__H
#define __SRC_GLOBAL__H

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <cstdlib>

#include <bagel_config.h>

namespace bagel {

extern void static_variables();

static void print_header() {
  std::cout << std::endl;
  std::cout << "  ===============================================================" << std::endl;
  std::cout << "    BAGEL - Freshly leavened quantum chemistry                   " << std::endl;
  std::cout << "  ===============================================================" << std::endl;
  std::cout << std::endl;
}

static void print_version() {
  std::cout << "    Version: " << VERSION << std::endl;
  std::cout << std::endl;
}

static void print_footer() {
  std::cout << std::endl;
  std::cout << "  " << std::endl;
  std::cout << "  ===============================================================" << std::endl;
  std::cout << std::endl;
}


template<typename T>
std::string getenv_multiple(const T& head) {
  char const* val = std::getenv(head);
  return val ? std::string(val) : "";
}


template<typename T, typename ...args>
std::string getenv_multiple(const T& head, const args&... tail) {
  char const* val = std::getenv(head);
  return val ? std::string(val) : getenv_multiple(tail...);
}

}

#endif

