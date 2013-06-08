//
// BAGEL - Parallel electron correlation program.
// Filename: global.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_GLOBAL__H
#define __SRC_GLOBAL__H

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <stdlib.h>
#include <boost/regex.hpp>

namespace bagel {

void static_variables(int, char**);

static void print_header() {
  std::cout << std::endl;
  std::cout << "  ===============================================================" << std::endl;
  std::cout << "    BAGEL - Freshly leavened quantum chemistry                   " << std::endl;
  std::cout << "  ===============================================================" << std::endl;
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
  char const* val = getenv(head);
  return val ? std::string(val) : "";
}


template<typename T, typename ...args>
std::string getenv_multiple(const T& head, const args&... tail) {
  char const* val = getenv(head);
  return val ? std::string(val) : getenv_multiple(tail...);
}

}

#endif

