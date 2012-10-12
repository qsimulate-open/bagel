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
#include <boost/regex.hpp>

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


static int count_string(const std::string inputfile, const std::string keyword) {
  try {
    std::ifstream ifs;
    ifs.open(inputfile.c_str());
    if (!ifs.is_open())
      throw std::runtime_error("input file cannot be opened.");

    boost::smatch what;
    boost::regex reg(keyword);
    int out = 0;
    while(true) {
      std::string sline;
      if (!getline(ifs, sline)) break;
      std::string::const_iterator start = sline.begin();
      std::string::const_iterator end   = sline.end();
      if (boost::regex_search(start, end, what, reg)) ++out;
    }

    ifs.close();
    return out;
  } catch (...) {
    throw;
  }
}

#endif

