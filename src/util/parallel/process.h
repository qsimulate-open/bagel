//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: process.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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

#ifndef __SRC_PARALLEL_PROCESS_H
#define __SRC_PARALLEL_PROCESS_H

#include <iostream>
#include <sstream>

namespace bagel {

class Process {
  protected:
    // original stream
    std::streambuf* cout_orig;
    // dummy stream
    std::stringstream ss_;

    int print_level_;
    bool muted_;

  public:
    Process();
    ~Process();

    void cout_on();
    void cout_off();

    int print_level() const { return print_level_; }
    void set_print_level(const int i) { print_level_ = i; }

};

}

#endif

