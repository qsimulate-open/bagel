//
// BAGEL - Parallel electron correlation program.
// Filename: timer.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#ifndef __BAGEL_UTIL_TIMER_H
#define __BAGEL_UTIL_TIMER_H

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <src/util/string.h>
#include <bagel_config.h>

namespace bagel {

class Timer {
  protected:
    std::chrono::high_resolution_clock::time_point tp_;
    const int level_;

  public:
    Timer(const int level = 0) : tp_(std::chrono::high_resolution_clock::now()), level_(level) {};

    // return duration in milliseconds
    double tick() {
      auto now = std::chrono::high_resolution_clock::now();
      double out = std::chrono::duration_cast<std::chrono::nanoseconds>(now - tp_).count()*1.0e-9;
      tp_ = now;
      return out;
    }

    // print out timing
    void tick_print(std::string title) {
      if (level_ == 0) {
        // top level printout
        std::cout << "       - " << std::left << std::setw(36) << title << std::right << std::setw(10) << std::setprecision(2) << tick() << std::endl;
      } else if (level_ == -1) {
        title = to_upper(title);
        std::cout << "    * " << std::left << std::setw(39) << title << std::right << std::setw(10) << std::setprecision(2) << tick() << std::endl;
#ifdef HAVE_MPI_H
//    } else if (level_ >= 1 && level_ <= resources__->proc()->print_level()) {
      } else if (level_ >= 1) {
        const std::string indent(13+2*level_, ' ');
        const std::string mark = (level_ == 1 ? "o" : (level_ == 2 ? "*" : "-"));
        std::cout << indent << std::left << mark << " " << std::setw(35) << title << std::right << std::setw(13) << std::setprecision(2) << tick() << std::endl;
#endif
      }
    }
};

}

#endif
