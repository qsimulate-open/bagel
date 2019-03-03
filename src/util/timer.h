//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: timer.h
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

#ifndef __BAGEL_UTIL_TIMER_H
#define __BAGEL_UTIL_TIMER_H

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <src/util/string_util.h>
#include <bagel_config.h>

namespace bagel {

class Timer {
  protected:
    std::chrono::high_resolution_clock::time_point tp_;
    const int level_;

  public:
    Timer(const int level = 0) : tp_(std::chrono::high_resolution_clock::now()), level_(level) {};

    // return duration in seconds
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
        std::cout << "       - " << std::left << std::setw(36) << title << std::right << std::setw(10) << std::fixed << std::setprecision(2) << tick() << std::endl;
      } else if (level_ == -1) {
        title = to_upper(title);
        std::cout << "    * " << std::left << std::setw(39) << title << std::right << std::setw(10) << std::fixed << std::setprecision(2) << tick() << std::endl;
#ifdef HAVE_MPI_H
      } else if (level_ >= 1 && level_ < 3) { // TODO for the time being suppressing the level 3 output
        const std::string indent(13+2*level_, ' ');
        const std::string mark = (level_ == 1 ? "o" : (level_ == 2 ? "*" : "-"));
        std::cout << indent << std::left << mark << " " << std::setw(35) << title << std::right << std::setw(13) << std::fixed << std::setprecision(2) << tick() << std::endl;
#endif
      }
    }
};

}

#endif
