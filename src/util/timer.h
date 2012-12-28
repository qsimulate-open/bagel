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

#ifndef __NEWINT_UTIL_TIMER_H
#define __NEWINT_UTIL_TIMER_H

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>

namespace bagel {

class Timer {
  protected:
    std::chrono::high_resolution_clock::time_point tp_; 

  public:
    Timer() : tp_(std::chrono::high_resolution_clock::now()) {};

    // return duration in milliseconds
    double tick() {
      auto now = std::chrono::high_resolution_clock::now();
      double out = std::chrono::duration_cast<std::chrono::milliseconds>(now - tp_).count()*0.001;
      tp_ = now;
      return out;
    }

    // print out timing
    void tick_print(const std::string& title, const int level = 0) {
      const std::string indent(15+2*level, ' ');
      const std::string mark = (level == 0 ? "o" : (level == 1 ? "*" : "-"));
      std::cout << indent << std::left << mark << " " << std::setw(35) << title << std::right << std::setprecision(2) << tick() << std::endl;
    }
};

}

#endif
