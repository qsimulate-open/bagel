//
// BAGEL - Parallel electron correlation program.
// Filename: muffle.h
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __SRC_UTIL_MUFFLE_H
#define __SRC_UTIL_MUFFLE_H

#include <iostream>

namespace bagel {

// Hides cout and restores it when the object is destroyed
class Muffle {
  private:
    std::stringstream trash_;
    std::streambuf* saved_;

  public:
    Muffle() {
      saved_ = std::cout.rdbuf();
      std::cout.rdbuf(trash_.rdbuf());
    }

    ~Muffle() {
      std::cout.rdbuf(saved_);
    }
};

}

#endif
