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
#include <fstream>

namespace bagel {

// Hides cout and restores it when the object is destroyed
class Muffle {
  private:
    std::shared_ptr<std::ostream> redirect_;
    std::streambuf* saved_;

  public:
    Muffle(std::string filename = "") {
      saved_ = std::cout.rdbuf();
      if (filename != "")
        redirect_ = std::make_shared<std::ofstream>(filename);
      else
        redirect_ = std::make_shared<std::ostringstream>();

      std::cout.rdbuf(redirect_->rdbuf());
    }

    ~Muffle() {
      std::cout.rdbuf(saved_);
    }
};

}

#endif
