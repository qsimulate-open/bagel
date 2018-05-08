//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: muffle.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __SRC_UTIL_MUFFLE_H
#define __SRC_UTIL_MUFFLE_H

#include <iostream>
#include <fstream>
#include <src/util/parallel/mpi_interface.h>

namespace bagel {

// Hides cout and restores it when the object is destroyed
// If split_nodes is set to true, then each node prints to a different file
class Muffle {
  private:
    std::shared_ptr<std::ostream> redirect_;
    std::streambuf* saved_;

  public:
    Muffle(std::string filename = "", const bool append = false, const bool split_nodes = false) {
      saved_ = std::cout.rdbuf();
      if (split_nodes)
        filename += ("_" + std::to_string(mpi__->world_rank()));

      if (filename != "" && (mpi__->world_rank() == 0 || split_nodes))
        redirect_ = append ? std::make_shared<std::ofstream>(filename, std::ios::app) : std::make_shared<std::ofstream>(filename);
      else
        redirect_ = std::make_shared<std::ostringstream>();

      std::cout.rdbuf(redirect_->rdbuf());
    }

    ~Muffle() {
      std::cout.rdbuf(saved_);
    }

    void unmute() { std::cout.rdbuf(saved_); }
    void mute()   { std::cout.rdbuf(redirect_->rdbuf()); }
};

}

#endif
