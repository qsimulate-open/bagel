//
// BAGEL - Parallel electron correlation program.
// Filename: vertex.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_PERIODIC_VERTEX_H
#define __SRC_PERIODIC_VERTEX_H

#include <array>
#include <bitset>
#include <iostream>
#include <src/util/constants.h>

namespace bagel {

class Vertex {
  protected:
    std::bitset<64> key_;
    std::array<double, 3> position_;

    void compute_multipoles();

  public:
    Vertex(std::bitset<64> key, std::array<double, 3> coord) : key_(key), position_(coord) { }
    ~Vertex() { }

    std::bitset<64> key() const { return key_; }
    std::bitset<3> node_key(const int i) const {
      std::bitset<3> out;
      out[0] = key_[i * 3    ];
      out[1] = key_[i * 3 + 1];
      out[2] = key_[i * 3 + 2];
      return out;
    }

    std::array<double, 3> position() const { return position_; }
    double position(const int i) const { return position_[i]; }
};

}
#endif
