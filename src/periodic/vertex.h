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

#include <src/util/constants.h>
#include <src/periodic/atomgroup.h>

namespace bagel {

class Vertex {
  protected:
    std::bitset<64> key_;
    std::shared_ptr<const AtomGroup> atomgroup_;
    std::array<double, 3> position_;

  public:
    Vertex(std::bitset<64> key, std::shared_ptr<const AtomGroup> atomgrp)
     : key_(key), atomgroup_(atomgrp), position_(atomgrp->position()) { }
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

    std::shared_ptr<const AtomGroup> atomgroup() const { return atomgroup_; }
    const std::vector<std::shared_ptr<const Atom>> atoms() const { return atomgroup_->atoms(); }
    std::shared_ptr<const Atom> atoms(const int i) const { return atomgroup_->atoms(i); }
    const std::vector<int> order_in_geom()  const { return atomgroup_->order_in_geom(); }
    int order_in_geom(const int i)  const { return atomgroup_->order_in_geom(i); }
    const std::vector<int> ishell()  const { return atomgroup_->ishell(); }
    int ishell(const int i)  const { return atomgroup_->ishell(i); }
    int nbasis() const { return atomgroup_->nbasis(); }
    int natom() const { return atomgroup_->natom(); }
    bool is_simple() const { return natom() == 1; }
};

}
#endif
