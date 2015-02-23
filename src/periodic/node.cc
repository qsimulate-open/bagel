//
// BAGEL - Parallel electron correlation program.
// Filename: node.cc
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


#include <src/periodic/node.h>

using namespace bagel;
using namespace std;

Node::Node(const std::bitset<nbit__> key, const int depth, std::shared_ptr<const Node> parent)
 : key_(key), depth_(depth), parent_(parent) {
  if (depth == 0)
    key_[0] = 1;

  is_complete_ = false;
  nbody_  = 0;
}


void Node::insert_vertex(shared_ptr<const Vertex> p) {

  assert(!is_complete_);
  bodies_.resize(nbody_ + 1);
  bodies_[nbody_] = p;
  ++nbody_;

  //cout << bodies_[nbody_-1]->position(0) << "  " << bodies_[nbody_-1]->position(1) << "  " << bodies_[nbody_-1]->position(2) << endl;
}



void Node::mark_complete() {

  is_complete_ = true;

  position_ = {{0.0, 0.0, 0.0}};
  double sum = 0.0;
  for (auto& body : bodies_) {
    position_[0] += body->atom()->atom_charge() * body->position(0);
    position_[1] += body->atom()->atom_charge() * body->position(1);
    position_[2] += body->atom()->atom_charge() * body->position(2);
    sum += body->atom()->atom_charge();
  }
  position_[0] /= sum;
  position_[1] /= sum;
  position_[2] /= sum;
}
