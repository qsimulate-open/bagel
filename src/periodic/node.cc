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

  nbody_  = 0;
}


void Node::insert_vertex(shared_ptr<const Vertex> p) {

  bodies_.resize(nbody_ + 1);
  bodies_[nbody_] = p;
  ++nbody_;

  //cout << bodies_[nbody_-1]->position(0) << "  " << bodies_[nbody_-1]->position(1) << "  " << bodies_[nbody_-1]->position(2) << endl;
}
