//
// BAGEL - Parallel electron correlation program.
// Filename: tree.h
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


#ifndef __SRC_PERIODIC_TREE_H
#define __SRC_PERIODIC_TREE_H

#include <vector>
#include <src/molecule/molecule.h>
#include <src/periodic/node.h>
#include <src/periodic/vertex.h>

namespace bagel {

class Node;
class Tree {
  protected:
    int max_height_;
    int nvertex_;
    std::vector<std::array<double, 3>> coordinates_;
    std::array<double, 3> position_;

    double box_length_;
    std::vector<std::bitset<nbit__>> particle_keys_;
    std::vector<std::shared_ptr<const Vertex>> leaves_;
    std::vector<std::shared_ptr<Node>> nodes_;
    int height_;


    void init();
    void build_tree();
    void get_particle_key(); // a place holder and (nbit__-1)/3 per coordinate
    void keysort();

  public:
    Tree(std::shared_ptr<const Molecule> mol, const int max_height = (nbit__ - 1)/3);
    ~Tree() { }

};

}
#endif
