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

#include <set>
#include <src/wfn/geometry.h>
#include <src/periodic/node.h>
#include <src/periodic/vertex.h>
#include <src/util/atommap.h>

namespace bagel {

class Node;
class Tree {
  protected:
    std::shared_ptr<const Geometry>geom_;
    int max_height_;
    bool do_contraction_;
    int nvertex_;
    int nbasis_;
    std::vector<std::array<double, 3>> coordinates_;
    std::array<double, 3> position_;

    std::vector<std::bitset<nbit__>> particle_keys_;
    std::vector<std::shared_ptr<const Vertex>> leaves_;
    std::vector<int> ordering_, shell_id_;
    int nnode_;
    std::vector<std::shared_ptr<Node>> nodes_;
    int height_;
    // to define well-separated distributions
    double thresh_;
    int ws_;

    // vertex contraction
    std::vector<std::shared_ptr<const AtomGroup>> atomgroup_;
    void contract_vertex();


    void init();
    void build_tree();
    void get_particle_key(); // a place holder and (nbit__-1)/3 per coordinate
    void keysort();

    std::shared_ptr<const ZMatrix> coulomb_;

  public:
    Tree(std::shared_ptr<const Geometry> geom, const int max_height = (nbit__ - 1)/3, const bool do_contract = false,
         const double thresh = PRIM_SCREEN_THRESH, const int ws = 1);
    ~Tree() { }

    void fmm(const int lmax, std::shared_ptr<const Matrix> density);
    std::shared_ptr<const ZMatrix> coulomb() { return coulomb_; }

    void print_tree_xyz() const;
};

}
#endif
