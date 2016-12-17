//
// BAGEL - Parallel electron correlation program.
// Filename: tree_sp.h
// Copyright (C) 2016 Toru Shiozaki
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


#ifndef __SRC_PERIODIC_TREE_SP_H
#define __SRC_PERIODIC_TREE_SP_H

#include <set>
#include <src/wfn/geometry.h>
#include <src/periodic/node_sp.h>
#include <src/periodic/vertex_sp.h>

namespace bagel {

class NodeSP;
class TreeSP {
  protected:
    std::shared_ptr<const Geometry>geom_;
    int max_height_;
    int lmax_;
    int nvertex_;
    int nbasis_;
    std::vector<std::array<double, 3>> coordinates_;
    std::array<double, 3> centre_;

    std::vector<std::bitset<nbit__>> particle_keys_;
    std::vector<std::shared_ptr<const VertexSP>> vertex_;
    std::vector<int> ordering_;
    int nnode_, nleaf_;
    std::vector<std::shared_ptr<NodeSP>> nodes_;
    int height_;
    std::vector<int> shell_in_geom_;
    // to define well-separated distributions
    double thresh_;
    int ws_;

    void init();
    void build_tree();
    void get_particle_key(); // a place holder and (nbit__-1)/3 per coordinate
    void keysort();

  public:
    TreeSP(std::shared_ptr<const Geometry> geom, const int max_height = (nbit__ - 1)/3, const int lmax = 10, const double thresh = PRIM_SCREEN_THRESH, const int ws = 2);

    ~TreeSP() { }

    void init_fmm(const bool dodf, const std::string auxfile) const;
    std::shared_ptr<const ZMatrix> fmm(std::shared_ptr<const Matrix> density = nullptr, const bool dodf = false, const double schwarz_thresh = 0.0) const;

    void print_tree_xyz() const;
    void print_leaves() const;
};

}
#endif
