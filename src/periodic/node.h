//
// BAGEL - Parallel electron correlation program.
// Filename: node.h
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


#ifndef __SRC_PERIODIC_NODE_H
#define __SRC_PERIODIC_NODE_H

#include <memory>
#include <vector>
#include <src/periodic/vertex.h>

namespace bagel {

class Node {
  friend class Tree;
  protected:
    std::bitset<nbit__> key_;
    int depth_;
    std::shared_ptr<const Node> parent_;
    std::vector<std::weak_ptr<const Node>> children_;

    std::array<double, 3> position_;
    bool is_complete_;
    bool is_leaf_;
    int nbody_, nchild_, nneighbour_, ninter_;
    double extent_;
    std::vector<std::shared_ptr<const Vertex>> bodies_;
    std::vector<std::shared_ptr<const Node>> interaction_list_;
    std::vector<std::shared_ptr<const Node>> neighbour_;

    void insert_vertex(std::shared_ptr<const Vertex>);
    void insert_child(std::shared_ptr<const Node> = NULL);
    void init();
    void get_interaction_list();
    void compute_position();
    void compute_extent(const double thresh = PRIM_SCREEN_THRESH);
    void insert_neighbour(std::shared_ptr<const Node> neigh, const bool is_neighbour = false, const int ws = 2);
    void make_interaction_list(const int ws = 2);

    int nbasis_;
    std::vector<std::shared_ptr<const ZMatrix>> multipoles_;
    std::shared_ptr<const ZMatrix> local_expansion_;
    void compute_multipoles(const int lmax = ANG_HRR_END);
    void compute_local_expansions(std::shared_ptr<const Matrix> density, const int lmax, std::vector<int> offset);
    std::shared_ptr<const ZMatrix> compute_Coulomb(std::shared_ptr<const Matrix> density, const int lmax, std::vector<int> offset);

  public:
    Node(const std::bitset<nbit__> key = 0, const int depth = 0, std::shared_ptr<const Node> parent = NULL);

    ~Node() { }

    std::bitset<nbit__> key() const { return key_; }
    int depth() const { return depth_; }

    std::shared_ptr<const Node> parent() const { return parent_; }
    std::array<double, 3> position() const { return position_; }
    double position(const int i) const { return position_[i]; }

    bool is_complete() const { return is_complete_; }
    bool is_leaf() const { return is_leaf_; }
    int nbody() const { return nbody_; }
    int nchild() const { return nchild_; }
    std::vector<std::shared_ptr<const Vertex>> bodies() const { return bodies_; }
    std::shared_ptr<const Vertex> bodies(const int i) const { return bodies_[i]; }
    std::shared_ptr<const Node> children(const int i) const { return children_[i].lock(); }

    int nbasis() const { return nbasis_; }
    double extent() const { return extent_; }
    int nneighbour() const { return nneighbour_; }
    std::vector<std::shared_ptr<const Node>> neighbour() const { return neighbour_; }
    std::vector<std::shared_ptr<const Node>> interaction_list() const { return interaction_list_; }

    std::vector<std::shared_ptr<const ZMatrix>> multipoles() const { return multipoles_; }
    std::shared_ptr<const ZMatrix> local_expansion() const { return local_expansion_; }
};

}
#endif
