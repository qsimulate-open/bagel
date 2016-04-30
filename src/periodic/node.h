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
#include <src/df/df.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/rys/naibatch.h>

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
    const double thresh_;
    std::shared_ptr<const DFDist> df_;

    void insert_vertex(std::shared_ptr<const Vertex>);
    void insert_child(std::shared_ptr<const Node> = NULL);
    void init();
    void get_interaction_list();
    void compute_extent(const double thresh = PRIM_SCREEN_THRESH);
    void insert_neighbour(std::shared_ptr<const Node> neigh, const bool is_neighbour = false, const int ws = 2);
    void make_interaction_list(const int ws);
    void form_df(const std::string auxfile);

    int nbasis_;
    bool is_same_as_parent_;
    int rank_;
    int iself_; // in neighbour_
    std::vector<std::shared_ptr<const ZMatrix>> multipoles_;
    std::vector<std::shared_ptr<const ZMatrix>> local_moment_;
    std::shared_ptr<const ZMatrix> local_expansion_;
    std::shared_ptr<const ZMatrix> exchange_;
    std::vector<std::shared_ptr<const ZMatrix>> child_local_expansion_;
    std::array<double, 3> compute_centre(std::array<std::shared_ptr<const Shell>, 2> shells);
    void compute_multipoles(const int lmax = ANG_HRR_END);
    std::shared_ptr<const ZMatrix> compute_NAI_far_field(const int lmax);
    void compute_local_expansions(std::shared_ptr<const Matrix> density, const int lmax, const std::vector<int> offsets);
    std::shared_ptr<const ZMatrix> compute_Coulomb(const int nbasis, std::shared_ptr<const Matrix> density, std::vector<int> offsets, const bool dodf = false, const std::vector<double> schwarz = std::vector<double>(), const double schwarz_thresh = 0.0);
    std::shared_ptr<const ZMatrix> compute_exact_Coulomb_FF(std::shared_ptr<const Matrix> density, std::vector<int> offsets);
    std::shared_ptr<const DFDist_ints<ERIBatch>> form_fit(const int nbas, const int naux, std::vector<std::shared_ptr<const Atom>> atoms, std::vector<std::shared_ptr<const Atom>> aux_atoms) const {
      return std::make_shared<const DFDist_ints<ERIBatch>>(nbas, naux, atoms, aux_atoms, thresh_, true /*J^-1/2*/, 0.0/*dum*/, false /*average*/, nullptr /*data2*/, true /*serial*/);
    }
    void sort_neighbours(std::vector<std::shared_ptr<const Node>> neighbours);

  public:
    Node(const std::bitset<nbit__> key = 0, const int depth = 0,
         std::shared_ptr<const Node> parent = NULL, const double thresh = PRIM_SCREEN_THRESH);

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

    bool is_same_as_parent() const { return is_same_as_parent_; }
    int rank() const { return rank_; }
    std::vector<std::shared_ptr<const ZMatrix>> multipoles() const { return multipoles_; }
    std::vector<std::shared_ptr<const ZMatrix>> local_moment() const { return local_moment_; }
    std::shared_ptr<const ZMatrix> local_moment(const int i) const { return local_moment_[i]; }
    std::shared_ptr<const ZMatrix> multipoles(const int i) const { return multipoles_[i]; }
    std::shared_ptr<const ZMatrix> local_expansion() const { return local_expansion_; }
    std::vector<std::shared_ptr<const ZMatrix>> child_local_expansion() const { return child_local_expansion_; }
    std::shared_ptr<const ZMatrix> child_local_expansion(const int i) const { return child_local_expansion_[i]; }
    void print_node() const;
};

}
#endif
