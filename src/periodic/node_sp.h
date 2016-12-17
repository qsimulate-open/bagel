//
// BAGEL - Parallel electron correlation program.
// Filename: node_sp.h
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


#ifndef __SRC_PERIODIC_NODE_SP_H
#define __SRC_PERIODIC_NODE_SP_H

#include <memory>
#include <vector>
#include <src/periodic/vertex_sp.h>
#include <src/df/df.h>
#include <src/integral/rys/eribatch.h>

namespace bagel {

class NodeSP {
  friend class TreeSP;
  protected:
    std::bitset<nbit__> key_;
    int depth_;
    std::shared_ptr<const NodeSP> parent_;
    std::vector<std::weak_ptr<const NodeSP>> child_;

    std::array<double, 3> centre_;
    bool is_complete_;
    bool is_leaf_;
    int nvertex_, nchild_, nneigh_, ninter_, nbasis0_, nbasis1_;
    double extent_;
    std::vector<std::shared_ptr<const VertexSP>> vertex_;
    std::vector<std::shared_ptr<const NodeSP>> interaction_list_;
    std::vector<std::shared_ptr<const NodeSP>> neigh_;
    double thresh_;
    int id_in_tree_;


    void insert_vertex(std::shared_ptr<const VertexSP>);
    void insert_child(std::shared_ptr<const NodeSP> = NULL);
    void init();
    void get_interaction_list();
    void compute_extent();
    bool is_neighbour(std::shared_ptr<const NodeSP> shells, const int ws);
    void insert_neigh(std::shared_ptr<const NodeSP> neigh, const bool is_neigh = false, const int ws = 2);
    void sort_neigh();
    void make_interaction_list(const int ws);

    bool is_same_as_parent_;
    int rank_;
    std::vector<std::shared_ptr<const ZMatrix>> multipole_;
    std::vector<std::shared_ptr<const ZMatrix>> local_moment_;
    std::shared_ptr<const ZMatrix> local_expansion_;
    std::vector<std::shared_ptr<const ZMatrix>> child_local_expansion_;
    std::array<double, 3> compute_centre(std::array<std::shared_ptr<const Shell>, 2> shells);
    void compute_multipoles(const int lmax = ANG_HRR_END);
    void compute_local_expansions(std::shared_ptr<const Matrix> density, const int lmax);
    std::shared_ptr<const ZMatrix> compute_Coulomb(const int dim, std::shared_ptr<const Matrix> density, std::vector<double> max_den, const bool dodf = false, const double schwarz_thresh = 0.0);

  public:
    NodeSP(const std::bitset<nbit__> key = 0, const int depth = 0,
         std::shared_ptr<const NodeSP> parent = NULL, const double thresh = PRIM_SCREEN_THRESH);

    ~NodeSP() { }

    std::bitset<nbit__> key() const { return key_; }
    int depth() const { return depth_; }

    std::shared_ptr<const NodeSP> parent() const { return parent_; }
    std::array<double, 3> centre() const { return centre_; }
    double centre(const int i) const { return centre_[i]; }
    double radius() const;

    bool is_complete() const { return is_complete_; }
    bool is_leaf() const { return is_leaf_; }
    int nvertex() const { return nvertex_; }
    int nchild() const { return nchild_; }
    std::vector<std::shared_ptr<const VertexSP>> vertex() const { return vertex_; }
    std::shared_ptr<const VertexSP> vertex(const int i) const { return vertex_[i]; }
    std::shared_ptr<const NodeSP> child(const int i) const { return child_[i].lock(); }

    double extent() const { return extent_; }
    int nneigh() const { return nneigh_; }
    int nbasis0() const { return nbasis0_; }
    int nbasis1() const { return nbasis1_; }
    std::vector<std::shared_ptr<const NodeSP>> neigh() const { return neigh_; }
    std::vector<std::shared_ptr<const NodeSP>> interaction_list() const { return interaction_list_; }

    bool is_same_as_parent() const { return is_same_as_parent_; }
    int rank() const { return rank_; }
    int id_in_tree() const { return id_in_tree_; }
    std::vector<std::shared_ptr<const ZMatrix>> multipole() const { return multipole_; }
    std::vector<std::shared_ptr<const ZMatrix>> local_moment() const { return local_moment_; }
    std::shared_ptr<const ZMatrix> local_moment(const int i) const { return local_moment_[i]; }
    std::shared_ptr<const ZMatrix> multipole(const int i) const { return multipole_[i]; }
    std::shared_ptr<const ZMatrix> local_expansion() const { return local_expansion_; }
    std::vector<std::shared_ptr<const ZMatrix>> child_local_expansion() const { return child_local_expansion_; }
    std::shared_ptr<const ZMatrix> child_local_expansion(const int i) const { return child_local_expansion_[i]; }
    void print_node() const;
};

}
#endif
