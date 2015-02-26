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

static const double pisq__ = pi__ * pi__;

Node::Node(const std::bitset<nbit__> key, const int depth, std::shared_ptr<const Node> parent)
 : key_(key), depth_(depth), parent_(parent) {
  if (depth == 0)
    key_[0] = 1;

  is_leaf_ = false;
  is_complete_ = false;
  nchild_  = 0;
  nbody_   = 0;
  nneighbour_ = 0;
}


void Node::insert_vertex(shared_ptr<const Vertex> p) {

  assert(!is_complete_);
  bodies_.resize(nbody_ + 1);
  bodies_[nbody_] = p;
  ++nbody_;
}


void Node::insert_child(shared_ptr<const Node> child) {

  if (child = NULL) {
    is_leaf_ = false;
  } else {
    assert(!is_leaf_);
    children_.resize(nchild_ + 1);
    children_[nchild_] = child;
    ++nchild_;
  }
}


void Node::init() {

  is_complete_ = true;
  if (nchild_ == 0) is_leaf_ = true;

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

  compute_extent();
}


void Node::compute_extent(const double thresh) {

  std::vector<std::shared_ptr<const Shell>> shells;
  for (auto& body : bodies_) shells.insert(shells.end(), body->atom()->shells().begin(), body->atom()->shells().end());
  const array<double, 3> centre = position_; // use charge centre for now

  extent_ = 0.0;
  for (auto& ish : shells)
    for (auto& jsh : shells) {
      const vector<double> exp0 = ish->exponents();
      const vector<double> exp1 = jsh->exponents();
      array<double, 3> AB;
      AB[0] = ish->position(0) - jsh->position(0);
      AB[1] = ish->position(1) - jsh->position(1);
      AB[2] = ish->position(2) - jsh->position(2);
      const double rsq = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
      const double lnthresh = log(thresh);

      for (auto& expi0 : exp0) {
        for (auto& expi1 : exp1) {
          const double cxp_inv = 1.0 / (expi0 + expi1);
          const double expi01 = expi0 * expi1;
          const double lda_kl = sqrt((- lnthresh - expi01 * rsq * cxp_inv + 0.75 * log(4.0 * expi01 / pisq__)) * cxp_inv);

          array<double, 3> tmp;
          tmp[0] = (ish->position(0) * expi0 + jsh->position(0) * expi1) * cxp_inv - centre[0];
          tmp[1] = (ish->position(1) * expi0 + jsh->position(1) * expi1) * cxp_inv - centre[1];
          tmp[2] = (ish->position(2) * expi0 + jsh->position(2) * expi1) * cxp_inv - centre[2];

          const double extent0 = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]) + lda_kl;
          if (extent0 > extent_) extent_ = extent0;
        }
      }
    }
}


void Node::insert_neighbour(shared_ptr<const Node> neigh, const bool is_neighbour, const int ws) {

  assert(neigh->depth() == depth_);
  if (!is_neighbour) {
    array<double, 3> r12;
    r12[0] = position_[0] - neigh->position(0);
    r12[1] = position_[1] - neigh->position(1);
    r12[2] = position_[2] - neigh->position(2);
    const double r = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
    if (r < (1.0 + ws) * (extent_ + neigh->extent())) {
      neighbour_.resize(nneighbour_ + 1);
      neighbour_[nneighbour_] = neigh;
      ++nneighbour_;
    }
  } else {
    neighbour_.resize(nneighbour_ + 1);
    neighbour_[nneighbour_] = neigh;
    ++nneighbour_;
  }
}
