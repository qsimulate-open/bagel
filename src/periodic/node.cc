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
#include <src/periodic/multipolebatch.h>
#include <src/periodic/localexpansion.h>
#include <src/integral/rys/eribatch.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;

Node::Node(const std::bitset<nbit__> key, const int depth, std::shared_ptr<const Node> parent, const double thresh)
 : key_(key), depth_(depth), parent_(parent), thresh_(thresh) {
  if (depth == 0)
    key_[0] = 1;

  is_leaf_ = false;
  is_complete_ = false;
  nchild_  = 0;
  nbody_   = 0;
  nneighbour_ = 0;
  ninter_ = 0;
}


void Node::insert_vertex(shared_ptr<const Vertex> p) {

  assert(!is_complete_);
  bodies_.resize(nbody_ + 1);
  bodies_[nbody_] = p;

  ++nbody_;
}


void Node::insert_child(shared_ptr<const Node> child) {

  if (child) {
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
  nbasis_ = 0;
  double sum = 0.0;
  for (auto& body : bodies_) {
    position_[0] += body->atom()->atom_charge() * body->position(0);
    position_[1] += body->atom()->atom_charge() * body->position(1);
    position_[2] += body->atom()->atom_charge() * body->position(2);
    sum += body->atom()->atom_charge();
    nbasis_ += body->atom()->nbasis();
  }
  position_[0] /= sum;
  position_[1] /= sum;
  position_[2] /= sum;

  compute_extent(thresh_);
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
    const double r = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
    const double tmp = pow((1.0 + ws) * (extent_ + neigh->extent()), 2);
    if (r < tmp) {
      neighbour_.resize(nneighbour_ + 1);
      neighbour_[nneighbour_] = neigh;
      ++nneighbour_;
    } else {
      interaction_list_.resize(ninter_ + 1);
      interaction_list_[ninter_] = neigh;
      ++ninter_;
    }
  } else {
    neighbour_.resize(nneighbour_ + 1);
    neighbour_[nneighbour_] = neigh;
    ++nneighbour_;
  }
}


void Node::make_interaction_list(const int ws) {

  ninter_ = 0;
  for (auto& pneighbour : parent_->neighbour()) {
    for (int n = 0; n != pneighbour->nchild(); ++n) {
      shared_ptr<const Node> child = pneighbour->children(n);
      array<double, 3> r12;
      r12[0] = position_[0] - child->position(0);
      r12[1] = position_[1] - child->position(1);
      r12[2] = position_[2] - child->position(2);
      const double r = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      if (r > (1.0 + ws) * (extent_ + child->extent())) {
        interaction_list_.resize(ninter_ + 1);
        interaction_list_[ninter_] = child;
        ++ninter_;
      }
    }
  }
}


void Node::compute_multipoles(const int lmax) {

  const int nmultipole = (lmax + 1) * (lmax + 1);
  multipoles_.resize(nmultipole);
  vector<ZMatrix> multipoles(nmultipole);
  for (int i = 0; i != nmultipole; ++i) {
    ZMatrix tmp(nbasis_, nbasis_);
    tmp.zero();
    multipoles[i] = tmp;
  }

  if (is_leaf_) { //compute multipole integrals
    size_t ob0 = 0;
    for (auto& a0 : bodies_) {
      for (auto& b0 : a0->atom()->shells()) {
        size_t ob1 = 0;
        for (auto& a1 : bodies_) {
          for (auto& b1 : a1->atom()->shells()) {
            MultipoleBatch mpole(array<shared_ptr<const Shell>, 2>{{b1, b0}}, position_, lmax);
            mpole.compute();
            for (int i = 0; i != nmultipole; ++i)
              multipoles[i].add_block(1.0, ob1, ob0, b1->nbasis(), b0->nbasis(), mpole.data(i));

            ob1 += b1->nbasis();
          }
        }
        ob0 += b0->nbasis();
      }
    }
  } else { //shift children's multipoles
    int offset = 0;
    for (int n = 0; n != nchild_; ++n) {
      shared_ptr<const Node> child = children_[n].lock();
      array<double, 3> r12;
      r12[0] = position_[0] - child->position(0);
      r12[1] = position_[1] - child->position(1);
      r12[2] = position_[2] - child->position(2);
      assert(nmultipole == child->multipoles().size());
      LocalExpansion shift(r12, child->multipoles(), lmax);
      vector<shared_ptr<const ZMatrix>> moment = shift.compute_shifted_moments();

      for (int i = 0; i != nmultipole; ++i)
        multipoles[i].add_block(1.0, offset, offset, child->nbasis(), child->nbasis(), moment[i]->data());

      offset += child->nbasis();
    }
  }

  for (int i = 0; i != nmultipole; ++i)
    multipoles_[i] = make_shared<const ZMatrix>(multipoles[i]);
}


void Node::compute_local_expansions(shared_ptr<const Matrix> density, const int lmax, vector<int> offset) {

  auto out = make_shared<ZMatrix>(nbasis_, nbasis_);
  out->zero();

  for (auto& distant_node : interaction_list_) {
    array<double, 3> r12;
    r12[0] = distant_node->position(0) - position_[0];
    r12[1] = distant_node->position(1) - position_[1];
    r12[2] = distant_node->position(2) - position_[2];
    LocalExpansion lx(r12, distant_node->multipoles(), lmax);
    vector<shared_ptr<const ZMatrix>> lmoments = lx.local_moments();

    const int dimb = distant_node->nbasis();
    // need to get a sub-density matrix corresponding to distant_node
    Matrix subden(dimb, dimb);
    subden.zero();
    size_t ob0 = 0;
    for (auto& body0 : distant_node->bodies()) {
      size_t ish0 = 0;
      for (auto& b0 : body0->atom()->shells()) {
        const int offset0 = offset[body0->ishell() + ish0];
        const size_t size0 = b0->nbasis();
        ++ish0;

        size_t ob1 = 0;
        for (auto& body1 : distant_node->bodies()) {
          size_t ish1 = 0;
          for (auto& b1 : body1->atom()->shells()) {
            const int offset1 = offset[body1->ishell() + ish1];
            const size_t size1 = b1->nbasis();
            ++ish1;

            shared_ptr<const Matrix> tmp = density->get_submatrix(offset1, offset0, size1, size0);
            subden.copy_block(ob1, ob0, size1, size0, tmp);
            ob1 += size1;
          }
        }
        ob0 += size0;
      }
    }

    const int nmultipole = (lmax + 1) * (lmax + 1);
    for (int i = 0; i != nmultipole; ++i) {
      complex<double> contract = 0.0;
      for (int j = 0; j != dimb; ++j)
        for (int k = 0; k != dimb; ++k)
          contract += lmoments[i]->element(k, j) * subden.element(k, j);
      *out += contract * *multipoles_[i];
    }
  }

  local_expansion_ = out;
}


shared_ptr<const ZMatrix> Node::compute_Coulomb(shared_ptr<const Matrix> density, const int lmax, vector<int> offset) {

  assert(is_leaf());
  auto out = make_shared<ZMatrix>(density->ndim(), density->mdim());
  out->zero();

  ZMatrix tmp(nbasis_, nbasis_);
  tmp.zero();
  const int nmultipole = (lmax + 1) * (lmax + 1);
  for (auto& distant_node : interaction_list_) {
    array<double, 3> r12;
    r12[0] = distant_node->position(0) - position_[0];
    r12[1] = distant_node->position(1) - position_[1];
    r12[2] = distant_node->position(2) - position_[2];
    LocalExpansion local(r12, distant_node->multipoles(), lmax);
    vector<shared_ptr<const ZMatrix>> lmoments = local.local_moments();

    const int dimb = distant_node->nbasis();
    Matrix subden(dimb, dimb);
    subden.zero();
    size_t ob0 = 0;
    for (auto& body0 : distant_node->bodies()) {
      size_t ish0 = 0;
      for (auto& b0 : body0->atom()->shells()) {
        const int offset0 = offset[body0->ishell() + ish0];
        const size_t size0 = b0->nbasis();
        ++ish0;

        size_t ob1 = 0;
        for (auto& body1 : distant_node->bodies()) {
          size_t ish1 = 0;
          for (auto& b1 : body1->atom()->shells()) {
            const int offset1 = offset[body1->ishell() + ish1];
            const size_t size1 = b1->nbasis();
            ++ish1;

            shared_ptr<const Matrix> tmp = density->get_submatrix(offset1, offset0, size1, size0);
            subden.copy_block(ob1, ob0, size1, size0, tmp);

            ob1 += size1;
          }
        }
        ob0 += size0;
      }
    }

    // translate local at box centre to particle positions
    for (auto& body : bodies_) {
      r12[0] = position_[0] - body->position(0);
      r12[1] = position_[1] - body->position(1);
      r12[2] = position_[2] - body->position(2);
      LocalExpansion local_shift(r12, lmoments, lmax);
      vector<shared_ptr<const ZMatrix>> new_lmoments = local_shift.compute_shifted_local_expansions();

      for (int i = 0; i != nmultipole; ++i) {
        complex<double> contract = 0.0;
        for (int j = 0; j != dimb; ++j)
          for (int k = 0; k != dimb; ++k)
            contract += new_lmoments[i]->element(k, j) * subden.element(k, j);

        tmp += contract * *multipoles_[i];
      }
    }
  }

  // add FF local expansions to coulomb matrix
  size_t offset0 = 0;
  for (auto& body : bodies_) {
    const size_t nbas = body->atom()->nbasis();
    ZMatrix sublocal = *(tmp.get_submatrix(offset0, offset0, nbas, nbas));

    const size_t offset1 =  offset[body->ishell()];
    out->add_block(1.0, offset1, offset1, nbas, nbas, sublocal.data());

    offset0 += nbas;
  }

  // compute near-field interactions using direct integration and add to far field
  vector<shared_ptr<const Shell>> basis;
  vector<int> new_offset;
  for (auto& close_node : neighbour_) {
    for (auto& close_body : close_node->bodies()) {
      const vector<shared_ptr<const Shell>> tmp = close_body->atom()->shells();
      basis.insert(basis.end(), tmp.begin(), tmp.end());
      const int nshell = close_body->atom()->nshell();
      vector<int> tmpoff(nshell);
      for (int i = 0; i != nshell; ++i)
        tmpoff[i] = offset[close_body->ishell() + i];

      new_offset.insert(new_offset.end(), tmpoff.begin(), tmpoff.end());
    }
  }

  const size_t size = basis.size();
  const double* density_data = density->data();

  for (int i0 = 0; i0 != size; ++i0) {
    const shared_ptr<const Shell>  b0 = basis[i0];
    const int b0offset = new_offset[i0];
    const int b0size = b0->nbasis();

    for (int i1 = 0; i1 != size; ++i1) {
      const shared_ptr<const Shell>  b1 = basis[i1];
      const int b1offset = new_offset[i1];
      const int b1size = b1->nbasis();

      for (int i2 = 0; i2 != size; ++i2) {
        const shared_ptr<const Shell>  b2 = basis[i2];
        const int b2offset = new_offset[i2];
        const int b2size = b2->nbasis();

        for (int i3 = 0; i3 != size; ++i3) {
          const shared_ptr<const Shell>  b3 = basis[i3];
          const int b3offset = new_offset[i3];
          const int b3size = b3->nbasis();

          if ((b0offset + b0size + b1offset + b1size) < (b2offset + b3offset)) continue;
          if ((b0offset + b0size + b2offset + b2size) < (b1offset + b3offset)) continue;

          array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
          ERIBatch eribatch(input, 0.0);
          eribatch.compute();
          const double* eridata = eribatch.data();

          for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
            const int j0n = j0 * density->ndim();
            for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
              const int j1n = j1 * density->ndim();
              for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
                const int j2n = j2 * density->ndim();
                for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                  const int j3n = j3 * density->ndim();

                  double eri = *eridata / nneighbour_;
                  if (j0 + j1 <  j2 + j3 || j0 + j2 <  j1 + j3) continue;
                  if (j0 + j1 == j2 + j3) eri *= 0.5;
                  if (j0 + j2 == j1 + j3) eri *= 0.5;

                  out->element(j1, j0) += density_data[j3n + j2] * eri;
                  out->element(j3, j2) += density_data[j1n + j0] * eri;
                  out->element(j0, j1) += density_data[j2n + j3] * eri;
                  out->element(j2, j3) += density_data[j0n + j1] * eri;
                }
              }
            }
          }

        }
      }
    }
  }

  return out;
}
