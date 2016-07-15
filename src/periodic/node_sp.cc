//
// BAGEL - Parallel electron correlation program.
// Filename: node_sp.cc
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


#include <src/util/f77.h>
#include <src/periodic/node_sp.h>
#include <src/periodic/multipolebatch.h>
#include <src/integral/os/overlapbatch.h>
#include <src/periodic/localexpansion.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;

NodeSP::NodeSP(const bitset<nbit__> key, const int depth, shared_ptr<const NodeSP> parent, const double thresh)
 : key_(key), depth_(depth), parent_(parent), thresh_(thresh) {
  if (depth == 0)
    key_[0] = 1;

  is_leaf_ = false;
  is_complete_ = false;
  nchild_  = 0;
  nvertex_   = 0;
  nneigh_ = 0;
  ninter_ = 0;
  is_same_as_parent_ = false;
  rank_ = 0;
}


void NodeSP::insert_vertex(shared_ptr<const VertexSP> p) {

  assert(!is_complete_);
  vertex_.resize(nvertex_ + 1);
  vertex_[nvertex_] = p;

  ++nvertex_;
}


void NodeSP::insert_child(shared_ptr<const NodeSP> child) {

  if (child) {
    assert(!is_leaf_);
    child_.resize(nchild_ + 1);
    child_[nchild_] = child;
    ++nchild_;
  }
}


void NodeSP::init() {

  is_complete_ = true;
  iself_ = -1;
  extent_ = -1.0;
  if (nchild_ == 0) is_leaf_ = true;

  if (!is_leaf_)
    child_local_expansion_.resize(nchild_);

  nbasis_ = 0;
  centre_ = {{0.0, 0.0, 0.0}};
  for (auto& v : vertex_) {
    centre_[0] += v->centre(0);
    centre_[1] += v->centre(1);
    centre_[2] += v->centre(2);
    nbasis_ += v->nbasis();
  }

  array<double, 3> r12;
  r12[0] = centre_[0] - parent_->centre(0);
  r12[1] = centre_[1] - parent_->centre(1);
  r12[2] = centre_[2] - parent_->centre(2);

  const double r = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
  if (r <= numerical_zero__) {
    is_same_as_parent_ = true;
    rank_ = parent_->rank();
  } else {
    rank_ = parent_->rank() + 1;
  }

  compute_extent();
}


bool NodeSP::is_neighbour(shared_ptr<const NodeSP> nd, const int ws) {

  array<double, 3> rvec = nd->centre();
  assert(extent_ > 0 && nd->extent() > 0);

  array<double, 3> v;
  v[0] = centre_[0] - rvec[0];
  v[1] = centre_[1] - rvec[1];
  v[2] = centre_[2] - rvec[2];
  const double r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return (r <= (1.0 + ws) * (extent_ + nd->extent()));
}


void NodeSP::compute_extent() {

  extent_ = 0.0;
  vector<shared_ptr<const Shell>> shells;
  for (auto& v : vertex_)
    if (v->extent() > extent_) extent_ = v->extent();
}


void NodeSP::insert_neigh(shared_ptr<const NodeSP> neigh, const bool is_neigh, const int ws) {

  assert(neigh->depth() == depth_);
  if (!is_neigh) {
    const bool check = is_neighbour(neigh, ws);
    array<double, 3> r12;
    r12[0] = centre_[0] - neigh->centre(0);
    r12[1] = centre_[1] - neigh->centre(1);
    r12[2] = centre_[2] - neigh->centre(2);
    const double r = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
    if (check) {
      neigh_.resize(nneigh_ + 1);
      neigh_[nneigh_] = neigh;
      if (r < numerical_zero__)
        iself_ = nneigh_;
      ++nneigh_;
    } else {
      interaction_list_.resize(ninter_ + 1);
      interaction_list_[ninter_] = neigh;
      ++ninter_;
    }
  } else {
    neigh_.resize(nneigh_ + 1);
    neigh_[nneigh_] = neigh;
    ++nneigh_;
  }
}


void NodeSP::make_interaction_list(const int ws) {

  ninter_ = 0;
  for (auto& pneigh : parent_->neigh()) {
    for (int n = 0; n != pneigh->nchild(); ++n) {
      shared_ptr<const NodeSP> child = pneigh->child(n);
      assert(child->depth() == depth_);
      bool is_close = false;
      for (int i = 0; i != nneigh_; ++i) {
        if (child->key() == neigh_[i]->key()) {
          is_close = true;
          break;
        }
      }
      if (!is_close) {
        interaction_list_.resize(ninter_ + 1);
        interaction_list_[ninter_] = child;
        ++ninter_;
      }
    }
  }
}


void NodeSP::compute_multipoles(const int lmax) {

  const int nmult = (lmax + 1) * (lmax + 1);
  multipole_.resize(nmult);
  vector<shared_ptr<ZMatrix>> multipole(nmult);
  for (int i = 0; i != nmult; ++i)
    multipole[i] = make_shared<ZMatrix>(nbasis_, nbasis_);

  bool skip = false;
  if (is_leaf_) { // shift sp's multipoles
    size_t offset = 0;
    for (auto& v : vertex_) {
      array<double, 3> r12;
      r12[0] = centre_[0] - v->centre(0);
      r12[1] = centre_[1] - v->centre(1);
      r12[2] = centre_[2] - v->centre(2);
      assert(nmult == v->nmult());
      LocalExpansion shift(r12, v->multipole(), lmax);
      vector<shared_ptr<const ZMatrix>> moment = shift.compute_shifted_multipoles();

      for (int i = 0; i != nmult; ++i)
        multipole[i]->copy_block(offset, offset, v->nbasis(), v->nbasis(), moment[i]->data());

      offset += v->nbasis();
    }
  } else { // shift children's multipoles
    shared_ptr<const NodeSP> child = child_[0].lock();
    if (child->is_same_as_parent())
      skip = true;

    if (!skip) {
      size_t offset = 0;
      for (int n = 0; n != nchild_; ++n) {
        shared_ptr<const NodeSP> child = child_[n].lock();
        array<double, 3> r12;
        r12[0] = centre_[0] - child->centre(0);
        r12[1] = centre_[1] - child->centre(1);
        r12[2] = centre_[2] - child->centre(2);
        assert(nmult == child->multipole().size());
        LocalExpansion shift(r12, child->multipole(), lmax);
        vector<shared_ptr<const ZMatrix>> moment = shift.compute_shifted_multipoles();

        for (int i = 0; i != nmult; ++i)
          multipole[i]->copy_block(offset, offset, child->nbasis(), child->nbasis(), moment[i]->data());

        offset += child->nbasis();
      }
    }
  }

  if (!skip) {
    for (int i = 0; i != nmult; ++i)
      multipole_[i] = multipole[i];
  } else {
    shared_ptr<const NodeSP> child = child_[0].lock();
    for (int i = 0; i != nmult; ++i)
      multipole_[i] = child->multipole(i);
  }
}


void NodeSP::compute_local_expansions(shared_ptr<const Matrix> density, const int lmax) {

  const int nmult = (lmax + 1) * (lmax + 1);

  vector<int> l_map(nmult);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l)
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      l_map[cnt] = l;

  auto lrs = make_shared<ZMatrix>(nbasis_, nbasis_);

  if (density) {
    for (auto& distant_node : interaction_list_) { // M2L
      array<double, 3> r12;
      r12[0] = centre_[0] - distant_node->centre(0);
      r12[1] = centre_[1] - distant_node->centre(1);
      r12[2] = centre_[2] - distant_node->centre(2);
      LocalExpansion lx(r12, distant_node->multipole(), lmax);
      vector<shared_ptr<const ZMatrix>> lmoments = lx.compute_local_moments();

      // get sub-density matrix D_tu (tu=far)
      const int dimb = distant_node->nbasis();
      Matrix den_tu(dimb, dimb);
      den_tu.zero();
      size_t ob0 = 0;
      size_t ob1 = 0;
      for (auto& v : distant_node->vertex()) {

        shared_ptr<const Shell> sh0 = v->shell0();
        shared_ptr<const Shell> sh1 = v->shell1();

        const int offset0 = v->offset(0);
        const int size0 = sh0->nbasis();
        const int offset1 = v->offset(1);
        const int size1 = sh1->nbasis();

        shared_ptr<const Matrix> tmp = density->get_submatrix(offset1, offset0, size1, size0);
        den_tu.copy_block(ob1, ob0, size1, size0, tmp);
        ob1 += size1;
        ob0 += size0;
      }

      for (int i = 0; i != nmult; ++i) {
        const double contract = lmoments[i]->get_real_part()->dot_product(den_tu);
        *lrs += pow(-1.0, l_map[i]) * contract * *multipole_[i];
      }

    }
  }

  local_expansion_ = make_shared<const ZMatrix>(*lrs);
}


shared_ptr<const ZMatrix> NodeSP::compute_Coulomb(const int dim, shared_ptr<const Matrix> density, const bool dodf, const vector<double> schwarz, const double schwarz_thresh) {

  assert(is_leaf());
  auto out = make_shared<ZMatrix>(dim, dim);
  out->zero();
  const size_t ndim = density->ndim();

  // add FF local expansions to coulomb matrix
  #if 0
  size_t ob0 = 0;
  size_t ob1 = 0;
  for (auto& v : vertex_) {
    shared_ptr<const Shell> sh0 = v->shell0();
    shared_ptr<const Shell> sh1 = v->shell1();
    const int size0 = sh0->nbasis();
    const int size1 = sh1->nbasis();

    ZMatrix sublocal = *(local_expansion_->get_submatrix(ob1, ob0, size1, size0));
    out->copy_block(v->offset(1), v->offset(0), size1, size0, sublocal.data());

    ob0 += size0;
    ob1 += size1;
  }
  #endif

  // compute near-field interactions using direct integration and add to far field

  if (!dodf && density) {
    const double* density_data = density->data();
    const int nshell = sqrt(schwarz.size());

    for (auto& v01 : vertex_) {
      shared_ptr<const Shell> b0 = v01->shell0();
      const int i0 = v01->shell_ind(0);
      const int b0offset = v01->offset(0);
      const int b0size = b0->nbasis();

      shared_ptr<const Shell> b1 = v01->shell1();
      const int i1 = v01->shell_ind(1);
      const int b1offset = v01->offset(1);
      const int b1size = b1->nbasis();

      for (auto& neigh : neigh_) {
        for (auto& v23 : neigh->vertex()) {
          shared_ptr<const Shell> b2 = v23->shell0();
          const int i2 = v23->shell_ind(0);
          const int b2offset = v23->offset(0);
          const int b2size = b2->nbasis();

          shared_ptr<const Shell> b3 = v23->shell1();
          const int i3 = v23->shell_ind(1);
          const int b3offset = v23->offset(1);
          const int b3size = b3->nbasis();

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
                  double eri = *eridata;
                  out->element(j0, j1) += density_data[j2n + j3] * eri;
                  out->element(j0, j2) -= density_data[j1n + j3] * eri * 0.5;
                }
              }
            }
          }

        }
      }
    }
  } else if (dodf && density) {
    throw runtime_error("No DF yet :-(");
  }

  return out;
}


void NodeSP::print_node() const {

  cout << "*** Node at " << setprecision(5) << centre_[0] << "   " << centre_[1] << "   " << centre_[2] << endl;
  cout << "has " << nvertex_ << " vertex: " << endl;
  for (auto& vertex : vertex_)
    cout << setw(5) << setprecision(5) << setw(10) << vertex->centre(0) << "   "
                                       << setw(10) << vertex->centre(1) << "   "
                                       << setw(10) << vertex->centre(2) << endl;
  cout << "and " << nneigh_ << " neighs" << endl;
  for (auto& neigh : neigh_)
    cout << setprecision(5) << setw(10) << neigh->centre(0) << "  "
                            << setw(10) << neigh->centre(1) << "  "
                            << setw(10) << neigh->centre(2) << endl;
}
