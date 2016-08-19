//
// BAGEL - Parallel electron correlation program.
// Filename: box.cc
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
#include <src/periodic/box.h>
#include <src/periodic/multipolebatch.h>
#include <src/integral/os/overlapbatch.h>
#include <src/periodic/localexpansion.h>
#include <src/integral/rys/eribatch.h>

using namespace bagel;
using namespace std;

static const double pisq__ = pi__ * pi__;

void Box::init() {

  centre_ = {{0, 0, 0}};
  nbasis0_ = 0;   nbasis1_ = 0;
  for (auto& i : sp_) {
    centre_[0] += i->centre(0);
    centre_[1] += i->centre(1);
    centre_[2] += i->centre(2);
    nbasis0_ += i->nbasis0();
    nbasis1_ += i->nbasis1();
  }
  centre_[0] /= nsp();
  centre_[1] /= nsp();
  centre_[2] /= nsp();

  extent_ = 0;
  for (auto& i : sp_) {
    if (i->schwarz() < 1e-20) continue;
    double tmp = 0;
    for (int j = 0; j != 3; ++j)
      tmp += pow(i->centre(j)-centre_[j], 2.0);
    const double ei = sqrt(tmp) + i->extent();
    assert(ei > 0);
    if (extent_ < ei) extent_ = ei;
  }
}


void Box::insert_sp(vector<shared_ptr<const ShellPair>> sp) {

  const int nsp = sp_.size();
  sp_.resize(nsp + sp.size());
  for (int i = 0; i != sp.size(); ++i)
    sp_[i + nsp] = sp[i];

}


void Box::insert_child(shared_ptr<const Box> child) {

  const int nchild = child_.size();
  if (child) {
    child_.resize(nchild + 1);
    child_[nchild] = child;
  }
}


void Box::get_neigh(vector<shared_ptr<Box>> box, const int ws) {

  neigh_.resize(box.size());
  int nn = 0;
  for (auto& b : box) {
    if (b->rank() == rank_ && is_neigh(b, ws)) {
      neigh_[nn] = b;
      ++nn;
    }
  }
  neigh_.resize(nn);
}


bool Box::is_neigh(shared_ptr<const Box> box, const int ws) const {

  double rr = 0;
  for (int i = 0; i != 3; ++i)
    rr += pow(centre_[i] - box->centre(i), 2);

  const bool out = (sqrt(rr) < (1+ws)*(extent_ + box->extent()));
  return out;
}


void Box::get_inter(vector<shared_ptr<Box>> box, const int ws) {

  inter_.resize(box.size());
  int ni = 0;
  for (auto& b : box) {
    if (b->rank() == rank_ && !is_neigh(b, ws) && (!parent_ || parent_->is_neigh(b->parent(), ws))) {
      inter_[ni] = b;
      ++ni;
    }
  }
  inter_.resize(ni);
}


void Box::compute_multipoles() {

  const int nmult = (lmax_ + 1) * (lmax_ + 1);

  multipole_.resize(nmult);
  vector<shared_ptr<ZMatrix>> multipole(nmult);
  for (int i = 0; i != nmult; ++i)
    multipole[i] = make_shared<ZMatrix>(nbasis1_, nbasis0_);

  if (nchild() == 0) { // leaf = shift sp's multipoles
    size_t offset0 = 0;
    size_t offset1 = 0;
    for (auto& v : sp_) {
      vector<shared_ptr<const ZMatrix> vmult(nmult) = v->multipoles(lmax_);
      array<double, 3> r12;
      r12[0] = centre_[0] - v->centre(0);
      r12[1] = centre_[1] - v->centre(1);
      r12[2] = centre_[2] - v->centre(2);
      LocalExpansion shift(r12, vmult, lmax_);
      vector<shared_ptr<const ZMatrix>> moment = shift.compute_shifted_multipoles();

      for (int i = 0; i != nmult; ++i)
        multipole[i]->copy_block(offset1, offset0, v->nbasis1(), v->nbasis0(), moment[i]->data());

      offset0 += v->nbasis0();
      offset1 += v->nbasis1();
    }
    assert(offset0 == nbasis0_);
    assert(offset1 == nbasis1_);
  } else { // shift children's multipoles
    size_t offset0 = 0;
    size_t offset1 = 0;
    for (int n = 0; n != nchild(); ++n) {
      shared_ptr<const Box> c = child(n);
      array<double, 3> r12;
      r12[0] = centre_[0] - c->centre(0);
      r12[1] = centre_[1] - c->centre(1);
      r12[2] = centre_[2] - c->centre(2);
      LocalExpansion shift(r12, c->multipole(), lmax_);
      vector<shared_ptr<const ZMatrix>> moment = shift.compute_shifted_multipoles();

      for (int i = 0; i != nmult; ++i)
        multipole[i]->copy_block(offset1, offset0, c->nbasis1(), c->nbasis0(), moment[i]->data());

      offset0 += c->nbasis0();
      offset1 += c->nbasis1();
    }
    assert(offset0 == nbasis0_);
    assert(offset1 == nbasis1_);
  }

  for (int i = 0; i != nmult; ++i)
    multipole_[i] = multipole[i];
}


void Box::compute_local_expansions(std::shared_ptr<const Matrix> density) {

}


shared_ptr<const ZMatrix> Box::compute_node_energy(shared_ptr<const Matrix> density, vector<double> max_den, const double schwarz_thresh) const {

  assert(nchild() == 0);
  auto out = make_shared<ZMatrix>(density->ndim(), density->ndim());
  out->zero();


  // FF: add local expansions

  int ninteg = 0;
  // NF: 4c integrals
  const int shift = sizeof(int) * 4;
  const double* density_data = density->data();
  const int nsh = sqrt(max_den.size());
  assert (nsh*nsh == max_den.size());

  for (auto& v01 : sp_) {
    shared_ptr<const Shell> b0 = v01->shell(1);
    const int i0 = v01->shell_ind(1);
    const int b0offset = v01->offset(1);
    const int b0size = b0->nbasis();

    shared_ptr<const Shell> b1 = v01->shell(0);
    const int i1 = v01->shell_ind(0);
    if (i1 < i0) continue;
    const int b1offset = v01->offset(0);
    const int b1size = b1->nbasis();

    const int i01 = i0 * nsh + i1;
    const double density_01 = max_den[i01] * 4.0;

    for (auto& neigh : neigh_) {
      for (auto& v23 : neigh->sp()) {
        shared_ptr<const Shell> b2 = v23->shell(1);
        const int i2 = v23->shell_ind(1);
        if (i2 < i0) continue;
        const int b2offset = v23->offset(1);
        const int b2size = b2->nbasis();

        shared_ptr<const Shell> b3 = v23->shell(0);
        const int i3 = v23->shell_ind(0);
        if (i3 < i2) continue;
        const int b3offset = v23->offset(0);
        const int b3size = b3->nbasis();

        const int i23 = i2 * nsh + i3;
        if (i23 < i01) continue;

        const double density_23 = max_den[i23] * 4.0;
        const double density_02 = max_den[i0 * nsh + i2];
        const double density_03 = max_den[i0 * nsh + i3];
        const double density_12 = max_den[i1 * nsh + i2];
        const double density_13 = max_den[i1 * nsh + i3];

        const double mulfactor = max(max(max(density_01, density_02),
                                         max(density_03, density_12)),
                                         max(density_13, density_23));
        const double integral_bound = mulfactor * v01->schwarz() * v23->schwarz();
        const bool skip_schwarz = integral_bound < schwarz_thresh;
        if (skip_schwarz) continue;

        ++ninteg;
        array<shared_ptr<const Shell>,4> input = {{b3, b2, b1, b0}};
        ERIBatch eribatch(input, mulfactor);
        eribatch.compute();
        const double* eridata = eribatch.data();

        for (int j0 = b0offset; j0 != b0offset + b0size; ++j0) {
          const int j0n = j0 * density->ndim();
          for (int j1 = b1offset; j1 != b1offset + b1size; ++j1) {
            if (j1 < j0) {
              eridata += b2size * b3size;
              continue;
            }
            const int j1n = j1 * density->ndim();
            const unsigned int nj01 = (j0 << shift) + j1;
            const double scale01 = (j0 == j1) ? 0.5 : 1.0;
            for (int j2 = b2offset; j2 != b2offset + b2size; ++j2) {
              const int j2n = j2 * density->ndim();
              for (int j3 = b3offset; j3 != b3offset + b3size; ++j3, ++eridata) {
                if (j3 < j2) continue;
                const int j3n = j3 * density->ndim();
                const unsigned int nj23 = (j2 << shift) + j3;
                if (nj23 < nj01 && i01 == i23) continue;
                const double scale23 = (j2 == j3) ? 0.5 : 1.0;
                const double scale = (nj01 == nj23) ? 0.25 : 0.5;

                const double eri = *eridata;
                const double intval = eri * scale * scale01 * scale23;
                const double intval4 = 4.0 * intval;
                out->element(j1, j0) += density_data[j2n + j3] * intval4;
                out->element(j3, j2) += density_data[j0n + j1] * intval4;
                out->element(max(j2, j0), min(j2,j0)) -= density_data[j1n + j3] * intval;
                out->element(j3, j0) -= density_data[j1n + j2] * intval;
                out->element(max(j1,j2), min(j1,j2)) -= density_data[j0n + j3] * intval;
                out->element(max(j1,j3), min(j1,j3)) -= density_data[j0n + j2] * intval;
              }
            }
          }
        }

      }
    }
  }
  //cout << "ninteg = " << ninteg << endl;

  return out;
}
