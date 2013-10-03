//
// BAGEL - Parallel electron correlation program.
// Filename: zharrison_compute.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/zfci/zharrison.h>
#include <src/math/davidson.h>
#include <src/util/taskqueue.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

/* Implementing the method as described by Harrison and Zarrabian */
shared_ptr<RelZDvec> ZHarrison::form_sigma(shared_ptr<const RelZDvec> ccvec, shared_ptr<const RelMOFile> jop, const vector<int>& conv) const {
  const int ij = norb_*norb_;
  auto sigmavec = make_shared<RelZDvec>(space_, nstate_);

#if 0
  shared_ptr<Determinants> base_det = space_->basedet();
  shared_ptr<Determinants> int_det = space_->finddet(nelea_-1,neleb_-1);

  /* d and e are only used in the alpha-beta case and exist in the (nalpha-1)(nbeta-1) spaces */
  auto d = make_shared<Dvec>(int_det, ij);
  auto e = make_shared<Dvec>(int_det, ij);

  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(2);
    if (conv[istate]) continue;
    shared_ptr<const Civec> cc = ccvec->data(istate);
    shared_ptr<Civec> sigma = sigmavec->data(istate);

    // (taskaa)
    sigma_aa(cc, sigma, jop);
    pdebug.tick_print("taskaa");

    // (taskbb)
    sigma_bb(cc, sigma, jop);
    pdebug.tick_print("taskbb");

    // (2ab) alpha-beta contributions
    /* Resembles more the Knowles & Handy FCI terms */
    d->zero();

    sigma_2ab_1(cc, d);
    pdebug.tick_print("task2ab-1");

    sigma_2ab_2(d, e, jop);
    pdebug.tick_print("task2ab-2");

    sigma_2ab_3(sigma, e);
    pdebug.tick_print("task2ab-3");
  }
#endif
throw logic_error("end of sigma");
  return sigmavec;
}

void ZHarrison::sigma_aa(shared_ptr<const RelZDvec> cc, shared_ptr<RelZDvec> sigma, shared_ptr<const RelMOFile> jop) const {
#if 0
  assert(cc->det() == sigma->det());

  shared_ptr<const Determinants> det = cc->det();
  const int lb = cc->lenb();

  const int norb = norb_;
  unique_ptr<double[]> h1(new double[norb*norb]);
  for (int i = 0, ij = 0; i < norb; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      h1[i + j*norb] = h1[j + i*norb] = jop->mo1e(ij);
    }
  }

  unique_ptr<double[]> h2(new double[norb*norb*norb*norb]);
  for (int i = 0, ijkl = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      for (int k = 0; k < norb; ++k) {
        for (int l = 0; l < norb; ++l, ++ijkl) {
          h2[ijkl] = jop->mo2e_hz(l,k,j,i) - jop->mo2e_hz(k,l,j,i);
        }
      }
    }
  }

  TaskQueue<HZTaskAA> tasks(det->lena());

  double* target = sigma->data();
  for (auto aiter = det->stringa().begin(); aiter != det->stringa().end(); ++aiter, target+=lb)
    tasks.emplace_back(cc, *aiter, target, h1.get(), h2.get());

  tasks.compute();
#endif
}

void ZHarrison::sigma_bb(shared_ptr<const RelZDvec> cc, shared_ptr<RelZDvec> sigma, shared_ptr<const RelMOFile> jop) const {
#if 0
  shared_ptr<const Civec> cc_trans = cc->transpose();
  auto sig_trans = make_shared<Civec>(cc_trans->det());

  sigma_aa(cc_trans, sig_trans, jop);

  sigma->ax_plus_y(1.0, *sig_trans->transpose(sigma->det()));
#endif
}

void ZHarrison::sigma_2ab_1(shared_ptr<const RelZDvec> cc, shared_ptr<RelZDvec> d) const {
#if 0
  const int norb = norb_;

  shared_ptr<const Determinants> bdet = cc->det(); // base
  shared_ptr<const Determinants> tdet = d->det();  // target

  const int lbt = tdet->lenb();
  const int lbs = bdet->lenb();
  const double* source_base = cc->data();

  TaskQueue<HZTaskAB1> tasks(norb*norb);

  for (int k = 0; k < norb; ++k) {
    for (int l = 0; l < norb; ++l) {
      double* target_base = d->data(k*norb + l)->data();
      tasks.emplace_back(tdet, lbs, source_base, target_base, k, l);
    }
  }

  tasks.compute();
#endif
}

void ZHarrison::sigma_2ab_2(shared_ptr<RelZDvec> d, shared_ptr<RelZDvec> e, shared_ptr<const RelMOFile> jop) const {
#if 0
  const int la = d->lena();
  const int lb = d->lenb();
  const int ij = d->ij();
  const int lenab = la*lb;
  dgemm_("n", "n", lenab, ij, ij, 1.0, d->data(), lenab, jop->mo2e_ptr(), ij, 0.0, e->data(), lenab);
#endif
}

void ZHarrison::sigma_2ab_3(shared_ptr<RelZDvec> sigma, shared_ptr<RelZDvec> e) const {
#if 0
  const shared_ptr<Determinants> base_det = space_->basedet();
  const shared_ptr<Determinants> int_det = space_->finddet(nelea_-1,neleb_-1);

  const int norb = norb_;
  const int lbt = base_det->lenb();
  const int lbs = int_det->lenb();
  double* target_base = sigma->data();

  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      const double* source_base = e->data(i*norb + j)->data();
      for (auto& aiter : int_det->phiupa(i)) {
        double *target = target_base + aiter.target*lbt;
        const double *source = source_base + aiter.source*lbs;
        for (auto& biter : int_det->phiupb(j)) {
          const double sign = aiter.sign * biter.sign;
          target[biter.target] += sign * source[biter.source];
        }
      }
    }
  }
#endif
}
