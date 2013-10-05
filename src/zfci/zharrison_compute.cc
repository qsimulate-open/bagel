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
#include <src/fci/hztasks.h>
#include <src/smith/prim_op.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

/* Implementing the method as described by Harrison and Zarrabian */
shared_ptr<RelZDvec> ZHarrison::form_sigma(shared_ptr<const RelZDvec> ccvec, shared_ptr<const RelMOFile> jop, const vector<int>& conv) const {
  const int ij = norb_*norb_;
  auto sigmavec = make_shared<RelZDvec>(space_, nstate_);

  // diagonal part
  for (auto& isp : space_->detmap()) { 
    shared_ptr<const Determinants> base_det = isp.second; 

    const bool noab = (base_det->nelea() == 0 || base_det->neleb() == 0);
    shared_ptr<const Determinants> int_det = noab ? shared_ptr<const Determinants>() : int_space_->finddet(base_det->nelea()-1, base_det->neleb()-1);

    /* d and e are only used in the alpha-beta case and exist in the (nalpha-1)(nbeta-1) spaces */
    shared_ptr<ZDvec> d, e;
    if (!noab) {
      d = make_shared<ZDvec>(int_det, ij);
      e = make_shared<ZDvec>(int_det, ij);
    }

    for (int istate = 0; istate != nstate_; ++istate) {
      Timer pdebug(2);
      if (conv[istate]) continue;
      shared_ptr<const ZCivec> cc = ccvec->find(base_det)->data(istate);
      shared_ptr<ZCivec> sigma = sigmavec->find(base_det)->data(istate);

      // (taskaa)
      sigma_aa(cc, sigma, jop);
      pdebug.tick_print("taskaa");

      // (taskbb)
      sigma_bb(cc, sigma, jop);
      pdebug.tick_print("taskbb");

      if (!noab) {
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
    }
  }
  return sigmavec;
}


void ZHarrison::sigma_aa(shared_ptr<const ZCivec> cc, shared_ptr<ZCivec> sigma, shared_ptr<const RelMOFile> jop, const bool beta) const {
  assert(cc->det() == sigma->det());

  shared_ptr<const Determinants> det = cc->det();

  bitset<2> bit2;
  bitset<4> bit4;
  if (beta) { bit2 = ~bit2; bit4 = ~bit4; }

  shared_ptr<const ZMatrix> h1 = jop->mo1e(bit2);
  auto h2 = make_shared<ZMatrix>(*jop->mo2e(bit4));
  SMITH::sort_indices<1,0,2,3,1,1,-1,1>(jop->mo2e(bit4)->data(), h2->data(), norb_, norb_, norb_, norb_);

  TaskQueue<HZTaskAA<complex<double>>> tasks(det->lena());

  const int lb = cc->lenb();
  complex<double>* target = sigma->data();
  for (auto aiter = det->stringa().begin(); aiter != det->stringa().end(); ++aiter, target+=lb)
    tasks.emplace_back(cc, *aiter, target, h1->data(), h2->data());

  tasks.compute();
}


void ZHarrison::sigma_bb(shared_ptr<const ZCivec> cc, shared_ptr<ZCivec> sigma, shared_ptr<const RelMOFile> jop) const {
  shared_ptr<const ZCivec> cc_trans = cc->transpose();
  auto sig_trans = make_shared<ZCivec>(cc_trans->det());

  sigma_aa(cc_trans, sig_trans, jop, /*beta*/true);

  sigma->ax_plus_y(1.0, *sig_trans->transpose(sigma->det()));
}


void ZHarrison::sigma_2ab_1(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  shared_ptr<const Determinants> bdet = cc->det(); // base
  shared_ptr<const Determinants> tdet = d->det();  // target

  const int lbt = tdet->lenb();
  const int lbs = bdet->lenb();
  const complex<double>* source_base = cc->data();

  TaskQueue<HZTaskAB1<complex<double>>> tasks(norb_*norb_);

  for (int k = 0; k < norb_; ++k) { // alpha
    for (int l = 0; l < norb_; ++l) { // beta
      complex<double>* target_base = d->data(l*norb_ + k)->data();
      tasks.emplace_back(tdet, lbs, source_base, target_base, k, l);
    }
  }

  tasks.compute();
}


void ZHarrison::sigma_2ab_2(shared_ptr<ZDvec> d, shared_ptr<ZDvec> e, shared_ptr<const RelMOFile> jop) const {
  const int ij = d->ij();
  const int lenab = d->lena()*d->lenb();
  ZMatrix tmp(*jop->mo2e(bitset<4>("1010")));
  SMITH::sort_indices<1,0,2,3,1,1,-1,1>(jop->mo2e(bitset<4>("0110"))->data(), tmp.data(), norb_, norb_, norb_, norb_);

  zgemm3m_("n", "t", lenab, ij, ij, 1.0, d->data(), lenab, tmp.data(), ij, 0.0, e->data(), lenab);
}


void ZHarrison::sigma_2ab_3(shared_ptr<ZCivec> sigma, shared_ptr<ZDvec> e) const {
  const shared_ptr<const Determinants> base_det = sigma->det();
  const shared_ptr<const Determinants> int_det = e->det();

  const int lbt = base_det->lenb();
  const int lbs = int_det->lenb();
  complex<double>* target_base = sigma->data();

  for (int i = 0; i < norb_; ++i) { // alpha
    for (int j = 0; j < norb_; ++j) { // beta
      const complex<double>* source_base = e->data(j*norb_ + i)->data();
      for (auto& aiter : int_det->phiupa(i)) {
        complex<double>* target = target_base + aiter.target*lbt;
        const complex<double>* source = source_base + aiter.source*lbs;
        for (auto& biter : int_det->phiupb(j)) {
          const double sign = aiter.sign * biter.sign;
          target[biter.target] += sign * source[biter.source];
        }
      }
    }
  }
}
