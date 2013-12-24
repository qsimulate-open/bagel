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
  auto sigmavec = make_shared<RelZDvec>(space_, nstate_);
  auto sigmavec_trans = sigmavec->clone(); // important note: the space stays the same after transposition

  for (auto& isp : space_->detmap()) {
    shared_ptr<const Determinants> base_det = isp.second;
    const int nelea = base_det->nelea();
    const int neleb = base_det->neleb();

    for (int istate = 0; istate != nstate_; ++istate) {
      if (conv[istate]) continue;

      shared_ptr<const ZCivec> cc = ccvec->find(nelea, neleb)->data(istate);
      sigma_one(cc, sigmavec, jop, istate, /*diag*/true, /*transpose*/false);

      shared_ptr<const ZCivec> cc_trans = cc->transpose();
      sigma_one(cc_trans, sigmavec_trans, jop, istate, /*diag*/false, /*transpose*/true);
    }
  }

  for (auto& isp : space_->detmap()) {
    const int nelea = isp.second->nelea();
    const int neleb = isp.second->neleb();
    shared_ptr<ZDvec> sigma = sigmavec->find(nelea, neleb);
    shared_ptr<ZDvec> sigma_trans = sigmavec_trans->find(neleb, nelea);

    // daxpy for each Civec
    for (int ist = 0; ist != nstate_; ++ist) {
      if (conv[ist]) continue;
      sigma->data(ist)->ax_plus_y(1.0, *sigma_trans->data(ist)->transpose());
    }
  }

  return sigmavec;
}


void ZHarrison::sigma_one(shared_ptr<const ZCivec> cc, shared_ptr<RelZDvec> sigmavec, shared_ptr<const RelMOFile> jop,
                          const int istate, const bool diag, const bool trans) const {
  Timer pdebug(2);

  const int ij = norb_*norb_;
  shared_ptr<const Determinants> base_det = cc->det();
  // diagonal output
  const int nelea = base_det->nelea();
  const int neleb = base_det->neleb();
  shared_ptr<ZCivec> sigma = sigmavec->find(nelea, neleb)->data(istate);

  sigma_aa(cc, sigma, jop, trans);
  pdebug.tick_print("taskaa");

  const bool noab = (base_det->nelea() == 0 || base_det->neleb() == 0);
  const bool noaa =  base_det->nelea() <= 1 || base_det->neleb()+1 > norb_;

  const bool output1 = base_det->nelea()-1 >= 0 && base_det->neleb()+1 <= norb_;
  const bool output2 = base_det->nelea()-2 >= 0 && base_det->neleb()+2 <= norb_;

  if (!noab && (diag || output1)) {
    shared_ptr<const Determinants> int_det = int_space_->finddet(nelea-1, neleb-1);
    auto d = make_shared<ZDvec>(int_det, ij);
    auto e = make_shared<ZDvec>(int_det, ij);

    // (2ab) alpha-beta contributions
    /* Resembles more the Knowles & Handy FCI terms */

    sigma_2e_annih_ab(cc, d);
    pdebug.tick_print("task2ab-1");

    if (diag) {
      // (a^+ b^+ b a) and (a^+ b^+ a b) contributions
      sigma_2e_h0101_h1001(d, e, jop);
      pdebug.tick_print("task2ab-2 (0)");

      sigma_2e_create_ab(sigma, e);
      pdebug.tick_print("task2ab-3 (0)");
    }

    if (output1) {
      // output area
      shared_ptr<ZCivec> sigma_1 = sigmavec->find(nelea-1, neleb+1)->data(istate);

      // (b^+ a) contribution
      sigma_1e_ab(cc, sigma_1, jop, trans);

      // (b^+b^+ b a) contribution
      sigma_2e_h<1,1,0,1>(d, e, jop, trans);
      pdebug.tick_print("task2ab-2 (+1)");
      sigma_2e_create_bb(sigma_1, e);
      pdebug.tick_print("task2ab-3 (+1)");
    }
  }

  if (!noaa) {
    shared_ptr<const Determinants> int_det = int_space_->finddet(nelea-2, neleb);
    auto d = make_shared<ZDvec>(int_det, ij);
    auto e = make_shared<ZDvec>(int_det, ij);

    sigma_2e_annih_aa(cc, d);
    pdebug.tick_print("task2aa-1 (2)");

    // (a^+ b^+ a a) contribution
    sigma_2e_h<0,1,0,0>(d, e, jop, trans);

    assert(neleb+1 <= norb_);
    // +1 sector
    shared_ptr<ZCivec> sigma_1 = sigmavec->find(nelea-1, neleb+1)->data(istate);
    sigma_2e_create_ab(sigma_1, e);

    // +2 sector
    if (base_det->neleb()+2 <= norb_) {
      shared_ptr<ZCivec> sigma_2 = sigmavec->find(nelea-2, neleb+2)->data(istate);
      // reusing
      sigma_2e_h<1,1,0,0>(d, e, jop, trans, 0.5);
      sigma_2e_create_bb(sigma_2, e);
    }
  }
}


void ZHarrison::sigma_aa(shared_ptr<const ZCivec> cc, shared_ptr<ZCivec> sigma, shared_ptr<const RelMOFile> jop, const bool trans) const {
  assert(cc->det()->nelea() == sigma->det()->nelea());
  assert(cc->det()->neleb() == sigma->det()->neleb());

  shared_ptr<const Determinants> det = cc->det();

  bitset<2> bit2;
  bitset<4> bit4;
  if (trans) { bit2 = ~bit2; bit4 = ~bit4; }

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


void ZHarrison::sigma_1e_ab(shared_ptr<const ZCivec> cc, shared_ptr<ZCivec> sigma, shared_ptr<const RelMOFile> jop, const bool trans) const {
  assert(cc->det()->nelea()-1 == sigma->det()->nelea());
  assert(cc->det()->neleb()+1 == sigma->det()->neleb());

  const int lbs = cc->lenb();
  const int lbt = sigma->lenb();

  bitset<2> bit2("10");
  if (trans) bit2 = ~bit2;
  shared_ptr<const ZMatrix> h1 = jop->mo1e(bit2);

  shared_ptr<const Determinants> ccdet = cc->det();
  shared_ptr<const Determinants> sigmadet = sigma->det();

  // One-electron part
  for (int i = 0; i != norb_; ++i) {
    for (auto& a : sigmadet->stringa()) {
      if (a[i]) continue;
      auto ca = a; ca.set(i);
      const complex<double>* source =    cc->data() + lbs * ccdet->lexical<0>(ca);
            complex<double>* target = sigma->data() + lbt * sigmadet->lexical<0>(a);
      const double asign = ccdet->sign<0>(ca, i);

      for (int j = 0; j != norb_; ++j) {
        for (auto& b : ccdet->stringb()) {
          if (b[j]) continue;
          auto cb = b; cb.set(j);
          const complex<double> fac = h1->element(j,i) * (sigmadet->sign<1>(b, j) * asign);
          target[sigmadet->lexical<1>(cb)] += fac * source[ccdet->lexical<1>(b)];
        }
      }

    }
  }
}


//////////////// functions for two electron annihilation ///////////////

void ZHarrison::sigma_2e_annih_aa(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  // compute |d> = a_i a_j|cc>
  const complex<double>* const source = cc->data();
  const size_t lb = cc->lenb();
  assert(lb == d->lenb());

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j != norb_; ++j) {
      if (i == j) continue;
      complex<double>* target = d->data(j+norb_*i)->data();
      for (auto& a : d->det()->stringa()) {
        if (a[i] || a[j]) continue;
        auto ca = a; ca.set(i); ca.set(j);
        const double factor = Determinants::sign(a, i, j) * (i < j ? -1.0 : 1.0);
        const size_t offas = cc->det()->lexical<0>(ca);
        const size_t offat = d->det()->lexical<0>(a);
        blas::ax_plus_y_n(factor, source+lb*offas, lb, target+lb*offat);
      }
    }
  }
}


void ZHarrison::sigma_2e_annih_ab(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  // compute |d> = b_l a_k|cc>
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


//////////////// functions for two electron creation ///////////////

void ZHarrison::sigma_2e_create_ab(shared_ptr<ZCivec> sigma, shared_ptr<const ZDvec> e) const {
  // compute |d> = a_i b_j|cc>
  const shared_ptr<const Determinants> int_det = e->det();
  const shared_ptr<const Determinants> base_det = sigma->det();

  const int lbt = sigma->det()->lenb();
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


void ZHarrison::sigma_2e_create_bb(shared_ptr<ZCivec> sigma, shared_ptr<const ZDvec> e) const {
  // compute |sigma> = b+_i b+_j|e_ij>
  const shared_ptr<const Determinants> base_det = sigma->det();
  const shared_ptr<const Determinants> int_det = e->det();

  const int lbt = base_det->lenb();
  const int lbs = int_det->lenb();
  const int la = base_det->lena();

  assert(base_det->lena() == int_det->lena());

  complex<double>* target_base = sigma->data();

  // TODO can be reduced by a factor of two by using symmetry
  // TODO not efficient code
  for (int i = 0; i < norb_; ++i) { // beta
    for (int j = 0; j < norb_; ++j) { // beta
      if (i == j) continue;
      const complex<double>* source_base = e->data(i+norb_*j)->data();

      for (size_t ia = 0; ia != la; ++ia) {
        complex<double>* target = target_base + ia*lbt;
        const complex<double>* source = source_base + ia*lbs;

        for (auto& b : int_det->stringb()) {
          if (b[i] || b[j]) continue;
          bitset<nbit__> cb = b;
          cb.set(i); cb.set(j);

          const double sign = (i < j ? 1.0 : -1.0) * Determinants::sign(b, i, j);

          target[base_det->lexical<1>(cb)] += sign * source[int_det->lexical<1>(b)];
        }

      }
    }
  }
}


//////////////// functions for multiplication of the Hamiltonian ///////////////

void ZHarrison::sigma_2e_h0101_h1001(shared_ptr<const ZDvec> d, shared_ptr<ZDvec> e, shared_ptr<const RelMOFile> jop) const {
  const int ij = d->ij();
  const int lenab = d->lena()*d->lenb();
  ZMatrix tmp(*jop->mo2e(bitset<4>("0101")));
  SMITH::sort_indices<1,0,2,3,1,1,-1,1>(jop->mo2e(bitset<4>("1001"))->data(), tmp.data(), norb_, norb_, norb_, norb_);

  zgemm3m_("n", "t", lenab, ij, ij, 1.0, d->data(), lenab, tmp.data(), ij, 0.0, e->data(), lenab);
}
