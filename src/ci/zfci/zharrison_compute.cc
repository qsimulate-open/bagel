//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zharrison_compute.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/ci/zfci/zharrison.h>
#include <src/ci/fci/hztasks.h>
#include <src/util/math/davidson.h>
#include <src/util/taskqueue.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;

/* Implementing the method as described by Harrison and Zarrabian */
#ifndef HAVE_MPI_H
shared_ptr<RelZDvec> ZHarrison::form_sigma(shared_ptr<const RelZDvec> ccvec, shared_ptr<const ZMOFile> jop, const vector<int>& conv) const {
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

#else

shared_ptr<RelZDvec> ZHarrison::form_sigma(shared_ptr<const RelZDvec> ccvec, shared_ptr<const ZMOFile> jop, const vector<int>& conv) const {
  auto sigmavec = make_shared<RelZDvec>(space_, nstate_);
  auto sigmavec_trans = sigmavec->clone(); // important note: the space stays the same after transposition

  int icnt = 0;
  for (auto& isp : space_->detmap()) {
    shared_ptr<const Determinants> base_det = isp.second;
    const int nelea = base_det->nelea();
    const int neleb = base_det->neleb();

    for (int istate = 0; istate != nstate_; ++istate) {
      if (conv[istate]) continue;

      shared_ptr<const ZCivec> cc = ccvec->find(nelea, neleb)->data(istate);
      icnt = sigma_one_parallel(icnt, cc, sigmavec, jop, istate, /*diag*/true, /*transpose*/false);

      shared_ptr<const ZCivec> cc_trans = cc->transpose();
      icnt = sigma_one_parallel(icnt, cc_trans, sigmavec_trans, jop, istate, /*diag*/false, /*transpose*/true);
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
      mpi__->allreduce(sigma->data(ist)->data(), sigma->data(ist)->size());
    }
  }

  return sigmavec;
}
#endif


void ZHarrison::sigma_one(shared_ptr<const ZCivec> cc, shared_ptr<RelZDvec> sigmavec, shared_ptr<const ZMOFile> jop,
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

      // (b^+b^+ b a) contribution
      sigma_2e_h<1,1,0,1>(d, e, jop, trans);
      pdebug.tick_print("task2ab-2 (+1)");
      sigma_2e_create_bb(sigma_1, e);
      pdebug.tick_print("task2ab-3 (+1)");
    }
  }

  if (output1) {
    // output area
    shared_ptr<ZCivec> sigma_1 = sigmavec->find(nelea-1, neleb+1)->data(istate);

    // (b^+ a) contribution
    sigma_1e_ab(cc, sigma_1, jop, trans);
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


void ZHarrison::sigma_aa(shared_ptr<const ZCivec> cc, shared_ptr<ZCivec> sigma, shared_ptr<const ZMOFile> jop, const bool trans) const {
  assert(cc->det()->nelea() == sigma->det()->nelea());
  assert(cc->det()->neleb() == sigma->det()->neleb());

  shared_ptr<const Determinants> det = cc->det();

  bitset<2> bit2;
  bitset<4> bit4;
  if (trans) { bit2 = ~bit2; bit4 = ~bit4; }

  shared_ptr<const ZMatrix> h1 = jop->mo1e(bit2);
  auto h2 = make_shared<ZMatrix>(*jop->mo2e(bit4));
  sort_indices<1,0,2,3,1,1,-1,1>(jop->mo2e(bit4)->data(), h2->data(), norb_, norb_, norb_, norb_);

  TaskQueue<HZTaskAA<complex<double>>> tasks(det->lena());

  const int lb = cc->lenb();
  complex<double>* target = sigma->data();
  for (auto& a : det->string_bits_a()) {
    tasks.emplace_back(cc, a, target, h1->data(), h2->data());
    target += lb;
  }

  tasks.compute();
}


void ZHarrison::sigma_1e_ab(shared_ptr<const ZCivec> cc, shared_ptr<ZCivec> sigma, shared_ptr<const ZMOFile> jop, const bool trans) const {
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
    size_t aindex = 0;
    for (auto& a : sigmadet->string_bits_a()) {
      if (a[i]) {
        ++aindex;
        continue;
      }
      auto ca = a; ca.set(i);
      const complex<double>* source =    cc->data() + lbs * ccdet->lexical<0>(ca);
            complex<double>* target = sigma->data() + lbt * aindex;
      const double asign = ccdet->sign<0>(ca, i);

      for (int j = 0; j != norb_; ++j) {
        size_t bindex = 0;
        for (auto& b : ccdet->string_bits_b()) {
          if (b[j]) {
            ++bindex;
            continue;
          }
          auto cb = b; cb.set(j);
          const complex<double> fac = h1->element(j,i) * (sigmadet->sign<1>(b, j) * asign);
          target[sigmadet->lexical<1>(cb)] += fac * source[bindex];
          ++bindex;
        }
      }
      ++aindex;
    }
  }
}


//////////////// functions for two electron annihilation ///////////////

void ZHarrison::sigma_2e_annih_aa(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  // compute |d> = a_i a_j|cc>
  const complex<double>* const source = cc->data();
  const size_t lb = cc->lenb();
  assert(lb == d->lenb());

  TaskQueue<function<void(void)>> tq(d->lena());
  size_t aindex = 0;
  for (auto& a : d->det()->string_bits_a()) {
    tq.emplace_back(
      [this, &d, &source, a, aindex, &cc, &lb] () {
        for (int i = 0; i != norb_; ++i) {
          for (int j = 0; j != norb_; ++j) {
            if (i == j) continue;
            if (a[i] || a[j]) continue;
            complex<double>* target = d->data(j+norb_*i)->data();
            auto ca = a; ca.set(i); ca.set(j);
            const double factor = Determinants::sign(a, i, j) * (i < j ? -1.0 : 1.0);
            const size_t offas = cc->det()->lexical<0>(ca);
            blas::ax_plus_y_n(factor, source+lb*offas, lb, target+lb*aindex);
          }
        }
      }
    );
    ++aindex;
  }
  tq.compute();
}


void ZHarrison::sigma_2e_annih_ab(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  // compute |d> = b_l a_k|cc>
  shared_ptr<const Determinants> bdet = cc->det(); // base
  shared_ptr<const Determinants> tdet = d->det();  // target

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

  TaskQueue<function<void(void)>> tq(lbt);
  size_t bindex = 0;
  for (auto& b : base_det->string_bits_b()) {
    tq.emplace_back(
      [this, &e, &int_det, &la, &lbs, &lbt, bindex, &target_base, b] () {
        for (int i = 0; i < norb_; ++i) { // beta
          for (int j = 0; j < i; ++j) { // beta
            if (!b[i] || !b[j]) { // looking for particles to annihilate
              continue;
            }
            const complex<double>* source_ij = e->data(i+norb_*j)->data();
            const complex<double>* source_ji = e->data(j+norb_*i)->data();
            bitset<nbit__> cb = b;
            cb.reset(i); cb.reset(j);
            const double sign = Determinants::sign(b, i, j);
            const size_t bdlex = int_det->lexical<1>(cb);

            zaxpy_(la, -sign, source_ij+bdlex, lbs, target_base+bindex, lbt);
            zaxpy_(la,  sign, source_ji+bdlex, lbs, target_base+bindex, lbt);
          }
        }
      }
    );
    ++bindex;
  }
  tq.compute();
}


//////////////// functions for multiplication of the Hamiltonian ///////////////

void ZHarrison::sigma_2e_h0101_h1001(shared_ptr<const ZDvec> d, shared_ptr<ZDvec> e, shared_ptr<const ZMOFile> jop) const {
  ZMatrix tmp(*jop->mo2e("0101"));
  sort_indices<1,0,2,3,1,1,-1,1>(jop->mo2e("1001")->data(), tmp.data(), norb_, norb_, norb_, norb_);
  contract(1.0, *d, {0,1,2}, tmp, {3,2}, 0.0, *e, {0,1,3});
}
