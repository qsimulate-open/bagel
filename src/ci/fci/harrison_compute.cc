//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zharrison_compute.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <src/ci/fci/harrison.h>
#include <src/ci/fci/hztasks.h>
#include <src/util/taskqueue.h>
#include <src/util/prim_op.h>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::HarrisonZarrabian)

using namespace std;
using namespace bagel;

// toggle for timing print out.
static const bool tprint = false;


/* Implementing the method as described by Harrison and Zarrabian */
shared_ptr<Dvec> HarrisonZarrabian::form_sigma(shared_ptr<const Dvec> ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const { // d and e are scratch area for D and E intermediates
  const int ij = norb_*norb_;
  assert(ccvec->ij() == nstate_);

  auto sigmavec = make_shared<Dvec>(ccvec->det(), nstate_);
  sigmavec->zero();

  shared_ptr<Determinants> base_det = space_->finddet(nelea_, neleb_);
  shared_ptr<Determinants> int_det = space_->finddet(nelea_-1,neleb_-1);

  /* d and e are only used in the alpha-beta case and exist in the (nalpha-1)(nbeta-1) spaces */
  auto d = make_shared<Dvec>(int_det, ij);
  auto e = make_shared<Dvec>(int_det, ij);

  for (int istate = 0; istate != nstate_; ++istate) {
    Timer pdebug(3);
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

  return sigmavec;
}


void HarrisonZarrabian::sigma_aa(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  assert(cc->det() == sigma->det());

  shared_ptr<const Determinants> det = cc->det();
  const int lb = cc->lenb();

  auto h1 = make_shared<Matrix>(norb_, norb_);
  for (int i = 0, ij = 0; i < norb_; ++i) {
    for (int j = 0; j <= i; ++j, ++ij) {
      h1->element(i, j) = h1->element(j, i) = jop->mo1e(ij);
    }
  }

  auto h2 = make_shared<Matrix>(*jop->mo2e());
  sort_indices<1,0,2,3,1,1,-1,1>(jop->mo2e()->data(), h2->data(), norb_, norb_, norb_, norb_);

  TaskQueue<HZTaskAA<double>> tasks(det->lena());

  double* target = sigma->data();
  for (auto& a : det->string_bits_a()) {
    tasks.emplace_back(cc, a, target, h1->data(), h2->data());
    target += lb;
  }

  tasks.compute();
}


void HarrisonZarrabian::sigma_bb(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  shared_ptr<const Civec> cc_trans = cc->transpose();
  auto sig_trans = make_shared<Civec>(cc_trans->det());

  sigma_aa(cc_trans, sig_trans, jop);

  sigma->ax_plus_y(1.0, *sig_trans->transpose(sigma->det()));
}


void HarrisonZarrabian::sigma_2ab_1(shared_ptr<const Civec> cc, shared_ptr<Dvec> d) const {
  const int norb = norb_;

  shared_ptr<const Determinants> bdet = cc->det(); // base
  shared_ptr<const Determinants> tdet = d->det();  // target

  const int lbs = bdet->lenb();
  const double* source_base = cc->data();

  TaskQueue<HZTaskAB1<double>> tasks(norb*norb);

  for (int k = 0; k < norb; ++k) {
    for (int l = 0; l < norb; ++l) {
      double* target_base = d->data(k*norb + l)->data();
      tasks.emplace_back(tdet, lbs, source_base, target_base, k, l);
    }
  }

  tasks.compute();
}


void HarrisonZarrabian::sigma_2ab_2(shared_ptr<Dvec> d, shared_ptr<Dvec> e, shared_ptr<const MOFile> jop) const {
  const int ij = d->ij();
  const int lenab = d->lena() * d->lenb();
  dgemm_("n", "n", lenab, ij, ij, 1.0, d->data(), lenab, jop->mo2e_ptr(), ij, 0.0, e->data(), lenab);
}


void HarrisonZarrabian::sigma_2ab_3(shared_ptr<Civec> sigma, shared_ptr<Dvec> e) const {
  const shared_ptr<const Determinants> base_det = sigma->det();
  const shared_ptr<const Determinants> int_det = e->det();

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
}
