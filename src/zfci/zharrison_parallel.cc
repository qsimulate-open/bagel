//
// BAGEL - Parallel electron correlation program.
// Filename: zharrison_parallel.cc
// Copyright (C) 2013 Jefferson E. Bates
//
// Author: Jefferson E. Bates <jefferson.bates@northwestern.edu>
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

using namespace std;
using namespace bagel;

#ifdef HAVE_MPI_H

int ZHarrison::sigma_one_parallel(const int icnt,   shared_ptr<const ZCivec> cc, shared_ptr<RelZDvec> sigmavec, shared_ptr<const RelMOFile> jop,
                                  const int istate, const bool diag, const bool trans) const {

  const int ij = norb_*norb_;
  shared_ptr<const Determinants> base_det = cc->det();
  // diagonal output
  const int nelea = base_det->nelea();
  const int neleb = base_det->neleb();

  const bool noab = (base_det->nelea() == 0 || base_det->neleb() == 0);
  const bool noaa =  base_det->nelea() <= 1 || base_det->neleb()+1 > norb_;
  const bool output1 = base_det->nelea()-1 >= 0 && base_det->neleb()+1 <= norb_;

  int cnt = icnt;
  auto comp = [&](const int c) { return c % mpi__->size() == mpi__->rank(); };

  if (comp(cnt++)) {
    shared_ptr<ZCivec> sigma = sigmavec->find(nelea, neleb)->data(istate);
    sigma_aa(cc, sigma, jop, trans);
  }

  if (!noab && diag) {
    if (comp(cnt++)) {
      shared_ptr<ZCivec> sigma = sigmavec->find(nelea, neleb)->data(istate);
      shared_ptr<const Determinants> int_det = int_space_->finddet(nelea-1, neleb-1);
      auto d = make_shared<ZDvec>(int_det, ij);
      auto e = make_shared<ZDvec>(int_det, ij);

      // (2ab) alpha-beta contributions
      /* Resembles more the Knowles & Handy FCI terms */

      sigma_2e_annih_ab(cc, d);

      // (a^+ b^+ b a) and (a^+ b^+ a b) contributions
      sigma_2e_h0101_h1001(d, e, jop);

      sigma_2e_create_ab(sigma, e);
    }
  }

  if (output1) {
    if (comp(cnt++)) {
      if (!noab) {
        shared_ptr<const Determinants> int_det = int_space_->finddet(nelea-1, neleb-1);
        auto d = make_shared<ZDvec>(int_det, ij);
        auto e = make_shared<ZDvec>(int_det, ij);

        sigma_2e_annih_ab(cc, d);

        // output area
        shared_ptr<ZCivec> sigma_1 = sigmavec->find(nelea-1, neleb+1)->data(istate);

        // (b^+b^+ b a) contribution
        sigma_2e_h<1,1,0,1>(d, e, jop, trans);
        sigma_2e_create_bb(sigma_1, e);
      }
      // output area
      shared_ptr<ZCivec> sigma_1 = sigmavec->find(nelea-1, neleb+1)->data(istate);

      // (b^+ a) contribution
      sigma_1e_ab(cc, sigma_1, jop, trans);
    }
  }

#if 0
  if (!noaa) {
    if (comp(cnt++)) {
      shared_ptr<const Determinants> int_det = int_space_->finddet(nelea-2, neleb);
      auto d = make_shared<ZDvec>(int_det, ij);
      auto e = make_shared<ZDvec>(int_det, ij);

      sigma_2e_annih_aa(cc, d);

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
#else
  if (!noaa) {
    if (comp(cnt++)) {
      shared_ptr<const Determinants> int_det = int_space_->finddet(nelea-2, neleb);
      auto d = make_shared<ZDvec>(int_det, ij);
      auto e = make_shared<ZDvec>(int_det, ij);

      sigma_2e_annih_aa(cc, d);

      // (a^+ b^+ a a) contribution
      sigma_2e_h<0,1,0,0>(d, e, jop, trans);

      assert(neleb+1 <= norb_);
      // +1 sector
      shared_ptr<ZCivec> sigma_1 = sigmavec->find(nelea-1, neleb+1)->data(istate);
      sigma_2e_create_ab(sigma_1, e);
    }
  }

  if (!noaa && base_det->neleb()+2 <= norb_) {
    if (comp(cnt++)) {
      shared_ptr<const Determinants> int_det = int_space_->finddet(nelea-2, neleb);
      auto d = make_shared<ZDvec>(int_det, ij);
      auto e = make_shared<ZDvec>(int_det, ij);

      sigma_2e_annih_aa(cc, d);

      // +2 sector
      shared_ptr<ZCivec> sigma_2 = sigmavec->find(nelea-2, neleb+2)->data(istate);
      sigma_2e_h<1,1,0,0>(d, e, jop, trans, 0.5);
      sigma_2e_create_bb(sigma_2, e);
    }
  }
#endif

  return cnt;
}

#endif
