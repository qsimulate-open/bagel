//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/caspt2/MSCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


MSCASPT2::MSCASPT2::MSCASPT2(const CASPT2::CASPT2& cas) {
  info_    = cas.info_;
  virt_    = cas.virt_;
  active_  = cas.active_;
  closed_  = cas.closed_;
  rvirt_   = cas.rvirt_;
  ractive_ = cas.ractive_;
  rclosed_ = cas.rclosed_;
  heff_    = cas.heff_;
  fockact_ = cas.fockact_;
  e0all_   = cas.e0all_;

  t2all_ = cas.t2all_;
  lall_  = cas.lall_;
  h1_ = cas.h1_;
  f1_ = cas.f1_;
  v2_ = cas.v2_;
  den1 = cas.h1_->clone();
  den2 = cas.h1_->clone();
  Den1 = cas.init_residual();

  rdm0all_ = cas.rdm0all_;
  rdm1all_ = cas.rdm1all_;
  rdm2all_ = cas.rdm2all_;
  rdm3all_ = cas.rdm3all_;
  rdm4all_ = cas.rdm4all_;
}


void MSCASPT2::MSCASPT2::solve_deriv() {
  Timer timer;
  const int nstates = info_->ciwfn()->nstates();
  const int target = info_->target();

  // first-order energy from the energy expression
  // TODO these can be computed more efficiently using mixed states
  {
    shared_ptr<Tensor> result = den1->clone();
    shared_ptr<Tensor> result2 = Den1->clone();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      const double jheff = (*heff_)(jst, target);
      for (int ist = 0; ist != nstates; ++ist) { // ket
        set_rdm(jst, ist);
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          const double isheff = (*heff_)(istate, target);
          l2 = t2all_[istate]->at(ist); // careful
          shared_ptr<Queue> queue = make_density1q(true, ist == jst);
          while (!queue->done())
            queue->next_compute();
          result->ax_plus_y(isheff*jheff, den1);

          shared_ptr<Queue> queue2 = make_density2q(true, ist == jst);
          while (!queue2->done())
            queue2->next_compute();
          result2->ax_plus_y(isheff*jheff, Den1);
        }
      }
    }
    den1_ = result->matrix();
    Den1_ = result2;
  }
  // second-order contribution from the lambda terms
  {
    // TODO probably not necessary
    shared_ptr<Tensor> result = den2->clone();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      for (int ist = 0; ist != nstates; ++ist) { // ket
        set_rdm(jst, ist);
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          l2 = lall_[istate]->at(ist);
          t2 = t2all_[istate]->at(jst);
          shared_ptr<Queue> queue = make_densityq(true, ist == jst);
          while (!queue->done())
            queue->next_compute();
          result->ax_plus_y(1.0, den2);
        }
      }
    }
    den2_ = result->matrix();
  }
  // first-order contribution from the lambda terms
  {
    // TODO probably not necessary
    shared_ptr<Tensor> result = den1->clone();
    shared_ptr<Tensor> result2 = Den1->clone();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      for (int ist = 0; ist != nstates; ++ist) { // ket
        set_rdm(jst, ist);

        l2 = lall_[jst]->at(ist);
        shared_ptr<Queue> queue = make_density1q(true, ist == jst);
        while (!queue->done())
          queue->next_compute();
        result->ax_plus_y(1.0, den1);

        shared_ptr<Queue> queue2 = make_density2q(true, ist == jst);
        while (!queue2->done())
          queue2->next_compute();
        result2->ax_plus_y(1.0, Den1);
      }
    }
    den1_->ax_plus_y(1.0, result->matrix());
    Den1_->ax_plus_y(1.0, result2);
  }
  // because of the convention...
  den1_->scale(0.5);
  Den1_->scale(0.5);
  timer.tick_print("Correlated density matrix evaluation");

  // CI derivative..
  ci_deriv_ = make_shared<Dvec>(info_->ref()->ciwfn()->det(), nstates);
  // outer most look over RDM derivatives
  for (int nst = 0; nst != nstates; ++nst) {
    const size_t cisize = ci_deriv_->data(nst)->size();
    const size_t cimax = 1000; // TODO interface to the input
    const size_t npass = (cisize-1)/cimax+1;
    const size_t chunk = (cisize-1)/npass+1;
    const double nheff = (*heff_)(nst, target);

    for (int ipass = 0; ipass != npass; ++ipass) {
      const size_t offset = ipass * chunk;
      const size_t size = min(chunk, cisize-offset);

      tie(ci_, rci_, rdm0deriv_, rdm1deriv_, rdm2deriv_, rdm3deriv_, rdm4deriv_)
        = SpinFreeMethod<double>::feed_rdm_deriv(info_, active_, fockact_, nst, offset, size);

      // output area
      deci = make_shared<Tensor>(vector<IndexRange>{ci_});
      deci->allocate();

      shared_ptr<Queue> dec;

      for (int lst = 0; lst != nstates; ++lst) {
        const double lheff = (*heff_)(lst, target);
        // derivative with respect to M
        for (int mst = 0; mst != nstates; ++mst) {
          const double mheff = (*heff_)(mst, target);

          l2 = t2all_[lst]->at(mst);
          dec = make_deci2q(/*zero*/true);
          while (!dec->done())
            dec->next_compute();
          blas::ax_plus_y_n(lheff*nheff, deci->vectorb()->data(), size, ci_deriv_->data(mst)->data()+offset);

          l2 = t2all_[lst]->at(nst);
          dec = make_deci3q(/*zero*/true);
          while (!dec->done())
            dec->next_compute();
          blas::ax_plus_y_n(lheff*mheff, deci->vectorb()->data(), size, ci_deriv_->data(mst)->data()+offset);
        }
      }

      // derivative with respect to M
      for (int mst = 0; mst != nstates; ++mst) {
        l2 = lall_[nst]->at(mst);
        dec = make_deci2q(/*zero*/true);
        while (!dec->done())
          dec->next_compute();
        l2 = lall_[mst]->at(nst);
        dec = make_deci3q(false);
        while (!dec->done())
          dec->next_compute();

        for (int lst = 0; lst != nstates; ++lst) {
          e0_ = e0all_[lst];
          l2 = lall_[lst]->at(mst);
          t2 = t2all_[lst]->at(nst);
          dec = make_deciq(false);
          while (!dec->done())
            dec->next_compute();

          l2 = t2all_[lst]->at(mst);
          t2 = lall_[lst]->at(nst);
          dec = make_deciq(false);
          while (!dec->done())
            dec->next_compute();
        }
        blas::ax_plus_y_n(1.0, deci->vectorb()->data(), size, ci_deriv_->data(mst)->data()+offset);
      }
      stringstream ss; ss << "CI derivative evaluation " << setw(5) << ipass+1 << " / " << npass
                          << "   (" << setw(2) << nst+1 << " /" << setw(2) << nstates << ")";
      timer.tick_print(ss.str());
    }
  }
}

#endif
