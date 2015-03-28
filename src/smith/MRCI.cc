//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/util/math/davidson.h>
#include <src/smith/extrap.h>
#include <src/smith/MRCI.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

MRCI::MRCI::MRCI(shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {
  eig_ = f1_->diag();
  nstates_ = ref->ciwfn()->nstates();

  for (int i = 0; i != nstates_; ++i) {
    auto tmp = make_shared<MultiTensor>(nstates_);
    for (auto& j : *tmp)
      j = init_amplitude();
    t2all_.push_back(tmp);
    sall_.push_back(tmp->clone());
    nall_.push_back(tmp->clone());
  }
}


void MRCI::MRCI::solve() {
  Timer timer;
  print_iteration();

  const double core_nuc = core_energy_ + ref_->geom()->nuclear_repulsion();

  // target state
  for (int istate = 0; istate != nstates_; ++istate) {
    const double refen = ref_->ciwfn()->energy(istate) - core_nuc;
    // takes care of ref coefficients
    t2all_[istate]->fac(istate) = 1.0;
    nall_[istate]->fac(istate)  = 1.0;
    sall_[istate]->fac(istate)  = refen;

    for (int jst = 0; jst != nstates_; ++jst) {
      set_rdm(jst, istate);
      s = sall_[istate]->at(jst);
      auto queue = make_sourceq(false);
      while (!queue->done())
        queue->next_compute();
    }
    for (int ist = 0; ist != nstates_; ++ist) {
      t2 = t2all_[istate]->at(ist);
      for (int jst = 0; jst != nstates_; ++jst) {
        set_rdm(jst, ist);
        n = nall_[istate]->at(jst);
        auto queue = make_normq(false);
        while (!queue->done())
          queue->next_compute();
      }
    }
  }

  DavidsonDiag_<Amplitude, Residual> davidson(nstates_, 10);

  // first iteration is trivial
  {
    vector<shared_ptr<const Amplitude>> a0;
    vector<shared_ptr<const Residual>> r0;
    for (int istate = 0; istate != nstates_; ++istate) {
      a0.push_back(make_shared<Amplitude>(t2all_[istate]->copy(), nall_[istate]->copy(), this));
      r0.push_back(make_shared<Residual>(sall_[istate]->copy(), this));
    }
    davidson.compute(a0, r0);
  }

  // set the result to t2
  {
    vector<shared_ptr<Residual>> res = davidson.residual();
    for (int i = 0; i != nstates_; ++i) {
      t2all_[i]->zero();
      update_amplitude(t2all_[i], res[i]->tensor());
    }
  }

  shared_ptr<MultiTensor> rtmp = t2all_[0]->clone();

  Timer mtimer;
  int iter = 0;
  for ( ; iter != ref_->maxiter(); ++iter) {

    // loop over state of interest
    vector<shared_ptr<const Amplitude>> a0;
    vector<shared_ptr<const Residual>> r0;
    for (int istate = 0; istate != nstates_; ++istate) {
      // first calculate left-hand-side vectors of t2 (named n)
      nall_[istate]->zero();
      for (int ist = 0; ist != nstates_; ++ist) {
        for (int jst = 0; jst != nstates_; ++jst) {
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          n  = nall_[istate]->at(jst);;
          auto queue = make_normq(false);
          while (!queue->done())
            queue->next_compute();
        }
      }

      // normalize t2 and n
      const double scal = 1.0 / sqrt(dot_product_transpose(nall_[istate], t2all_[istate]));
      nall_[istate]->scale(scal);
      t2all_[istate]->scale(scal);

      a0.push_back(make_shared<Amplitude>(t2all_[istate]->copy(), nall_[istate]->copy(), this));

      // compute residuals (named r)
      rtmp->zero();
      for (int ist = 0; ist != nstates_; ++ist) { // ket sector
        for (int jst = 0; jst != nstates_; ++jst) { // bra sector
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          r = rtmp->at(jst);
          auto queue = make_residualq(false);
          while (!queue->done())
            queue->next_compute();
        }
      }

      // <ab/ij| T |0_ist> Eref_ist.
      {
        shared_ptr<MultiTensor> m = t2all_[istate]->copy();
        // First weighted T2 amplitude
        for (int ist = 0; ist != nstates_; ++ist)
          m->at(ist)->scale(ref_->ciwfn()->energy(ist) - core_nuc);
        // then add it to residual
        for (int ist = 0; ist != nstates_; ++ist) {
          for (int jst = 0; jst != nstates_; ++jst) {
            set_rdm(jst, ist);
            t2 = m->at(ist);
            n  = rtmp->at(jst);;
            auto queue = make_normq(false);
            while (!queue->done())
              queue->next_compute();
          }
        }
      }

      {
        shared_ptr<MultiTensor> m = rtmp->copy();
        int j = 0;
        for (auto& i : *sall_[istate]) {
          double tmp = 0.0;
          for (auto& j : *t2all_[istate])
            tmp += dot_product_transpose(i, j);
          m->fac(j++) += tmp;
        }
        r0.push_back(make_shared<Residual>(m, this));
      }
    }

    energy_ = davidson.compute(a0, r0);

    // find new trial vectors
    vector<shared_ptr<Residual>> res = davidson.residual();
    vector<bool> conv(nstates_, false);
    for (int i = 0; i != nstates_; ++i) {
      const double err = res[i]->tensor()->rms();
      print_iteration(iter, energy_[i]+core_nuc, err, mtimer.tick(), i);

      t2all_[i]->zero();
      update_amplitude(t2all_[i], res[i]->tensor());

      if (err < ref_->thresh()) conv[i] = true;
    }
    if (nstates_ > 1) cout << endl;

    if (all_of(conv.begin(), conv.end(), [](bool i){ return i;})) break;
  }
  print_iteration(iter == ref_->maxiter());
  timer.tick_print("MRCI energy evaluation");
}


void MRCI::MRCI::solve_deriv() {
  throw std::logic_error("Nuclear gradients not implemented for MRCI");
}


#endif
