//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI.cc
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


#include <src/util/math/davidson.h>
#include <src/smith/extrap.h>
#include <src/smith/mrci/MRCI.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

MRCI::MRCI::MRCI(shared_ptr<const SMITH_Info<double>> ref) : SpinFreeMethod(ref) {
  eig_ = f1_->diag();
  nstates_ = ref->ciwfn()->nstates();

  for (int i = 0; i != nstates_; ++i) {
    auto tmp = make_shared<MultiTensor>(nstates_);
    for (auto& j : *tmp)
      j = init_amplitude();
    t2all_.push_back(tmp);

    auto tmp2 = make_shared<MultiTensor>(nstates_);
    for (auto& j : *tmp2)
      j = init_residual();
    sall_.push_back(tmp2);
    nall_.push_back(tmp2->copy());
  }
}


void MRCI::MRCI::solve() {
  Timer timer;
  print_iteration();

  const double core_nuc = core_energy_ + info_->geom()->nuclear_repulsion();

  // target state
  for (int istate = 0; istate != nstates_; ++istate) {
    const double refen = info_->ciwfn()->energy(istate) - core_nuc;
    // takes care of ref coefficients
    t2all_[istate]->fac(istate) = 1.0;
    nall_[istate]->fac(istate)  = 1.0;
    sall_[istate]->fac(istate)  = refen;

    for (int jst = 0; jst != nstates_; ++jst) {
      set_rdm(jst, istate);
      s = sall_[istate]->at(jst);
      auto queue = make_sourceq(false, jst == istate);
      while (!queue->done())
        queue->next_compute();
    }
  }

  DavidsonDiag_<Amplitude<double>, Residual<double>, Matrix> davidson(nstates_, info_->davidson_subspace());

  // first iteration is trivial
  {
    vector<shared_ptr<const Amplitude<double>>> a0;
    vector<shared_ptr<const Residual<double>>> r0;
    for (int istate = 0; istate != nstates_; ++istate) {
      a0.push_back(make_shared<Amplitude<double>>(t2all_[istate]->copy(), nall_[istate]->copy(), this));
      r0.push_back(make_shared<Residual<double>>(sall_[istate]->copy(), this));
    }
    energy_ = davidson.compute(a0, r0);
    for (int istate = 0; istate != nstates_; ++istate)
      assert(fabs(energy_[istate]+core_nuc - info_->ciwfn()->energy(istate)) < 1.0e-8);
  }

  // set the result to t2
  {
    vector<shared_ptr<Residual<double>>> res = davidson.residual();
    for (int i = 0; i != nstates_; ++i) {
      t2all_[i]->zero();
      e0_ = e0all_[i];
      update_amplitude(t2all_[i], res[i]->tensor());
    }
  }

  shared_ptr<MultiTensor> rtmp = nall_[0]->copy();
  rtmp->zero();

  Timer mtimer;
  vector<bool> conv(nstates_, false);
  for (int iter = 0; iter != info_->maxiter(); ++iter) {

    // loop over state of interest
    vector<shared_ptr<const Amplitude<double>>> a0;
    vector<shared_ptr<const Residual<double>>> r0;
    for (int istate = 0; istate != nstates_; ++istate) {
      if (conv[istate]) {
        a0.push_back(nullptr);
        r0.push_back(nullptr);
        continue;
      }
      // first calculate left-hand-side vectors of t2 (named n)
      nall_[istate]->zero();
      for (int ist = 0; ist != nstates_; ++ist) {
        for (int jst = 0; jst != nstates_; ++jst) {
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          n  = nall_[istate]->at(jst);
          auto queue = make_normq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();
        }
      }

      // normalize t2 and n
      const double scal = 1.0 / sqrt(detail::real(dot_product_transpose(nall_[istate], t2all_[istate])));
      nall_[istate]->scale(scal);
      t2all_[istate]->scale(scal);

      a0.push_back(make_shared<Amplitude<double>>(t2all_[istate]->copy(), nall_[istate]->copy(), this));

      // compute residuals (named r)
      rtmp->zero();
      for (int ist = 0; ist != nstates_; ++ist) { // ket sector
        for (int jst = 0; jst != nstates_; ++jst) { // bra sector
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          r = rtmp->at(jst);
          auto queue = make_residualq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();
          diagonal(r, t2);
        }
      }

      // <ab/ij| T |0_ist> Eref_ist.
      {
        shared_ptr<MultiTensor> m = t2all_[istate]->copy();
        for (int ist = 0; ist != nstates_; ++ist) {
          // First weighted T2 amplitude
          m->at(ist)->scale(info_->ciwfn()->energy(ist) - core_nuc);
          // then add it to residual
          for (int jst = 0; jst != nstates_; ++jst) {
            set_rdm(jst, ist);
            t2 = m->at(ist);
            n  = rtmp->at(jst);
            auto queue = make_normq(false, jst == ist);
            while (!queue->done())
              queue->next_compute();
          }
        }
      }

      {
        shared_ptr<MultiTensor> m = rtmp->copy();
        for (int ist = 0; ist != nstates_; ++ist)
          m->fac(ist) = dot_product_transpose(sall_[ist], t2all_[istate]);
        r0.push_back(make_shared<Residual<double>>(m, this));
      }
    }

    energy_ = davidson.compute(a0, r0);

    // find new trial vectors
    vector<shared_ptr<Residual<double>>> res = davidson.residual();
    for (int i = 0; i != nstates_; ++i) {
      const double err = res[i]->tensor()->norm() / pow(res[i]->tensor()->size(), 0.25);
      print_iteration(iter, energy_[i]+core_nuc, err, mtimer.tick(), i);

      t2all_[i]->zero();
      conv[i] = err < info_->thresh();
      if (!conv[i]) {
        e0_ = e0all_[i];
        update_amplitude(t2all_[i], res[i]->tensor());
      }
    }
    if (nstates_ > 1) cout << endl;

    if (all_of(conv.begin(), conv.end(), [](bool i){ return i;})) break;
  }
  print_iteration(!all_of(conv.begin(), conv.end(), [](bool i){ return i;}));
  timer.tick_print("MRCI energy evaluation");

  // Davidson corrections...
  {
    cout << endl;
    vector<double> energy_q(nstates_);
    vector<shared_ptr<Amplitude<double>>> ci = davidson.civec();
    for (int i = 0; i != nstates_; ++i) {
      double c = 0.0;
      double eref = 0.0;
      for (int j = 0; j != nstates_; ++j) {
        const double cn = norm(ci[i]->tensor()->fac(j));
        c += cn;
        eref += cn * info_->ciwfn()->energy(j);
      }
      eref /= c;
      const double eq = energy_[i]+core_nuc + (energy_[i]+core_nuc-eref)*(1.0-c)/c;
      print_iteration(0, eq, 0.0, 0.0, i);
    }
    cout << endl;
  }
  timer.tick_print("MRCI+Q energy evaluation");
}


void MRCI::MRCI::solve_gradient(const int targetJ, const int targetI, shared_ptr<const NacmType> nacmtype, const bool nocider) {
  throw std::logic_error("Nuclear gradients not implemented for MRCI");
}


#endif
