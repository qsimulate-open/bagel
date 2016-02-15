//
// BAGEL - Parallel electron correlation program.
// Filename: caspt2grad_util.cc
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

#include <bagel_config.h>
#include <src/smith/caspt2grad.h>

using namespace std;
using namespace bagel;

shared_ptr<DFFullDist> CASPT2Grad::contract_D1(shared_ptr<const DFFullDist> full) const {
#ifdef COMPILE_SMITH
  const int nclosed = ref_->nclosed();
  const int nocc = ref_->nocc();
  const int nall = nocc + ref_->nvirt();

  // TODO D1 must be parallelised as it is very big.
  // construct D1 to be used in Y4 and Y5

  auto dm2 = d2_->matrix2();
  auto D1 = make_shared<btas::Tensor4<double>>(nocc,nall,nocc,nall);
  fill(D1->begin(), D1->end(), 0.0);
  {
    auto is_cc = [&](const int& i) { return i < nclosed; };
    auto is_closed = [&](const int& i) { return i < nclosed && i >= ncore_; };
    auto is_act = [&](const int& i) { return i >= nclosed && i < nocc; };
    auto is_virt = [&](const int& i) { return i >= nocc; };

    // resizing dm2_(le,kf) to dm2_(lt,ks). no resort necessary.
    for (int s = 0; s != nall; ++s) // extend
      for (int k = 0; k != nocc; ++k)
        for (int t = 0; t != nall; ++t) // extend
          for (int l = 0; l != nocc; ++l) {
            // TODO ugly code - there should be a smart way of doing this! just need a logic to test if it has to be symmetrized
            if (k >=  ncore_ && l >= ncore_) {
              // ccaa, cxaa, xxaa
              if (is_virt(t) && is_virt(s)) {
                (*D1)(l, t, k, s) = dm2->element(l-ncore_+(nocc-ncore_)*(t-nclosed), k-ncore_+(nocc-ncore_)*(s-nclosed));
                // cxaa
                if (is_act(k) ^ is_act(l))
                  (*D1)(l, t, k, s) += dm2->element(k-ncore_+(nocc-ncore_)*(s-nclosed), l-ncore_+(nocc-ncore_)*(t-nclosed));
              // ccxa, ccxx
              } else if (t >= nclosed && s >= nclosed && is_closed(k) && is_closed(l)) {
                (*D1)(l, t, k, s) = dm2->element(l-ncore_+(nocc-ncore_)*(t-nclosed), k-ncore_+(nocc-ncore_)*(s-nclosed));
                // ccxa
                if ((t < nocc) ^ (s < nocc))
                  (*D1)(l, t, k, s) += dm2->element(k-ncore_+(nocc-ncore_)*(s-nclosed), l-ncore_+(nocc-ncore_)*(t-nclosed));
              // cxxa, xcxa
              } else if ((is_act(k) ^ is_act(l)) && ((t < nocc) ^ (s < nocc)) && (t >= nclosed && s >= nclosed)) {
                (*D1)(l, t, k, s)  = dm2->element(l-ncore_+(nocc-ncore_)*(t-nclosed), k-ncore_+(nocc-ncore_)*(s-nclosed));
                (*D1)(l, t, k, s) += dm2->element(k-ncore_+(nocc-ncore_)*(s-nclosed), l-ncore_+(nocc-ncore_)*(t-nclosed));
              // cxxx
              } else if ((is_act(k) ^ is_act(l)) && is_act(t) && is_act(s)) {
                (*D1)(l, t, k, s)  = dm2->element(l-ncore_+(nocc-ncore_)*(t-nclosed), k-ncore_+(nocc-ncore_)*(s-nclosed));
                (*D1)(l, t, k, s) += dm2->element(k-ncore_+(nocc-ncore_)*(s-nclosed), l-ncore_+(nocc-ncore_)*(t-nclosed));
              // xxxa
              } else if (((k >= nclosed) && (l >= nclosed)) && ((t < nocc) ^ (s < nocc)) && (t >= nclosed && s >= nclosed)) {
                (*D1)(l, t, k, s)  = dm2->element(l-ncore_+(nocc-ncore_)*(t-nclosed), k-ncore_+(nocc-ncore_)*(s-nclosed));
                (*D1)(l, t, k, s) += dm2->element(k-ncore_+(nocc-ncore_)*(s-nclosed), l-ncore_+(nocc-ncore_)*(t-nclosed));
              }
            }

            // c(cc)x, c(cc)a, x(cc)a // TODO maybe better to work with inactive Fock to deal with this contribution
            if ((is_cc(k) && is_cc(l) && (is_act(s) ^ is_act(t))) // c(cc)x
            || ((is_act(k) ^ is_act(l)) && ((is_cc(s) && is_virt(t)) || (is_virt(s) && is_cc(t)))) // x(cc)a
            ||  (is_cc(k) && is_cc(l) && ((is_cc(s) && is_virt(t)) || (is_virt(s) && is_cc(t))))) { // c(cc)a
              if (s == k)
                (*D1)(l, t, k, s)  += 2.0*d11_->element(l, t);
              if (s == l)
                (*D1)(l, t, k, s)  -= d11_->element(k, t);
              if (t == l)
                (*D1)(l, t, k, s)  += 2.0*d11_->element(k, s);
              if (t == k)
                (*D1)(l, t, k, s)  -= d11_->element(l, s);
            }
          }
  }
  return full->apply_2rdm(*D1);
#else
  return full->clone(); // dummy
#endif
}
