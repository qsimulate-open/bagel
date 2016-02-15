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
  const int n = nocc-ncore_;

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
    auto is_exc = [&](const int& i, const int& j) { return is_virt(i) || (is_act(i) && is_cc(j)); };

    // resizing dm2_(le,kf) to dm2_(lt,ks). no resort necessary.
    for (int s = 0; s != nall-nclosed; ++s) // extend
      for (int k = 0; k != nocc-ncore_; ++k)
        for (int t = 0; t != nall-nclosed; ++t) // extend
          for (int l = 0; l != nocc-ncore_; ++l) {
            const int ss = s + nclosed;
            const int tt = t + nclosed;
            const int kk = k + ncore_;
            const int ll = l + ncore_;
            // ccaa, cxaa, xxaa
            if (is_virt(tt) && is_virt(ss)) {
              (*D1)(ll,tt,kk,ss) = dm2->element(l+n*t, k+n*s);
              // cxaa
              if (is_act(kk) ^ is_act(ll))
                (*D1)(ll,tt,kk,ss) += dm2->element(k+n*s, l+n*t);
            // ccxa, ccxx
            } else if (is_closed(kk) && is_closed(ll)) {
              (*D1)(ll,tt,kk,ss) = dm2->element(l+n*t, k+n*s);
              // ccxa
              if (is_act(tt) ^ is_act(ss))
                (*D1)(ll,tt,kk,ss) += dm2->element(k+n*s, l+n*t);
            // cxxa, xcxa
            } else if ((is_act(kk) ^ is_act(ll)) && (is_act(tt) ^ is_act(ss))) {
              (*D1)(ll,tt,kk,ss)  = dm2->element(l+n*t, k+n*s);
              (*D1)(ll,tt,kk,ss) += dm2->element(k+n*s, l+n*t);
            // cxxx
            } else if ((is_act(kk) ^ is_act(ll)) && is_act(tt) && is_act(ss)) {
              (*D1)(ll,tt,kk,ss)  = dm2->element(l+n*t, k+n*s);
              (*D1)(ll,tt,kk,ss) += dm2->element(k+n*s, l+n*t);
            // xxxa
            } else if (is_act(kk) && is_act(ll) && (is_act(tt) ^ is_act(ss))) {
              (*D1)(ll,tt,kk,ss)  = dm2->element(l+n*t, k+n*s);
              (*D1)(ll,tt,kk,ss) += dm2->element(k+n*s, l+n*t);
            }
          }

    for (int s = 0; s != nclosed; ++s)
      for (int t = 0; t != nall; ++t)
        for (int l = 0; l != nocc; ++l)
          // c(cc)x, c(cc)a, x(cc)a
          if (is_exc(t, l)) {
            (*D1)(l, t, s, s)  += 2.0*d11_->element(l, t); // (g|lt)dlt -> (g|ss)
            (*D1)(s, t, l, s)  -=     d11_->element(l, t); // (g|st)dlt -> (g|ls)
            (*D1)(s, s, l, t)  += 2.0*d11_->element(l, t); // (g|ss)dlt -> (g|lt)
            (*D1)(l, s, s, t)  -=     d11_->element(l, t); // (g|ls)dlt -> (g|st)
          }

  }
  return full->apply_2rdm(*D1);
#else
  return full->clone(); // dummy
#endif
}
