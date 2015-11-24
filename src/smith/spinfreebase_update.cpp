//
// BAGEL - Parallel electron correlation program.
// Filename: spinfreebase_denom.cpp
// Copyright (C) 2015 Toru Shiozaki
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

#ifdef SPINFREEMETHOD_DETAIL

template<typename DataType>
shared_ptr<TATensor<DataType,4>> SpinFreeMethod<DataType>::update_amplitude(shared_ptr<TATensor<DataType,4>> t, shared_ptr<const TATensor<DataType,4>> r) const {
  auto tt = make_shared<MultiTATensor<DataType,4>>(vector<DataType>{0.0}, vector<shared_ptr<TATensor<DataType,4>>>{t});
  // ... not good code...
  auto r0 = const_pointer_cast<TATensor<DataType,4>>(r);
  auto rr = make_shared<MultiTATensor<DataType,4>>(vector<DataType>{0.0}, vector<shared_ptr<TATensor<DataType,4>>>{r0});

  shared_ptr<MultiTATensor<DataType,4>> tnew = update_amplitude(tt, rr);
  return tnew->at(0);
}


template<typename DataType>
shared_ptr<MultiTATensor<DataType,4>> SpinFreeMethod<DataType>::update_amplitude(shared_ptr<MultiTATensor<DataType,4>> t, shared_ptr<const MultiTATensor<DataType,4>> r) const {
  if (t->nref() != r->nref())
    throw logic_error("something is wrong. SpinFreeMethod::update_amplitude");

  const int fac2 = is_same<DataType,double>::value ? 1.0 : 2.0;
  const int ncore = info_->ncore() * fac2;
  const int nocc  = (info_->nclosed() + info_->nact()) * fac2;

  for (int ist = 0; ist != t->nref(); ++ist) {
    auto rist = r->at(ist);
    t->fac(ist) = 0.0;

    // not the best structure, but I am assuming that this does not take too much time...
    for (int jst = 0; jst != t->nref(); ++jst) {
      auto tjst = t->at(jst);

      if (ist == jst) { // AACC
        LazyTATensor<DataType,4,DenomAACC> d({virt_, virt_, closed_, closed_}, DenomAACC(eig_, nocc, nocc, ncore, ncore));
        if (is_same<DataType,double>::value)
          (*tjst)("c0,a1,c2,a3") -= ((*rist)("c0,a1,c2,a3") * (1.0/6.0) + (*rist)("c0,a3,c2,a1") * (1.0/12.0)) * d("a1,a3,c0,c2");
        else
          (*tjst)("c0,a1,c2,a3") -= (*rist)("c0,a1,c2,a3") * 0.25 * d("a1,a3,c0,c2");
      }
      { // AAXX
        LazyTATensor<DataType,3,DenomGen2<1,1>> d({ortho2_, virt_, virt_}, DenomGen2<1,1>(e0_, denom_->denom_xx(), eig_, nocc, nocc));
        shared_ptr<const TATensor<DataType,3>> s = denom_->tashalf_xx({ortho2_, active_, active_});
        (*tjst)("x00,a1,x20,a3") -= (*rist)("x0,a1,x2,a3") * (*s)("o4,x0,x2") * d("o4,a1,a3") * 0.5 * (*s)("o4,x00,x20");
      }
      { // AACX
        LazyTATensor<DataType,4,DenomGen1<1,1,-1>> d({ortho1_, virt_, virt_, closed_}, DenomGen1<1,1,-1>(e0_, denom_->denom_x(), eig_, nocc, nocc, ncore));
        shared_ptr<const TATensor<DataType,2>> s = denom_->tashalf_x({ortho1_, active_});
        if (is_same<DataType,double>::value)
          (*tjst)("x00,a1,c2,a3") -= ((*rist)("c2,a3,x0,a1")*(2.0/3.0) + (*rist)("c2,a1,x0,a3")*(1.0/3.0)) * (*s)("o4,x0") * d("o4,a1,a3,c2") * (*s)("o4,x00");
        else
          (*tjst)("x00,a1,c2,a3") -= (*rist)("c2,a3,x0,a1") * 0.5 * (*s)("o4,x0") * d("o4,a1,a3,c2") * (*s)("o4,x00");
      }
      { // AXCC
        LazyTATensor<DataType,4,DenomGen1<1,-1,-1>> d({ortho1_, virt_, closed_, closed_}, DenomGen1<1,-1,-1>(e0_, denom_->denom_h(), eig_, nocc, ncore, ncore));
        shared_ptr<const TATensor<DataType,2>> s = denom_->tashalf_h({ortho1_, active_});
        if (is_same<DataType,double>::value)
          (*tjst)("c0,a1,c2,x30") -= ((*rist)("c2,x3,c0,a1")*(2.0/3.0) + (*rist)("c0,x3,c2,a1")*(1.0/3.0)) * (*s)("o4,x3") * d("o4,a1,c0,c2") * (*s)("o4,x30");
        else
          (*tjst)("c0,a1,c2,x30") -= (*rist)("c2,x3,c0,a1") * 0.5 * (*s)("o4,x3") * d("o4,a1,c0,c2") * (*s)("o4,x3");
      }
      { // XXCC
        LazyTATensor<DataType,3,DenomGen2<-1,-1>> d({ortho2_, closed_, closed_}, DenomGen2<-1,-1>(e0_, denom_->denom_hh(), eig_, ncore, ncore));
        shared_ptr<const TATensor<DataType,3>> s = denom_->tashalf_hh({ortho2_, active_, active_});
        (*tjst)("c0,x10,c2,x30") -= (*rist)("c0,x1,c2,x3") * (*s)("o4,x1,x3") * d("o4,c0,c2") * 0.5 * (*s)("o4,x10,x30");
      }
      { // XXCX
        LazyTATensor<DataType,2,DenomGen3<-1>> d({ortho3_, closed_}, DenomGen3<-1>(e0_, denom_->denom_xxh(), eig_, ncore));
        shared_ptr<const TATensor<DataType,4>> s = denom_->tashalf_xxh({ortho3_, active_, active_, active_});
        (*tjst)("x00,x10,c2,x30") -= (*rist)("c2,x3,x0,x1") * (*s)("o4,x0,x1,x3") * d("o4,c2") * (*s)("o4,x00,x10,x30");
      }
      { // AXXX
        LazyTATensor<DataType,2,DenomGen3<1>> d({ortho3_, virt_}, DenomGen3<1>(e0_, denom_->denom_xhh(), eig_, nocc));
        shared_ptr<const TATensor<DataType,4>> s = denom_->tashalf_xhh({ortho3_, active_, active_, active_});
        (*tjst)("x00,a1,x20,x30") -= (*rist)("x2,x3,x0,a1") * (*s)("o4,x0,x2,x3") * d("o4,a1") * (*s)("o4,x00,x20,x30");
      }
      {
        LazyTATensor<DataType,3,DenomGen2<1,-1>> d({ortho2t_, virt_, closed_}, DenomGen2<1,-1>(e0_, denom_->denom_xh(), eig_, nocc, ncore));
        if (is_same<DataType,double>::value) {
          shared_ptr<const TATensor<DataType,3>> s0 = denom_->template tashalf_xh_nonrel<0>({ortho2t_, active_, active_});
          shared_ptr<const TATensor<DataType,3>> s1 = denom_->template tashalf_xh_nonrel<1>({ortho2t_, active_, active_});

          // TODO not the most efficient code
          (*tjst)("c0,a1,x20,x30") -= (*rist)("x2,x3,c0,a1") * (*s0)("o4,x2,x3") * d("o4,a1,c0") * (*s0)("o4,x20,x30")
                                    + (*rist)("c0,x3,x2,a1") * (*s1)("o4,x2,x3") * d("o4,a1,c0") * (*s0)("o4,x20,x30");
          (*tjst)("x20,a1,c0,x30") -= (*rist)("x2,x3,c0,a1") * (*s0)("o4,x2,x3") * d("o4,a1,c0") * (*s1)("o4,x20,x30")
                                    + (*rist)("c0,x3,x2,a1") * (*s1)("o4,x2,x3") * d("o4,a1,c0") * (*s1)("o4,x20,x30");

        } else {
          shared_ptr<const TATensor<DataType,3>> s = denom_->tashalf_xh_rel({ortho2t_, active_, active_});
          (*tjst)("c0,a1,x20,x30") -= (*rist)("x2,x3,c0,a1") * (*s)("o4,x2,x3") * d("o4,a1,c0") * (*s)("o4,x20,x30");
        }
      }
    } // jst loop
  } // ist loop

  return t;
}

#endif
