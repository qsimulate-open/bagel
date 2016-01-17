//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spinfreebase_denom.cpp
// Copyright (C) 2015 Toru Shiozaki
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

  const double e0 = e0_;
  const VecView eig(eig_);

  for (int ist = 0; ist != t->nref(); ++ist) {
    auto rist = r->at(ist);
    t->fac(ist) = 0.0;

    // not the best structure, but I am assuming that this does not take too much time...
    for (int jst = 0; jst != t->nref(); ++jst) {
      auto tjst = t->at(jst);

      if (ist == jst) { // AACC
        TATensor<DataType,4> i0({closed_, virt_, closed_, virt_}, true);
        if (is_same<DataType,double>::value)
          i0("c0,a1,c2,a3") = (*rist)("c0,a1,c2,a3") * (1.0/6.0) + (*rist)("c0,a3,c2,a1") * (1.0/12.0);
        else
          i0("c0,a1,c2,a3") = (*rist)("c0,a1,c2,a3") * 0.25;
        foreach_inplace(i0, [=](typename TATensor<DataType,4>::value_type& tile) {
          auto range = tile.range();
          auto lo = range.lobound();
          auto up = range.upbound();
          size_t n = 0;
          for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
            for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
              for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
                for (size_t i3 = lo[3]; i3 != up[3]; ++i3)
                  tile[n++] /= - eig(i3+ncore) + eig(i2+nocc) - eig(i1+ncore) + eig(i0+nocc);
        });
        (*tjst)("c0,a1,c2,a3") -= i0("c0,a1,c2,a3");
      }
      { // AAXX
        shared_ptr<const TATensor<DataType,3>> sist = denom_->tashalf_xx(ist);
        shared_ptr<const TATensor<DataType,3>> sjst = denom_->tashalf_xx(jst);
        TATensor<DataType,3> i0({ortho2_, virt_, virt_}, true);
        i0("o4,a1,a3") = (*rist)("x0,a1,x2,a3") * (*sist)("o4,x0,x2");
        const VecView denom = denom_->denom_xx();
        foreach_inplace(i0, [=](typename TATensor<DataType,3>::value_type& tile) {
          auto range = tile.range();
          auto lo = range.lobound();
          auto up = range.upbound();
          size_t n = 0;
          for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
            for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
              for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
                tile[n++] /= denom(i2) + eig(i1+nocc) + eig(i0+nocc) - e0;
        });
        (*tjst)("x0,a1,x2,a3") -= i0("o4,a1,a3") * 0.5 * (*sjst)("o4,x0,x2").conj();
      }
      { // AACX
        shared_ptr<const TATensor<DataType,2>> sist = denom_->tashalf_x(ist);
        shared_ptr<const TATensor<DataType,2>> sjst = denom_->tashalf_x(jst);
        TATensor<DataType,4> i0({ortho1_, virt_, virt_, closed_}, true);
        if (is_same<DataType,double>::value)
          i0("o4,a1,a3,c2") = ((*rist)("c2,a3,x0,a1") * 2.0 + (*rist)("c2,a1,x0,a3")) * ((*sist)("o4,x0") * (1.0/3.0));
        else
          i0("o4,a1,a3,c2") = (*rist)("c2,a3,x0,a1") * ((*sist)("o4,x0") * 0.5);
        const VecView denom = denom_->denom_x();
        foreach_inplace(i0, [=](typename TATensor<DataType,4>::value_type& tile) {
          auto range = tile.range();
          auto lo = range.lobound();
          auto up = range.upbound();
          size_t n = 0;
          for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
            for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
              for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
                for (size_t i3 = lo[3]; i3 != up[3]; ++i3)
                  tile[n++] /= denom(i3) + eig(i2+nocc) + eig(i1+nocc) - eig(i0+ncore) - e0;
        });
        (*tjst)("x0,a1,c2,a3") -=  i0("o4,a1,a3,c2") * (*sjst)("o4,x0").conj();
      }
      { // AXCC
        shared_ptr<const TATensor<DataType,2>> sist = denom_->tashalf_h(ist);
        shared_ptr<const TATensor<DataType,2>> sjst = denom_->tashalf_h(jst);
        TATensor<DataType,4> i0({ortho1_, virt_, closed_, closed_}, true);
        if (is_same<DataType,double>::value)
          i0("o4,a1,c0,c2") = ((*rist)("c2,x3,c0,a1") * 2.0 + (*rist)("c0,x3,c2,a1")) * ((*sist)("o4,x3") * (1.0/3.0));
        else
          i0("o4,a1,c0,c2") = (*rist)("c2,x3,c0,a1") * ((*sist)("o4,x3") * 0.5);
        const VecView denom = denom_->denom_h();
        foreach_inplace(i0, [=](typename TATensor<DataType,4>::value_type& tile) {
          auto range = tile.range();
          auto lo = range.lobound();
          auto up = range.upbound();
          size_t n = 0;
          for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
            for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
              for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
                for (size_t i3 = lo[3]; i3 != up[3]; ++i3)
                  tile[n++] /= denom(i3) + eig(i2+nocc) - eig(i1+ncore) - eig(i0+ncore) - e0;
        });
        (*tjst)("c0,a1,c2,x3") -= i0("o4,a1,c0,c2") * (*sjst)("o4,x3").conj();
      }
      { // XXCC
        shared_ptr<const TATensor<DataType,3>> sist = denom_->tashalf_hh(ist);
        shared_ptr<const TATensor<DataType,3>> sjst = denom_->tashalf_hh(jst);
        TATensor<DataType,3> i0({ortho2_, closed_, closed_}, true);
        i0("o4,c0,c2") = (*rist)("c0,x1,c2,x3") * ((*sist)("o4,x1,x3") * 0.5);
        const VecView denom = denom_->denom_hh();
        foreach_inplace(i0, [=](typename TATensor<DataType,3>::value_type& tile) {
          auto range = tile.range();
          auto lo = range.lobound();
          auto up = range.upbound();
          size_t n = 0;
          for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
            for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
              for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
                tile[n++] /= denom(i2) - eig(i1+ncore) - eig(i0+ncore) - e0;
        });
        (*tjst)("c0,x1,c2,x3") -= i0("o4,c0,c2") * (*sjst)("o4,x1,x3").conj();
      }
      { // XXCX
        shared_ptr<const TATensor<DataType,4>> sist = denom_->tashalf_xhh(ist);
        shared_ptr<const TATensor<DataType,4>> sjst = denom_->tashalf_xhh(jst);
        TATensor<DataType,2> i0({ortho3_, closed_}, true);
        i0("o4,c2") = (*rist)("c2,x3,x0,x1") * (*sist)("o4,x0,x1,x3");
        const VecView denom = denom_->denom_xhh();
        foreach_inplace(i0, [=](typename TATensor<DataType,2>::value_type& tile) {
          auto range = tile.range();
          auto lo = range.lobound();
          auto up = range.upbound();
          size_t n = 0;
          for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
            for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
              tile[n++] /= denom(i1) - eig(i0+ncore) - e0;
        });
        (*tjst)("x0,x1,c2,x3") -= i0("o4,c2") * (*sjst)("o4,x0,x1,x3").conj();
      }
      { // AXXX
        shared_ptr<const TATensor<DataType,4>> sist = denom_->tashalf_xxh(ist);
        shared_ptr<const TATensor<DataType,4>> sjst = denom_->tashalf_xxh(jst);
        TATensor<DataType,2> i0({ortho3_, virt_}, true);
        i0("o4,a1") = (*rist)("x2,x3,x0,a1") * (*sist)("o4,x0,x2,x3");
        const VecView denom = denom_->denom_xxh();
        foreach_inplace(i0, [=](typename TATensor<DataType,2>::value_type& tile) {
          auto range = tile.range();
          auto lo = range.lobound();
          auto up = range.upbound();
          size_t n = 0;
          for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
            for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
              tile[n++] /= denom(i1) + eig(i0+nocc) - e0;
        });
        (*tjst)("x0,a1,x2,x3") -= i0("o4,a1") * (*sjst)("o4,x0,x2,x3").conj();
      }
      { // XXCA
        const VecView denom = denom_->denom_xh();
        TATensor<DataType,3> i0({ortho2t_, virt_, closed_}, true);
        auto lambda = [=](typename TATensor<DataType,3>::value_type& tile) {
          auto range = tile.range();
          auto lo = range.lobound();
          auto up = range.upbound();
          size_t n = 0;
          for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
            for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
              for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
                tile[n++] /= denom(i2) + eig(i1+nocc) - eig(i0+ncore) - e0;
        };
        if (is_same<DataType,double>::value) {
          shared_ptr<const TATensor<DataType,3>> s0ist = denom_->template tashalf_xh<true>(ist);
          shared_ptr<const TATensor<DataType,3>> s0jst = denom_->template tashalf_xh<true>(jst);
          shared_ptr<const TATensor<DataType,3>> s1ist = denom_->template tashalf_xh<false>(ist);
          shared_ptr<const TATensor<DataType,3>> s1jst = denom_->template tashalf_xh<false>(jst);
          i0("o4,a1,c0") = (*rist)("x2,x3,c0,a1") * (*s0ist)("o4,x2,x3") + (*rist)("c0,x3,x2,a1") * (*s1ist)("o4,x2,x3");
          foreach_inplace(i0, lambda);
          (*tjst)("c0,a1,x2,x3") -= i0("o4,a1,c0") * (*s0jst)("o4,x2,x3");
          (*tjst)("x2,a1,c0,x3") -= i0("o4,a1,c0") * (*s1jst)("o4,x2,x3");
        } else {
          shared_ptr<const TATensor<DataType,3>> sist = denom_->template tashalf_xh<true>(ist);
          shared_ptr<const TATensor<DataType,3>> sjst = denom_->template tashalf_xh<true>(jst);
          i0("o4,a1,c0") = (*rist)("x2,x3,c0,a1") * (*sist)("o4,x2,x3");
          foreach_inplace(i0, lambda);
          (*tjst)("c0,a1,x2,x3") -= i0("o4,a1,c0") * (*sjst)("o4,x2,x3").conj();
        }
      }
    } // jst loop
  } // ist loop

  return t;
}

#endif
