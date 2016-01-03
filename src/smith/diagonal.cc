//
// BAGEL - Parallel electron correlation program.
// Filename: diagonal.cc
// Copyright (C) 2015 Shiozaki group
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

#include <src/smith/mrci/MRCI.h>
#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relcaspt2/RelCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void CASPT2::CASPT2::diagonal(shared_ptr<TATensor<double,4>> r, shared_ptr<const TATensor<double,4>> t) const {
  const int ncore = info_->ncore();
  const int nocc  = info_->nclosed() + info_->nact();
  const VecView eig = eig_;
  TATensor<double,4> i0({closed_, virt_, closed_, virt_}, true);
  i0("c2,a3,c0,a1") = (*t)("c0,a1,c2,a3")*8.0 - (*t)("c0,a3,c2,a1")*4.0;
  foreach_inplace(i0, [&](typename TATensor<double,4>::value_type& tile) {
    auto range = tile.range();
    auto lo = range.lobound();
    auto up = range.upbound();
    size_t n = 0;
    for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
      for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
        for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
          for (size_t i3 = lo[3]; i3 != up[3]; ++i3)
            tile[n++] *= - eig(i3+ncore) + eig(i2+nocc) - eig(i1+ncore) + eig(i0+nocc);
  });
  (*r)("c2,a3,c0,a1") += i0("c2,a3,c0,a1");
}


void RelCASPT2::RelCASPT2::diagonal(shared_ptr<TATensor<complex<double>,4>> r, shared_ptr<const TATensor<complex<double>,4>> t) const {
  const int ncore = (info_->ncore())*2;
  const int nocc  = (info_->nclosed() + info_->nact())*2;
  const VecView eig = eig_;
  TATensor<complex<double>,4> i0({closed_, virt_, closed_, virt_}, true);
  i0("c2,a3,c0,a1") = (*t)("c0,a1,c2,a3");
  foreach_inplace(i0, [&](typename TATensor<complex<double>,4>::value_type& tile) {
    auto range = tile.range();
    auto lo = range.lobound();
    auto up = range.upbound();
    size_t n = 0;
    for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
      for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
        for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
          for (size_t i3 = lo[3]; i3 != up[3]; ++i3)
            tile[n++] *= (- eig(i3+ncore) + eig(i2+nocc) - eig(i1+ncore) + eig(i0+nocc)) * 4.0;
  });
  (*r)("c2,a3,c0,a1") += i0("c2,a3,c0,a1");
}


// this function takes care of 4-external.
void MRCI::MRCI::diagonal(shared_ptr<TATensor<double,4>> r, shared_ptr<const TATensor<double,4>> t) const {
  const bool diag = (*rdm0_)("") == 1.0;
  if (diag)
    (*r)("c3,a4,c1,a2") += (*v2_)("a6,a2,a5,a4") * ((*t)("c1,a6,c3,a5")*8.0 - (*t)("c1,a5,c3,a6")*4.0);
  (*r)("c2,a3,x0,a1") += (*v2_)("a5,a1,a4,a3") * (((*t)("x1,a5,c2,a4")*2.0 - (*t)("x1,a4,c2,a5")) * (*rdm1_)("x1,x0"));
  (*r)("x1,a2,x0,a1") += (*v2_)("a4,a1,a3,a2") * ((*t)("x3,a3,x2,a4") * ((*rdm2_)("x3,x1,x2,x0") * 2.0));
}


// this function takes care of 4-external.
void RelMRCI::RelMRCI::diagonal(shared_ptr<TATensor<std::complex<double>,4>> r, shared_ptr<const TATensor<std::complex<double>,4>> t) const {
  const bool diag = (*rdm0_)("") == 1.0;
  if (diag)
    (*r)("c3,a4,c1,a2") += (*v2_)("a6,a2,a5,a4") * ((*t)("c1,a6,c3,a5") * 4.0);
  (*r)("c2,a3,x0,a1") += (*v2_)("a5,a1,a4,a3") * ((*t)("x1,a5,c2,a4") * ((*rdm1_)("x1,x0") * 2.0));
  (*r)("x1,a2,x0,a1") += (*v2_)("a4,a1,a3,a2") * ((*t)("x3,a3,x2,a4") * ((*rdm2_)("x3,x1,x2,x0") * 2.0));
}

#endif
