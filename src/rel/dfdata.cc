//
// BAGEL - Parallel electron correlation program.
// Filename: dfhalfcomplex.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <src/rel/dfdata.h>

using namespace std;
using namespace bagel;

DFData::DFData(shared_ptr<const DFDist> df, pair<int, int> coord, const int alpha) : RelDFBase(coord, alpha), dfdata_(df), swap_(false) {

}

DFData::DFData(const DFData& o, bool bas, bool coo) : RelDFBase(o), dfdata_(o.df()), swap_(o.swap_) {

  if (bas) {
    basis_.first ^= 1;
    basis_.second ^= 1;
  }

  if (coo) {
    std::swap(coord_.first, coord_.second); 
    swap_ ^= true;
    std::swap(sigma1_, sigma2_);
    std::swap(spinor_[0], spinor_[1]);
    fac_ = conj(fac_);
  }

}


//swap basis
shared_ptr<const DFData> DFData::opp() {
  return shared_ptr<const DFData>(new DFData(*this, true, false));
}


//swap coord
shared_ptr<const DFData> DFData::swap() {
  return shared_ptr<const DFData>(new DFData(*this, false, true));
}


//swap both
shared_ptr<const DFData> DFData::opp_and_swap() {
  return shared_ptr<const DFData>(new DFData(*this, true, true));
}


const tuple<int, int, int, int> DFData::compute_index_Jop() const {
  // 4x4 ZMatrix either starting at 0,0 (large) 2n,0 (large,small) or 0,2n (small,large) or 2n,2n (small)
  return make_tuple(basis(0), basis(1), basis(0)^1, basis(1)^1);
}

