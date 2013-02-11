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

DFData::DFData(shared_ptr<const DFDist> df, pair<int, int> coord, const std::vector<int> alpha) : RelDFBase(coord, alpha), dfdata_(df), swap_(false) {
  common_init();
}

DFData::DFData(const DFData& o, bool coo) : RelDFBase(o), dfdata_(o.df()), swap_(o.swap_) {
  common_init();

  if (coo) {
    swap_ ^= true;
    for (auto& i : basis_) i->swap();
    std::swap(coord_.first, coord_.second); 
    std::swap(sigma1_, sigma2_);
  }

}


//swap coord
shared_ptr<const DFData> DFData::swap() const {
  return shared_ptr<const DFData>(new DFData(*this, true));
}


vector<shared_ptr<DFHalfComplex>> DFData::compute_half_transform(array<shared_ptr<const Matrix>,4> rc, array<shared_ptr<const Matrix>,4> ic) const {
  vector<shared_ptr<DFHalfComplex>> out;
  for (auto& i : basis())
    out.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(shared_from_this(), i, rc, ic)));
  return out;
}
