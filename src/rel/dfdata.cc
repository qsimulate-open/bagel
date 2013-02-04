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

DFData::DFData(shared_ptr<const DFDist> df, pair<int, int> coord) : dfdata_(df), coord_(coord), swap_(false) {

  if ((coord_.first == Comp::Z) ^ (coord_.second == Comp::Z))
    basis_ = make_pair(Basis::a, Basis::b);
  else
    basis_ = make_pair(Basis::a, Basis::a);

}

DFData::DFData(const DFData& o, bool bas, bool coo) : dfdata_(o.df()), coord_(o.coord()), basis_(o.basis_), swap_(o.swap_) {
  if (bas) {
    basis_.first ^= 1;
    basis_.second ^= 1;
  }
  if (coo) {
    std::swap(coord_.first, coord_.second); 
    swap_ ^= true;
  }
}


double DFData::cross_coeff() const {
  assert(cross() == true);
  return ((coord_.first == Comp::Y) ^ (coord_.second == Comp::Y)) ? -1.0 : 1.0;
}


shared_ptr<const DFData> DFData::opp() {
  return shared_ptr<const DFData>(new DFData(*this, true, false));
}


shared_ptr<const DFData> DFData::swap() {
  return shared_ptr<const DFData>(new DFData(*this, false, true));
}


shared_ptr<const DFData> DFData::opp_and_swap() {
  return shared_ptr<const DFData>(new DFData(*this, true, true));
}


const int DFData::coeff_index() const {
  return coord_.first == Comp::L ? 0 : 2;
}
