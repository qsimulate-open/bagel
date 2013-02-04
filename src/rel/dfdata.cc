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

DFData::DFData(shared_ptr<const DFDist> df, pair<int, int> coord) : dfdata_(df), coord_(coord) {

  cross_ = true;
  if (coord_.first == coord_.second) cross_ = false;

  if (coord_.first != coord_.second && (coord_.first == 2 || coord_.second == 2))
    basis_ = make_pair(0,0);
  else
    basis_ = make_pair(0,1);

}

DFData::DFData(shared_ptr<const DFData> df, bool basis, bool coord) : dfdata_(df->df()), coord_(df->coord()), basis_(df->basis_) {
  cross_ = true;
  if (coord_.first == coord_.second) cross_ = false;
  if (basis) basis_ = df->opp();
  if (coord) coord_ = df->swap();
}

double DFData::cross_coeff() const {
  assert(cross_ == true);
  if (coord_.first == 1 || coord_.second == 1) {
    return -1.0;
  } else {
    return 1.0;
  }

}

shared_ptr<const DFData> DFData::opp_basis() {
  return shared_ptr<const DFData>(new DFData(shared_from_this(), true, false));
}

shared_ptr<const DFData> DFData::swap_df() {
  return shared_ptr<const DFData>(new DFData(shared_from_this(), false, true));
}

shared_ptr<const DFData> DFData::opp_and_swap() {
  return shared_ptr<const DFData>(new DFData(shared_from_this(), true, true));
}

pair<int, int> DFData::opp() const {
  pair<int, int> basis = basis_;
  basis.first  = 1 - basis.first;
  basis.second = 1 - basis.second;
  return basis;
}

pair<int, int> DFData::swap() const {
  pair<int, int> coord = coord_;
  int interm = coord.first;
  coord.first = coord.second;
  coord.second = interm;
  return coord;
}

const int DFData::coeff_index() const {
  if (coord_.first == 3)
    return 0;
  else
    return 2;
}
