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

  coeff_ = calc_coeff(coord_, basis_);

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

  coeff_ = calc_coeff(coord_, basis_);
}


double DFData::cross_coeff() const {
  assert(cross() == true);
  return ((coord_.first == Comp::Y) ^ (coord_.second == Comp::Y)) ? -1.0 : 1.0;
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


int DFData::coeff_index() const {
  return coord_.first == Comp::L ? 0 : 2;
}


const tuple<int, int, int, int> DFData::compute_index_Jop() const {
  // 4x4 ZMatrix either starting at 0,0 (large) or 2n,2n (small)
  const int start = coord_.first == DFData::Comp::L ? 0 : 2;
  // put transposed Matrices in submatrix opposite original
  const int opp1 =  1^basis_.first;
  const int opp2 =  1^basis_.second;

  const int index1 = start + basis_.first;
  const int index2 = start + basis_.second;
  const int index3 = start + opp1;
  const int index4 = start + opp2;

  return make_tuple(index1, index2, index3, index4);
}


const tuple<int, int, int, int> DFData::compute_index_mixed_Jop() const {
  // 4x4 ZMatrix either starting at 2n,0 (large,small) or 0,2n (small,large)
  const int start1 = coord_.first == DFData::Comp::L ? 0 : 2;
  const int start2 = coord_.second == DFData::Comp::L ? 0 : 2;
  // put transposed Matrices in submatrix opposite original
  const int opp1 =  1^basis_.first;
  const int opp2 =  1^basis_.second;

  const int index1 = start1 + basis_.first;
  const int index2 = start2 + basis_.second;
  const int index3 = start1 + opp1;
  const int index4 = start2 + opp2;

  return make_tuple(index1, index2, index3, index4);
}

pair<complex<double>, complex<double>> DFData::calc_coeff(pair<const int, const int> coord, pair<const int, const int> basis) {
#if 1
  // Xa Xb, Ya Yb, Za Zb respectively
  const double tc = 1.0 / (2.0 * c__);
  const complex<double> rcoeff (1.0, 0.0);
  const complex<double> icoeff (0.0, 1.0);
  // TODO check again the y coeff
  complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {icoeff, -icoeff}, {rcoeff, -rcoeff}};
  pair<complex<double>, complex<double>> out = make_pair(rcoeff,rcoeff);
  double power = 0.0;

  if (coord.first != DFData::Comp::L) {
    power += 1.0;
    out.first *= coeff[coord.first][basis.first];
    out.first = out.first * ::pow(tc, power);
  }

  power = 0.0;

  if (coord.second != DFData::Comp::L) {
    power += 1.0;
    out.second *= coeff[coord.second][basis.second];
    out.second = out.second * ::pow(tc, power);
  }
  return out;
#endif


}


