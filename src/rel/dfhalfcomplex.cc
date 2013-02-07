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


#include <src/rel/dfhalfcomplex.h>

using namespace std;
using namespace bagel;

DFHalfComplex::DFHalfComplex(shared_ptr<const DFData> df, shared_ptr<const Matrix> rcoeff, shared_ptr<const Matrix> icoeff)
                              : coord_(df->coord()), basis_(df->basis()) {

  coeff_ = calc_coeff(coord_, basis_);
  shared_ptr<DFHalfDist> rhalfbj;
  shared_ptr<DFHalfDist> ihalfbj;

  if (df->swapped()) {
    rhalfbj = df->df()->compute_half_transform_swap(rcoeff);
    ihalfbj = df->df()->compute_half_transform_swap(icoeff); 
  } else {
    rhalfbj = df->df()->compute_half_transform(rcoeff);
    ihalfbj = df->df()->compute_half_transform(icoeff); 
  }

  dfhalf_[0] = rhalfbj->apply_J();
  dfhalf_[1] = ihalfbj->apply_J();

}


void DFHalfComplex::set_sum_diff() {
  df2_[0] = dfhalf_[0]->copy();
  df2_[0]->daxpy(1.0, dfhalf_[1]);
  df2_[1] = dfhalf_[0]->copy();
  df2_[1]->daxpy(-1.0, dfhalf_[1]);
}


void DFHalfComplex::zaxpy(std::complex<double> a, std::shared_ptr<const DFHalfComplex> o) {
  if (imag(a) == 0.0) {
    const double fac = real(a);
    dfhalf_[0]->daxpy(fac, o->dfhalf_[0]);
    dfhalf_[1]->daxpy(fac, o->dfhalf_[1]);
  } else if (real(a) == 0.0) {
    const double fac = imag(a);
    dfhalf_[0]->daxpy(-fac, o->dfhalf_[1]);
    dfhalf_[1]->daxpy( fac, o->dfhalf_[0]);
  } else {
    throw logic_error("DFHalfComplex::zaxpy can be called by real or imaginary coeff (and not complex)");
  }
}

const tuple<int, int> DFHalfComplex::compute_index_Exop(shared_ptr<const DFHalfComplex> o) const {
  const pair<const int, const int> basis2 = o->basis();
  const pair<const int, const int> coord2 = o->coord();

  // 4x4 ZMatrix starting at 0,0 (large, large) or 0,2n (large, small) or 2n,0 (small, large) or 2n,2n (small)
  const int start1 = coord_.first == DFData::Comp::L ? 0 : 2;
  const int start2 = coord2.first == DFData::Comp::L ? 0 : 2;
  //go from small large to large small or vice versa
  const int index1 = start1 + basis_.second;
  const int index2 = start2 + basis2.second;

  return make_tuple(index1, index2);
}


int DFHalfComplex::coeff_matrix() const {
  return coord_.first == DFData::Comp::L ? basis_.second : basis_.second + 2;
}


bool DFHalfComplex::matches(shared_ptr<DFHalfComplex> o) const {
  return coord_.second == o->coord().second && basis_.second == o->basis().second;
}


complex<double> DFHalfComplex::factor(shared_ptr<const DFHalfComplex> o) const {
  pair<const int, const int> coord2 = o->coord();
  pair<const int, const int> basis2 = o->basis();

  const complex<double> rcoeff (1.0, 0.0);
  const complex<double> icoeff (0.0, 1.0);

  complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {icoeff, -icoeff}, {rcoeff, -rcoeff}};
//complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {-icoeff, icoeff}, {rcoeff, -rcoeff}};

  complex<double> prod = coeff[coord2.first][basis2.first] / coeff[coord_.first][basis_.first];

  return prod;
}

pair<complex<double>, complex<double>> DFHalfComplex::calc_coeff(pair<const int, const int> coord, pair<const int, const int> basis) {
#if 1
  // Xa Xb, Ya Yb, Za Zb respectively
  const double tc = 1.0 / (2.0* c__);
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



