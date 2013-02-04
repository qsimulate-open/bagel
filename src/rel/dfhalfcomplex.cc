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

  dim_ = rcoeff->ndim();
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


complex<double> DFHalfComplex::compute_coeff(pair<const int, const int> basis2, pair<const int, const int> coord2) {
  const int large = 3;
  const double tc = 1.0 / (2.0* c__);
  double power = 0.0;
  const complex<double> rcoeff (1.0, 0.0);
  const complex<double> icoeff (0.0, 1.0);

  // Xa Xb, Ya Yb, Za Zb respectively
  complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {icoeff, -icoeff}, {rcoeff, -rcoeff}};

  complex<double> coeff1 = rcoeff, coeff2 = rcoeff, coeff3 = rcoeff, coeff4 = rcoeff;

  // |coord_ coord_*)(coord2* coord2)
  if (coord_.first != large) {
    power+=2.0;
    coeff1 = coeff[coord_.first][basis_.first];
    coeff2 = (coord_.second == 1) ? -coeff[coord_.second][basis_.second] : coeff[coord_.second][basis_.second]; 
  }
  if (coord2.first != large) {
    power+=2.0;
    coeff3 = (coord2.first == 1) ? -coeff[coord2.first][basis2.first] : coeff[coord2.first][basis2.first];
    coeff4 = coeff[coord2.second][basis2.second];
  }

  complex<double> out = coeff1 * coeff2 * coeff3 * coeff4 * ::pow(tc, power);
  return out;
}


const tuple<int, int> DFHalfComplex::compute_index_Exop(pair<const int, const int> basis2, pair<const int, const int> coord2) {

  const int large = 3;

  // 4x4 ZMatrix starting at 0,0 (large, large) or 0,2n (large, small) or 2n,0 (small, large) or 2n,2n (small)
  const int start1 = (coord_.first == large ? 0 : 2);
  const int start2 = (coord2.first == large ? 0 : 2);
  //go from small large to large small or vice versa
  const int index1 = start1 + basis_.second;
  const int index2 = start2 + basis2.second;

  return make_tuple(index1, index2);
}


const tuple<int, int, int, int> DFHalfComplex::compute_index_Jop(pair<const int, const int> basis, pair<const int, const int> coord) {
  const int large = 3;
  // 4x4 ZMatrix either starting at 0,0 (large) or 2n,2n (small)
  int start = (coord.first == large ? 0 : 2);
  // put transposed Matrices in submatrix opposite original
  const int opp1 =  1^basis.first;
  const int opp2 =  1^basis.second;

  const int index1 = start + basis.first;
  const int index2 = start + basis.second;
  const int index3 = start + opp1;
  const int index4 = start + opp2; 

  return make_tuple(index1, index2, index3, index4);
}


const int DFHalfComplex::coeff_matrix() const {
  const int large = 3;
  return coord_.first == large ? basis_.second : basis_.second + 2;
}
