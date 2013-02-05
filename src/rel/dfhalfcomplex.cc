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


complex<double> DFHalfComplex::compute_coeff(shared_ptr<const DFData> o) const {
  const pair<const int, const int> basis2 = o->basis();
  const pair<const int, const int> coord2 = o->coord();
  return compute_coeff(basis2, coord2);
}


complex<double> DFHalfComplex::compute_coeff(shared_ptr<const DFHalfComplex> o) const {
  const pair<const int, const int> basis2 = o->basis();
  const pair<const int, const int> coord2 = o->coord();
  return compute_coeff(basis2, coord2);
}


complex<double> DFHalfComplex::compute_coeff(pair<const int, const int> basis2, pair<const int, const int> coord2) const {
  const double tc = 1.0 / (2.0* c__);

  // Xa Xb, Ya Yb, Za Zb respectively
  const complex<double> rcoeff (1.0, 0.0);
  const complex<double> icoeff (0.0, 1.0);
  // TODO check again the y coeff
  complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {icoeff, -icoeff}, {rcoeff, -rcoeff}};
//complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {-icoeff, icoeff}, {rcoeff, -rcoeff}};

  complex<double> prod = 1.0;
  double power = 0.0;

  // |coord_ coord_*)(coord2* coord2)
  if (coord_.first != DFData::Comp::L) {
    power += 2.0;
    prod *= coeff[coord_.first][basis_.first];
    // if the bra vector is y, we need -1 for taking the conjugate
    prod *= conj(coeff[coord_.second][basis_.second]);
  }
  if (coord2.first != DFData::Comp::L) {
    power += 2.0;
    // if the bra vector is y, we need -1 for taking the conjugate
    prod *= conj(coeff[coord2.first][basis2.first]);
    prod *= coeff[coord2.second][basis2.second];
  }

  return prod * ::pow(tc, power);
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


const tuple<int, int, int, int> DFHalfComplex::compute_index_Jop(shared_ptr<const DFData> o) const {
  const pair<const int, const int> basis = o->basis();
  const pair<const int, const int> coord = o->coord();
  // 4x4 ZMatrix either starting at 0,0 (large) or 2n,2n (small)
  const int start = coord.first == DFData::Comp::L ? 0 : 2;
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
