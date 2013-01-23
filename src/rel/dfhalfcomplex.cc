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

DFHalfComplex::DFHalfComplex(const shared_ptr<const DFDist> df, shared_ptr<const Matrix> rcoeff, shared_ptr<const Matrix> icoeff, 
                             const bool swap, pair<const int, const int> coord, pair<const int, const int> basis) : coord_(coord), basis_(basis) {

  dim_ = rcoeff->ndim();
  shared_ptr<DFHalfDist> rhalfbj;
  shared_ptr<DFHalfDist> ihalfbj;

  if (swap == true) {
    rhalfbj = df->compute_half_transform_swap(rcoeff);
    ihalfbj = df->compute_half_transform_swap(icoeff); 
  } else {
    rhalfbj = df->compute_half_transform(rcoeff);
    ihalfbj = df->compute_half_transform(icoeff); 
  }

  dfdata_[0] = rhalfbj->apply_J();
  dfdata_[1] = ihalfbj->apply_J();
}

#if 1
array<shared_ptr<Matrix>, 2> DFHalfComplex::form_2index_complex(shared_ptr<DFHalfComplex> dfc) {

//  pair<const int, const int> index = make_pair(basis_.second, dfc->basis().second);

  array<shared_ptr<Matrix>, 2> data;
  data[0] = shared_ptr<Matrix>(new Matrix(dim_, dim_));
  data[1] = shared_ptr<Matrix>(new Matrix(dim_, dim_));
  data[0]->zero();
  data[1]->zero();

  pair<const int, const double> coeff = compute_coeff(dfc->basis(), dfc->coord());

  *data[coeff.first] += *dfdata_[0]->form_2index(dfc->get_real(), 1.0) * coeff.second;
  *data[coeff.first] += *dfdata_[1]->form_2index(dfc->get_imag(), -1.0) * coeff.second;
  *data[1-coeff.first] += *dfdata_[0]->form_2index(dfc->get_imag(), 1.0) * coeff.second * (coeff.first == 1 ? -1.0 : 1.0);
  *data[1-coeff.first] += *dfdata_[1]->form_2index(dfc->get_real(), 1.0) * coeff.second * (coeff.first == 1 ? -1.0 : 1.0);

  return data;
}
#endif

pair<const int, const double> DFHalfComplex::compute_coeff(pair<const int, const int> basis2, pair<const int, const int> coord2) {
  const complex<double> rcoeff (1.0, 0.0);
  const complex<double> icoeff (0.0, 1.0);

  // Xa Xb, Ya Yb, Za Zb respectively
  complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {icoeff, -icoeff}, {rcoeff, -rcoeff}};

  // |coord_ coord_*)(coord2* coord2)
  complex<double> coeff1 = coeff[coord_.first][basis_.first];
  complex<double> coeff2 = (coord_.second == 1) ? -coeff[coord_.second][basis_.second] : coeff[coord_.second][basis_.second]; 
  complex<double> coeff3 = (coord2.first == 1) ? -coeff[coord2.first][basis2.first] : coeff[coord2.first][basis2.first];
  complex<double> coeff4 = coeff[coord2.second][basis2.second];

  complex<double> out = coeff1*coeff2*coeff3*coeff4;
  return (out.imag() == 0.0) ? make_pair(0, out.real()) : make_pair(1, out.imag());
}
