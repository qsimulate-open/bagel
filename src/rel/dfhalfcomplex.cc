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

#if 0
array<shared_ptr<Matrix>, 2> DFHalfComplex::form_2index_large_large(shared_ptr<DFHalfComplex> dfc) {

  initialize_data_();

  *data_[0] -= *dfdata_[0]->form_2index(dfc->get_real(), 1.0);
  *data_[0] += *dfdata_[1]->form_2index(dfc->get_imag(), -1.0);
  *data_[1] -= *dfdata_[1]->form_2index(dfc->get_real(), 1.0);
  *data_[1] -= *dfdata_[0]->form_2index(dfc->get_imag(), -1.0);

  return data_;

}
array<shared_ptr<Matrix>, 2> DFHalfComplex::compute_large_Jop(shared_ptr<const DFDist> df, shared_ptr<const Matrix> trcoeff, shared_ptr<const Matrix> ticoeff) {

  initialize_data_();

  *data_[0] += *df->compute_Jop(dfdata_[0], trcoeff, true);
  *data_[0] += *df->compute_Jop(dfdata_[1], ticoeff, true);
  *data_[1] -= *df->compute_Jop(dfdata_[0], ticoeff, true);
  *data_[1] += *df->compute_Jop(dfdata_[1], trcoeff, true);

  return data_;
}
#endif

array<shared_ptr<Matrix>, 2> DFHalfComplex::form_2index(shared_ptr<DFHalfComplex> dfc) {

  initialize_data_();

  pair<const int, const double> coeff = compute_coeff(dfc->basis(), dfc->coord());

  *data_[coeff.first] -= *dfdata_[0]->form_2index(dfc->get_real(), 1.0) * coeff.second;
  *data_[coeff.first] += *dfdata_[1]->form_2index(dfc->get_imag(), -1.0) * coeff.second;
  *data_[1-coeff.first] -= *dfdata_[0]->form_2index(dfc->get_imag(), -1.0) * coeff.second * (coeff.first == 1 ? -1.0 : 1.0);
  *data_[1-coeff.first] -= *dfdata_[1]->form_2index(dfc->get_real(), 1.0) * coeff.second * (coeff.first == 1 ? -1.0 : 1.0);

  return data_;
}

array<shared_ptr<Matrix>, 2> DFHalfComplex::compute_Jop(shared_ptr<const DFDist> df, shared_ptr<const Matrix> trocoeff, 
                 shared_ptr<const Matrix> tiocoeff, pair<const int, const int> basis, pair<const int, const int> coord) {

  initialize_data_();

  pair<const int, const double> coeff = compute_coeff(basis, coord);

  *data_[coeff.first] += *df->compute_Jop(dfdata_[0], trocoeff, true) * coeff.second;
  *data_[coeff.first] += *df->compute_Jop(dfdata_[1], tiocoeff, true) * coeff.second;
  *data_[1-coeff.first] -= *df->compute_Jop(dfdata_[0], tiocoeff, true) * coeff.second * (coeff.first == 1 ? -1.0 : 1.0);
  *data_[1-coeff.first] += *df->compute_Jop(dfdata_[1], trocoeff, true) * coeff.second * (coeff.first == 1 ? -1.0 : 1.0);

  return data_;
}

#if 0
array<shared_ptr<Matrix>, 2> DFHalfComplex::compute_small_Jop(shared_ptr<const DFDist> dfs, array<shared_ptr<const Matrix>, 4> trocoeff, 
                 array<shared_ptr<const Matrix>, 4> tiocoeff, pair<const int, const int> basis, pair<const int, const int> coord) {

  if(coord.first == -1)
    throw logic_error("Can only call compute_small_Jop with small DFDist quantity");

  initialize_data_();

  pair<const int, const double> coeff = compute_coeff(basis, coord);

  // Match trocoeff index to index in dfock
  const int coeff_basis = basis_.second + 3;

  *data_[coeff.first] += *dfs->compute_Jop(dfdata_[0], trocoeff[coeff_basis], true) * coeff.second;
  *data_[coeff.first] -= *dfs->compute_Jop(dfdata_[1], tiocoeff[coeff_basis], true) * coeff.second;
  *data_[1-coeff.first] += *dfs->compute_Jop(dfdata_[0], tiocoeff[coeff_basis], true) * coeff.second * (coeff.first == 1 ? -1.0 : 1.0);
  *data_[1-coeff.first] += *dfs->compute_Jop(dfdata_[1], trocoeff[coeff_basis], true) * coeff.second * (coeff.first == 1 ? -1.0 : 1.0);

  return data_;
}
#endif

pair<const int, const double> DFHalfComplex::compute_coeff(pair<const int, const int> basis2, pair<const int, const int> coord2) {
  const complex<double> rcoeff (1.0, 0.0);
  const complex<double> icoeff (0.0, 1.0);

  // Xa Xb, Ya Yb, Za Zb respectively
  complex<double> coeff[3][2] = {{rcoeff, rcoeff}, {icoeff, -icoeff}, {rcoeff, -rcoeff}};

  complex<double> coeff1 = rcoeff, coeff2 = rcoeff, coeff3 = rcoeff, coeff4 = rcoeff;

  // |coord_ coord_*)(coord2* coord2)
  if (coord_.first != -1) {
    coeff1 = coeff[coord_.first][basis_.first];
    coeff2 = (coord_.second == 1) ? -coeff[coord_.second][basis_.second] : coeff[coord_.second][basis_.second]; 
  }
  if (coord2.first != -1) {
    coeff3 = (coord2.first == 1) ? -coeff[coord2.first][basis2.first] : coeff[coord2.first][basis2.first];
    coeff4 = coeff[coord2.second][basis2.second];
  }

  
  complex<double> out = coeff1*coeff2*coeff3*coeff4;
  return (out.imag() == 0.0) ? make_pair(0, out.real()) : make_pair(1, out.imag());
}

void DFHalfComplex::initialize_data_() {
  data_[0] = shared_ptr<Matrix>(new Matrix(dim_, dim_));
  data_[1] = shared_ptr<Matrix>(new Matrix(dim_, dim_));
  data_[0]->zero();
  data_[1]->zero();
}
  
