//
// BAGEL - Parallel electron correlation program.
// Filename: coeff.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <algorithm>
#include <src/scf/coeff.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;
using namespace btas;

BOOST_CLASS_EXPORT_IMPLEMENT(Coeff)

Coeff::Coeff(const Matrix& inp) : Matrix(inp.ndim(), inp.mdim()) {
  copy_n(inp.data(), size(), data());
}


Coeff::Coeff(Matrix&& inp) : Matrix(move(inp)) {
}


Coeff::Coeff(vector<shared_ptr<const Coeff>> coeff_vec) : Matrix(num_basis(coeff_vec), num_basis(coeff_vec)) {

  double* cdata = data();
  for(auto icoeff = coeff_vec.begin(); icoeff != coeff_vec.end(); ++icoeff) {
    const double* cur_data = (*icoeff)->data();

    int cur_nstart = 0;
    for(auto iz0 = coeff_vec.begin(); iz0 != icoeff; ++iz0) { cur_nstart += (*iz0)->ndim(); }

    int cur_nbasis = (*icoeff)->ndim();

    int cur_nend = ndim() - (cur_nstart + cur_nbasis);

    int cur_mdim = (*icoeff)->mdim();

    /* The matrix is initialized to zero, so don't need to fill in zeros, just step past them. */
    for(int mm = 0; mm != cur_mdim; ++mm) {
      /* Step past first zeros */
      cdata += cur_nstart;
      /* Copy elements from current Coeff */
      cdata = copy_n(cur_data, cur_nbasis, cdata);
      cur_data += cur_nbasis;
      /* Step past empty elements at the end */
      cdata += cur_nend;
    }
  }
}


int Coeff::num_basis(vector<shared_ptr<const Coeff>> coeff_vec) const {
  return accumulate(coeff_vec.begin(), coeff_vec.end(), 0, [](const int& a, shared_ptr<const Coeff>& b) { return a+b->ndim(); });
}


shared_ptr<Matrix> Coeff::form_weighted_density_rhf(const int n, const vector<double>& e) const {
  auto out = make_shared<Matrix>(ndim(), ndim());
  for (int i = 0; i != n; ++i) {
    auto sl = slice(i, i+1);
    contract(2.0*e[i], *sl, {0,1}, *sl, {2,1}, 1.0, *out, {0,2});
  }
  return out;
}


pair<shared_ptr<Matrix>, shared_ptr<Matrix>> Coeff::split(const int nrow1, const int nrow2) const {
  auto out1 = make_shared<Matrix>(nrow1, mdim());
  auto out2 = make_shared<Matrix>(nrow2, mdim());

  assert(nrow1+nrow2 == ndim());

  const double* source = data();
  double* data1 = out1->data();
  double* data2 = out2->data();

  for (int m = 0; m != mdim(); ++m, data1+=out1->ndim(), data2+=out2->ndim(), source+=ndim()) {
    copy_n(source,       nrow1, data1);
    copy_n(source+nrow1, nrow2, data2);
  }

  return {out1, out2};
}
