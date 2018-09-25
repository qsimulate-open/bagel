//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: coeff.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <algorithm>
#include <cassert>
#include <src/util/f77.h>
#include <src/wfn/coeff.h>

using namespace std;
using namespace bagel;
using namespace btas;


template <typename MatType, class Enable>
Coeff_<MatType, Enable>::Coeff_(const MatType& inp) : MatType(inp.ndim(), inp.mdim()) {
  copy_n(inp.data(), size(), data());
}


template <typename MatType, class Enable>
Coeff_<MatType, Enable>::Coeff_(MatType&& inp) : MatType(move(inp)) {
}


template <typename MatType, class Enable>
Coeff_<MatType, Enable>::Coeff_(vector<shared_ptr<const Coeff_<MatType>>> coeff_vec) : MatType(num_basis(coeff_vec), num_basis(coeff_vec)) {

  DataType* cdata = data();
  for(auto icoeff = coeff_vec.begin(); icoeff != coeff_vec.end(); ++icoeff) {
    const DataType* cur_data = (*icoeff)->data();

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


template <typename MatType, class Enable>
int Coeff_<MatType, Enable>::num_basis(vector<shared_ptr<const Coeff_<MatType>>> coeff_vec) const {
  return accumulate(coeff_vec.begin(), coeff_vec.end(), 0, [](const int& a, shared_ptr<const Coeff_<MatType>>& b) { return a+b->ndim(); });
}


template <typename MatType, class Enable>
shared_ptr<MatType> Coeff_<MatType, Enable>::form_weighted_density_rhf(const int n, const VecView e) const {
  auto out = make_shared<MatType>(ndim(), ndim());
  for (int i = 0; i != n; ++i) {
    auto sl = slice(i, i+1);
    contract(2.0*e(i), sl, {0,1}, sl, {2,1}, 1.0, *out, {0,2}, false, true);
  }
  return out;
}


template <typename MatType, class Enable>
pair<shared_ptr<MatType>, shared_ptr<MatType>> Coeff_<MatType, Enable>::split(const int nrow1, const int nrow2) const {
  auto out1 = make_shared<MatType>(nrow1, mdim());
  auto out2 = make_shared<MatType>(nrow2, mdim());

  assert(nrow1+nrow2 == ndim());

  const DataType* source = data();
  DataType* data1 = out1->data();
  DataType* data2 = out2->data();

  for (int m = 0; m != mdim(); ++m, data1+=out1->ndim(), data2+=out2->ndim(), source+=ndim()) {
    copy_n(source,       nrow1, data1);
    copy_n(source+nrow1, nrow2, data2);
  }

  return {out1, out2};
}

template class bagel::Coeff_<Matrix>;
template class bagel::Coeff_<ZMatrix>;

BOOST_CLASS_EXPORT_IMPLEMENT(Coeff)
BOOST_CLASS_EXPORT_IMPLEMENT(ZCoeff)

