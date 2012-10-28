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
// The BAGEL package is free software; you can redistribute it and\/or modify
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

#include <algorithm>
#include <src/scf/coeff.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;

Coeff::Coeff(const Matrix1e& inp) : Matrix1e(inp.geom()) {

  ndim_ = inp.ndim();
  mdim_ = inp.mdim();
  copy_n(inp.data(), nbasis_*nbasis_, data());

}

Coeff::Coeff(vector<shared_ptr<const Coeff> > coeff_vec) : Matrix1e(supergeom(coeff_vec)) {
  /* Force ndim_ = nbasis_ to avoid any unexpected behavior */
  ndim_ = nbasis_; mdim_ = 0;
  for(auto icoeff = coeff_vec.begin(); icoeff != coeff_vec.end(); ++icoeff) {
    mdim_ += (*icoeff)->mdim();
  }

  double* cdata = data();
  for(auto icoeff = coeff_vec.begin(); icoeff != coeff_vec.end(); ++icoeff) {
    double* cur_data = (*icoeff)->data();

    int cur_nstart = 0;
    for(auto iz0 = coeff_vec.begin(); iz0 != icoeff; ++iz0) { cur_nstart += (*iz0)->nbasis(); }

    int cur_nbasis = (*icoeff)->nbasis();

    int cur_nend = nbasis_ - (cur_nstart + cur_nbasis);

    int cur_mdim = (*icoeff)->mdim();

    /* The matrix is initialized to zero, so don't need to fill in zeros, just step past them. */
    for(int mm = 0; mm != cur_mdim; ++mm) {
      /* Step past first zeros */
      cdata += cur_nstart;
      /* Copy elements from current Coeff */
      cdata = copy(cur_data, cur_data + cur_nbasis, cdata);
      cur_data += cur_nbasis;
      /* Step past empty elements at the end */
      cdata += cur_nend;
    }
  }
}

Coeff::~Coeff() {

}

shared_ptr<const Geometry> Coeff::supergeom(vector<shared_ptr<const Coeff> > coeff_vec) {
  vector<shared_ptr<const Geometry> > geovec;
  for(auto icoeff = coeff_vec.begin(); icoeff != coeff_vec.end(); ++icoeff) {
    geovec.push_back((*icoeff)->geom());
  }

  shared_ptr<const Geometry> out(new const Geometry(geovec));

  return out;
}

shared_ptr<Matrix1e> Coeff::form_density_rhf(const int n, const int offset) const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_));
  double* out_data = out->data() + offset*nbasis_;

  dgemm_("N", "T", nbasis_, nbasis_, n, 2.0, data(), nbasis_, data(), nbasis_, 0.0, out_data, nbasis_);

  return out;
}


shared_ptr<Matrix1e> Coeff::form_weighted_density_rhf(const int n, const vector<double>& e, const int offset) const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_));
  double* out_data = out->data() + offset*nbasis_;
  double* cdata = data();
  for (int i = 0; i != n; ++i, cdata += nbasis_) {
    dgemm_("N", "T", nbasis_, nbasis_, 1, 2.0*e[i], cdata, nbasis_, cdata, nbasis_, 1.0, out_data, nbasis_);
  }
  return out;
}


pair<shared_ptr<Coeff>, shared_ptr<Coeff> > Coeff::split(const int nrow1, const int nrow2) const {
  shared_ptr<Coeff> out1(new Coeff(geom_, nrow1, mdim_));
  shared_ptr<Coeff> out2(new Coeff(geom_, nrow2, mdim_));

  assert(nrow1+nrow2 == ndim_);

  const double* source = data();
  double* data1 = out1->data();
  double* data2 = out2->data();

  for (int m = 0; m != mdim_; ++m, data1+=out1->nbasis(), data2+=out2->nbasis(), source+=nbasis_) {
    copy(source,       source+nrow1,       data1);
    copy(source+nrow1, source+nrow1+nrow2, data2);
  }

  return make_pair(out1, out2);
}
