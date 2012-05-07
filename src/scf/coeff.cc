//
// Newint - Parallel electron correlation program.
// Filename: coeff.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <algorithm>
#include <src/scf/coeff.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>

using namespace std;

Coeff::Coeff(const Matrix1e& inp) : Matrix1e(inp.geom()) {

  const int ndim_ = inp.ndim();
  const int mdim_ = inp.mdim();
  dcopy_(nbasis_*nbasis_, inp.data(), 1, data(), 1); 

}


Coeff::~Coeff() {

}


shared_ptr<Matrix1e> Coeff::form_density_rhf(const int n, const int offset) const {
  shared_ptr<Matrix1e> out(new Matrix1e(geom_));
  double* out_data = out->data() + offset*nbasis_;

  dgemm_("N", "T", nbasis_, nbasis_, n, 2.0, data(), nbasis_, data(), nbasis_, 0.0, out_data, nbasis_); 

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
