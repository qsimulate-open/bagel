//
// BAGEL - Parallel electron correlation program.
// Filename: reldffull.cc
// Copyright (C) 2013 Toru Shiozaki
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


#include <src/rel/reldffull.h>

using namespace std;
using namespace bagel;

RelDFFull::RelDFFull(shared_ptr<const RelDFHalf> df, array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff) : RelDFBase(*df) {

  basis_ = df->basis();
  if (basis_.size() != 1)
    throw logic_error("RelDFFull should be called with basis_.size() == 1");

  const int index = basis_.front()->basis(1);

  // TODO this could be cheaper by using a zgemm3m-type algorithm
  shared_ptr<DFFullDist> rfullbj = df->get_real()->compute_second_transform(rcoeff[index]);
              rfullbj->daxpy(-1.0, df->get_imag()->compute_second_transform(icoeff[index]));

  shared_ptr<DFFullDist> ifullbj = df->get_imag()->compute_second_transform(rcoeff[index]);
              ifullbj->daxpy( 1.0, df->get_real()->compute_second_transform(icoeff[index]));

  dffull_[0] = rfullbj->apply_J();
  dffull_[1] = ifullbj->apply_J();
}


RelDFFull::RelDFFull(const RelDFFull& o) : RelDFBase(o.cartesian_) {
  basis_ = o.basis_;
  dffull_[0] = o.dffull_[0]->copy();
  dffull_[1] = o.dffull_[1]->copy();
}


void RelDFFull::zaxpy(std::complex<double> a, std::shared_ptr<const RelDFFull> o) {
  if (imag(a) == 0.0) {
    const double fac = real(a);
    dffull_[0]->daxpy(fac, o->dffull_[0]);
    dffull_[1]->daxpy(fac, o->dffull_[1]);
  } else if (real(a) == 0.0) {
    const double fac = imag(a);
    dffull_[0]->daxpy(-fac, o->dffull_[1]);
    dffull_[1]->daxpy( fac, o->dffull_[0]);
  } else {
    const double rfac = real(a);
    dffull_[0]->daxpy(rfac, o->dffull_[0]);
    dffull_[1]->daxpy(rfac, o->dffull_[1]);
    const double ifac = imag(a);
    dffull_[0]->daxpy(-ifac, o->dffull_[1]);
    dffull_[1]->daxpy( ifac, o->dffull_[0]);
  }
}


list<shared_ptr<RelDFHalfB>> RelDFFull::back_transform(array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff) const {
  list<shared_ptr<RelDFHalfB>> out;
  assert(basis_.size() == 1);
  const int alpha = basis_[0]->alpha_comp();

  for (int i = 0; i != 4; ++i) {
    // Note that icoeff should be scaled by -1.0

    shared_ptr<DFHalfDist> real = dffull_[0]->back_transform(rcoeff[i]);
    real->daxpy( 1.0, dffull_[1]->back_transform(icoeff[i]));

    shared_ptr<DFHalfDist> imag = dffull_[0]->back_transform(icoeff[i]);
    imag->daxpy(-1.0, dffull_[0]->back_transform(icoeff[i]));
    
    out.push_back(make_shared<RelDFHalfB>(array<shared_ptr<DFHalfDist>,2>{{real, imag}}, i, alpha));
  }
  return out;
}
