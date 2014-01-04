//
// BAGEL - Parallel electron correlation program.
// Filename: sohcore.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/scf/sohcore.h>

using namespace std;
using namespace bagel;

SOHcore::SOHcore(const shared_ptr<const Geometry> geom, const shared_ptr<const SOHcore_base> h)
            : Matrix(2 * geom->nbasis(), 2 * geom->nbasis()), geom_(geom), hcore_(h) {
  form_sohcore();
}

// TODO: so1, so2, ecp
shared_ptr<const Matrix> SOHcore::so1() {
  return make_shared<Matrix>(hcore_->ndim(), hcore_->mdim());
}

shared_ptr<const Matrix> SOHcore::so2() {
  return make_shared<Matrix>(hcore_->ndim(), hcore_->mdim());
}

shared_ptr<const Matrix> SOHcore::ecp() {
  shared_ptr<Matrix> out = make_shared<Matrix>(hcore_->ndim(), hcore_->mdim());
  *out += *hcore_;
  return out;
}

void SOHcore::form_sohcore() {
  shared_ptr<const Matrix> ecpmat = ecp();

  add_block(1.0, 0, 0, hcore_->ndim(), hcore_->mdim(), ecpmat);
  add_block(1.0, hcore_->ndim(), hcore_->mdim(), hcore_->ndim(), hcore_->mdim(), ecpmat);

  add_block(1.0, 0, hcore_->mdim(), hcore_->ndim(), hcore_->mdim(), so1());
  add_block(1.0, hcore_->ndim(), 0, hcore_->ndim(), hcore_->mdim(), so2());
}
