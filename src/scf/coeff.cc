//
// Newint - Parallel electron correlation program.
// Filename: coeff.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#include <src/scf/coeff.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>


Coeff::Coeff(const Matrix1e& inp) : Matrix1e(inp.geom()) {

  const int ndim_ = inp.ndim();
  const int mdim_ = inp.mdim();
  dcopy_(nbasis_*nbasis_, inp.data(), 1, data(), 1); 

}


Coeff::~Coeff() {

}


Matrix1e Coeff::form_density_rhf() const {
  const int nocc = geom_->nocc() / 2;
  assert(geom_->nocc() % 2 == 0);

  Matrix1e out(geom_);
  double* out_data = out.data();

  dgemm_("N", "T", nbasis_, nbasis_, nocc, 1.0, data(), nbasis_, data(), nbasis_, 0.0, out_data, nbasis_); 

  return out;
}
 

Matrix1e Coeff::form_core_density_rhf() const {
  const int nocc = geom_->nfrc() / 2;
  assert(geom_->nfrc() % 2 == 0);

  Matrix1e out(geom_);
  double* out_data = out.data();

  dgemm_("N", "T", nbasis_, nbasis_, nocc, 1.0, data(), nbasis_, data(), nbasis_, 0.0, out_data, nbasis_); 

  return out;
}
 

