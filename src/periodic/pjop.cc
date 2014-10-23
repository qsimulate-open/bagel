//
// BAGEL - Parallel electron correlation program.
// Filename: pjop.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/periodic/pdfdist_ints.h>

using namespace std;
using namespace bagel;

shared_ptr<Matrix> PDFDist_ints::pcompute_Jop_from_coeff(shared_ptr<const VectorB> coeff) const {
  if (block_.size() != 2) throw logic_error("pcompute_Jop in periodic SCF needs block_.size() == 2");
  shared_ptr<Matrix> out = block_[0]->form_mat(coeff->slice(block_[0]->astart(), block_[0]->astart()+block_[0]->asize()));

  if (!serial_)
    out->allreduce();

  return out;
}


shared_ptr<VectorB> PDFDist_ints::pcompute_coeff(const shared_ptr<const Matrix> density) const {
  if (block_.size() != 2) throw logic_error("pcompute_Jop in periodic SCF needs block_.size() == 2");

  if (!data2_) throw logic_error("PDFDist_ints::pcompute_coeff was called without 2-index integrals");

  // get modified 3-index
  Matrix tmp = *projector_ * *data2_;
  contract(1.0, *block_[1], {3, 1, 2}, tmp, {0, 3}, 0.0, *block_[1], {0, 1, 2});
  // get 3-index coeff from the charged and chargeless part
  blas::ax_plus_y_n(1.0, coeffC_->data(), block_[1]->size(), block_[1]->data());

  auto out = make_shared<VectorB>(naux_);
  // contract with density
  shared_ptr<VectorB> coeff = block_[1]->form_vec(density);
  copy_n(coeff->data(), block_[1]->asize(), out->data()+block_[1]->astart());

  if (!serial_)
    out->allreduce();

  return out;
}


shared_ptr<Matrix> PDFDist_ints::pcompute_Jop(const shared_ptr<const Matrix> density) const {
  shared_ptr<const VectorB> coeff = pcompute_coeff(density);

  return pcompute_Jop_from_coeff(coeff);
}
