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

#include <src/periodic/pdfdist.h>

using namespace std;
using namespace bagel;

shared_ptr<PData> PDFDist::pcompute_Jop_from_coeff(shared_ptr<const VectorB> coeff) const {

  auto out = make_shared<PData>(nbasis_, ncell());
  for (int i = 0; i != ncell(); ++i) {
    // lattice sum with NAI
    auto jmat = make_shared<Matrix>(nbasis_, nbasis_);
    jmat->zero();
    for (int j = 0; j != ncell(); ++j) {
      shared_ptr<DFBlock> data3 = dfdist_[i]->data3_in_cell(j);
      // contract with coeff
      shared_ptr<Matrix> tmp = data3->form_mat(coeff->slice(data3->astart(), data3->astart() + data3->asize()));
      *jmat += *tmp;
      // add NAI contribution
      *jmat += 2.0 * *dfdist_[i]->nai_in_cell(j);
    }
    (*out)[i] = make_shared<ZMatrix>(*jmat , complex<double>(1.0, 0.0));
  }

  if (!serial_)
    out->allreduce();

  return out;
}


shared_ptr<VectorB> PDFDist::pcompute_coeff(const shared_ptr<const PData> density) const {

  if (!data2_) throw logic_error("PDFDist::pcompute_coeff was called without 2-index integrals");

  auto coeff1 = make_shared<VectorB>(naux_);
  auto coeff2 = make_shared<VectorB>(naux_);

  for (int i = 0; i != ncell(); ++i) {
    // get charged coeff by contracting with density
    auto tmp1 = make_shared<VectorB>(naux_);
    shared_ptr<btas::Tensor3<double>> coeffC = dfdist_[i]->coeffC();
    contract(1.0, group(*coeffC, 1, 3), {0, 1}, group(*(density->pdata(i)->get_real_part()), 0, 2), {1}, 0.0, *tmp1, {0});
    *coeff1 += *tmp1;

    // contract 3-index with density
    shared_ptr<DFBlock> data3 = dfdist_[i]->block(0);
    auto tmp2 = make_shared<VectorB>(naux_);
    contract(1.0, group(*data3, 1, 3), {0, 1}, group(*(density->pdata(i)->get_real_part()), 0, 2), {1}, 0.0, *tmp2, {0});
    *coeff2 += *tmp2;
  }

  // contract coeff1 with 2-index
  auto tmp = make_shared<VectorB>(naux_);
  contract(1.0, *eta_, {0, 1}, *coeff1, {1}, 0.0, *tmp, {0});
  *coeff2 -= *tmp;

  // get chargeless coeff
  *coeff2 = *data2_ * *coeff2;

  // get coeff from the charged and chargeless part
  auto out = make_shared<VectorB>(naux_);
  *out = *coeff1 + *coeff2;

  if (!serial_)
    out->allreduce();

  return out;
}


shared_ptr<PData> PDFDist::pcompute_Jop(const shared_ptr<const PData> density) const {
  shared_ptr<const VectorB> coeff = pcompute_coeff(density);

  return pcompute_Jop_from_coeff(coeff);
}
