//
// BAGEL - Parallel electron correlation program.
// Filename: qvec.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/casscf/qvec.h>

using namespace std;
using namespace bagel;


Qvec::Qvec(const int n, const int m, shared_ptr<const Matrix> coeff, const size_t nclosed, shared_ptr<const FCI> fci, shared_ptr<const RDM<2>> rdm)
 : Matrix(n,m) {

  assert(n == coeff->mdim());

  // one index transformed integrals (active)
  shared_ptr<const DFHalfDist> half = fci->jop()->mo2e_1ext();

  // J^{-1}(D|xy)
  shared_ptr<const DFFullDist> full = half->compute_second_transform(coeff->slice(nclosed, nclosed+m))->apply_JJ();

  // [D|tu] = (D|xy)Gamma_xy,tu
  shared_ptr<const DFFullDist> prdm = full->apply_2rdm(rdm->data());

  // (r,u) = (rt|D)[D|tu]
  shared_ptr<const Matrix> tmp = half->form_2index(prdm, 1.0);

  // MO transformation of the first index
  *this = *coeff % *tmp;

}
