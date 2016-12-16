//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: qvec.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/multi/casscf/qvec.h>

using namespace std;
using namespace bagel;


Qvec::Qvec(const int n, const int m, shared_ptr<const Matrix> coeff, const size_t nclosed,
           shared_ptr<const DFHalfDist> half, shared_ptr<const RDM<2>> rdm)
 : Matrix(n,m) {

  assert(n == coeff->mdim());

  // J^{-1}(D|xy)
  // TODO : DFDistT needs to be modified to handle cases where number of nodes is larger than half->nocc() * cdata.mdim()
  shared_ptr<const DFFullDist> full;
  if (half->nocc() * coeff->mdim() > mpi__->size()) {
    full = half->apply_JJ()->compute_second_transform(coeff->slice(nclosed, nclosed+m));
  } else {
    full = half->compute_second_transform(coeff->slice(nclosed, nclosed+m))->apply_JJ();
  }


  // [D|tu] = (D|xy)Gamma_xy,tu
  shared_ptr<const DFFullDist> prdm = full->apply_2rdm(*rdm);

  // (r,u) = (rt|D)[D|tu]
  shared_ptr<const Matrix> tmp = half->form_2index(prdm, 1.0);

  // MO transformation of the first index
  *this = *coeff % *tmp;

}
