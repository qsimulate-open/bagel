//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcasscf_coeff.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/multi/zcasscf/zcasscf.h>

using namespace std;
using namespace bagel;


void ZCASSCF::kramers_adapt(shared_ptr<ZMatrix> o, const int nvirt) const {
  assert(o->ndim() == o->mdim() && (nclosed_ + nact_ + nvirt)*2 == o->ndim());

  auto kramers_adapt_block = [this,&nvirt,&o](const int jst, const int ist, const int joff, const int ioff) {
    for (int i = 0; i != ist; ++i)
      for (int j = 0; j != jst; ++j) {
        o->element(joff+j, ioff+i) = 0.5*(o->element(joff+j, ioff+i) + conj(o->element(joff+j+jst, ioff+i+ist)));
        o->element(joff+j+jst,ioff+i+ist) = conj(o->element(joff+j,ioff+i));

        o->element(joff+jst+j, ioff+i) = 0.5*(o->element(joff+jst+j, ioff+i) - conj(o->element(joff+j, ioff+ist+i)));
        o->element(joff+j, ioff+ist+i) = - conj(o->element(joff+jst+j, ioff+i));
      }
  };

  const array<int,3> stride{{nclosed_, nact_, nvirt}};
  const array<int,3> offset{{0, 2*nclosed_, 2*nocc_}};
  for (int ii = 0; ii != 3; ++ii)
    for (int jj = 0; jj != 3; ++jj)
      kramers_adapt_block(stride[jj], stride[ii], offset[jj], offset[ii]);
}


void ZCASSCF::kramers_adapt(shared_ptr<ZRotFile> o, const int nclosed, const int nact, const int nvirt) {
  for (int i = 0; i != nclosed; ++i) {
    for (int j = 0; j != nvirt; ++j) {
      o->ele_vc(j, i) = (o->ele_vc(j, i) + conj(o->ele_vc(j+nvirt, i+nclosed))) * 0.5;
      o->ele_vc(j+nvirt, i+nclosed) = conj(o->ele_vc(j, i));

      o->ele_vc(j+nvirt, i) = (o->ele_vc(j+nvirt, i) - conj(o->ele_vc(j, i+nclosed))) * 0.5;
      o->ele_vc(j, i+nclosed) = - conj(o->ele_vc(j+nvirt, i));
    }
  }
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nvirt; ++j) {
      o->ele_va(j, i) = (o->ele_va(j, i) + conj(o->ele_va(j+nvirt, i+nact))) * 0.5;
      o->ele_va(j+nvirt, i+nact) = conj(o->ele_va(j, i));

      o->ele_va(j+nvirt, i) = (o->ele_va(j+nvirt, i) - conj(o->ele_va(j, i+nact))) * 0.5;
      o->ele_va(j, i+nact) = - conj(o->ele_va(j+nvirt, i));
    }
  }
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nclosed; ++j) {
      o->ele_ca(j, i) = (o->ele_ca(j, i) + conj(o->ele_ca(j+nclosed, i+nact))) * 0.5;
      o->ele_ca(j+nclosed, i+nact) = conj(o->ele_ca(j, i));

      o->ele_ca(j+nclosed, i) = (o->ele_ca(j+nclosed, i) - conj(o->ele_ca(j, i+nact))) * 0.5;
      o->ele_ca(j, i+nact) = - conj(o->ele_ca(j+nclosed, i));
    }
  }
}


shared_ptr<ZCoeff_Kramers> ZCASSCF::nonrel_to_relcoeff(shared_ptr<const Matrix> nr_coeff) const {
  // constructs a relativistic coefficient for electronic components from a non-rel coefficient
  const int n = nr_coeff->ndim();
  const int m = nr_coeff->mdim();
  assert(nvirt_ - nneg_/2 == nvirtnr_);

  // compute T^(-1/2)
  shared_ptr<ZMatrix> t12 = overlap_->get_submatrix(n*2, n*2, n, n);
  t12 = t12->tildex(thresh_overlap_);

  // compute S^(1/2)
  shared_ptr<ZMatrix> shalf = overlap_->get_submatrix(0, 0, n, n);
  shalf = shalf->tildex(thresh_overlap_);

  if (t12->mdim() != shalf->mdim())
    throw runtime_error("Different linear dependency for the overlap and kinetic matrices in conversion to relativistic coefficients.");

  // compute positronic orbital coefficients
  auto tcoeff = make_shared<ZMatrix>(n, m);
  tcoeff->add_real_block(1.0, 0, 0, n, m, *nr_coeff);
  *tcoeff = *t12 * (*shalf % *tcoeff);

  // build output coefficient matrix
  auto out = make_shared<ZCoeff_Kramers>(4*n, nr_coeff->localized(), nclosed_, nact_, nvirtnr_, nneg_);
  assert(out->mdim() == 4*tcoeff->mdim());
  out->copy_real_block(1.0, 0, 0, n, m, *nr_coeff);
  out->copy_real_block(1.0, n, 2*m, n, m, *nr_coeff);
  out->copy_block(2*n, m, n, m, *tcoeff);
  out->copy_block(3*n, 3*m, n, m, *tcoeff);
  return out;
}


// Eliminates the positronic entries for the given rot file
void ZCASSCF::zero_positronic_elements(shared_ptr<ZRotFile> rot) {
  int nr_nvirt = nvirt_ - nneg_/2;
  for (int i = 0; i != nclosed_*2; ++i) {
    for (int j = 0; j != nneg_/2; ++j) {
      rot->ele_vc(j + nr_nvirt, i) = 0.0;
      rot->ele_vc(j + nr_nvirt + nvirt_, i) = 0.0;
    }
  }
  for (int i = 0; i != nact_*2; ++i) {
    for (int j = 0; j != nneg_/2; ++j) {
      rot->ele_va(j + nr_nvirt, i) = 0.0;
      rot->ele_va(j + nr_nvirt + nvirt_, i) = 0.0;
    }
  }
}
