//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf_coeff.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/multi/zcasscf/zcasscf.h>

using namespace std;
using namespace bagel;


void ZCASSCF::kramers_adapt(shared_ptr<ZMatrix> o, const int nvirt) const {
  // function to enforce time-reversal symmetry
  //    for a complex matrix o, that is SYMMETRIC under time reversal

  auto kramers_adapt_block = [this](shared_ptr<ZMatrix> o, unsigned int tfac, int nq1, int nq2, int joff, int ioff, int nvirt) {
    // function to enforce time-reversal symmetry for a given block of a complex matrix o.
    // tfac                 symmetry factor under time-reversal (symmetric : t=1, antisymm : t=-1)
    // nq1, nq2             nclosed_, nact_, nvirt
    // ioff, joff           column and row offsets, respectively
    assert( o->ndim() == o->mdim() && (nclosed_ + nact_ + nvirt)*2 == o->ndim() );
    assert( tfac == 1 || tfac == -1);
    const double t = tfac == 1 ? 1.0 : -1.0;
    for (int i = 0; i != nq2; ++i) {
      for (int j = 0; j != nq1; ++j) {
        // Diagonal contributions : "A" matrices in notation of T. Suae thesis
        o->element(joff+j, ioff+i) = ( o->element(joff+j, ioff+i) + conj(o->element(joff+j+nq1, ioff+i+nq2)) ) * 0.5;
        o->element(joff+j+nq1,ioff+i+nq2) = t * conj(o->element(joff+j,ioff+i));

        // Diagonal contributions : "B" matrices in notation of T. Suae thesis
        o->element(joff+nq1+j, ioff+i) = ( o->element(joff+nq1+j, ioff+i) - conj(o->element(joff+j, ioff+nq2+i)) ) * 0.5;
        o->element(joff+j, ioff+nq2+i) = -t * conj(o->element(joff+nq1+j,ioff+i));
      }
    }
  };

  assert(o->ndim() == o->mdim() && (nclosed_ + nact_ + nvirt)*2 == o->ndim());
  const array<int,3> a0 {{nclosed_, nact_, nvirt}};
  const array<int,3> a1 {{0, 2*nclosed_, 2*nocc_}};
  for (int ii = 0; ii !=3; ++ii) {
    for (int jj = 0; jj !=3; ++jj) {
      kramers_adapt_block(o,1,a0[jj],a0[ii],a1[jj],a1[ii],nvirt);
    }
  }
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


// TODO rewrite this so it can tolerate linear dependency
shared_ptr<RelCoeff_Kramers> ZCASSCF::nonrel_to_relcoeff(shared_ptr<const Matrix> nr_coeff) const {
  // constructs a relativistic coefficient for electronic components from a non-rel coefficient
  const int n = nr_coeff->ndim();
  const int m = nr_coeff->mdim();
  assert(nvirt_ - nneg_/2 == nvirtnr_);

  // compute T^(-1/2)
  shared_ptr<ZMatrix> t12 = overlap_->get_submatrix(n*2, n*2, n, n);
  assert(t12->inverse_half(1.0e-10));

  // compute S^(1/2)
  shared_ptr<ZMatrix> shalf = overlap_->get_submatrix(0, 0, n, n);
  VectorB eig2(shalf->mdim());
  shalf->diagonalize(eig2);
  for (int k = 0; k != shalf->mdim(); ++k) {
    if (real(eig2[k]) >= 0.0)
      blas::scale_n(sqrt(sqrt(eig2(k))), shalf->element_ptr(0, k), shalf->ndim());
  }
  *shalf = *shalf ^ *shalf;

  // compute positronic orbital coefficients
  auto tcoeff = make_shared<ZMatrix>(n, m);
  tcoeff->add_real_block(1.0, 0, 0, n, n, *nr_coeff);
  *tcoeff = *t12 * *shalf * *tcoeff;

  // build output coefficient matrix
  auto out = make_shared<RelCoeff_Kramers>(4*n, nr_coeff->localized(), nclosed_, nact_, nvirtnr_, nneg_);
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
      rot->ele_vc(j + nr_nvirt, i) =  complex<double> (0.0,0.0);
      rot->ele_vc(j + nr_nvirt + nvirt_, i) =  complex<double> (0.0,0.0);
    }
  }
  for (int i = 0; i != nact_*2; ++i) {
    for (int j = 0; j != nneg_/2; ++j) {
      rot->ele_va(j + nr_nvirt, i) =  complex<double> (0.0,0.0);
      rot->ele_va(j + nr_nvirt + nvirt_, i) =  complex<double> (0.0,0.0);
    }
  }
}
