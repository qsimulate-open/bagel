//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd/orbital/asd_bfgs_denom.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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

#include <src/asd/orbital/asd_bfgs.h>

using namespace std;
using namespace bagel;

shared_ptr<const ASD_RotFile> ASD_BFGS::compute_denom(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> mcfock) const {
  auto out = make_shared<ASD_RotFile>(nclosed_, nact_, nvirt_, rasA_, rasB_);

  auto cfockd = make_shared<Matrix>(*cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * *rdm1);
  // ia part (4.7a)
  if (nvirt_ && nclosed_) {
    double* target = out->ptr_vc();
    for (int i = 0; i != nclosed_; ++i)
      for (int j = 0; j != nvirt_; ++j)
        *target++ = 4.0*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_)) - 4.0*(cfock->element(i,i)+afock->element(i,i));
  }
  // ra part (4.7b)
  if (nvirt_ && nact_) {
    double* target = out->ptr_va();
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nvirt_; ++j)
        *target++ = 2.0*rdm1->element(i,i)*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_))
                  - 2.0*cfockd->element(i,i) - 2.0*qxr->element(i+nclosed_, i);
  }
  // it part (4.7c)
  if (nclosed_ && nact_) {
    double* target = out->ptr_ca();
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nclosed_; ++j)
        *target++ = 4.0*(cfock->element(i+nclosed_, i+nclosed_)+afock->element(i+nclosed_, i+nclosed_) - cfock->element(j,j) - afock->element(j,j))
                  + 2.0*rdm1->element(i,i)*(cfock->element(j,j)+afock->element(j,j)) - 2.0*cfockd->element(i,i) - 2.0*qxr->element(i+nclosed_, i);
  }

  // tu part
  if (nact_) {
    double* target = out->ptr_aa();
    for (int i = 0; i != nactA_; ++i)  //A
      for (int j = nactA_; j != nact_; ++j)  //B
        *target++ = 2.0*(rdm1->element(i,i)*(cfock->element(j+nclosed_,j+nclosed_) + afock->element(j+nclosed_,j+nclosed_)) - mcfock->element(i,i)
                       + rdm1->element(j,j)*(cfock->element(i+nclosed_,i+nclosed_) + afock->element(i+nclosed_,i+nclosed_)) - mcfock->element(j,j)
                       - 2.0*rdm1->element(i,j)*cfock->element(i+nclosed_,j+nclosed_));
  }
  if (nactA_) {
    if (rasA_[0]) {//21
      double* target = out->ptr_aa21A();
      for (int i = 0; i != rasA_[0]; ++i) //RAS1(A)
        for (int j = rasA_[0]; j != rasA_[0]+rasA_[1]; ++j) //RAS2(A)
          *target++ = 2.0*(rdm1->element(i,i)*(cfock->element(j+nclosed_,j+nclosed_) + afock->element(j+nclosed_,j+nclosed_)) - mcfock->element(i,i)
                         + rdm1->element(j,j)*(cfock->element(i+nclosed_,i+nclosed_) + afock->element(i+nclosed_,i+nclosed_)) - mcfock->element(j,j)
                         - 2.0*rdm1->element(i,j)*cfock->element(i+nclosed_,j+nclosed_));
    }
    if (rasA_[0] && rasA_[2]) {//31
      double* target = out->ptr_aa31A();
      for (int i = 0; i != rasA_[0]; ++i) //RAS1(A)
        for (int j = rasA_[0]+rasA_[1]; j != nactA_; ++j)  //RAS3(A)
          *target++ = 2.0*(rdm1->element(i,i)*(cfock->element(j+nclosed_,j+nclosed_) + afock->element(j+nclosed_,j+nclosed_)) - mcfock->element(i,i)
                         + rdm1->element(j,j)*(cfock->element(i+nclosed_,i+nclosed_) + afock->element(i+nclosed_,i+nclosed_)) - mcfock->element(j,j)
                         - 2.0*rdm1->element(i,j)*cfock->element(i+nclosed_,j+nclosed_));
    }
    if (rasA_[2]) {//32
      double* target = out->ptr_aa32A();
      for (int i = rasA_[0]; i != rasA_[0]+rasA_[1]; ++i)  //RAS2(A)
        for (int j = rasA_[0]+rasA_[1]; j != nactA_; ++j) //RAS3(A)
          *target++ = 2.0*(rdm1->element(i,i)*(cfock->element(j+nclosed_,j+nclosed_) + afock->element(j+nclosed_,j+nclosed_)) - mcfock->element(i,i)
                         + rdm1->element(j,j)*(cfock->element(i+nclosed_,i+nclosed_) + afock->element(i+nclosed_,i+nclosed_)) - mcfock->element(j,j)
                         - 2.0*rdm1->element(i,j)*cfock->element(i+nclosed_,j+nclosed_));
    }
  }
  if (nactB_) {
    if (rasB_[0]) {//21
      double* target = out->ptr_aa21B();
      for (int i = nactA_; i != nactA_+rasB_[0]; ++i) //RAS1(B)
        for (int j = nactA_+rasB_[0]; j != nactA_+rasB_[0]+rasB_[1]; ++j) //RAS2(B)
          *target++ = 2.0*(rdm1->element(i,i)*(cfock->element(j+nclosed_,j+nclosed_) + afock->element(j+nclosed_,j+nclosed_)) - mcfock->element(i,i)
                         + rdm1->element(j,j)*(cfock->element(i+nclosed_,i+nclosed_) + afock->element(i+nclosed_,i+nclosed_)) - mcfock->element(j,j)
                         - 2.0*rdm1->element(i,j)*cfock->element(i+nclosed_,j+nclosed_));
    }
    if (rasB_[0] && rasB_[2]) {//31
      double* target = out->ptr_aa31B();
      for (int i = nactA_; i != nactA_+rasB_[0]; ++i) //RAS1(B)
        for (int j = nactA_+rasB_[0]+rasB_[1]; j != nact_; ++j) //RAS3(B)
          *target++ = 2.0*(rdm1->element(i,i)*(cfock->element(j+nclosed_,j+nclosed_) + afock->element(j+nclosed_,j+nclosed_)) - mcfock->element(i,i)
                         + rdm1->element(j,j)*(cfock->element(i+nclosed_,i+nclosed_) + afock->element(i+nclosed_,i+nclosed_)) - mcfock->element(j,j)
                         - 2.0*rdm1->element(i,j)*cfock->element(i+nclosed_,j+nclosed_));
    }
    if (rasB_[2]) {//32
      double* target = out->ptr_aa32B();
      for (int i = nactA_+rasB_[0]; i != nactA_+rasB_[0]+rasB_[1]; ++i) //RAS2(B)
        for (int j = nactA_+rasB_[0]+rasB_[1]; j != nact_; ++j) //RAS3(B)
          *target++ = 2.0*(rdm1->element(i,i)*(cfock->element(j+nclosed_,j+nclosed_) + afock->element(j+nclosed_,j+nclosed_)) - mcfock->element(i,i)
                         + rdm1->element(j,j)*(cfock->element(i+nclosed_,i+nclosed_) + afock->element(i+nclosed_,i+nclosed_)) - mcfock->element(j,j)
                         - 2.0*rdm1->element(i,j)*cfock->element(i+nclosed_,j+nclosed_));
    }

  }


  const double thresh = 1.0e-8;
  //examine approximate diagonal hessian elements
  for (int i = 0; i != out->size(); ++i) {
    if (out->data(i) < 0.0) {
      cout << setw(10) << setprecision(6) << out->data(i) << endl;
      throw runtime_error("Element of initial diagonal hessian < 0");
    }
    if (fabs(out->data(i)) < thresh) {
      out->data(i) = 1.0e10;
    }
  }
  return out;
}
