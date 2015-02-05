//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbopt/ras/bfgs_denom.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim: <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <src/asd/orbopt/ras/bfgs.h>

using namespace std;
using namespace bagel;

shared_ptr<const ASDRASRotFile> ASDRASBFGS::compute_denom(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> mcfock) const {
  auto out = make_shared<ASDRASRotFile>(nclosed_, nact_, nvirt_, nactA_, nactB_, rasA_, rasB_);
//const double tiny = 1.0e-15;

  shared_ptr<Matrix> cfockd;
  if (nact_) {
    cfockd = make_shared<Matrix>(*cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * *rdm1);
    //TODO check symmetric??
  }

  // ia part (4.7a)
  if (nvirt_ && nclosed_) {
    double* target = out->ptr_vc();
    for (int i = 0; i != nclosed_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        *target++ = 4.0*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_)) - 4.0*(cfock->element(i,i)+afock->element(i,i));
      }
    }
  }
  // ra part (4.7b)
  if (nvirt_ && nact_) {
    double* target = out->ptr_va();
    for (int i = 0; i != nact_; ++i) {
    //if (occup_[i] < tiny) continue;
      for (int j = 0; j != nvirt_; ++j) {
    //  *target++ = 2.0*occup_[i]*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_))
    //            - 2.0*occup_[i]*cfock->element(i+nclosed_, i+nclosed_) - 2.0*qxr->element(i+nclosed_, i);
        *target++ = 2.0*rdm1->element(i,i)*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_))
                  - 2.0*cfockd->element(i,i) - 2.0*qxr->element(i+nclosed_, i);
      }
    }
  }
  // it part (4.7c)
  if (nclosed_ && nact_) {
    double* target = out->ptr_ca();
    for (int i = 0; i != nact_; ++i) {
//    if (occup_[i] < tiny) continue;
      for (int j = 0; j != nclosed_; ++j) {
//      *target++ = 4.0*(cfock->element(i+nclosed_, i+nclosed_)+afock->element(i+nclosed_, i+nclosed_) - cfock->element(j,j) - afock->element(j,j))
//                + 2.0*occup_[i]*(cfock->element(j,j)+afock->element(j,j)) - 2.0*occup_[i]*cfock->element(i+nclosed_, i+nclosed_) - 2.0*qxr->element(i+nclosed_, i);
        *target++ = 4.0*(cfock->element(i+nclosed_, i+nclosed_)+afock->element(i+nclosed_, i+nclosed_) - cfock->element(j,j) - afock->element(j,j))
                  + 2.0*rdm1->element(i,i)*(cfock->element(j,j)+afock->element(j,j)) - 2.0*cfockd->element(i,i) - 2.0*qxr->element(i+nclosed_, i);
      }
    }
  }
  // tu part
  if (nact_) {
    double* target = out->ptr_aa();
    for (int i = nactA_; i != nact_; ++i) { //B
      for (int j = 0; j != nactA_; ++j) { //A
      //*target++ = 0.01;
        *target++ = 2.0*( 
                    - mcfock->element(j,j) - mcfock->element(i,i) 
                    + rdm1->element(j,j)*cfock->element(i+nclosed_,i+nclosed_) + rdm1->element(i,i)*cfock->element(j+nclosed_,j+nclosed_) 
                    - rdm1->element(i,j)*cfock->element(j+nclosed_,i+nclosed_) - rdm1->element(j,i)*cfock->element(i+nclosed_,j+nclosed_)
                  //+ 1.0*(afock->element(i+nclosed_,i+nclosed_) + afock->element(j+nclosed_,j+nclosed_))
                    + 2.0*(afock->element(i+nclosed_,i+nclosed_) + afock->element(j+nclosed_,j+nclosed_))
                    );
      }
    }
  }

  const double thresh = 1.0e-8;
  for (int i = 0; i != out->size(); ++i)
    if (fabs(out->data(i)) < thresh) {
      out->data(i) = 1.0e10;
    }
  return out;
}


