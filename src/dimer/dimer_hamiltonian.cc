//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_hamiltonian.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <tuple>

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <vector>

#include <src/scf/geometry.h>
#include <src/scf/coeff.h>
#include <src/scf/fock.h>
#include <src/scf/scf.h>
#include <src/fci/harrison.h>
#include <src/fci/space.h>
#include <src/dimer/dimer.h>
#include <src/util/matrix.h>

using namespace std;
using namespace bagel;

/************************************************************************************
*  Dimer::Dimer(shared_ptr<Geometry> A, vector<double> displacement)                *
*                                                                                   *
************************************************************************************/
void Dimer::hamiltonian() {
  const int nstatesA = ccvecs_.first->ij();
  const int nstatesB = ccvecs_.second->ij();

  nstates_ = make_pair(nstatesA, nstatesB);

  dimerstates_ = nstatesA*nstatesB;

  hamiltonian_ = shared_ptr<Matrix>(new Matrix(dimerstates_, dimerstates_));

  // create jop_ object (which effectively includes all closed-closed interactions)
  jop_ = shared_ptr<MOFile>(new Jop(sref_, ncore_, ncore_ + norb_, scoeff_, string("HZ")));
  // For (temporary?) coding simplicity, I will make jopA_ and jopB_, since with those I can directly use already written functions
  jopA_ = shared_ptr<MOFile>(new Jop(refs_.first, ncoreA_, ncoreA_ + norbA_, coeffs_.front(), string("HZ")));
  if (symmetric_) jopB_ = jopA_;
  else jopB_ = shared_ptr<MOFile>(new Jop(refs_.second, ncoreB_, ncoreB_ + norbB_, coeffs_.front(), string("HZ")));

  // compute close-close
  cout << "  o Computing closed-closed interactions" << endl;
  hamiltonian_ += compute_closeclose();

  // closed-active interactions
  cout << "  o Computing closed-active interactions" << endl;
  hamiltonian_ += compute_closeactive();

  // active-active interactions
  cout << "  o Computing active-active interactions" << endl;
  hamiltonian_ += compute_activeactive();
}

shared_ptr<Matrix> Dimer::compute_closeclose() {
  double core = sgeom_->nuclear_repulsion() + jop_->core_energy();

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));
  out->add_diag(core);

  return out;
}

shared_ptr<Matrix> Dimer::compute_closeactive() {
  // First, compute 4-index quantities from DF
  const int ncore = ncore_.first + ncore_.second;
  const int nact  = nact_.first + nact_.second;

  const int nactA = nact_.first;
  const int nactB = nact_.second;

  double *cdata_core = scoeff_->data();
  double *cdata_act  = scoeff_->data() + nact*dimerbasis_;

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  const double *outdata = out->data();

  // Half transforms
  shared_ptr<DF_Half> half_ir = sgeom_->df()->compute_half_transform(cdata_core, ncore);
  shared_ptr<DF_Half> half_ar = sgeom_->df()->compute_half_transform(cdata_act, nact);

  // Full transforms
  shared_ptr<DF_Full> buf_ii = half_ir->compute_second_transform(cdata_core, ncore)->apply_J();
  shared_ptr<DF_Full> buf_ia = half_ir->compute_second_transform(cdata_act, nact)->apply_J();
  shared_ptr<DF_Full> buf_aa = half_ar->compute_second_transform(cdata_act,nact)->apply_J();

  // Form 4-index quantities of interest
  unique_ptr<double[]> Kop_aaii = buf_aa->form_4index(buf_ii);
  unique_ptr<double[]> Kop_aiia = buf_ia->form_4index();

  // I think I will try to preserve the aaii and aiia ordering of the 4-index quantities and just alter the above sections to make sure it works
  // Tab = \sum_i [ 2*(ab|ii) - (ia|ib) ]
  unique_ptr<double[]> Tab(new double[nactA*nactA + nactB*nactB]);
  
  double *tdata = Tab.get();
  for(int a = 0; a < nactA; ++a) {
    for(int b = 0; b < nactA; ++b) {
      double temp = 0.0;
      for(int i = 0; i < ncore; ++i) {
        temp += 2.0*Kop_aaii[act<0>(a) + nact*act<0>(b) + nact*nact*i + nact*nact*ncore*i] - Kop_aiia[act<0>(a) + nact*i + nact*ncore*i + nact*ncore*nact*act<0>(b)];
      }
      *tdata++ = temp;
    }

  for(int a = 0; a < nactB; ++a) {
    for(int b = 0; b < nactB; ++b) {
      double temp = 0.0;
      for(int i = 0; i < ncore; ++i) {
        temp += 2.0*Kop_aaii[act<1>(a) + nact*act<1>(b) + nact*nact*i + nact*nact*ncore*i] - Kop_aiia[act<1>(a) + nact*i + nact*ncore*i + nact*ncore*nact*act<1>(b)];
      }
      *tdata++ = temp;
    }
  }

  shared_ptr<Determinants> neut_det = space_->finddet(0,0);

  const int lena = neut_det->lena();
  const int lenb = neut_det->lenb();

  // Okay, I have Tab. Now for the actual matrix elements.
  for(int stateA = 0; stateA < nstatesA_; ++stateA) {
    const double *ccA = ccvecs_.first->data(stateA)->data();
    for(int stateB = 0; stateB < nstatesB_; ++stateB) {
      const double *ccB = ccvecs_.second->data(stateB)->data();
      const int stateAB = stateA + stateB*nstatesA_;

      // ab belongs to unit A
      for(int stateAp = 0; stateAp < nstatesA_; ++stateAp) {
        const int stateABp = stateAp + stateB*nstatesA_;
        double element = 0;

        const double *ccAp = ccvecs_.first->data(stateAp)->data();

        const int abA = nactA*nactA;
        for(int ab = 0; ab < abA; ++ab) {
          const double t = tdata[ab];

          // alpha-alpha
          for(auto& iter : neut_det->phia(ab)) {
            double sign = static_cast<double>(get<1>(iter));
            element += t*sign*ddot_(lenb, ccA + get<2>(iter)*lenb, 1, ccAp+get<0>(iter)*lenb, 1);
          }

          // beta-beta
          for(auto& iter : neut_det->phib(ab)) {
            double sign = static_cast<double>(get<1>(iter));
            element += t*sign*ddot_(lena, ccA + get<2>(iter), lenb, ccAp + get<0>(iter), lenb);
          }
        }

        outdata[stateAB + ndstates_*stateABp] += element;
      }

      // ab belongs to unit B
      for(int stateBp = 0; stateBp < nstatesB_; ++stateBp) {
        const int state ABp = stateA + stateBp*nstatesA_;
        double element = 0;

        const double *ccBp = ccvecs_.second->data(stateBp)->data();

        const int abB = nactB*nactB;
        for(int ab = 0; ab < abB; ++ab) {
          const double t = tdata[ab + nactA*nactA];

          // alpha-alpha
          for(auto iter = neut_det->phia(ab)->begin(); iter != neut_det->phia(ab)->end(); ++iter) {
            double sign = static_cast<double>(get<1>(*iter));
            element += t*sign*ddot_(lenb, ccB + get<2>(*iter)*lenb, 1, ccBp+get<0>(*iter)*lenb, 1);
          }

          // beta-beta
          for(auto iter = neut_det->phib(ab)->begin(); iter != neut_det->phib(ab)->end(): ++iter) {
            double sign = static_cast<double>(get<1>(*iter));
            element += t*sign*ddot_(lena, ccB + get<2>(*iter), lenb, ccBp + get<0>(*iter), lenb);
          }
        }

        outdata[stateAB + ndstates_*stateABp] += element;
      }
    }
  }

  //Done?
  return out;
}

shared_ptr<Matrix> Dimer::compute_intra_activeactive() {
  shared_ptr<Determinants> neut_det = space_->finddet(0,0);
  const int lenab = neut_det->lena() * neut_det->lenb();

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));
  const double *odata = out->data();

  // first do H^{AA}_{AA} case
  shared_ptr<Dvec> sigmavecAA = form_sigma(ccvecs_.first, jops_.first);

  for(int stateA = 0; stateA < nstatesA_; ++stateA) {
    for(int stateB = 0; stateB < nstatesB_; ++stateB) {
      const int stateAB = stateA + stateB*nstatesA_;
      const double *sdataA = sigmavecAA->data(stateA)->data();
      for(int stateAp = 0; stateAp < nstatesA_; ++stateAp) {
        const int stateABp = stateAp + stateB*nstatesA_;
        const double *cdataAp = ccvecs_.first->data(stateAp)->data();
        double element = ddot_(lenab, sdataA, 1, cdataAp, 1);
        
        odata[stateAB + ndstates_*stateABp] += element;
      }
    }
  }
  
  // now do H^{BB}_{BB} case
  shared_ptr<Dvec> sigmavecBB;
  sigmavecBB = symmetric_ ? sigmavecAA : form_sigma(ccvecs_.second, jops_.second);

  for(int stateA = 0; stateA < nstatesA_; ++stateA) {
    for(int stateB = 0; stateB < nstatesB_; ++stateB) {
      const int stateAB = stateA + stateB*nstatesA_;
      const double *sdataB = sigmavecBB->data(stateB)->data();
      for(int stateBp = 0; stateBp < nstatesB_; ++stateBp) {
        const int stateABp = stateA + stateBp*nstatesA_;
        const double *cdataBp = ccvecs_.second->data(stateBp)->data();
        double element = ddot_(lenab, sdataB, 1, cdataBp, 1);
        
        odata[stateAB + ndstates_*stateABp] += element;
      }
    }
  }

  return out;
}


//And now for the hard part
shared_ptr<Matrix> Dimer::compute_inter_activeactive() {
  // ijA and ijB need to match the compression of the det
  const int ijA = nact_.first*nact_.first; 
  const int ijB = nact_.second*nact_.second; 

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  shared_ptr<Matrix> E_ac = form_EFmatrices(ccvecs_.first, nstates_.first, ijA);

  // are the monomers identical?
  shared_ptr<Matrix> F_bd;
  if(symmetric_) F_bd = E_ac;
  else F_bd = form_EFmatrices(ccvecs_.second, nstates_.second, ijB);

  // build <ab|cd> - <ab|dc> matrix
  shared_ptr<Matrix> H_abcd = form_Hmatrix(jop_);

  shared_ptr<Matrix> out(new Matrix( (*E_ac) * (*H_abcd) ^ (*F_bd) ));

  // reorder, currently (AA',BB'), want (AB, A'B')
  shared_ptr<Matrix> out2(new Matrix(out));

  double* odata = out->data();
  double* odata2 = out2->data();
  int index = 0;

  const int nstatesA = nstates_.first;
  const int nstatesB = nstates_.second;
  for(int A = 0; A < nstatesA; ++A) {
    for(int B = 0; B < nstatesB; ++B) {
      for(int Ap = 0; Ap < nstatesA; ++Ap) {
        for(int Bp = 0; Bp < nstatesB; ++Bp, ++odata2) {
          *odata2 = odata[A + nstatesA*Ap + (B + nstatesB*Bp)*nstatesA*nstatesA];
        }
      }
    }
  }

  return out2;
}

shared_ptr<Matrix> form_EFmatrices(shared_ptr<const Dvec> > ccvec, const int ij, const int nstates) {
  shared_ptr<Matrix> out(new Matrix(nstates*nstates, ij));

  double *edata = out.data();

  shared_ptr<NewDvec> c(new NewDvec(ccvec->det(), ij));

  const int la = ccvec->det()->lena();
  const int lb = ccvec->det()->lenb();

  for(int state = 0; state < nstates; ++state) {
    const double* source_base = ccvec->data(state)->data();

    c->zero();

    for(int ac = 0; ac < ij; ++ac) {
      double *target_base = c->data(ac)->data();

      // alpha-alpha
      for(auto& iter : ccvec->det()->phia(ac) ) {
        const double sign = static_cast<double>(get<1>(iter));
        const double* target_array = target_base + get<2>(iter)*lb;
        daxpy_(lb, sign, source_base + get<0>(iter)*lb, 1, target_array, 1);
      }

      // beta-beta
      for(auto& iter : ccvec->det()->phib(ac)) {
        const double sign = static_cast<double>(get<1>(iter));
        const double* target_array = target_base + get<2>(iter);
        daxpy_(la, sign, source_base + get<0>(iter), lb, target_array, lb);
      }
    }

    // | C > ^A_ac is done
    const int lab = la*lb;
    for(int statep = 0; statep < nstates; ++statep) {
      const double *adata = ccvec->data(state)->data();
      for(int ac = 0; ac < ij; ++ac) {
        const double *cdata = c->data(ac)->data();
        *edata++ = ddot_(lab, adata, 1, cdata, 1);
      }
    }
  }
}


shared_ptr<Matrix> form_Hmatrix(shared_ptr<MOFile> jop_, const int ijA, const int ijB) {
  shared_ptr<Matrix> out(new Matrix(ijA, ijB));
  double* odata = out->data();
  
  for(int b = 0; b < nbasis_.second; ++b) {
    for(int d = 0; d < nbasis_.second; ++d) {
      for(int a = 0; a < nbasis_.first; ++a) {
        for(int c = 0; c < nbasis_.first; ++c) {
          *odata++ = jop_->mo2e(act<0>(a), act<1>(b), act<0>(c), act<1>(d)) - jop_->mo2e(act<0>(a), act<1>(b), act<1>(d), act<0>(c))
        }
      }
    }
  }

  return out;
}
