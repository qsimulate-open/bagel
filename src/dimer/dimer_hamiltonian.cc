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
*  General note: For the moment, everything is written so it will be easiest to     *
* debug, using as much code as possible from other functions already written. This  *
* will eventually be fixed.                                                         *
************************************************************************************/

void Dimer::hamiltonian() {
  const int nstatesA = ccvecs_.first->ij();
  const int nstatesB = ccvecs_.second->ij();

  const int ncore = ncore_.first + ncore_.second;
  const int nact  = nact_.first + nact_.second;

  nstates_ = make_pair(nstatesA, nstatesB);

  dimerstates_ = nstatesA*nstatesB;

  hamiltonian_ = shared_ptr<Matrix>(new Matrix(dimerstates_, dimerstates_));

  // create jop_ object (which effectively includes all closed-closed interactions)
  jop_ = shared_ptr<MOFile>(new Jop(sref_, ncore, ncore + nact, scoeff_, string("HZ")));
  // For (temporary?) coding simplicity, I will make jopA_ and jopB_, since with those I can directly use already written functions
  jops_.first = shared_ptr<MOFile>(new Jop(refs_.first, ncore_.first, ncore_.first + nact_.first, coeffs_.first, string("HZ")));
  if (symmetric_) jops_.second = jops_.first;
  else jops_.second = shared_ptr<MOFile>(new Jop(refs_.second, ncore_.second, ncore_.second + nact_.second, coeffs_.first, string("HZ")));

  // compute close-close
  cout << "  o Computing closed-closed interactions" << endl;
  *hamiltonian_ += *compute_closeclose();

  // closed-active interactions
  cout << "  o Computing closed-active interactions" << endl;
  *hamiltonian_ += *compute_closeactive();

  // intramolecular active-active interactions
  cout << "  o Computing intramolecular active-active interactions" << endl;
  *hamiltonian_ += *compute_intra_activeactive();

  // intermolecular active-active interactions
  cout << "  o Computing intermolecular active-active interactions" << endl;
  *hamiltonian_ += *compute_inter_activeactive();
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

  const int nstatesA = nstates_.first;
  const int nstatesB = nstates_.second;

  double *cdata_core = scoeff_->data();
  double *cdata_act  = scoeff_->data() + nact*dimerbasis_;

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  double *outdata = out->data();

  const int ijA = (nactA*(nactA+1))/2;
  const int ijB = (nactB*(nactB+1))/2;

  unique_ptr<double[]> Hab(new double[ijA + ijB]);
  double *hdata = Hab.get();
  // Sort 1e data so only the intramolecular terms are included  
  for (int i = 0; i < nactA; ++i) {
    for (int j = 0; j <= i; ++j) {
      *hdata++ = jop_->mo1e(j,i);
    }
  }

  for (int i =0; i < nactB; ++i) {
    for (int j = 0; j <= i; ++j) {
      *hdata++ = jop_->mo1e(j+nactA,i+nactA);
    }
  }

  shared_ptr<Determinants> neut_det = space_->finddet(0,0);

  const int lena = neut_det->lena();
  const int lenb = neut_det->lenb();
  const int lenab = lena*lenb;

  shared_ptr<Dvec> sigmavecA = form_sigma_1e(ccvecs_.first, Hab.get(), ijA);

  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = stateA + stateB*nstatesA;
      const double *sdataA = sigmavecA->data(stateA)->data();
      for(int stateAp = 0; stateAp < nstatesA; ++stateAp) {
        const int stateABp = stateAp + stateB*nstatesA;
        const double *cdataAp = ccvecs_.first->data(stateAp)->data();
        double element = ddot_(lenab, sdataA, 1, cdataAp, 1); 
        outdata[stateAB + dimerstates_*stateABp] += element;
      }
    }
  }

  shared_ptr<Dvec> sigmavecB;
  sigmavecB = symmetric_ ? sigmavecA : form_sigma_1e(ccvecs_.second, Hab.get() + ijA, ijB);

  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = stateA + stateB*nstatesA;
      const double *sdataB = sigmavecB->data(stateB)->data();
      for(int stateBp = 0; stateBp < nstatesB; ++stateBp) {
        const int stateABp = stateA + stateBp*nstatesA;
        const double *cdataBp = ccvecs_.second->data(stateBp)->data();
        double element = ddot_(lenab, sdataB, 1, cdataBp, 1);

        outdata[stateAB + dimerstates_*stateABp] += element;
      }
    }
  }

  //Done?
  return out;
}

shared_ptr<Matrix> Dimer::compute_intra_activeactive() {
  shared_ptr<Determinants> neut_det = space_->finddet(0,0);
  const int lenab = neut_det->lena() * neut_det->lenb();

  const int nstatesA = nstates_.first;
  const int nstatesB = nstates_.second;

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));
  double *odata = out->data();

  // first do H^{AA}_{AA} case
  shared_ptr<Dvec> sigmavecAA = form_sigma_2e(ccvecs_.first, jops_.first, nact_.first);

  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = stateA + stateB*nstatesA;
      const double *sdataA = sigmavecAA->data(stateA)->data();
      for(int stateAp = 0; stateAp < nstatesA; ++stateAp) {
        const int stateABp = stateAp + stateB*nstatesA;
        const double *cdataAp = ccvecs_.first->data(stateAp)->data();
        double element = ddot_(lenab, sdataA, 1, cdataAp, 1);
        
        odata[stateAB + dimerstates_*stateABp] += element;
      }
    }
  }
  
  // now do H^{BB}_{BB} case
  shared_ptr<Dvec> sigmavecBB;
  sigmavecBB = symmetric_ ? sigmavecAA : form_sigma_2e(ccvecs_.second, jops_.second, nact_.second);

  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = stateA + stateB*nstatesA;
      const double *sdataB = sigmavecBB->data(stateB)->data();
      for(int stateBp = 0; stateBp < nstatesB; ++stateBp) {
        const int stateABp = stateA + stateBp*nstatesA;
        const double *cdataBp = ccvecs_.second->data(stateBp)->data();
        double element = ddot_(lenab, sdataB, 1, cdataBp, 1);
        
        odata[stateAB + dimerstates_*stateABp] += element;
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

  // build <ab|cd> - <ab|dc> matrix
  shared_ptr<Matrix> H_abcd = form_Hmatrix(ijA, ijB);

  // alpha-alpha
  shared_ptr<Matrix> E_ac_alpha = form_EFmatrices_alpha(ccvecs_.first, nstates_.first, ijA);

  // are the monomers identical?
  shared_ptr<Matrix> F_bd_alpha;
  if(symmetric_) F_bd_alpha = E_ac_alpha;
  else F_bd_alpha = form_EFmatrices_alpha(ccvecs_.second, nstates_.second, ijB);

  shared_ptr<Matrix> out(new Matrix( (*E_ac_alpha) * (*H_abcd) ^ (*F_bd_alpha) ));

  // beta-beta
  shared_ptr<Matrix> E_ac_beta = form_EFmatrices_beta(ccvecs_.first, nstates_.first, ijA);

  // are the monomers identical?
  shared_ptr<Matrix> F_bd_beta;
  if(symmetric_) F_bd_beta = E_ac_beta;
  else F_bd_beta = form_EFmatrices_beta(ccvecs_.second, nstates_.second, ijB);

  *out += (*E_ac_beta) * (*H_abcd) ^ (*F_bd_beta);

  // reorder, currently (AA',BB'), want (AB, A'B')
  shared_ptr<Matrix> out2(new Matrix(*out));
  out2->zero();

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

shared_ptr<Matrix> Dimer::form_EFmatrices_alpha(shared_ptr<const Dvec> ccvec, const int ij, const int nstates) const {
  shared_ptr<Matrix> out(new Matrix(nstates*nstates, ij));

  double *edata = out->data();

  shared_ptr<Dvec> c(new Dvec(ccvec->det(), ij));

  const int la = ccvec->det()->lena();
  const int lb = ccvec->det()->lenb();

  for(int state = 0; state < nstates; ++state) {
    const double* source_base = ccvec->data(state)->data();

    c->zero();

    for(int ac = 0; ac < ij; ++ac) {
      double* target_base = c->data(ac)->data();

      // alpha-alpha
      for(auto& iter : ccvec->det()->phia(ac) ) {
        const double sign = static_cast<double>(get<1>(iter));
        double* target_array = target_base + get<2>(iter)*lb;
        const double* source = source_base + get<0>(iter)*lb;
        daxpy_(lb, sign, source, 1, target_array, 1);
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

  return out;
}

shared_ptr<Matrix> Dimer::form_EFmatrices_beta(shared_ptr<const Dvec> ccvec, const int ij, const int nstates) const{
  shared_ptr<Matrix> out(new Matrix(nstates*nstates, ij));

  double *edata = out->data();

  shared_ptr<Dvec> c(new Dvec(ccvec->det(), ij));

  const int la = ccvec->det()->lena();
  const int lb = ccvec->det()->lenb();

  for(int state = 0; state < nstates; ++state) {
    const double* source_base = ccvec->data(state)->data();

    c->zero();

    for(int ac = 0; ac < ij; ++ac) {
      double *target_base = c->data(ac)->data();

      // beta-beta
      for(auto& iter : ccvec->det()->phib(ac)) {
        const double sign = static_cast<double>(get<1>(iter));
        double* target_array = target_base + get<2>(iter);
        const double* source = source_base + get<0>(iter);
        daxpy_(la, sign, source, lb, target_array, lb);
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

  return out;
}

shared_ptr<Matrix> Dimer::form_Hmatrix(const int ijA, const int ijB) const {
  shared_ptr<Matrix> out(new Matrix(ijA, ijB));
  double* odata = out->data();
  const int nbasisA = nbasis_.first;
  const int nbasisB = nbasis_.second;
  
  for(int b = 0; b < nbasisB; ++b) {
    for(int d = 0; d < nbasisB; ++d) {
      for(int a = 0; a < nbasisA; ++a) {
        for(int c = 0; c < nbasisA; ++c) {
          *odata++ = jop_->mo2e_hz(act<0>(a), act<1>(b), act<0>(c), act<1>(d)) - jop_->mo2e_hz(act<0>(a), act<1>(b), act<1>(d), act<0>(c));
        }
      }
    }
  }

  return out;
}
