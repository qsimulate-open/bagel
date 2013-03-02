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

#include <src/wfn/geometry.h>
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

  cout << " ===== Starting construction of dimer Hamiltonian ===== " << endl;

  // TODO not compatible with having different active spaces
  space_ = shared_ptr<Space>(new Space(ccvecs_.first->det(), 1));
  hamiltonian_ = shared_ptr<Matrix>(new Matrix(dimerstates_, dimerstates_));

  // create jop_ object (which effectively includes all closed-closed interactions)
  jop_ = shared_ptr<DimerJop>(new DimerJop(sref_, ncore, ncore + nact_.first, ncore + nact, scoeff_));

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

  hamiltonian_->print("Dimer Hamiltonian", dimerstates_);

  vector<double> eigvec(dimerstates_, 0.0);
  hamiltonian_->diagonalize(eigvec.data());

  cout << endl << " ===== Adiabatic state energies ===== " << endl;

  int istate = 0;
  for( auto& ieig : eigvec ) {
    cout << "   state " << setw(3) << istate << ": " << setprecision(12) << setw(16) << ieig << endl;
     ++istate;
  }
}

shared_ptr<Matrix> Dimer::compute_closeclose() {
  double core = sgeom_->nuclear_repulsion() + jop_->core_energy();

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));
  out->add_diag(core);

  return out;
}

shared_ptr<Matrix> Dimer::compute_closeactive() {
  const int ncore = ncore_.first + ncore_.second;
  const int nact  = nact_.first + nact_.second;

  const int nactA = nact_.first;
  const int nactB = nact_.second;

  const int nstatesA = nstates_.first;
  const int nstatesB = nstates_.second;

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  const int ijA = nactA*nactA;
  const int ijB = nactB*nactB;

  shared_ptr<Determinants> neut_det = space_->finddet(0,0);
  const int lenab = neut_det->lena() * neut_det->lenb();

  shared_ptr<Dvec> sigmavecA = form_sigma_1e(ccvecs_.first, jop_->mo1e_first(), ijA);

  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = dimerstate(stateA, stateB);
      const double *sdataA = sigmavecA->data(stateA)->data();
      for(int stateAp = 0; stateAp < nstatesA; ++stateAp) {
        const int stateABp = dimerstate(stateAp, stateB);
        const double *cdataAp = ccvecs_.first->data(stateAp)->data();
        double element = ddot_(lenab, sdataA, 1, cdataAp, 1); 
        out->element(stateAB, stateABp) += element;
      }
    }
  }

  shared_ptr<Dvec> sigmavecB = form_sigma_1e(ccvecs_.second, jop_->mo1e_second(), ijB);

  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = dimerstate(stateA, stateB);
      const double *sdataB = sigmavecB->data(stateB)->data();
      for(int stateBp = 0; stateBp < nstatesB; ++stateBp) {
        const int stateABp = dimerstate(stateA, stateBp);
        const double *cdataBp = ccvecs_.second->data(stateBp)->data();
        double element = ddot_(lenab, sdataB, 1, cdataBp, 1);
        out->element(stateAB, stateABp) += element;
      }
    }
  }

  //Done
  return out;
}

shared_ptr<Matrix> Dimer::compute_intra_activeactive() {
  shared_ptr<Determinants> neut_det = space_->finddet(0,0);
  const int lenab = neut_det->lena() * neut_det->lenb();

  const int nstatesA = nstates_.first;
  const int nstatesB = nstates_.second;

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  // first do H^{AA}_{AA} case
  shared_ptr<Dvec> sigmavecAA = form_sigma_2e(ccvecs_.first, jop_->mo2e_first(), nact_.first);

  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = dimerstate(stateA, stateB);
      const double *sdataA = sigmavecAA->data(stateA)->data();
      for(int stateAp = 0; stateAp < nstatesA; ++stateAp) {
        const int stateABp = dimerstate(stateAp, stateB);
        const double *cdataAp = ccvecs_.first->data(stateAp)->data();
        double element = ddot_(lenab, sdataA, 1, cdataAp, 1);
        
        out->element(stateAB, stateABp) += element;
      }
    }
  }
  
  // now do H^{BB}_{BB} case
  shared_ptr<Dvec> sigmavecBB = form_sigma_2e(ccvecs_.second, jop_->mo2e_second(), nact_.second);

  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = dimerstate(stateA, stateB);
      const double *sdataB = sigmavecBB->data(stateB)->data();
      for(int stateBp = 0; stateBp < nstatesB; ++stateBp) {
        const int stateABp = dimerstate(stateA, stateBp);
        const double *cdataBp = ccvecs_.second->data(stateBp)->data();
        double element = ddot_(lenab, sdataB, 1, cdataBp, 1);
        
        out->element(stateAB, stateABp) += element;
      }
    }
  }

  return out;
}


//And now for the hard part
shared_ptr<Matrix> Dimer::compute_inter_activeactive() {
  // ijA and ijB need to match the compression of the det
  const int nactA = nact_.first;
  const int ijA = nactA*nactA;
  const int nstatesA = nstates_.first;

  const int nactB = nact_.second;
  const int ijB = nactB*nactB;
  const int nstatesB = nstates_.second;

  // alpha-alpha
  shared_ptr<Matrix> E_ac_alpha = form_EFmatrices_alpha(ccvecs_.first, nstatesA, ijA);
  shared_ptr<Matrix> F_bd_alpha = form_EFmatrices_alpha(ccvecs_.second, nstatesB, ijB);

  // beta-beta
  shared_ptr<Matrix> E_ac_beta = form_EFmatrices_beta(ccvecs_.first, nstatesA, ijA);
  shared_ptr<Matrix> F_bd_beta = form_EFmatrices_beta(ccvecs_.second, nstatesB, ijB);

  // build JK and J matrices
  #define SPLITJK
  #ifndef SPLITJK
  shared_ptr<Matrix> JK_abcd = form_JKmatrix(ijA, ijB);
  shared_ptr<Matrix> J_abcd = form_Jmatrix(ijA, ijB);

  shared_ptr<Matrix> tmp(new Matrix(nstatesA*nstatesA, nstatesB*nstatesB));

  *tmp += (*E_ac_alpha) * (*JK_abcd) ^ (*F_bd_alpha);
  *tmp += (*E_ac_beta) * (*JK_abcd) ^ (*F_bd_beta);

  *tmp += (*E_ac_alpha) * (*J_abcd) ^ (*F_bd_beta);
  *tmp += (*E_ac_beta) * (*J_abcd) ^ (*F_bd_alpha);

  // reorder, currently (AA',BB'), want (AB, A'B')
  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  for(int B = 0; B < nstatesB; ++B) {
    for(int A = 0; A < nstatesA; ++A) {
      const int AB = dimerstate(A,B);
      for(int Bp = 0; Bp < nstatesB; ++Bp) {
        const int BBp = Bp + B*nstatesB;
        for(int Ap = 0; Ap < nstatesA; ++Ap) {
          const int ABp = dimerstate(Ap,Bp);
          const int AAp = Ap + A*nstatesA;
          out->element(AB,ABp) = tmp->element(AAp,BBp);
        }
      }
    }
  }
  #else
  shared_ptr<Matrix> J_abcd = form_Jmatrix(ijA, ijB);
  shared_ptr<Matrix> K_abcd = form_JKmatrix(ijA, ijB);

  shared_ptr<Matrix> Jtmp(new Matrix(nstatesA*nstatesA, nstatesB*nstatesB));
  
  *Jtmp += (*E_ac_alpha) * (*J_abcd) ^ (*F_bd_alpha);
  *Jtmp += (*E_ac_beta) * (*J_abcd) ^ (*F_bd_beta);

  *Jtmp += (*E_ac_alpha) * (*J_abcd) ^ (*F_bd_beta);
  *Jtmp += (*E_ac_beta) * (*J_abcd) ^ (*F_bd_alpha);

  shared_ptr<Matrix> Ktmp(new Matrix(nstatesA*nstatesA, nstatesB*nstatesB));

  *Ktmp += (*E_ac_alpha) * (*K_abcd) ^ (*F_bd_alpha);
  *Ktmp += (*E_ac_beta) * (*K_abcd) ^ (*F_bd_beta);

  *Ktmp *= -1.0;
  
  shared_ptr<Matrix> Jout(new Matrix(dimerstates_, dimerstates_));
  shared_ptr<Matrix> Kout(new Matrix(dimerstates_, dimerstates_));
  for(int B = 0; B < nstatesB; ++B) {
    for(int A = 0; A < nstatesA; ++A) {
      const int AB = dimerstate(A,B);
      for(int Bp = 0; Bp < nstatesB; ++Bp) {
        const int BBp = Bp + B*nstatesB;
        for(int Ap = 0; Ap < nstatesA; ++Ap) {
          const int ABp = dimerstate(Ap, Bp);
          const int AAp = Ap + A*nstatesA;
          Jout->element(AB,ABp) = Jtmp->element(AAp,BBp);
          Kout->element(AB,ABp) = Ktmp->element(AAp,BBp);
        }
      }
    }
  }

  Jout->print( " J ", dimerstates_);
  Kout->print( " K ", dimerstates_);

  shared_ptr<Matrix> out(new Matrix(*Jout + *Kout));
  #endif

  return out;
}

shared_ptr<Matrix> Dimer::form_EFmatrices_alpha(shared_ptr<const Dvec> ccvec, const int nstates, const int ij) const {
  shared_ptr<Matrix> tmp(new Matrix(ij, nstates*nstates));

  double *edata = tmp->data();
  shared_ptr<const Determinants> det = space_->finddet(0,0);


  shared_ptr<Dvec> c(new Dvec(det, ij));

  const int la = det->lena();
  const int lb = det->lenb();
  const int lab = la*lb;

  for(int state = 0; state < nstates; ++state) {
    const double* source_base = ccvec->data(state)->data();

    c->zero();

    for(int ac = 0; ac < ij; ++ac) {
      double* target_base = c->data(ac)->data();

      // alpha-alpha
      for(auto& iter : det->phia(ac) ) {
        const double sign = static_cast<double>(iter.sign);
        double* target = target_base + iter.target*lb;
        const double* source = source_base + iter.source*lb;
        daxpy_(lb, sign, source, 1, target, 1);
      }
    }

    // | C > ^A_ac is done
    for(int statep = 0; statep < nstates; ++statep) {
      const double *adata = ccvec->data(statep)->data();
      for(int ac = 0; ac < ij; ++ac, ++edata) {
        const double *cdata = c->data(ac)->data();
        *edata = ddot_(lab, adata, 1, cdata, 1);

      }
    }
  }

  return tmp->transpose();
}

shared_ptr<Matrix> Dimer::form_EFmatrices_beta(shared_ptr<const Dvec> ccvec, const int nstates, const int ij) const{
  shared_ptr<Matrix> tmp(new Matrix(ij, nstates*nstates));

  double *edata = tmp->data();
  shared_ptr<const Determinants> det = space_->finddet(0,0);

  shared_ptr<Dvec> c(new Dvec(det, ij));

  const int la = det->lena();
  const int lb = det->lenb();
  const int lab = la*lb;

  for(int state = 0; state < nstates; ++state) {
    const double* source_base = ccvec->data(state)->data();

    c->zero();

    for(int ac = 0; ac < ij; ++ac) {
      double *target_base = c->data(ac)->data();

      // beta-beta
      for(auto& iter : det->phib(ac)) {
        const double sign = static_cast<double>(iter.sign);
        double* target = target_base + iter.target;
        const double* source = source_base + iter.source;
        daxpy_(la, sign, source, lb, target, lb);
      }
    }

    // | C > ^A_ac is done
    for(int statep = 0; statep < nstates; ++statep) {
      const double *adata = ccvec->data(statep)->data();
      for(int ac = 0; ac < ij; ++ac, ++edata) {
        const double *cdata = c->data(ac)->data();
        *edata = ddot_(lab, adata, 1, cdata, 1);
      }
    }
  }

  return tmp->transpose();
}

shared_ptr<Matrix> Dimer::form_JKmatrix(const int ijA, const int ijB) const {
  shared_ptr<Matrix> out(new Matrix(ijA, ijB));
  double* odata = out->data();
  const int nactA = nact_.first;
  const int nactB = nact_.second;

  for(int b = 0; b < nactB; ++b) {
    for(int d = 0; d < nactB; ++d) {
      for(int a = 0; a < nactA; ++a) {
        for(int c = 0; c < nactA; ++c, ++odata) {
          #ifndef SPLITJK
          *odata = jop_->mo2e_hz(act<0>(a), act<1>(b), act<0>(c), act<1>(d)) - jop_->mo2e_hz(act<0>(a), act<1>(b), act<1>(d), act<0>(c));
          #else
          *odata = jop_->mo2e_hz(act<0>(a), act<1>(b), act<1>(d), act<0>(c));
          #endif
        }
      }
    }
  }

  return out;
}

shared_ptr<Matrix> Dimer::form_Jmatrix(const int ijA, const int ijB) const {
  shared_ptr<Matrix> out(new Matrix(ijA, ijB));
  double* odata = out->data();
  const int nactA = nact_.first;
  const int nactB = nact_.second;

  for(int b = 0; b < nactB; ++b) {
    for(int d = 0; d < nactB; ++d) {
      for(int a = 0; a < nactA; ++a) {
        for(int c = 0; c < nactA; ++c, ++odata) {
          *odata = jop_->mo2e_hz(act<0>(a), act<1>(b), act<0>(c), act<1>(d));
        }
      }
    }
  }

  return out;
}
