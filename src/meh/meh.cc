//
// BAGEL - Parallel electron correlation program.
// Filename: meh.cc
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
MultiExcitonHamiltonian::MultiExcitonHamiltonian(shared_ptr<DimerCISpace> cispace) :
  ref_(cispace->ref()), coeff_(cispace->ref()), cispace_(cispace), 
  dimerbasis_(cispace->dimerbasis()), dimerclosed_(cispace->dimerclosed()), dimeractive_(cispace->dimeractive()),
  nact_(cispace->nact()), nbasis_(cispace->nbasis()), nstates_(cispace->nstates())
{
  common_init();
}

void MultiExcitonHamiltonian::common_init() {
  dimerstates_ = nstates_.first * nstates_.second;

  jop_ = shared_ptr<DimerJop>(new DimerJop(ref_, dimerclosed_, dimerclosed_ + nact_.first, dimerclosed_ + nact, coeff_));
}

void MultiExcitonHamiltonian::compute() {
  cout << " ===== Starting construction of dimer Hamiltonian ===== " << endl;

  hamiltonian_ = shared_ptr<Matrix>(new Matrix(dimerstates_, dimerstates_));

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
  double core = ref_->geom()->nuclear_repulsion() + jop_->core_energy();

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));
  out->add_diag(core);

  return out;
}

shared_ptr<Matrix> Dimer::compute_closeactive() {
  const int nclosed = dimerclosed_;
  const int nact = dimeractive_;

  const int nactA = nact_.first;
  const int nactB = nact_.second;

  const int nstatesA = nstates_.first;
  const int nstatesB = nstates_.second;

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  {
    shared_ptr<const Determinants> detA = cispace_->finddet<0>(0,0);
    shared_ptr<const Dvec> ccvecA = cispace_->ccvec<0>();

    const int lenab = detA->lena() * detA->lenb();

    shared_ptr<Dvec> sigmavecA = form_sigma_1e(ccvecA, jop_->mo1e_first(), nactA*nactA);

    for(int stateA = 0; stateA < nstatesA; ++stateA) {
      for(int stateB = 0; stateB < nstatesB; ++stateB) {
        const int stateAB = dimerstate(stateA, stateB);
        const double *sdataA = sigmavecA->data(stateA)->data();
        for(int stateAp = 0; stateAp < nstatesA; ++stateAp) {
          const int stateABp = dimerstate(stateAp, stateB);
          const double *cdataAp = ccvecA->data(stateAp)->data();
          double value = ddot_(lenab, sdataA, 1, cdataAp, 1); 
          out->element(stateAB, stateABp) += value;
        }
      }
    }
  }

  {
    shared_ptr<const Determinants> detB = cispace_->finddet<1>(0,0);
    shared_ptr<const Dvec> ccvecB = cispace_->ccvec<1>();

    const int lenab = detB->lena() * detB->lenb();

    shared_ptr<Dvec> sigmavecB = form_sigma_1e(ccvecB, jop_->mo1e_second(), nstatesB*nstatesB);

    for(int stateA = 0; stateA < nstatesA; ++stateA) {
      for(int stateB = 0; stateB < nstatesB; ++stateB) {
        const int stateAB = dimerstate(stateA, stateB);
        const double *sdataB = sigmavecB->data(stateB)->data();
        for(int stateBp = 0; stateBp < nstatesB; ++stateBp) {
          const int stateABp = dimerstate(stateA, stateBp);
          const double *cdataBp = ccvecB->data(stateBp)->data();
          double value = ddot_(lenab, sdataB, 1, cdataBp, 1);
          out->element(stateAB, stateABp) += value;
        }
      }
    }
  }

  return out;
}

shared_ptr<Matrix> Dimer::compute_intra_activeactive() {
  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  const int nstatesA = nstates_.first;
  const int nstatesB = nstates_.second;

  const int nactA = nact_.first;
  const int nactB = nact_.second;

  // first H^{AA}_{AA}
  {
    shared_ptr<const Determinants> detA = cispace_->finddet<0>(0,0);
    shared_ptr<const Dvec> ccvecA = cispace_->ccvec<0>();

    const int lenab = detA->lena() * detA->lenb();

    shared_ptr<Dvec> sigmavecAA = form_sigma_2e(ccvecA, jop_->mo2e_first(), nactA);

    for(int stateA = 0; stateA < nstatesA; ++stateA) {
      for(int stateB = 0; stateB < nstatesB; ++stateB) {
        const int stateAB = dimerstate(stateA, stateB);
        const double *sdataA = sigmavecAA->data(stateA)->data();
        for(int stateAp = 0; stateAp < nstatesA; ++stateAp) {
          const int stateABp = dimerstate(stateAp, stateB);
          const double *cdataAp = ccvecA->data(stateAp)->data();
          double value; = ddot_(lenab, sdataA, 1, cdataAp, 1);
          
          out->element(stateAB, stateABp) += value;;
        }
      }
    }
  }
  
  // now do H^{BB}_{BB} case
  {
    shared_ptr<const Determinants> detB = cispace_->finddet<1>(0,0);
    shared_ptr<const Dvec> ccvecB = cispace_->ccvec<1>();

    const int lenab = detB->lena() * detB->lenb();
    shared_ptr<Dvec> sigmavecBB = form_sigma_2e(ccvecB, jop_->mo2e_second(), nactB);

    for(int stateA = 0; stateA < nstatesA; ++stateA) {
      for(int stateB = 0; stateB < nstatesB; ++stateB) {
        const int stateAB = dimerstate(stateA, stateB);
        const double *sdataB = sigmavecBB->data(stateB)->data();
        for(int stateBp = 0; stateBp < nstatesB; ++stateBp) {
          const int stateABp = dimerstate(stateA, stateBp);
          const double *cdataBp = ccvecB->data(stateBp)->data();
          double value; = ddot_(lenab, sdataB, 1, cdataBp, 1);
          
          out->element(stateAB, stateABp) += value;;
        }
      }
    }
  }

  return out;
}

shared_ptr<Matrix> Dimer::compute_inter_activeactive() {
  shared_ptr<const Dvec> ccvecA = cispace_->ccvec<0>();
  shared_ptr<const Dvec> ccvecB = cispace_->ccvec<1>();

  const int nactA = nact_.first;
  const int ijA = nactA*nactA;
  const int nstatesA = nstates_.first;

  const int nactB = nact_.second;
  const int ijB = nactB*nactB;
  const int nstatesB = nstates_.second;

  // alpha-alpha
  shared_ptr<Matrix> gamma_AA_alpha = form_gamma_alpha(ccvecA);
  shared_ptr<Matrix> gamma_BB_alpha = form_gamma_alpha(ccvecB);

  // beta-beta
  Matrix gamma_AA_beta = *form_gamma_beta(ccvecA);
  Matrix gamma_BB_beta = *form_gamma_beta(ccvecB);

  // build J and K matrices
  shared_ptr<Matrix> Jmatrix, Kmatrix;
  tie(Jmatrix, Kmatrix) = form_JKmatrices();

  Matrix tmp(nstatesA*nstatesA, nstatesB*nstatesB);

  tmp += gamma_AA_alpha * (*Kmatrix) ^ gamma_BB_alpha;
  tmp += gamma_AA_beta * (*Kmatrix) ^ gamma_BB_beta;

  tmp += (gamma_AA_alpha + gamma_AA_beta) * (*Jmatrix) ^ (gamma_BB_alpha + gamma_BB_beta);

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
          out->element(AB,ABp) = tmp(AAp,BBp);
        }
      }
    }
  }

  return out;
}

shared_ptr<Matrix> Dimer::form_gamma_alpha(shared_ptr<const Dvec> ccvec) const {
  const int nstates = ccvec->ij();

  shared_ptr<const Determinants> det = ccvec->det();
  const int norb = det->norb();
  const int ij = norb * norb;

  Matrix tmp(ij, nstates*nstates);

  double *edata = tmp.data();

  Dvec c(det, ij);

  const int la = det->lena();
  const int lb = det->lenb();
  const int lab = la*lb;

  for(int state = 0; state < nstates; ++state) {
    const double* source_base = ccvec->data(state)->data();

    c.zero();

    for(int ac = 0; ac < ij; ++ac) {
      double* target_base = c.data(ac)->data();

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
        const double *cdata = c.data(ac)->data();
        *edata = ddot_(lab, adata, 1, cdata, 1);

      }
    }
  }

  return tmp.transpose();
}

shared_ptr<Matrix> Dimer::form_gamma_beta(shared_ptr<const Dvec> ccvec) const{
  const int nstates = ccvec->ij();

  shared_ptr<const Determinants> det = ccvec->det();
  const int norb = det->norb();
  const int ij = norb * norb;

  Matrix tmp(ij, nstates*nstates);

  double *edata = tmp.data();

  Dvec c(det, ij);

  const int la = det->lena();
  const int lb = det->lenb();
  const int lab = la*lb;

  for(int state = 0; state < nstates; ++state) {
    const double* source_base = ccvec->data(state)->data();

    c.zero();

    for(int ac = 0; ac < ij; ++ac) {
      double *target_base = c.data(ac)->data();

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
        const double *cdata = c.data(ac)->data();
        *edata = ddot_(lab, adata, 1, cdata, 1);
      }
    }
  }

  return tmp.transpose();
}

template<int A, int B, int C, int D>
pair<shared_ptr<Matrix>, shared_ptr<Matrix>> Dimer::form_Jmatrices() const {
  const int nactA = nact_.first;
  const int nactB = nact_.second;

  int ijA = 1;
  int unitA = 4 - (A + B + C + D);
  for ( int i = 0; i < unitA; ++i ) ijA *= nactA;

  int ijB = 1;
  unitB = A + B + C + D;
  for ( int i = 0; i < unitB; ++i ) ijB *= nactB;

  shared_ptr<Matrix> Jout(new Matrix(ijA, ijB));
  shared_ptr<Matrix> Kout(new Matrix(ijA, ijB));

  // Because of the templating, all of the index mess SHOULD be done at compile time
  for(int d = 0; d < (D == 0 ? nactA : nactB); ++d) {
    for(int c = 0; c < (C == 0 ? nactA : nactB); ++c) {
      for(int b = 0; b < (B == 0 ? nactA : nactB); ++b) {
        for(int a = 0; a < (A == 0 ? nactA : nactB); ++a) {
          int iA, jB;
          std::tie(iA, jB) = index<A,B,C,D>(a,b,c,d);

          Jout->element(iA,jB) = jop_->mo2e_hz(active<A>(a), active<B>(b), active<C>(c), active<D>(d));
          Kout->element(iA,jB) = jop_->mo2e_hz(active<A>(a), active<B>(b), active<D>(d), active<C>(c));
        }
      }
    }
  }

  return std::make_pair(Jout,Kout);
}
