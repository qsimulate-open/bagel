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
#include <src/dimer/dimer_cispace.h>
#include <src/dimer/dimer.h>
#include <src/util/matrix.h>
#include <src/meh/meh.h>

using namespace std;
using namespace bagel;

/************************************************************************************
*  General note: For the moment, everything is written so it will be easiest to     *
* debug, using as much code as possible from other functions already written. This  *
* will eventually be fixed.                                                         *
************************************************************************************/
MultiExcitonHamiltonian::MultiExcitonHamiltonian(shared_ptr<Dimer> dimer, shared_ptr<DimerCISpace> cispace) :
  ref_(dimer->sref()), coeff_(dimer->scoeff()), cispace_(cispace), 
  dimerbasis_(dimer->dimerbasis()), dimerclosed_(dimer->sref()->nclosed()), dimeractive_(dimer->sref()->nact()),
  nact_(dimer->nact()), nbasis_(dimer->nbasis()), nstates_(cispace->nstates())
{
  common_init();
}

void MultiExcitonHamiltonian::common_init() {
  dimerstates_ = nstates_.first * nstates_.second;

  jop_ = shared_ptr<DimerJop>(new DimerJop(ref_, dimerclosed_, dimerclosed_ + nact_.first, dimerclosed_ + dimeractive_, coeff_));
  cispace_->complete();

  // Process DimerCISpace to form and organize needed Civecs
  vector<shared_ptr<Civec>> vec_m0_q0_A;
  vector<shared_ptr<Civec>> vec_m0_q0_B;

  vector<shared_ptr<Civec>> vec_m1_qC_A;
  vector<shared_ptr<Civec>> vec_m_1_qC_A;

  vector<shared_ptr<Civec>> vec_m1_qC_B;
  vector<shared_ptr<Civec>> vec_m_1_qC_B;

  vector<shared_ptr<Civec>> vec_m1_qA_A;
  vector<shared_ptr<Civec>> vec_m_1_qA_A;

  vector<shared_ptr<Civec>> vec_m1_qA_B;
  vector<shared_ptr<Civec>> vec_m_1_qA_B;

  // Start by processing the singlet states
  {
    shared_ptr<const Dvec> tmpvec = cispace_->ccvec<0>(0,0);
    vec_m0_q0_A.insert(vec_m0_q0_A.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());
  }

  {
    shared_ptr<const Dvec> tmpvec = cispace_->ccvec<1>(0,0);
    vec_m0_q0_B.insert(vec_m0_q0_B.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());
  }

  // Now add in cation states
  if (cispace_->cations()) {
    {
      shared_ptr<const Dvec> tmpvec = cispace_->ccvec<0>(0,-1);
      vec_m1_qC_A.insert(vec_m1_qC_A.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());

      shared_ptr<Determinants> det = cispace_->add_det<0>(-1,0);
      tmpvec = tmpvec->spinflip(det);
      vec_m_1_qC_A.insert(vec_m_1_qC_A.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());
    }

    {
      shared_ptr<const Dvec> tmpvec = cispace_->ccvec<1>(0,-1);
      vec_m1_qC_B.insert(vec_m1_qC_B.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());

      shared_ptr<Determinants> det = cispace_->add_det<1>(-1,0);
      tmpvec = tmpvec->spinflip(det);
      vec_m_1_qC_B.insert(vec_m_1_qC_B.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());
    }
  }

  // Now for the anion states
  if(cispace_->anions()) {
    {
      shared_ptr<const Dvec> tmpvec = cispace_->ccvec<0>(1,0);
      vec_m1_qA_A.insert(vec_m1_qA_A.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());

      shared_ptr<Determinants> det = cispace_->add_det<0>(0,1);
      tmpvec = tmpvec->spinflip(det);
      vec_m_1_qA_A.insert(vec_m_1_qA_A.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());
    }

    {
      shared_ptr<const Dvec> tmpvec = cispace_->ccvec<1>(1,0);
      vec_m1_qA_B.insert(vec_m1_qA_B.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());

      shared_ptr<Determinants> det = cispace_->add_det<1>(0,1);
      tmpvec = tmpvec->spinflip(det);
      vec_m_1_qA_B.insert(vec_m_1_qA_B.end(), tmpvec->dvec().begin(), tmpvec->dvec().end());
    }
  }

  // Package like civecs into Dvecs
  dimerstates_ = 0;

  // First, singlets
  shared_ptr<Dvec> S_A(new Dvec(vec_m0_q0_A));
  shared_ptr<Dvec> S_B(new Dvec(vec_m0_q0_B));

  subspaces_.push_back(DimerSubspace(ChargeSpin::SS, dimerstates_, make_pair(S_A,S_B)));
  dimerstates_ += subspaces_.back().dimerstates();

  // Now, AC, if we've got em
  if (cispace_->anions() && cispace_->cations()) {
    shared_ptr<Dvec> Aalpha_A(new Dvec(vec_m1_qA_A));
    shared_ptr<Dvec> Abeta_A(new Dvec(vec_m_1_qA_A));
    shared_ptr<Dvec> Aalpha_B(new Dvec(vec_m1_qA_B));
    shared_ptr<Dvec> Abeta_B(new Dvec(vec_m_1_qA_B));

    shared_ptr<Dvec> Calpha_A(new Dvec(vec_m1_qC_A));
    shared_ptr<Dvec> Cbeta_A(new Dvec(vec_m_1_qC_A));
    shared_ptr<Dvec> Calpha_B(new Dvec(vec_m1_qC_B));
    shared_ptr<Dvec> Cbeta_B(new Dvec(vec_m_1_qC_B));

    subspaces_.push_back(DimerSubspace(ChargeSpin::AaCb, dimerstates_, make_pair(Aalpha_A,Cbeta_B)));
    dimerstates_ += subspaces_.back().dimerstates();

    subspaces_.push_back(DimerSubspace(ChargeSpin::AbCa, dimerstates_, make_pair(Abeta_A,Calpha_B)));
    dimerstates_ += subspaces_.back().dimerstates();
    
    subspaces_.push_back(DimerSubspace(ChargeSpin::CaAb, dimerstates_, make_pair(Calpha_A,Abeta_B)));
    dimerstates_ += subspaces_.back().dimerstates();

    subspaces_.push_back(DimerSubspace(ChargeSpin::CbAa, dimerstates_, make_pair(Cbeta_A,Aalpha_B)));
    dimerstates_ += subspaces_.back().dimerstates();
  }
}

Coupling MultiExcitonHamiltonian::coupling_type(DimerSubspace& AB, DimerSubspace& ApBp) {
  pair<int,int> neleaAB = make_pair(AB.ci<0>()->det()->nelea(), AB.ci<1>()->det()->nelea());
  pair<int,int> nelebAB = make_pair(AB.ci<0>()->det()->neleb(), AB.ci<1>()->det()->neleb());

  pair<int,int> neleaApBp = make_pair(ApBp.ci<0>()->det()->nelea(), ApBp.ci<1>()->det()->nelea());
  pair<int,int> nelebApBp = make_pair(ApBp.ci<0>()->det()->neleb(), ApBp.ci<1>()->det()->neleb());

  // AlphaTransfer and BetaTransfer
  pair<int,int> AT = make_pair(neleaApBp.first - neleaAB.first, neleaApBp.second - neleaAB.second);
  pair<int,int> BT = make_pair(nelebApBp.first - nelebAB.first, nelebApBp.second - nelebAB.second);

  /************************************************************
  *  BT\AT  | ( 0, 0) | (+1,-1) | (-1,+1) | (+2,-2) | (-2,+2) *
  *-----------------------------------------------------------*
  * ( 0, 0) |  diag   |  aET    |  -aET   |  aaET   | -aaET   *
  * (+1,-1) |  bET    |  dABT   |  ABflp  |         |         *
  * (-1,+1) | -bET    | BAflp   | -dABT   |         |         *
  * (+2,-2) |  bbET   |         |         |         |         *
  * (-2,+2) | -bbET   |         |         |         |         *
  ************************************************************/

  const int icouple = coupling_index(AT, BT);

  if      ( icouple == coupling_index( 0, 0, 0, 0) ) return Coupling::diagonal;
  else if ( icouple == coupling_index( 0, 0,+1,-1) ) return Coupling::bET;
  else if ( icouple == coupling_index( 0, 0,-1,+1) ) return Coupling::inv_bET;
  else if ( icouple == coupling_index(+1,-1, 0, 0) ) return Coupling::aET;
  else if ( icouple == coupling_index(+1,-1,+1,-1) ) return Coupling::abET;
  else if ( icouple == coupling_index(+1,-1,-1,+1) ) return Coupling::baFlip;
  else if ( icouple == coupling_index(-1,+1, 0, 0) ) return Coupling::inv_aET;
  else if ( icouple == coupling_index(-1,+1,+1,-1) ) return Coupling::abFlip;
  else if ( icouple == coupling_index(-1,+1,-1,+1) ) return Coupling::inv_abET;
  else if ( icouple == coupling_index(+2,-2, 0, 0) ) return Coupling::aaET;
  else if ( icouple == coupling_index(-2,+2, 0, 0) ) return Coupling::inv_aaET;
  else if ( icouple == coupling_index( 0, 0,+2,-2) ) return Coupling::bbET;
  else if ( icouple == coupling_index( 0, 0,-2,+2) ) return Coupling::inv_bbET;
  else                                               return Coupling::none;
}

void MultiExcitonHamiltonian::compute() {
  cout << endl << " ===== Starting construction of dimer Hamiltonian ===== " << endl;

  hamiltonian_ = shared_ptr<Matrix>(new Matrix(dimerstates_, dimerstates_));

  cout << "  o Computing diagonal blocks" << endl; 
  for (auto& subspace : subspaces_) {
    hamiltonian_->add_block(subspace.offset(), subspace.offset(), compute_diagonal_block(subspace));
  }

  cout << "  o Computing off-diagonal blocks" << endl;
  const int spaceij = subspaces_.size();
  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

      shared_ptr<Matrix> block = couple_blocks(*iAB, *jAB);

      hamiltonian_->add_block(ioff, joff, block);
      hamiltonian_->add_block(joff, ioff, block->transpose());
    }
  }

  adiabats_ = shared_ptr<Matrix>(new Matrix(*hamiltonian_));

  energies_ = vector<double>(dimerstates_, 0.0);
  adiabats_->diagonalize(energies_.data());
}


std::shared_ptr<Matrix> MultiExcitonHamiltonian::form_gamma(std::shared_ptr<const Dvec> ccvecA, std::shared_ptr<const Dvec> ccvecAp, std::shared_ptr<Quantization> action) const { 
  const int nstatesA = ccvecA->ij();
  const int nstatesAp = ccvecAp->ij();

  std::shared_ptr<const Determinants> detA = ccvecA->det();
  const int norb = detA->norb();
  const int ij = action->ij(norb);

  Matrix tmp(ij, nstatesA*nstatesAp);

  double *edata = tmp.data();

  for(int state = 0; state < nstatesA; ++state) {
    std::shared_ptr<Dvec> c = action->compute(ccvecA->data(state));

    // | C > ^A_ac is done
    for(int statep = 0; statep < nstatesAp; ++statep) {
      for(int ac = 0; ac < ij; ++ac, ++edata) {
        *edata = c->data(ac)->ddot(*ccvecAp->data(statep)); 
      }   
    }   
  }

  return tmp.transpose();
}

void MultiExcitonHamiltonian::reorder_matrix(const double* source, double* target, 
  const int nA, const int nAp, const int nB, const int nBp) const
{
  const int nstatesAB = nA * nB;
  const int nstatesAA = nA * nAp;

  for (int Bp = 0; Bp < nBp; ++Bp) {
    for (int Ap = 0; Ap < nAp; ++Ap) {
      const int ABp = Ap + nAp * Bp;
      for (int B = 0; B < nB; ++B) {
        const int BBp = B + nB * Bp;
        for (int A = 0; A < nA; ++A) {
          const int AAp = A + nA * Ap;
          const int AB = A + nA * B;
          target[AB + nstatesAB * ABp] = source[AAp + nstatesAA * BBp];
        }
      }
    }
  }
}

void MultiExcitonHamiltonian::print_hamiltonian(const string title, const int nstates) {
  hamiltonian_->print(title, nstates);
}

void MultiExcitonHamiltonian::print_energies(const string title, const int nstates) {
  const int end = std::min(nstates, dimerstates_);
  cout << endl << " ===== " << title << " =====" << endl;
  for (int istate = 0; istate < end; ++istate)
    cout << "   state  " << setw(3) << istate << ": " << setprecision(12) << setw(16) << energies_.at(istate) << endl;
}

void MultiExcitonHamiltonian::print_adiabats(const string title, const int nstates) {
  adiabats_->print(title, nstates);
}

void MultiExcitonHamiltonian::print(const int nstates) {
  print_hamiltonian("ME Hamiltonian", nstates);
  cout << endl;
  print_energies("Adiabatic Energies", nstates);
  cout << endl;
  print_adiabats("Adiabatic States", nstates);
}
