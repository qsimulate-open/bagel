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
#include <src/util/lexical_cast.h>
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

  // Process DimerCISpace to form and organize needed Civecs
  vector<vector<shared_ptr<Civec>>> collection_A(8);
  vector<vector<shared_ptr<Civec>>> collection_B(8);

  // Start by processing the singlet states
  state_inserter(collection_A, CS::S, cispace_->ccvec<0>(0,0));
  state_inserter(collection_B, CS::S, cispace_->ccvec<1>(0,0));

  // Now add in cation states
  if (cispace_->cations()) {
    state_inserter(collection_A, CS::Ca, cispace_->ccvec<0>(0,-1));
    auto flippedA = cispace_->ccvec<0>(0,-1)->spinflip(cispace_->add_det<0>(-1,0));
    cispace_->insert<0>(flippedA);
    state_inserter(collection_A, CS::Cb, flippedA);

    state_inserter(collection_B, CS::Ca, cispace_->ccvec<1>(0,-1));
    auto flippedB = cispace_->ccvec<1>(0,-1)->spinflip(cispace_->add_det<1>(-1,0));
    cispace_->insert<1>(flippedB);
    state_inserter(collection_B, CS::Cb, flippedB);
  }

  // Now for the anion states
  if (cispace_->anions()) {
    state_inserter(collection_A, CS::Aa, cispace_->ccvec<0>(1,0));
    auto flippedA = cispace_->ccvec<0>(1,0)->spinflip(cispace_->add_det<0>(0,1));
    cispace_->insert<0>(flippedA);
    state_inserter(collection_A, CS::Ab, flippedA);

    state_inserter(collection_B, CS::Aa, cispace_->ccvec<1>(1,0));
    auto flippedB = cispace_->ccvec<1>(1,0)->spinflip(cispace_->add_det<1>(0,1));
    cispace_->insert<1>(flippedB);
    state_inserter(collection_B, CS::Ab, flippedB);
  }

  // And finally triplet states
  if (cispace_->triplets()) {
    {
      state_inserter(collection_A, CS::Ta, cispace_->ccvec<0>(1,-1));

      auto T0A = cispace_->ccvec<0>(1,-1)->spin_lower(cispace_->add_det<0>(0,0));
      cispace_->insert<0>(T0A);
      state_inserter(collection_A, CS::T0, T0A);

      auto TbA = cispace_->ccvec<0>(1,-1)->spinflip(cispace_->add_det<0>(-1,1));
      cispace_->insert<0>(TbA);
      state_inserter(collection_A, CS::Tb, TbA);
    }
    {
      state_inserter(collection_B, CS::Ta, cispace_->ccvec<1>(1,-1));

      auto T0B = cispace_->ccvec<1>(1,-1)->spin_lower(cispace_->add_det<1>(0,0));
      cispace_->insert<1>(T0B);
      state_inserter(collection_B, CS::T0, T0B);

      auto TbB = cispace_->ccvec<1>(1,-1)->spinflip(cispace_->add_det<1>(-1,1));
      cispace_->insert<1>(TbB);
      state_inserter(collection_B, CS::Tb, TbB);
    }
  }

  // Package like civecs into Dvecs
  dimerstates_ = 0;

  // First, singlets
  shared_ptr<Dvec> SA(new Dvec(collection_A.at(static_cast<int>(CS::S))));
  shared_ptr<Dvec> SB(new Dvec(collection_B.at(static_cast<int>(CS::S))));

  subspaces_.push_back(DimerSubspace(dimerstates_, " S", " S", make_pair(SA,SB)));

  // Now, AC, if we've got em
  if (cispace_->anions() && cispace_->cations()) {
    shared_ptr<Dvec> Aa_A(new Dvec(collection_A.at(static_cast<int>(CS::Aa))));
    shared_ptr<Dvec> Ab_A(new Dvec(collection_A.at(static_cast<int>(CS::Ab))));
    shared_ptr<Dvec> Aa_B(new Dvec(collection_B.at(static_cast<int>(CS::Aa))));
    shared_ptr<Dvec> Ab_B(new Dvec(collection_B.at(static_cast<int>(CS::Ab))));

    shared_ptr<Dvec> Ca_A(new Dvec(collection_A.at(static_cast<int>(CS::Ca))));
    shared_ptr<Dvec> Cb_A(new Dvec(collection_A.at(static_cast<int>(CS::Cb))));
    shared_ptr<Dvec> Ca_B(new Dvec(collection_B.at(static_cast<int>(CS::Ca))));
    shared_ptr<Dvec> Cb_B(new Dvec(collection_B.at(static_cast<int>(CS::Cb))));

    subspaces_.emplace_back(dimerstates_, "Aa", "Cb", make_pair(Aa_A, Cb_B));
    subspaces_.emplace_back(dimerstates_, "Ab", "Ca", make_pair(Ab_A, Ca_B));
    subspaces_.emplace_back(dimerstates_, "Ca", "Ab", make_pair(Ca_A, Ab_B));
    subspaces_.emplace_back(dimerstates_, "Cb", "Aa", make_pair(Cb_A, Aa_B));
  }

  if (cispace_->triplets()) {
    shared_ptr<Dvec> Ta_A(new Dvec(collection_A.at(static_cast<int>(CS::Ta))));
    shared_ptr<Dvec> T0_A(new Dvec(collection_A.at(static_cast<int>(CS::T0))));
    shared_ptr<Dvec> Tb_A(new Dvec(collection_A.at(static_cast<int>(CS::Tb))));

    shared_ptr<Dvec> Ta_B(new Dvec(collection_B.at(static_cast<int>(CS::Ta))));
    shared_ptr<Dvec> T0_B(new Dvec(collection_B.at(static_cast<int>(CS::T0))));
    shared_ptr<Dvec> Tb_B(new Dvec(collection_B.at(static_cast<int>(CS::Tb))));

    subspaces_.emplace_back(dimerstates_, "Ta", "Tb", make_pair(Ta_A, Tb_B));
    subspaces_.emplace_back(dimerstates_, "T_", "T_", make_pair(T0_A, T0_B));
    subspaces_.emplace_back(dimerstates_, "Tb", "Ta", make_pair(Tb_A, Ta_B));
  }

  cispace_->complete();
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
        const int BBp = Bp + nBp * B;
        for (int A = 0; A < nA; ++A) {
          const int AAp = Ap + nAp * A;
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

void MultiExcitonHamiltonian::print_adiabats(const double thresh, const string title, const int nstates) {
  const int end = std::min(nstates, dimerstates_);
  cout << endl << " ===== " << title << " =====" << endl;
  for (int istate = 0; istate < end; ++istate) {
    cout << "   state  " << setw(3) << istate << ": " << setprecision(12) << setw(16) << energies_.at(istate) << endl;
    double *eigendata = adiabats_->element_ptr(0,istate);
    for(auto& subspace : subspaces_) {
      const int nA = subspace.nstates<0>();
      const int nB = subspace.nstates<1>();
      for(int i = 0; i < nA; ++i) {
        for(int j = 0; j < nB; ++j, ++eigendata) {
          if ( (*eigendata)*(*eigendata) > thresh ) {
            cout << "      " << subspace.string(i,j) << setprecision(12) << setw(20) << *eigendata << endl;
          }
        }
      }
    }
    cout << endl;
  }
}

void MultiExcitonHamiltonian::print(const int nstates, const double thresh) {
  print_adiabats(thresh, "Adiabatic States", nstates);
}
