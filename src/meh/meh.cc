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
#include <src/dimer/dimer_prop.h>
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
  vector<vector<shared_ptr<Civec>>> collection_A(static_cast<int>(CS::MAX));
  vector<vector<shared_ptr<Civec>>> collection_B(static_cast<int>(CS::MAX));

  // Start by processing the singlet states
  state_inserter<0>(collection_A, CS::S, 0, 0);
  state_inserter<1>(collection_B, CS::S, 0, 0);

  // Now add in cation states
  if (cispace_->cations()) {
    state_inserter<0>(collection_A, CS::Ca, 0, -1);
    state_inserter<1>(collection_B, CS::Ca, 0, -1);
  }
  // dication states
  if (cispace_->dications()) {
    state_inserter<0>(collection_A, CS::dCa, 0, -2);
    state_inserter<1>(collection_B, CS::dCa, 0, -2);
  }

  // Now for the anion states
  if (cispace_->anions()) {
    state_inserter<0>(collection_A, CS::Aa, 1, 0);
    state_inserter<1>(collection_B, CS::Aa, 1, 0);
  }

  // dianion states
  if (cispace_->dianions()) {
    state_inserter<0>(collection_A, CS::dAa, 2, 0);
    state_inserter<1>(collection_B, CS::dAa, 2, 0);
  }

  // triplet states
  if (cispace_->triplets()) {
    state_inserter<0>(collection_A, CS::Ta, 1, -1);
    state_inserter<1>(collection_B, CS::Ta, 1, -1);
  }

  // quintets states
  if (cispace_->quintets()) {
    state_inserter<0>(collection_A, CS::Qaa, 2, -2);
    state_inserter<1>(collection_B, CS::Qaa, 2, -2);
  }

  // Package like civecs into Dvecs
  dimerstates_ = 0;

  vector<shared_ptr<Dvec>> dvecs_A(static_cast<int>(CS::MAX));
  vector<shared_ptr<Dvec>> dvecs_B(static_cast<int>(CS::MAX));

  for(int i = 0; i < static_cast<int>(CS::MAX); ++i) {
    if (!collection_A.at(i).empty()) dvecs_A.at(i) = shared_ptr<Dvec>(new Dvec(collection_A.at(i)));
    if (!collection_B.at(i).empty()) dvecs_B.at(i) = shared_ptr<Dvec>(new Dvec(collection_B.at(i)));
  }

  // First, singlets
  int ss = static_cast<int>(CS::S);
  subspaces_.emplace_back(dimerstates_, "  S", "  S", make_pair(dvecs_A.at(ss),dvecs_B.at(ss)));

  // Now, AC, if we've got em
  if (cispace_->anions() && cispace_->cations()) {
    int aa = static_cast<int>(CS::Aa);
    int cb = static_cast<int>(CS::Cb);
    subspaces_.emplace_back(dimerstates_, " Aa", " Cb", make_pair(dvecs_A.at(aa), dvecs_B.at(cb)));
    subspaces_.emplace_back(dimerstates_, " Ab", " Ca", make_pair(dvecs_A.at(aa+1), dvecs_B.at(cb-1)));
    subspaces_.emplace_back(dimerstates_, " Ca", " Ab", make_pair(dvecs_A.at(cb-1), dvecs_B.at(aa+1)));
    subspaces_.emplace_back(dimerstates_, " Cb", " Aa", make_pair(dvecs_A.at(cb), dvecs_B.at(aa)));
  }

  // Now, dAC, if we've got em
  if (cispace_->dianions() && cispace_->dications()) {
    int aa = static_cast<int>(CS::dAa);
    int cb = static_cast<int>(CS::dCb);
    subspaces_.emplace_back(dimerstates_, "dAa", "dCb", make_pair(dvecs_A.at(aa), dvecs_B.at(cb)));
    subspaces_.emplace_back(dimerstates_, "dA0", "dC0", make_pair(dvecs_A.at(aa+1), dvecs_B.at(cb-1)));
    subspaces_.emplace_back(dimerstates_, "dAb", "dCa", make_pair(dvecs_A.at(aa+2), dvecs_B.at(cb-2)));

    subspaces_.emplace_back(dimerstates_, "dCa", "dAb", make_pair(dvecs_A.at(cb-2), dvecs_B.at(aa+2)));
    subspaces_.emplace_back(dimerstates_, "dC0", "dA0", make_pair(dvecs_A.at(cb-1), dvecs_B.at(aa+1)));
    subspaces_.emplace_back(dimerstates_, "dCb", "dAa", make_pair(dvecs_A.at(cb), dvecs_B.at(aa)));
  }

  // Triplets
  if (cispace_->triplets()) {
    int ta = static_cast<int>(CS::Ta);
    int tb = static_cast<int>(CS::Tb);
    subspaces_.emplace_back(dimerstates_, " Ta", " Tb", make_pair(dvecs_A.at(ta), dvecs_B.at(tb)));
    subspaces_.emplace_back(dimerstates_, " T_", " T_", make_pair(dvecs_A.at(ta+1), dvecs_B.at(tb-1)));
    subspaces_.emplace_back(dimerstates_, " Tb", " Ta", make_pair(dvecs_A.at(ta+2), dvecs_B.at(tb-2)));
  }

  // Quintets
  if (cispace_->quintets()) {
    int qa = static_cast<int>(CS::Qaa);
    int qb = static_cast<int>(CS::Qbb);
    subspaces_.emplace_back(dimerstates_, "Qaa", "Qbb", make_pair(dvecs_A.at(qa), dvecs_B.at(qb)));
    subspaces_.emplace_back(dimerstates_, "Qa ", "Qb ", make_pair(dvecs_A.at(qa+1), dvecs_B.at(qb-1)));
    subspaces_.emplace_back(dimerstates_, "Q0 ", "Q0 ", make_pair(dvecs_A.at(qa+2), dvecs_B.at(qb-2)));
    subspaces_.emplace_back(dimerstates_, "Qb ", "Qa ", make_pair(dvecs_A.at(qa+3), dvecs_B.at(qb-3)));
    subspaces_.emplace_back(dimerstates_, "Qbb", "Qaa", make_pair(dvecs_A.at(qa+4), dvecs_B.at(qb-4)));
  }

  cispace_->complete();
}

const Coupling MultiExcitonHamiltonian::coupling_type(const DimerSubspace& AB, const DimerSubspace& ApBp) const {
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

  cout << "  o Computing diagonal blocks." << endl;
  for (auto& subspace : subspaces_) {
    hamiltonian_->add_block(subspace.offset(), subspace.offset(), compute_diagonal_block(subspace));
  }

  cout << "  o Computing off-diagonal blocks." << endl;
  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

      shared_ptr<Matrix> block = couple_blocks(*iAB, *jAB);

      hamiltonian_->add_block(ioff, joff, block);
      hamiltonian_->add_block(joff, ioff, block->transpose());
    }
  }

  cout << "  o Diagonalizing ME Hamiltonian" << endl;
  adiabats_ = shared_ptr<Matrix>(new Matrix(*hamiltonian_));
  energies_ = vector<double>(dimerstates_, 0.0);
  adiabats_->diagonalize(energies_.data());

  cout << "  o Constructing spin matrix." << endl;

  spin_ = shared_ptr<Matrix>(new Matrix(dimerstates_, dimerstates_));
  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

      shared_ptr<Matrix> spin_block = spin_couple_blocks(*iAB, *jAB);

      spin_->add_block(ioff, joff, spin_block);
      spin_->add_block(joff, ioff, spin_block->transpose());
    }
    spin_->add_block(ioff, ioff, compute_diagonal_spin_block(*iAB));
  }


  spinadiabats_ = shared_ptr<Matrix>(new Matrix( (*adiabats_) % (*spin_) * (*adiabats_) ));

  cout << "  o Computing properties" << endl;
  DimerDipole dipole = DimerDipole(ref_, dimerclosed_, dimerclosed_ + nact_.first, dimerclosed_ + dimeractive_, coeff_);
  array<string,3> mu_labels = {{"x", "y", "z"}};
  for (int i = 0; i < 3; ++i) {
    string label("mu_");
    label += mu_labels[i];
    shared_ptr<Matrix> tmp = compute_1e_prop(dipole.dipoles<0>(i), dipole.dipoles<1>(i), dipole.cross_dipole(i), dipole.core_dipole(i));
    shared_ptr<Matrix> prop(new Matrix( (*adiabats_) % (*tmp) * (*adiabats_) ));
    properties_.push_back(make_pair(label, prop));
  }
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_1e_prop(shared_ptr<const Matrix> hAA, shared_ptr<const Matrix> hBB, shared_ptr<const Matrix> hAB, const double core) const {

  shared_ptr<Matrix> out(new Matrix(dimerstates_, dimerstates_));

  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

      shared_ptr<Matrix> out_block = compute_offdiagonal_1e(*iAB, *jAB, hAB);

      out->add_block(ioff, joff, out_block);
      out->add_block(joff, ioff, out_block->transpose());
    }
    out->add_block(ioff, ioff, compute_diagonal_1e(*iAB, hAA->data(), hBB->data(), core));
  }

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_diagonal_1e(const DimerSubspace& AB, const double* hAA, const double* hBB, const double diag) const {
  shared_ptr<const Dvec> ccvecA = AB.ci<0>();
  shared_ptr<const Dvec> ccvecB = AB.ci<1>();
  shared_ptr<const Dvec> sigmavecA = form_sigma_1e(ccvecA, hAA);
  shared_ptr<const Dvec> sigmavecB = form_sigma_1e(ccvecB, hBB);

  const int nstatesA = AB.nstates<0>();
  const int nstatesB = AB.nstates<1>();
  const int dimerstates = AB.dimerstates();

  shared_ptr<Matrix> out(new Matrix(dimerstates, dimerstates));

  for(int stateAp = 0; stateAp < nstatesA; ++stateAp) {
    for(int stateBp = 0; stateBp < nstatesB; ++stateBp) {
      const int stateApBp = AB.dimerindex(stateAp,stateBp);
      // hAA
      for(int stateA = 0; stateA < nstatesA; ++stateA) {
        // stateB = stateBp
        const int stateAB = AB.dimerindex(stateA,stateBp);
        out->element(stateAB, stateApBp) += sigmavecA->data(stateA)->ddot(*ccvecA->data(stateAp));
      }

      // hBB
      for(int stateB = 0; stateB < nstatesB; ++stateB) {
        // stateA = stateAp
        const int stateAB = AB.dimerindex(stateAp, stateB);
        out->element(stateAB, stateApBp) += sigmavecB->data(stateB)->ddot(*ccvecB->data(stateBp));
      }
    }
  }

  out->add_diag(diag);

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_offdiagonal_1e(const DimerSubspace& AB, const DimerSubspace& ApBp, shared_ptr<const Matrix> hAB) const {
  Coupling term_type = coupling_type(AB,ApBp);

  shared_ptr<Quantization> operatorA;
  shared_ptr<Quantization> operatorB;
  int neleA = AB.ci<0>()->det()->nelea() + AB.ci<0>()->det()->neleb();

  switch(term_type) {
    case Coupling::aET :
      operatorA = shared_ptr<Quantization>(new OneBody<SQ::CreateAlpha>());
      operatorB = shared_ptr<Quantization>(new OneBody<SQ::AnnihilateAlpha>());
      break;
    case Coupling::inv_aET :
      operatorA = shared_ptr<Quantization>(new OneBody<SQ::AnnihilateAlpha>());
      operatorB = shared_ptr<Quantization>(new OneBody<SQ::CreateAlpha>());
      --neleA;
      break;
    case Coupling::bET :
      operatorA = shared_ptr<Quantization>(new OneBody<SQ::CreateBeta>());
      operatorB = shared_ptr<Quantization>(new OneBody<SQ::AnnihilateBeta>());
      break;
    case Coupling::inv_bET :
      operatorA = shared_ptr<Quantization>(new OneBody<SQ::AnnihilateBeta>());
      operatorB = shared_ptr<Quantization>(new OneBody<SQ::CreateBeta>());
      --neleA;
      break;
    default :
      return shared_ptr<Matrix>(new Matrix(AB.dimerstates(), ApBp.dimerstates()));
  }

  Matrix gamma_A = *form_gamma(AB.ci<0>(), ApBp.ci<0>(), operatorA);
  Matrix gamma_B = *form_gamma(AB.ci<1>(), ApBp.ci<1>(), operatorB);
  Matrix tmp = gamma_A * (*hAB) ^ gamma_B;

  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));
  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<0>());

  if ( (neleA % 2) == 1 ) out->scale(-1.0);

  return out;
}


shared_ptr<Matrix> MultiExcitonHamiltonian::form_gamma(shared_ptr<const Dvec> ccvecA, shared_ptr<const Dvec> ccvecAp, shared_ptr<Quantization> action) const { 
  const int nstatesA = ccvecA->ij();
  const int nstatesAp = ccvecAp->ij();

  shared_ptr<const Determinants> detA = ccvecA->det();
  const int norb = detA->norb();
  const int ij = action->ij(norb);

  Matrix tmp(ij, nstatesA*nstatesAp);

  double *edata = tmp.data();

  for(int state = 0; state < nstatesA; ++state) {
    shared_ptr<Dvec> c = action->compute(ccvecA->data(state));

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

void MultiExcitonHamiltonian::print_hamiltonian(const string title, const int nstates) const {
  hamiltonian_->print(title, nstates);
}

void MultiExcitonHamiltonian::print_adiabats(const double thresh, const string title, const int nstates) const {
  const int end = min(nstates, dimerstates_);
  cout << endl << " ===== " << title << " =====" << endl;
  for (int istate = 0; istate < end; ++istate) {
    cout << "   state  " << setw(3) << istate << ": " << setprecision(12) << setw(16) << energies_.at(istate) << ", <S^2> = " << setprecision(4) << setw(6) << spinadiabats_->element(istate,istate) << endl;
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

void MultiExcitonHamiltonian::print_property(const string label, shared_ptr<const Matrix> property , const int nstates) const {
  const string indent("   ");
  const int nprint = min(nstates, property->ndim());

  cout << indent << " " << label << "    |0>";
  for (int istate = 1; istate < nprint; ++istate) cout << "         |" << istate << ">";
  cout << endl;
  for (int istate = 0; istate < nprint; ++istate) {
    cout << indent << "<" << istate << "|";
    for (int jstate = 0; jstate < nprint; ++jstate) {
      cout << setw(12) << setprecision(6) << property->element(jstate, istate);
    }   
    cout << endl;
  }
  cout << endl;
}

void MultiExcitonHamiltonian::print(const int nstates, const double thresh) const {
  print_adiabats(thresh, "Adiabatic States", nstates);
  for (auto& prop : properties_) print_property(prop.first, prop.second, nstates);
}
