//
// BAGEL - Parallel electron correlation program.
// Filename: meh.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __meh_meh_h
#define __meh_meh_h

#include <utility>

#include <src/scf/coeff.h>
#include <src/wfn/reference.h>
#include <src/dimer/dimer.h>
#include <src/dimer/dimer_cispace.h>
#include <src/dimer/dimer_jop.h>
#include <src/fci/dvec.h>
#include <src/util/matrix.h>
#include <src/meh/quantization.h>

namespace bagel {

/************************************************************************************
*  This class computes the Hamiltonian matrix for multiexciton states described by  *
* CAS calculations                                                                  *
************************************************************************************/

using MatrixPtr = std::shared_ptr<Matrix>;

enum class Coupling {
  none = 0,
  diagonal = 1,
  aET = 2,
  inv_aET = -2,
  bET = 3,
  inv_bET = -3,
  abFlip = 4,
  baFlip = -4,
  abET = 5,
  inv_abET = -5,
  aaET = 6,
  inv_aaET = -6,
  bbET = 7,
  inv_bbET = -7
};

enum class CS {
  S, Aa, Ab, dAa, dA0, dAb, Ca, Cb, dCa, dC0, dCb, Ta, T0, Tb, Qaa, Qa, Q0, Qb, Qbb, MAX
};

class DimerSubspace {
  protected:
    const int offset_;
    const int nstatesA_;
    const int nstatesB_;
    const std::string stringA_;
    const std::string stringB_;

    std::pair<std::shared_ptr<Dvec>, std::shared_ptr<Dvec>> ci_;

  public:
    DimerSubspace(int& _offset, const std::string _stringA, const std::string _stringB, std::pair<std::shared_ptr<Dvec>, std::shared_ptr<Dvec>> _ci) : 
      offset_(_offset), nstatesA_(_ci.first->ij()), nstatesB_(_ci.second->ij()), stringA_(_stringA), stringB_(_stringB),
       ci_(_ci) { _offset += dimerstates(); }

    const int offset() const { return offset_; }
    const int dimerstates() const { return nstatesA_ * nstatesB_; }
    const int dimerindex(const int iA, const int iB) const { return (iA + iB*nstatesA_); }
    const std::string string(const int i, const int j) const {
      std::string out = stringA_ + lexical_cast<std::string>(i) + std::string(" ") + stringB_ + lexical_cast<std::string>(j);
      return out;
    }

    template <int unit> const int nstates() const { return ( unit == 0 ? nstatesA_ : nstatesB_ ); }
    template <int unit> std::shared_ptr<const Dvec> ci() const { return ( unit == 0 ? ci_.first : ci_.second ); }
    
};

class MultiExcitonHamiltonian {
   protected:
      std::shared_ptr<const Dimer> dimer_;
      std::shared_ptr<const Reference> ref_;
      std::shared_ptr<const Coeff> coeff_;

      std::shared_ptr<DimerJop> jop_;
      std::shared_ptr<DimerCISpace> cispace_;

      std::vector<DimerSubspace> subspaces_;

      MatrixPtr hamiltonian_;
      MatrixPtr adiabats_; // Eigenvectors of adiabatic states
      MatrixPtr spin_; //S^2 matrix
      MatrixPtr spinadiabats_; // S^2 matrix transformed into adiabatic basis
      std::vector<std::pair<std::string, MatrixPtr>> properties_;

      std::vector<double> energies_; // Adiabatic energies

      // Total system quantities
      const int dimerbasis_;
      const int dimerclosed_;
      const int dimeractive_;
      int dimerstates_;

      // Localized quantities
      std::pair<const int, const int> nact_;
      std::pair<const int, const int> nbasis_;
      std::pair<const int, const int> nstates_;

   public:
      MultiExcitonHamiltonian(std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerCISpace> cispace);

      int dimerstate(const int A, const int B) const { return (A + B*nstates_.first); };

      void compute();

      void print_hamiltonian(const std::string title = "MultiExciton Hamiltonian", const int nstates = 10) const;
      void print_adiabats(const double thresh = 0.05, const std::string title = "Adiabats", const int nstates = 10) const;
      void print_property(const std::string label, std::shared_ptr<const Matrix>, const int size = 10) const ;
      void print(const int nstates = 10, const double thresh = 0.0001) const;

      const Coupling coupling_type(const DimerSubspace& AB, const DimerSubspace& ApBp) const;

   private:
      void common_init();
      void reorder_matrix(const double* source, double* target, const int nA, const int nAp, const int nB, const int nBp) const;

      template <int A, int B, int C, int D> std::pair<int, int> index(int a, int b, int c, int d) const {
        int iA = 0, jB = 0;
        if (A == 0) iA += a; else jB += a;
        if (B == 0) iA = b + iA*nact_.first; else jB = b + jB*nact_.second;
        if (C == 0) iA = c + iA*nact_.first; else jB = c + jB*nact_.second;
        if (D == 0) iA = d + iA*nact_.first; else jB = d + jB*nact_.second;
        return std::make_pair(iA,jB);
      }

      template <int unit> int active(int a) const { return (a + unit*nact_.first); }

      int coupling_index(std::pair<int,int> AT, std::pair<int,int> BT) const { return coupling_index(AT.first,AT.second,BT.first,BT.second); }
      int coupling_index(const int a, const int b, const int c, const int d) const {
        return (a + b*large__ + c*large__*large__ + d*large__*large__*large__);
      }

      MatrixPtr compute_1e_prop(std::shared_ptr<const Matrix> hAA, std::shared_ptr<const Matrix> hBB, std::shared_ptr<const Matrix> hAB, const double core) const;
      MatrixPtr compute_offdiagonal_1e(const DimerSubspace& AB, const DimerSubspace& ApBp, std::shared_ptr<const Matrix> hAB) const;
      MatrixPtr compute_diagonal_1e(const DimerSubspace& subspace, const double* hAA, const double* hBB, const double diag) const;

      // Diagonal block stuff
      MatrixPtr compute_diagonal_block(DimerSubspace& subspace);
      MatrixPtr compute_diagonal_spin_block(DimerSubspace& subspace);

      MatrixPtr compute_intra_2e(DimerSubspace& subspace);

      std::shared_ptr<Dvec> form_sigma_1e(std::shared_ptr<const Dvec> ccvec, const double* hdata) const;
      std::shared_ptr<Dvec> form_sigma_2e(std::shared_ptr<const Dvec> ccvec, const double* mo2e_ptr) const;

      void sigma_2aa(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, const double* mo2e_ptr, const int nact) const;
      void sigma_2bb(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, const double* mo2e_ptr, const int nact) const;
      void sigma_2ab_1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d, const int nact) const;
      void sigma_2ab_2(std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e, const double* mo2e_ptr) const;
      void sigma_2ab_3(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e, const int nact) const;
      
      // Helper functions
      template<int unit>
      void state_inserter(std::vector<std::vector<std::shared_ptr<Civec>>>& ccvec, const CS cs, const int qa, const int qb);

      template<int A, int B, int C, int D> std::pair<MatrixPtr, MatrixPtr> form_JKmatrices() const;
      MatrixPtr form_gamma(std::shared_ptr<const Dvec> ccvecA, std::shared_ptr<const Dvec> ccvecAp, std::shared_ptr<Quantization> action) const;

      // Off-diagonal stuff
      MatrixPtr couple_blocks(DimerSubspace& AB, DimerSubspace& ApBp); // Off-diagonal driver for H
      MatrixPtr spin_couple_blocks(DimerSubspace& AB, DimerSubspace& ApBp); // Off-diagonal driver for S^2

      MatrixPtr compute_inter_2e(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_aET(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_bET(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_abFlip(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_abET(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_aaET(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_bbET(DimerSubspace& AB, DimerSubspace& ApBp);
};

template<int unit>
void MultiExcitonHamiltonian::state_inserter(std::vector<std::vector<std::shared_ptr<Civec>>>& ccvec, const CS cs, const int qa, const int qb) {
  const double thresh = 1.0e-3;
  std::shared_ptr<Dvec> currentvec = cispace_->ccvec<unit>(qa,qb);

  const int nspin = currentvec->det()->nspin();
  const int mult = nspin + 1;

  // assuming high spin, expectation value should be
  double expectation = static_cast<double>(nspin) * 0.5;
  expectation = expectation * (expectation + 1.0);

  const int cs_int = static_cast<int>(cs);
  const int size = currentvec->ij();
  for(int i = 0; i < size; ++i) {
    if ( std::abs(expectation - currentvec->data(i)->spin_expectation()) < thresh ) {
      ccvec.at(cs_int).push_back(currentvec->data(i));
    }
    else {std::cout << "Removing one state of type (qa,qb) = (" << qa << "," << qb << ") due to spin contamination." << std::endl;}
  }


  for(int i = 1; i < mult; ++i) {
    currentvec = currentvec->spin_lower(cispace_->add_det<unit>(qa-i,qb+i));
    cispace_->insert<unit>(currentvec);
    const int cs_position = static_cast<int>(cs) + i;
    for(int j = 0; j < size; ++j) {
      const double norm = currentvec->data(j)->norm();
      if ( norm < numerical_zero__ ) throw std::runtime_error("Spin lowering operator yielded no state.");
      currentvec->data(j)->scale(1.0/norm);
      if ( std::abs(expectation - currentvec->data(j)->spin_expectation()) < thresh ) {
        ccvec.at(cs_position).push_back(currentvec->data(j));
      }
    }
  }
}

template<int A, int B, int C, int D>
std::pair<MatrixPtr, MatrixPtr> MultiExcitonHamiltonian::form_JKmatrices() const {
  const int nactA = nact_.first;
  const int nactB = nact_.second;

  int ijA = 1;
  int unitA = 4 - (A + B + C + D);
  for ( int i = 0; i < unitA; ++i ) ijA *= nactA;

  int ijB = 1;
  int unitB = A + B + C + D;
  for ( int i = 0; i < unitB; ++i ) ijB *= nactB;

  MatrixPtr Jout(new Matrix(ijA, ijB));
  MatrixPtr Kout(new Matrix(ijA, ijB));

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

}

#endif
