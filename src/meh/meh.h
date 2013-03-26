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

namespace bagel {

/************************************************************************************
*  This class computes the Hamiltonian matrix for multiexciton states described by  *
* CAS calculations                                                                  *
************************************************************************************/

typedef std::shared_ptr<Matrix> MatrixPtr;

enum class Coupling {
  none = 0,
  diagonal = 1, // Probably won't be used
  alphaET = 2,
  inv_alphaET = -2,
  betaET = 3,
  inv_betaET = -3,
  ABflip = 4,
  BAflip = -4,
  alphabetaET = 5,
  inv_alphabetaET = -5
};

enum class ChargeSpin {
  SS = 0,
  T0T0 = 0,
  AaCb = 1,
  AbCa = 2,
  CaAb = 3,
  CbAa = 4,
  TaTb = 5,
  TbTa = 6
};

// What started off as a simple structure is now becoming a bonafide helper class
class DimerSubspace {
  protected:
    const ChargeSpin chargespin_;

    const int offset_;
    const int nstatesA_;
    const int nstatesB_;

    std::pair<std::shared_ptr<Dvec>, std::shared_ptr<Dvec>> ci_;

  public:
    DimerSubspace(const ChargeSpin _cs, const int _offset, std::pair<std::shared_ptr<Dvec>, std::shared_ptr<Dvec>> _ci) : 
      chargespin_(_cs), offset_(_offset), nstatesA_(_ci.first->ij()), nstatesB_(_ci.second->ij()), ci_(_ci) {}

    const ChargeSpin chargespin() const { return chargespin_; }
    const int offset() const { return offset_; }
    const int dimerstates() const { return nstatesA_ * nstatesB_; }
    const int dimerindex(const int iA, const int iB) { return (iA + iB*nstatesA_); }

    template <int unit> const int nstates() { return ( unit == 0 ? nstatesA_ : nstatesB_ ); }
    template <int unit> std::shared_ptr<const Dvec> ci() { return ( unit == 0 ? ci_.first : ci_.second ); }
    
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

      std::vector<double> energies_; // Adiabatic energies

      // Total system quantities
      const int dimerbasis_;
      const int dimerclosed_;
      const int dimeractive_;
      int dimerstates_;

      // Static variables (for a little extra clarity)
      static const int Alpha = 0;
      static const int Beta = 1;
      static const int Annihilate = 0;
      static const int Create = 1;

      // Localized quantities
      std::pair<const int, const int> nact_;
      std::pair<const int, const int> nbasis_;
      std::pair<const int, const int> nstates_;

   public:
      MultiExcitonHamiltonian(std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerCISpace> cispace);

      int dimerstate(const int A, const int B) const { return (A + B*nstates_.first); };

      void compute();

      void print_hamiltonian(const std::string title = "MultiExciton Hamiltonian", const int nstates = 10);
      void print_energies(const std::string title = "Adiabatic state energies", const int nstates = 10);
      void print_adiabats(const std::string title = "Adiabats", const int nstates = 10);
      void print(const int nstates = 10);

      Coupling coupling_type(DimerSubspace& AB, DimerSubspace& ApBp);

   private:
      void common_init();
      void reorder_matrix(const double* source, double* target, const int nA, const int nAp, const int nB, const int nBp) const;
      void reorder_matrix(const double* source, double* target, const int nA, const int nB) { reorder_matrix(source,target,nA,nA,nB,nB); }

      template <int A, int B, int C, int D> std::pair<int, int> index(int a, int b, int c, int d) const {
        int iA = 0, jB = 0;
        if (A == 0) iA += a; else jB += a;
        if (B == 0) iA = b + iA*nact_.first; else jB = b + jB*nact_.second;
        if (C == 0) iA = c + iA*nact_.first; else jB = c + jB*nact_.second;
        if (D == 0) iA = d + iA*nact_.first; else jB = d + jB*nact_.second;
        return std::make_pair(iA,jB);
      }

      template <int unit> int active(int a) const { return (a + unit*nact_.first); }

      const int coupling_index(std::pair<int,int> AT, std::pair<int,int> BT) const {
        return coupling_index(AT.first,AT.second,BT.first,BT.second);
      }
      const int coupling_index(const int a, const int b, const int c, const int d) const {
        return (a + b*large__ + c*large__*large__ + d*large__*large__*large__);
      }

      // Diagonal block stuff
      MatrixPtr compute_diagonal_block(DimerSubspace& subspace);

      MatrixPtr compute_closeclose(DimerSubspace& subspace);
      MatrixPtr compute_closeactive(DimerSubspace& subspace);
      MatrixPtr compute_intra_activeactive(DimerSubspace& subspace);
      MatrixPtr compute_inter_activeactive(DimerSubspace& subspace);

      std::shared_ptr<Dvec> form_sigma_1e(std::shared_ptr<const Dvec> ccvec, double* hdata) const;
      std::shared_ptr<Dvec> form_sigma_2e(std::shared_ptr<const Dvec> ccvec, double* mo2e_ptr) const;

      void sigma_2aa(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, double* mo2e_ptr, const int nact) const;
      void sigma_2bb(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, double* mo2e_ptr, const int nact) const;
      void sigma_2ab_1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d, const int nact) const;
      void sigma_2ab_2(std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e, double* mo2e_ptr) const;
      void sigma_2ab_3(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e, const int nact) const;
      
      // gamma = < A' | a^\dagger c | A >
      template<int spin> MatrixPtr form_gamma(std::shared_ptr<const Dvec> ccvec) const { return form_gamma<spin,spin>(ccvec,ccvec); }
      template<int spin1, int spin2> MatrixPtr form_gamma(std::shared_ptr<const Dvec> ccvecA, std::shared_ptr<const Dvec> ccvecAp) const;

      // Off-diagonal stuff
      MatrixPtr couple_blocks(DimerSubspace& AB, DimerSubspace& ApBp); // Off-diagonal driver

      MatrixPtr compute_alphaET(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_betaET(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_ABflip(DimerSubspace& AB, DimerSubspace& ApBp);
      MatrixPtr compute_alphabetaET(DimerSubspace& AB, DimerSubspace& ApBp);

      // Helper functions
      template<int A, int B, int C, int D> std::pair<MatrixPtr, MatrixPtr> form_JKmatrices() const;

      // Just these for now, add them as I need them
      template<int spin, int oper> std::shared_ptr<Dvec> operator_q(std::shared_ptr<const Civec>) const;
      template<int spin1, int spin2> std::shared_ptr<Dvec> operator_ca(std::shared_ptr<const Civec>) const;
};

template<int spin1, int spin2>
std::shared_ptr<Matrix> MultiExcitonHamiltonian::form_gamma(std::shared_ptr<const Dvec> ccvecA, std::shared_ptr<const Dvec> ccvecAp) const {
  const int nstatesA = ccvecA->ij();
  const int nstatesAp = ccvecAp->ij();

  std::shared_ptr<const Determinants> detA = ccvecA->det();
  const int norb = detA->norb();
  const int ij = norb * norb;

  Matrix tmp(ij, nstatesA*nstatesAp);

  double *edata = tmp.data();

  for(int state = 0; state < nstatesA; ++state) {
    std::shared_ptr<Dvec> c = operator_ca<spin1,spin2>(ccvecA->data(state));

    // | C > ^A_ac is done
    for(int statep = 0; statep < nstatesAp; ++statep) {
      for(int ac = 0; ac < ij; ++ac, ++edata) {
        *edata = c->data(ac)->ddot(*ccvecAp->data(statep)); 
      }   
    }   
  }

  return tmp.transpose();
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

template<int spin, int oper> std::shared_ptr<Dvec> MultiExcitonHamiltonian::operator_q(std::shared_ptr<const Civec> ccvec) const {
  std::shared_ptr<const Determinants> source_det = ccvec->det();
  std::shared_ptr<const Determinants> target_det = ( (spin == Alpha) ? 
    (oper == Annihilate ? ccvec->det()->remalpha() : ccvec->det()->addalpha() ) : 
    (oper == Annihilate ? ccvec->det()->rembeta() : ccvec->det()->addbeta() ));

  const int norb = target_det->norb();

  const int source_lena = source_det->lena();
  const int source_lenb = source_det->lenb();

  const int target_lena = target_det->lena();
  const int target_lenb = target_det->lenb();

  const int length = (spin == Alpha ? target_lenb : target_lena);
  const int source_start = ( spin == Alpha ? source_lenb : 1 );
  const int target_start = ( spin == Alpha ? target_lenb : 1 );
  const int source_stride = ( spin == Alpha ? 1 : source_lenb );
  const int target_stride = ( spin == Alpha ? 1 : target_lenb );

  std::shared_ptr<Dvec> out(new Dvec(target_det, norb));

  const double* source_base = ccvec->data();

  for (int i = 0; i < norb; ++i) {
    double* target_base = out->data(i)->data();

    for(auto& iter : (spin == Alpha ? (oper == Annihilate ? source_det->phidowna(i) : source_det->phiupa(i)) : (oper == Annihilate ? source_det->phidownb(i) : source_det->phiupb(i)) ) ) {
      const double sign = static_cast<double>( iter.sign );
      double* target = target_base + target_start * iter.target;
      const double* source = source_base + source_start * iter.source;
      daxpy_(length, sign, source, source_stride, target, target_stride);
    }
  }

  return out;
}

template<int spin1, int spin2> std::shared_ptr<Dvec> MultiExcitonHamiltonian::operator_ca(std::shared_ptr<const Civec> ccvec) const {
  if (spin1 == spin2) {
    std::shared_ptr<const Determinants> det = ccvec->det();
    const int norb = det->norb();

    const int lena = det->lena();
    const int lenb = det->lenb();

    const int length = (spin1 == Alpha ? lenb : lena);
    const int start = ( spin1 == Alpha ? lenb : 1 );
    const int stride = ( spin1 == Alpha ? 1 : lenb );

    const int sizeij = norb * norb;
    std::shared_ptr<Dvec> out(new Dvec(det, sizeij));

    const double* source_base = ccvec->data();

    for (int ij = 0; ij < sizeij; ++ij) {
      double* target_base = out->data(ij)->data();
      for(auto& iter : (spin1==Alpha ? det->phia(ij) : det->phib(ij))) {
        const double sign = static_cast<double>(iter.sign);
        double* target = target_base + start * iter.target;
        const double* source = source_base + start * iter.source;
        daxpy_(length, sign, source, stride, target, stride);
      }
    }

    return out;
  }
  else {
    const int norb = ccvec->det()->norb();

    std::shared_ptr<Dvec> intermediate = operator_q<spin2,Annihilate>(ccvec);

    std::vector<std::shared_ptr<Dvec>> tmp_vec;
    for (int i = 0; i < norb; ++i) {
      tmp_vec.push_back(operator_q<spin1,Create>(intermediate->data(i)));
    }

    std::shared_ptr<Dvec> out(new Dvec(tmp_vec.front()->det(), norb*norb));
    int ij = 0;
    for (int i = 0; i < norb; ++i) {
      for (int j = 0; j < norb; ++j, ++ij) {
        out->data(ij) = tmp_vec.at(j)->data(i);
      }
    }

    return out;
  }
}

//template<int spin> std::shared_ptr<Dvec> operator_aa(std::shared_ptr<const Civec>);
//template<int spin> std::shared_ptr<Dvec> operator_ac(std::shared_ptr<const Civec>);
//template<int spin> std::shared_ptr<Dvec> operator_cc(std::shared_ptr<const Civec>);

}

#endif
