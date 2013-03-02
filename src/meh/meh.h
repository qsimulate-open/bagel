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
#include <src/fci/dvec.h>
#include <src/util/matrix.h>

namespace bagel {

/************************************************************************************
*  This class computes the Hamiltonian matrix for multiexciton states described by  *
* CAS calculations                                                                  *
************************************************************************************/

typedef std::shared_ptr<Matrix> MatrixPtr;

class MultiExcitonHamiltonian {
   protected:
      std::shared_ptr<const Dimer> dimer_;
      std::shared_ptr<const Reference> ref_;
      std::shared_ptr<const Coeff> coeff_;

      std::shared_ptr<DimerJop> jop_;
      std::shared_ptr<DimerCISpace> cispace_;

      MatrixPtr hamiltonian_;

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

   private:
      void common_init();

      template <int A, int B, int C, int D> std::pair<int, int> index(int a, int b, int c, int d) const {
        int iA = 0, jB = 0;
        if (A == 0) iA += a; else jB += a;
        if (B == 0) iA = b + iA*nact_.first; else jB = b + jB*nact_.second;
        if (C == 0) iA = c + iA*nact_.first; else jB = c + jB*nact_.second;
        if (D == 0) iA = d + iA*nact_.first; else jB = d + jB*nact_.second;
        return std::make_pair(iA,jB);
      }

      template <int unit> int active(int a) const { return (a + unit*nact_.first); }

      MatrixPtr compute_closeclose();
      MatrixPtr compute_closeactive();
      MatrixPtr compute_intra_activeactive();
      MatrixPtr compute_inter_activeactive();

      std::shared_ptr<Dvec> form_sigma_1e(std::shared_ptr<const Dvec> ccvec, double* hdata, const int ij) const;
      std::shared_ptr<Dvec> form_sigma_2e(std::shared_ptr<const Dvec> ccvec, double* mo2e_ptr, const int nact) const;

      void sigma_2aa(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, double* mo2e_ptr, const int nact) const;
      void sigma_2bb(std::shared_ptr<const Civec> cc, std::shared_ptr<Civec> sigma, double* mo2e_ptr, const int nact) const;
      void sigma_2ab_1(std::shared_ptr<const Civec> cc, std::shared_ptr<Dvec> d, const int nact) const;
      void sigma_2ab_2(std::shared_ptr<Dvec> d, std::shared_ptr<Dvec> e, double* mo2e_ptr) const;
      void sigma_2ab_3(std::shared_ptr<Civec> sigma, std::shared_ptr<Dvec> e, const int nact) const;
      
      // gamma = < A' | a^\dagger c | A >
      MatrixPtr form_gamma_alpha(std::shared_ptr<const Dvec> ccvec) const;
      MatrixPtr form_gamma_beta(std::shared_ptr<const Dvec> ccvec) const;

      template<int A, int B, int C, int D> std::pair<MatrixPtr, MatrixPtr> form_JKmatrices() const;
};

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

}

#endif
