//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dimer_jop.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//



#ifndef __BAGEL_DIMER_JOP_H
#define __BAGEL_DIMER_JOP_H

#include <src/ci/fci/mofile.h>

namespace bagel {

class DimerJop : public Jop {
  protected:
    std::pair<std::shared_ptr<MOFile>, std::shared_ptr<MOFile>> jops_;

    // Array is big enough to store all possible coulomb matrices just for simplicity
    std::array<std::shared_ptr<const Matrix>, 16> matrices_;
    std::shared_ptr<const Matrix> cross_mo1e_;

    std::pair<int, int> nact_;

  public:
    DimerJop(const std::shared_ptr<const Reference> ref, const int nstart, const int nfenceA, const int nfenceB,
      std::shared_ptr<const Coeff> coeff); // note that in DimerJop, I'm forcing a HZ Jop
    DimerJop(const int nactA, const int nactB, std::shared_ptr<CSymMatrix> mo1e, std::shared_ptr<Matrix> mo2e);

    std::shared_ptr<const Matrix> cross_mo1e() const { return cross_mo1e_; }

    template<int unit> std::shared_ptr<MOFile> monomer_jop() const { return ( unit == 0 ? jops_.first : jops_.second ); }

    /**
      \brief Generates a subset of the two-electron integrals \f$\langle ab|cd\rangle = (ac|bd)\f$
              where each of \f$a,b,c,d\f$ can belong to either monomer \f$A\f$ or \f$B\f$.
      \details Returns a Matrix(\f$N_A,N_B\f$) where \f$N_A=(\mbox{norb}_A)^d\f$,
               \f$N_B=(\mbox{norb}_B)^{(4-d)}\f$ and \f$d\f$ is the number of \f$A\f$ operators.

      \tparam A specifies whether the first index belongs to monomer \f$A\f$ (0) or \f$B\f$ (1)
      \tparam B specifies whether the second index belongs to monomer \f$A\f$ (0) or \f$B\f$ (1)
      \tparam C specifies whether the third index belongs to monomer \f$A\f$ (0) or \f$B\f$ (1)
      \tparam D specifies whether the fourth index belongs to monomer \f$A\f$ (0) or \f$B\f$ (1)
    */
    template<int A, int B, int C, int D>
    std::shared_ptr<const Matrix> coulomb_matrix();
    /// Same as non-const version but will not create the matrix if it doesn't already exist
    template<int A, int B, int C, int D>
    std::shared_ptr<const Matrix> coulomb_matrix() const;

  private:
    void common_init(const int norbA, const int norbB);

    template <int unit> int nact() const { return ( unit == 0 ? nact_.first : nact_.second ); }
    template <int unit> int active(const int a) const { return a + unit * nact_.first; }

    template <int A, int B, int C, int D>
    std::pair<int, int> index(const int a, const int b, const int c, const int d) const;
};

template<int A, int B, int C, int D>
std::shared_ptr<const Matrix> DimerJop::coulomb_matrix() {
  // First check to see if it's already stored
  const int cindex = A + 2*B + 4*C + 8*D;
  if (matrices_[cindex]) {
    return matrices_[cindex];
  }
  else {
    const int nactA = nact_.first;
    const int nactB = nact_.second;

    int ijA = 1;
    int unitA = 4 - (A + B + C + D);
    for ( int i = 0; i < unitA; ++i ) ijA *= nactA;

    int ijB = 1;
    int unitB = A + B + C + D;
    for ( int i = 0; i < unitB; ++i ) ijB *= nactB;

    auto out = std::make_shared<Matrix>(ijA, ijB);

    for(int d = 0; d < nact<D>(); ++d) {
      for(int c = 0; c < nact<C>(); ++c) {
        for(int b = 0; b < nact<B>(); ++b) {
          for(int a = 0; a < nact<A>(); ++a) {
            int iA, jB;
            std::tie(iA, jB) = index<A,B,C,D>(a,b,c,d);

            out->element(iA,jB) = mo2e_hz(active<A>(a), active<B>(b), active<C>(c), active<D>(d));
          }
        }
      }
    }

    out->localize();
    matrices_[cindex] = out;
    return out;
  }
}

template<int A, int B, int C, int D>
std::shared_ptr<const Matrix> DimerJop::coulomb_matrix() const { return matrices_[A + 2*B + 4*C + 8*D]; }

template <int A, int B, int C, int D> std::pair<int, int> DimerJop::index(int a, int b, int c, int d) const {
  int iA = 0, jB = 0;
  if (D == 0) iA = d; else jB = d;
  if (C == 0) iA = c + iA*nact_.first; else jB = c + jB*nact_.second;
  if (B == 0) iA = b + iA*nact_.first; else jB = b + jB*nact_.second;
  if (A == 0) iA = a + iA*nact_.first; else jB = a + jB*nact_.second;
  return {iA,jB};
}

}

#endif
