//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_spin_coupling.hpp
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifdef ASD_HEADERS

#ifndef BAGEL_ASD_SPIN_COUPLING_H
#define BAGEL_ASD_SPIN_COUPLING_H

template <class VecType>
void ASDSpinMap<VecType>::couple_blocks(const DimerSubspace<VecType>& AB, const DimerSubspace<VecType>& ApBp) {
  const Coupling term_type = coupling_type(AB, ApBp);

  auto spin_block = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  if ( (term_type == Coupling::abFlip) || (term_type == Coupling::baFlip) ) {
    std::shared_ptr<VecType> SAp, SBp;

    std::shared_ptr<const VecType> A = AB.template ci<0>();
    std::shared_ptr<const VecType> B = AB.template ci<1>();
    switch (term_type) {
      case Coupling::abFlip :
        SAp = ApBp.template ci<0>()->spin_lower(A->det());
        SBp = ApBp.template ci<1>()->spin_raise(B->det());
        break;
      case Coupling::baFlip :
        SAp = ApBp.template ci<0>()->spin_raise(A->det());
        SBp = ApBp.template ci<1>()->spin_lower(B->det());
        break;
      default: assert(false); break; // Control should never be able to reach here...
    }

    const int nA = AB.template nstates<0>();
    const int nB = AB.template nstates<1>();
    const int nAp = ApBp.template nstates<0>();
    const int nBp = ApBp.template nstates<1>();

    std::vector<double> AdotAp, BdotBp;

    for (int iAp = 0; iAp < nAp; ++iAp) {
      for (int iA = 0; iA < nA; ++iA) {
        AdotAp.push_back(SAp->data(iAp)->dot_product(*A->data(iA)));
      }
    }

    for (int iBp = 0; iBp < nBp; ++iBp) {
      for (int iB = 0; iB < nB; ++iB) {
        BdotBp.push_back(SBp->data(iBp)->dot_product(*B->data(iB)));
      }
    }


    for (int iBp = 0; iBp < nBp; ++iBp) {
      for (int iAp = 0; iAp < nAp; ++iAp) {
        const int iABp = ApBp.dimerindex(iAp,iBp);
        int iAB = 0;
        for (int iB = 0; iB < nB; ++iB) {
          for (int iA = 0; iA < nA; ++iA, ++iAB) {
            spin_block->element(iAB,iABp) += AdotAp[iA+ nA*iAp] * BdotBp[iB + nB*iBp];
          }
        }
      }
    }

    const int n = spin_block->ndim();
    const int m = spin_block->mdim();

    const int ioff = AB.offset();
    const int joff = ApBp.offset();

    for (int ispin = 0; ispin < n; ++ispin) {
      for (int jspin = 0; jspin < m; ++jspin) {
        const double ele = spin_block->element(ispin, jspin);
        if ( std::fabs(ele) > 1.0e-4)  {
          this->emplace(std::make_pair(ispin + ioff, jspin + joff), ele);
          this->emplace(std::make_pair(jspin + joff, ispin + ioff), ele);
        }
      }
    }
  }
}


// diagonal parts:
//   S^2 = [ (S^A)^2 + (S^B)^2 ] + (S^2 - (S^A)^2 - (S^B)^2)
//       = [ (S^A)^2 + (S^B)^2 ] + 2 S_z^A S_z^B
template <class VecType>
void ASDSpinMap<VecType>::diagonal_block(const DimerSubspace<VecType>& subspace) {
  std::shared_ptr<const VecType> Ap = subspace.template ci<0>();
  std::shared_ptr<const VecType> Bp = subspace.template ci<1>();
  std::shared_ptr<const VecType> SA = Ap->spin();
  std::shared_ptr<const VecType> SB = Bp->spin();

  const int nA = subspace.template nstates<0>();
  const int nB = subspace.template nstates<1>();

  const double sz_AB = 0.5 * static_cast<double>(Ap->det()->nspin() * Bp->det()->nspin());

  auto spin_block = std::make_shared<Matrix>(nA*nB, nA*nB);

  std::vector<double> AdotAp;
  std::vector<double> BdotBp;

  for (int iAp = 0; iAp < nA; ++iAp) {
    for (int iA = 0; iA < nA; ++iA)
      AdotAp.push_back(SA->data(iA)->dot_product(*Ap->data(iAp)));
  }

  for (int iBp = 0; iBp < nB; ++iBp) {
    for (int iB = 0; iB < nB; ++iB)
      BdotBp.push_back(SB->data(iB)->dot_product(*Bp->data(iBp)));
  }

  for (int iBp = 0; iBp < nB; ++iBp) {
    for (int iAp = 0; iAp < nA; ++iAp) {
      const int iABp = subspace.dimerindex(iAp, iBp);
      for (int iA = 0; iA < nA; ++iA) // iB = iBp
        spin_block->element(subspace.dimerindex(iA,iBp),iABp) += AdotAp[iA+ nA*iAp];

      for (int iB = 0; iB < nB; ++iB) // iA = iAp
        spin_block->element(subspace.dimerindex(iAp,iB),iABp) += BdotBp[iB+ nB*iBp];

      spin_block->element(iABp, iABp) += sz_AB;
    }
  }

  const int n = spin_block->ndim();
  const int offset = subspace.offset();

  for (int ispin = 0; ispin < n; ++ispin) {
    for (int jspin = 0; jspin < ispin; ++jspin) {
      const double ele = spin_block->element(ispin,jspin);
      if (std::fabs(ele) > 1.0e-4) {
        this->emplace(std::make_pair(ispin + offset, jspin + offset), ele);
        this->emplace(std::make_pair(jspin + offset, ispin + offset), ele);
      }
    }
    const double ele = spin_block->element(ispin, ispin);
    if (std::fabs(ele) > 1.0e-4) this->emplace(std::make_pair(ispin + offset, ispin + offset), ele);
  }
}

#endif

#endif
