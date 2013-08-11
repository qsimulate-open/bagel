//
// BAGEL - Parallel electron correlation program.
// Filename: meh_spin_coupling.cc
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <src/meh/meh.h>

using namespace bagel;
using namespace std;

void MultiExcitonHamiltonian::spin_couple_blocks(DimerSubspace& AB, DimerSubspace& ApBp) {
  const Coupling term_type = coupling_type(AB, ApBp);

  auto spin_block = make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  if ( (term_type == Coupling::abFlip) || (term_type == Coupling::baFlip) ) {
    shared_ptr<Dvec> SA, SB;

    shared_ptr<const Dvec> Ap = ApBp.ci<0>();
    shared_ptr<const Dvec> Bp = ApBp.ci<1>();
    switch (term_type) {
      case Coupling::abFlip :
        SA = AB.ci<0>()->spin_lower(Ap->det());
        SB = AB.ci<1>()->spin_raise(Bp->det());
        break;
      case Coupling::baFlip :
        SA = AB.ci<0>()->spin_raise(Ap->det());
        SB = AB.ci<1>()->spin_lower(Bp->det());
        break;
      default: assert(false); break; // Control should never be able to reach here...
    }

    const int nA = AB.nstates<0>();
    const int nB = AB.nstates<1>();
    const int nAp = ApBp.nstates<0>();
    const int nBp = ApBp.nstates<1>();

    vector<double> AdotAp, BdotBp;

    for (int iAp = 0; iAp < nAp; ++iAp) {
      for (int iA = 0; iA < nA; ++iA) {
        AdotAp.push_back(SA->data(iA)->dot_product(*Ap->data(iAp)));
      }
    }

    for (int iBp = 0; iBp < nBp; ++iBp) {
      for (int iB = 0; iB < nB; ++iB) {
        BdotBp.push_back(SB->data(iB)->dot_product(*Bp->data(iBp)));
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
        if ( fabs(ele) > 1.0e-4) spin_->insert(ispin + ioff, jspin + joff, ele);
      }
    }
  }
}


void MultiExcitonHamiltonian::compute_diagonal_spin_block(DimerSubspace& subspace) {
  shared_ptr<const Dvec> Ap = subspace.ci<0>();
  shared_ptr<const Dvec> Bp = subspace.ci<1>();
  shared_ptr<const Dvec> SA = Ap->spin();
  shared_ptr<const Dvec> SB = Bp->spin();

  const int nA = subspace.nstates<0>();
  const int nB = subspace.nstates<1>();

  const double sz_AB = 0.5 * static_cast<double>(Ap->det()->nspin() * Bp->det()->nspin());

  auto spin_block = make_shared<Matrix>(nA*nB, nA*nB);

  vector<double> AdotAp;
  vector<double> BdotBp;

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
      if (fabs(ele) > 1.0e-4) spin_->insert(ispin + offset, jspin + offset, ele);
    }
    spin_->diagonal(ispin + offset) = spin_block->element(ispin, ispin);
  }
}
