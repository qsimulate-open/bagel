//
// BAGEL - Parallel electron correlation program.
// Filename: asd_init.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifdef ASD_HEADERS

#ifndef BAGEL_ASD_INIT_H
#define BAGEL_ASD_INIT_H

template <class VecType>
ASD<VecType>::ASD(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerCISpace_base<VecType>> cispace)
 : ASD_base(input, dimer), cispace_(cispace) {

  Timer timer;

  cispace_->complete();
  std::cout << "Computing intermediates..." << std::endl;
  cispace_->super_complete();
  std::cout << "  o completing CI spin space: " << timer.tick() << std::endl;

  // Organize subspaces
  dimerstates_ = 0;
  int maxspin = 0;
  // Process DimerCISpace to form and organize needed Civecs and figure out maxspin
  for ( auto& aiter : cispace_->template cispace<0>() ) {
    SpaceKey akey = aiter.first;
    // values of S_b that could give the proper spin:
    //   S_a-nspin_, S_a-nspin+2,...,S_a+nspin
    for (int Sb = std::abs(akey.S-std::abs(nspin_)); Sb <= akey.S+std::abs(nspin_); Sb+=2) {
      SpaceKey bkey( Sb, nspin_ - akey.m_s, charge_ - akey.q );
      std::shared_ptr<const VecType> bspace = cispace_->template ccvec<1>(bkey);
      if ( bspace ) {
        subspaces_.emplace_back(dimerstates_, akey, bkey, make_pair(aiter.second, bspace));
        maxspin = std::max(akey.S+Sb, maxspin);
      }
    }
  }
  max_spin_ = maxspin + 1;

  std::cout << "!@# number of elements in 'subspaces_' = " << subspaces_.size() << std::endl;//TODO:delete
}

template <class VecType>
std::shared_ptr<Matrix> ASD<VecType>::compute_1e_prop(std::shared_ptr<const Matrix> hAA, std::shared_ptr<const Matrix> hBB, std::shared_ptr<const Matrix> hAB, const double core) const {

  auto out = std::make_shared<Matrix>(dimerstates_, dimerstates_);

  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

// TODO remove this comment once the gammaforst issue has been fixed (bra and ket have been exchanged)
      std::array<MonomerKey,4> keys {{ jAB->template monomerkey<0>(), jAB->template monomerkey<1>(),
                                       iAB->template monomerkey<0>(), iAB->template monomerkey<1>() }};
      std::shared_ptr<Matrix> out_block = compute_offdiagonal_1e_H(keys, hAB);

      out->add_block(1.0, joff, ioff, out_block->ndim(), out_block->mdim(), out_block);
      out->add_block(1.0, ioff, joff, out_block->mdim(), out_block->ndim(), out_block->transpose());
    }
    std::shared_ptr<const Matrix> tmp = compute_diagonal_1e(*iAB, hAA->data(), hBB->data(), core);
    out->add_block(1.0, ioff, ioff, tmp->ndim(), tmp->mdim(), tmp);
  }

  return out;
}



#endif

#endif
