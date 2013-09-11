//
// BAGEL - Parallel electron correlation program.
// Filename: meh_gamma_coupling.hpp
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

#ifdef MEH_HEADERS

#ifndef BAGEL_MEH_GAMMA_COUPLING
#define BAGEL_MEH_GAMMA_COUPLING

template <class VecType>
void MultiExcitonHamiltonian<VecType>::gamma_couple_blocks(DSubSpace& AB, DSubSpace& ApBp) {
  Coupling term_type = coupling_type(AB, ApBp);

  DSubSpace* space1 = &AB;
  DSubSpace* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);

  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }

  // Throughout, I'm going to consider space1 to be the ket state, space2 to be the bra state
  if (term_type != Coupling::none) {
    const int ioffset = space1->offset();
    const int joffset = space2->offset();

    std::shared_ptr<const VecType> ket_A = space1->template ci<0>();
    std::shared_ptr<const VecType> ket_B = space1->template ci<1>();

    std::shared_ptr<const VecType> bra_A = space2->template ci<0>();
    std::shared_ptr<const VecType> bra_B = space2->template ci<1>();

    switch(term_type) {
      case Coupling::none :
        assert(false); // Control should never be able to reach here
        break;
      case Coupling::diagonal :
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha);
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta);
        break;
      case Coupling::aET : // Alpha ET
        // One-body aET
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::CreateAlpha);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateAlpha);

        //Two-body aET, type 1
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::CreateAlpha);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta);

        //Two-body aET, type 2
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha);
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateAlpha);
        break;
      case Coupling::bET : // Beta ET
        // One-body bET
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::CreateBeta);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateBeta);
        //Two-body bET, type 1
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::CreateBeta);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta);

        //Two-body aET, type 2
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta);
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateBeta);
        break;
      case Coupling::abFlip :
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha);
        break;
      case Coupling::abET :
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::CreateBeta, GammaSQ::CreateAlpha);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha);
        break;
      case Coupling::aaET :
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha);
        break;
      case Coupling::bbET :
        gammaforest_->template insert<0>(ket_A, ioffset, bra_A, joffset, GammaSQ::CreateBeta, GammaSQ::CreateBeta);
        gammaforest_->template insert<1>(ket_B, ioffset, bra_B, joffset, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta);
        break;
      default :
        assert(false); break; // Control should never reach here
    }
  }
}

#endif

#endif
