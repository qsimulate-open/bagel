//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_coupling.hpp
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

#ifndef __SRC_ASD_GAMMA_COUPLING_HPP
#define __SRC_ASD_GAMMA_COUPLING_HPP

template <class VecType, int N>
void GammaForest<VecType,N>::couple_blocks(const DimerSubspace<VecType>& AB, const DimerSubspace<VecType>& ApBp) {
  static_assert(N == 2, "the following assumes N == 2");
  Coupling term_type = coupling_type(AB, ApBp);

  auto* space1 = &AB;
  auto* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);

  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }

  // Throughout, I'm going to consider space1 to be the ket state, space2 to be the bra state
  if (term_type != Coupling::none) {
    const std::array<int,2> ioffset {{space1->template tag<0>(), space1->template tag<1>()}};
    const std::array<int,2> joffset {{space2->template tag<0>(), space2->template tag<1>()}};

    std::shared_ptr<const VecType> ket_A = space1->template ci<0>();
    std::shared_ptr<const VecType> ket_B = space1->template ci<1>();

    std::shared_ptr<const VecType> bra_A = space2->template ci<0>();
    std::shared_ptr<const VecType> bra_B = space2->template ci<1>();

    switch(term_type) {
      case Coupling::none :
        assert(false); // Control should never be able to reach here
        break;
      case Coupling::diagonal :
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha});
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta});
        break;
      case Coupling::aET : // Alpha ET
        // One-body aET
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::CreateAlpha});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateAlpha});

        //Two-body aET, type 1
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::CreateAlpha});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta});

        //Two-body aET, type 2
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha});
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateAlpha});
        break;
      case Coupling::bET : // Beta ET
        // One-body bET
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::CreateBeta});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateBeta});
        //Two-body bET, type 1
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::CreateBeta});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta});

        //Two-body aET, type 2
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta});
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateBeta});
        break;
      case Coupling::abFlip :
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha});
        break;
      case Coupling::abET :
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::CreateBeta, GammaSQ::CreateAlpha});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});
        break;
      case Coupling::aaET :
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});
        break;
      case Coupling::bbET :
        this->insert<0>(ket_A, ioffset[0], bra_A, joffset[0], {GammaSQ::CreateBeta, GammaSQ::CreateBeta});
        this->insert<1>(ket_B, ioffset[1], bra_B, joffset[1], {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});
        break;
      default :
        assert(false); break; // Control should never reach here
    }
  }
}

#endif
#endif
