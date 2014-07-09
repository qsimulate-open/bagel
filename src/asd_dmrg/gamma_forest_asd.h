//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest_one.h
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __MEH_GAMMA_FOREST_ASD_H
#define __MEH_GAMMA_FOREST_ASD_H

#include <vector>
#include <list>

#include <src/meh/gamma_forest.h>

namespace bagel {

template <class VecType>
class GammaForestASD : public GammaForest<VecType, 1> {
  public:
    GammaForestASD(std::map<SpaceKey, std::shared_ptr<const VecType>> monomer_states) {
      std::vector<std::list<GammaSQ>> possible_couplings = {
        {GammaSQ::AnnihilateAlpha},
        {GammaSQ::AnnihilateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},
        {GammaSQ::AnnihilateAlpha,  GammaSQ::AnnihilateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha},
      };

      for (auto& i : monomer_states) {
        for (auto& coupling : possible_couplings) {
          int new_charge = i.first.q;
          int new_ms = i.first.m_s;
          for (GammaSQ& operation : coupling) {
            const int dq = (operation == GammaSQ::AnnihilateAlpha || operation == GammaSQ::AnnihilateBeta) ? 1 : -1;
            const int sign = (operation == GammaSQ::AnnihilateBeta || operation == GammaSQ::CreateBeta) ? 1 : -1;
            new_charge += dq;
            new_ms += dq*sign;
          }

          for (auto& j : monomer_states) {
            if (new_charge == j.first.q && new_ms == j.first.m_s) {
#ifdef DEBUG
              std::cout << "inserting: <" << j.first.q << ", " << j.first.m_s << "|";
              for (auto opiter = coupling.rbegin(); opiter != coupling.rend(); ++opiter) {
                GammaSQ operation = *opiter;
                if ( operation == GammaSQ::AnnihilateAlpha || operation == GammaSQ::CreateAlpha )
                  std::cout << "(A)";
                else
                  std::cout << "(B)";
                if ( operation == GammaSQ::CreateAlpha || operation == GammaSQ::CreateBeta )
                  std::cout << "^t";
              }
              std::cout << "|" << i.first.q << ", " << i.first.m_s << ">" << std::endl;
#endif
              this->template insert<0>(i.second, i.first.tag(), j.second, j.first.tag(), coupling);
            }
          }
        }
      }

      this->compute();

#ifdef DEBUG
      this->for_each_branch([] (std::shared_ptr<GammaBranch<VecType>> b) { for (auto& i : b->gammas()) i.second->print(); });
#endif
    }

};

}

#endif
