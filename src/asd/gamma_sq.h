//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_forest.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __BAGEL_ASD_GAMMA_SQ_H
#define __BAGEL_ASD_GAMMA_SQ_H

#include <bitset>

namespace bagel {

enum class GammaSQ {
  CreateAlpha = 0,
  AnnihilateAlpha = 1,
  CreateBeta = 2,
  AnnihilateBeta = 3
};

namespace {
  bool is_alpha(const GammaSQ& a) { return !std::bitset<2>(static_cast<int>(a))[1]; }
  bool is_beta(const GammaSQ& a) { return std::bitset<2>(static_cast<int>(a))[1]; }
  bool is_creation(const GammaSQ& a) { return !std::bitset<2>(static_cast<int>(a))[0]; }
  bool is_annihilation(const GammaSQ& a) { return std::bitset<2>(static_cast<int>(a))[0]; }
  GammaSQ conjugate(const GammaSQ& a) { return GammaSQ(static_cast<int>(a) ^ 1); }
}

inline std::ostream& operator<<(std::ostream& out, const GammaSQ value){
  static std::map<GammaSQ, std::string> strings;
  if (strings.size() == 0) {
#define INSERT_ELEMENT(p) strings[p] = #p
    INSERT_ELEMENT(GammaSQ::CreateAlpha); INSERT_ELEMENT(GammaSQ::AnnihilateAlpha); INSERT_ELEMENT(GammaSQ::CreateBeta); INSERT_ELEMENT(GammaSQ::AnnihilateBeta);
#undef INSERT_ELEMENT
  }
  return out << std::setw(25) << std::left << strings[value] << std::right;
}

}

#endif
