//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest.h
// Copyright (C) 2013 Shane Parker
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
