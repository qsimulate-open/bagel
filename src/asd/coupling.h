//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: coupling.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_ASD_COUPLING_H
#define __SRC_ASD_COUPLING_H

#include <src/asd/dimersubspace.h>

namespace bagel {

/// Enumeration of the possible couplings between dimer blocks
enum class Coupling {
  none = 0,       ///< no coupling
  diagonal = 1,   ///< no change in occupation patterns
  aET = 2,        ///< alpha transfer (A --> B)
  inv_aET = -2,   ///< inverse alpha transfer (B --> A)
  bET = 3,        ///< beta transfer (A --> B)
  inv_bET = -3,   ///< inverse beta transfer (B --> A)
  abFlip = 4,     ///< alpha-beta flip
  baFlip = -4,    ///< beta-alpha flip
  abET = 5,       ///< alpha+beta transfer (A --> B)
  inv_abET = -5,  ///< inverse alpha+beta (B --> A)
  aaET = 6,       ///< alpha+alpha transfer (A --> B)
  inv_aaET = -6,  ///< inverse alpha+alpha transfer (B --> A)
  bbET = 7,       ///< beta+beta transfer (A --> B)
  inv_bbET = -7   ///< inverse beta+beta transfer (B --> A)
};

inline std::ostream& operator<<(std::ostream& out, const Coupling value){
  static std::map<Coupling, std::string> strings;
  if (strings.size() == 0) {
#define INSERT_ELEMENT(p) strings[p] = #p
    INSERT_ELEMENT(Coupling::none);
    INSERT_ELEMENT(Coupling::diagonal);
    INSERT_ELEMENT(Coupling::aET);
    INSERT_ELEMENT(Coupling::inv_aET);
    INSERT_ELEMENT(Coupling::bET);
    INSERT_ELEMENT(Coupling::inv_bET);
    INSERT_ELEMENT(Coupling::abFlip);
    INSERT_ELEMENT(Coupling::baFlip);
    INSERT_ELEMENT(Coupling::abET);
    INSERT_ELEMENT(Coupling::inv_abET);
    INSERT_ELEMENT(Coupling::aaET);
    INSERT_ELEMENT(Coupling::inv_aaET);
    INSERT_ELEMENT(Coupling::bbET);
    INSERT_ELEMENT(Coupling::inv_bbET);
#undef INSERT_ELEMENT
  }
  return out << std::setw(25) << std::left << strings[value] << std::right;
}

Coupling coupling_type(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp);
Coupling coupling_type(const std::array<MonomerKey,4>& keys);

}

#endif
