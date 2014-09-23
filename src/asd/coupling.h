//
// BAGEL - Parallel electron correlation program.
// Filename: coupling.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_ASD_COUPLING_H
#define __SRC_ASD_COUPLING_H

#include <src/asd/dimersubspace.h>

namespace bagel {

/// Enumeration of the possible couplings between dimer blocks
enum class Coupling {
  none = 0,       ///< no coupling
  diagonal = 1,   ///< no change in occupation patterns
  aET = 2,        ///< alpha transfer (A <-- B)
  inv_aET = -2,   ///< inverse alpha transfer (A --> B)
  bET = 3,        ///< beta transfer (A <-- B)
  inv_bET = -3,   ///< inverse beta transfer (A --> B)
  abFlip = 4,     ///< alpha-beta flip
  baFlip = -4,    ///< beta-alpha flip
  abET = 5,       ///< alpha+beta transfer (A <-- B)
  inv_abET = -5,  ///< inverse alpha+beta (A --> B)
  aaET = 6,       ///< alpha+alpha transfer (A <-- B)
  inv_aaET = -6,  ///< inverse alpha+alpha transfer (A --> B)
  bbET = 7,       ///< beta+beta transfer (A <-- B)
  inv_bbET = -7   ///< inverse beta+beta transfer (A --> B)
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
