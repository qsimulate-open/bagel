//
// BAGEL - Parallel electron correlation program.
// Filename: coupling.cc
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

#include <src/asd/coupling.h>

using namespace std;
using namespace bagel;

Coupling bagel::coupling_type(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) {
  array<MonomerKey,4> keys {{ AB.monomerkey<0>(), AB.monomerkey<1>(), ApBp.monomerkey<0>(), ApBp.monomerkey<1>()}};
  return coupling_type(keys);
}


Coupling bagel::coupling_type(const array<MonomerKey,4>& keys) {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  // <AB|
  pair<int,int> neleaAB {A.nelea(), B.nelea()};
  pair<int,int> nelebAB {A.neleb(), B.neleb()};

  // |A'B'>
  pair<int,int> neleaApBp {Ap.nelea(), Bp.nelea()};
  pair<int,int> nelebApBp {Ap.neleb(), Bp.neleb()};

  // AlphaTransfer(AT) and BetaTransfer(BT)
  pair<int,int> AT {neleaAB.first - neleaApBp.first, neleaAB.second - neleaApBp.second};
  pair<int,int> BT {nelebAB.first - nelebApBp.first, nelebAB.second - nelebApBp.second};

  /************************************************************
  *  BT\AT  | ( 0, 0) | (+1,-1) | (-1,+1) | (+2,-2) | (-2,+2) *
  *-----------------------------------------------------------*
  * ( 0, 0) |  diag   |  aET    |  -aET   |  aaET   | -aaET   *
  * (+1,-1) |  bET    |  dABT   |  ABflp  |         |         *
  * (-1,+1) | -bET    | BAflp   | -dABT   |         |         *
  * (+2,-2) |  bbET   |         |         |         |         *
  * (-2,+2) | -bbET   |         |         |         |         *
  ************************************************************/

  const auto icouple = make_tuple(AT.first, AT.second, BT.first, BT.second);

  if      (icouple == make_tuple( 0, 0, 0, 0)) return Coupling::diagonal;
  else if (icouple == make_tuple( 0, 0,+1,-1)) return Coupling::bET;
  else if (icouple == make_tuple( 0, 0,-1,+1)) return Coupling::inv_bET;
  else if (icouple == make_tuple(+1,-1, 0, 0)) return Coupling::aET;
  else if (icouple == make_tuple(+1,-1,+1,-1)) return Coupling::abET;
  else if (icouple == make_tuple(+1,-1,-1,+1)) return Coupling::baFlip;
  else if (icouple == make_tuple(-1,+1, 0, 0)) return Coupling::inv_aET;
  else if (icouple == make_tuple(-1,+1,+1,-1)) return Coupling::abFlip;
  else if (icouple == make_tuple(-1,+1,-1,+1)) return Coupling::inv_abET;
  else if (icouple == make_tuple(+2,-2, 0, 0)) return Coupling::aaET;
  else if (icouple == make_tuple(-2,+2, 0, 0)) return Coupling::inv_aaET;
  else if (icouple == make_tuple( 0, 0,+2,-2)) return Coupling::bbET;
  else if (icouple == make_tuple( 0, 0,-2,+2)) return Coupling::inv_bbET;
  else                                         return Coupling::none;
}
