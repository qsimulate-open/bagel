//
// BAGEL - Parallel electron correlation program.
// Filename: relcoeff.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/wfn/relcoeff.h>
#include <cassert>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace bagel;

RelCoeff::RelCoeff(const ZMatrix& _coeff, const int _nocc, const int _nact, const int _nvirt, const int _nneg, const bool move_neg)
 : ZMatrix(_coeff.ndim(), _coeff.mdim(), _coeff.localized()), nbasis_(ndim()/4), nocc_(_nocc), nact_(_nact), nvirt_(_nvirt), nneg_(_nneg) {
  assert(ndim()%4 == 0);
  assert(2 * (nocc_ + nact_ + nvirt_) + nneg_ == mdim());

  if (!move_neg) {
    // copy input matrix directly
    copy_block(0, 0, ndim(), mdim(), _coeff);
  } else {
    // move positronic orbitals to end of virtual space
    copy_block(0, 0,      ndim(), npos(), _coeff.slice(nneg_, nneg_+npos()));
    copy_block(0, npos(), ndim(), nneg_,  _coeff.slice(0,     nneg_));
  }
}


#if 0
BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff)
BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff_Striped)
BOOST_CLASS_EXPORT_IMPLEMENT(RelCoeff_Block)
#endif
