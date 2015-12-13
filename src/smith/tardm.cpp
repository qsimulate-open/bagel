//
// BAGEL - Parallel electron correlation program.
// Filename: tardm.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifdef SPINFREEMETHOD_DETAIL

template<>
void SpinFreeMethod<complex<double>>::feed_rdm_ta() {

  // first set the CI vectors to TA tensors.
  shared_ptr<const RelCIWfn> ciwfn = info_->ciwfn();
  shared_ptr<const RelZDvec> reldvec = ciwfn->civectors();
  map<pair<int, int>, shared_ptr<const ZDvec>> dvecs = reldvec->dvecs();

  const int nstates = info_->ciwfn()->nstates();
  for (int n = 0; n != nstates; ++n) {

  }

}

template<>
void SpinFreeMethod<double>::feed_rdm_ta() {
  assert(false);
}

#endif
