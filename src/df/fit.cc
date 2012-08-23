//
// Newint - Parallel electron correlation program.
// Filename: fit.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/df/fit.h>
#include <src/rysint/libint.h>

using namespace std;

pair<const double*, shared_ptr<RysInt> > ERIFit::compute_batch(array<shared_ptr<const Shell>,4>& input) { 
#ifdef LIBINT_INTERFACE
  shared_ptr<Libint> eribatch(new Libint(input));
#else
  shared_ptr<ERIBatch> eribatch(new ERIBatch(input, 2.0));
#endif
  eribatch->compute();
  return make_pair(eribatch->data(), eribatch);
}


pair<const double*, shared_ptr<RysInt> > YukawaFit::compute_batch(array<shared_ptr<const Shell>,4>& input) {
  shared_ptr<SlaterBatch> slaterbatch(new SlaterBatch(input, 0.0, gamma_, true)); // TODO true meas it computes Yukawa and Slater together, but Slater is discarded
  slaterbatch->compute();
  return make_pair(slaterbatch->data2(), slaterbatch);
}


pair<const double*, shared_ptr<RysInt> > SlaterFit::compute_batch(array<shared_ptr<const Shell>,4>& input) {
  shared_ptr<SlaterBatch> slaterbatch(new SlaterBatch(input, 0.0, gamma_, false));
  slaterbatch->compute();
  return make_pair(slaterbatch->data(), slaterbatch);
}

