//
// BAGEL - Parallel electron correlation program.
// Filename: complexparalleldf.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/df/complexparalleldf.h>

using namespace std;
using namespace bagel;


ComplexParallelDF::ComplexParallelDF(const size_t naux, const size_t nb1, const size_t nb2, const int nbl, shared_ptr<const ComplexParallelDF> df, shared_ptr<Matrix> dat)
 : nblock_(nbl), naux_(naux), nindex1_(nb1), nindex2_(nb2), df_(df), data2_(dat),
   comp_(df ? df->swap_ : true), swap_(df ? df->swap_ : false), serial_(df ? df->serial_ : false) { }


// TODO Probably will have to generalize this for complex<double> a
void ComplexParallelDF::ax_plus_y(const double a, const shared_ptr<const ComplexParallelDF> o) {
  throw runtime_error("ComplexParallelDF::ax_plus_y should be right for real a, but has not been verified");
  assert(nblock_ == o->nblock_);

  {
    auto j = o->dfdata_[0]->block_.begin();
    for (auto& i : dfdata_[0]->block_)
      i->ax_plus_y(a, *j++);
  }
  {
    auto j = o->dfdata_[1]->block_.begin();
    for (auto& i : dfdata_[1]->block_)
      i->ax_plus_y(a, *j++);
  }
}


// TODO Probably will have to generalize this for complex<double> a
void ComplexParallelDF::scale(const double a) {
  throw runtime_error("ComplexParallelDF::scale should be right for real a, but has not been verified");
  for (auto& i : dfdata_[0]->block_)
    i->scale(a);
  for (auto& i : dfdata_[1]->block_)
    i->scale(a);
}
