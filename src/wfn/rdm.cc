//
// BAGEL - Parallel electron correlation program.
// Filename: rdm.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <src/wfn/rdm.h>

using namespace bagel;
using namespace std;

RDM_base::RDM_base(const int n, const int rank) : norb_(n), rank_(rank) {
  assert(rank > 0);
  dim_ = 1lu;
  for (int i = 0; i != rank; ++i) dim_ *= n;
  data_ = unique_ptr<double[]>(new double[dim_*dim_]);
}


RDM_base::RDM_base(const RDM_base& o) : norb_(o.norb_), dim_(o.dim_), rank_(o.rank_) {
  data_ = unique_ptr<double[]>(new double[dim_*dim_]);
  copy_n(o.data(), dim_*dim_, data());
}
