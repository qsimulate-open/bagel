//
// BAGEL - Parallel electron correlation program.
// Filename: smallnaibatch.cc
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


#include <src/rel/smallnaibatch.h>
#include <src/osint/momentbatch.h>

using namespace std;
using namespace bagel;


SmallNAIBatch::SmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Geometry> geom)
  : shells_(info), aux_{{shells_[0]->kinetic_balance_uncont(), shells_[1]->kinetic_balance_uncont()}}, nai_(new NAIBatch(aux_, geom)),
    size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()) {

  // TODO do you need 4? or 3?
  data_ = unique_ptr<double[]>(new double[size_block_*4]);
}


void SmallNAIBatch::compute() {
  // first compute uncontracted NAI with auxiliary basis (cartesian)
  nai_->compute();

  // then we need to have momentum integrals
  MomentBatch coeff0(array<shared_ptr<const Shell>,2>{{aux_[0], shells_[0]->cartesian_shell()}}); 
  MomentBatch coeff1(array<shared_ptr<const Shell>,2>{{aux_[0], shells_[0]->cartesian_shell()}});

assert(false);
}
