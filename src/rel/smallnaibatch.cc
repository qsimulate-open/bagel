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
  : geom_(geom), shells_(info), aux_{{shells_[0]->kinetic_balance_uncont(), shells_[1]->kinetic_balance_uncont()}},
    size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()), stack_(resources__->get()) {

  // TODO do you need 4? or 3?
  data_ = stack_->get(size_block_*4);
}


SmallNAIBatch::~SmallNAIBatch() {
  stack_->release(size_block_*4, data_);
  resources__->release(stack_);
}

void SmallNAIBatch::compute() {
  // first compute uncontracted NAI with auxiliary basis (cartesian)
  NAIBatch nai(aux_, geom_, stack_);
  nai.compute();

  // then we need to have momentum integrals
  {
    // first half transformation
    MomentBatch coeff0(array<shared_ptr<const Shell>,2>{{aux_[0], shells_[0]->cartesian_shell()}}, stack_); 
  }
  {
    // second half transformation
    MomentBatch coeff1(array<shared_ptr<const Shell>,2>{{aux_[0], shells_[0]->cartesian_shell()}}, stack_);
  }
  {
    // optionally transformation to spherical
  }

}
