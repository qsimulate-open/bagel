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
#include <src/rysint/carsphlist.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


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
  const size_t s0size = shells_[0]->nbasis();
  const size_t s1size = shells_[1]->nbasis();
  const size_t a0size = aux_[0]->nbasis();
  const size_t a1size = aux_[1]->nbasis();
  double* ints = stack_->get(3 * s0size * a1size);
  {
    // first half transformation
    // momentum integrals (x,y,z)
    MomentBatch coeff0(array<shared_ptr<const Shell>,2>{{shells_[0]->cartesian_shell(), aux_[0]}}, stack_); 

    // shell[0] runs faster 
    double* tmparea = stack_->get(s0size * a0size);
    for (int i = 0; i != 3; ++i) {
      if (shells_[0]->spherical()) {
        const int carsphindex = shells_[0]->angular_number() * ANG_HRR_END + 0; // only transform shell
        const int nloop = shells_[0]->num_contracted() * a0size; 
        carsphlist.carsphfunc_call(carsphindex, nloop, coeff0.data(), tmparea);
      } else {
        assert(coeff0.size_block() == a0size*s0size);
        copy(coeff0.data(), coeff0.data()+coeff0.size_block(), tmparea); 
      }
      dgemm_("N", "N", s0size, a1size, a0size, 1.0, tmparea, s0size, nai.data(), a0size, 0.0, ints+i*s0size*a1size, s0size);
    }
    stack_->release(shells_[0]->nbasis()*aux_[0]->nbasis(), tmparea);
  }
  {
    // second half transformation
    MomentBatch coeff1(array<shared_ptr<const Shell>,2>{{shells_[1]->cartesian_shell(), aux_[1]}}, stack_);
  }

  stack_->release(3*shells_[0]->nbasis()*aux_[1]->nbasis(), ints);

}
