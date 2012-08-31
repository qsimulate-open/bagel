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
#include <src/osint/overlapbatch.h>
#include <src/rysint/carsphlist.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


SmallNAIBatch::SmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Geometry> geom)
  : geom_(geom), shells_(info), aux_{{shells_[0]->kinetic_balance_uncont(), shells_[1]->kinetic_balance_uncont()}},
    size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()), stack_(resources__->get()) {

  // TODO this can be done outside
  ovl0_ = stack_->get(aux_[0]->nbasis()*aux_[0]->nbasis());
  ovl1_ = stack_->get(aux_[1]->nbasis()*aux_[1]->nbasis());
  {
    OverlapBatch ovl0(array<shared_ptr<const Shell>,2>{{aux_[0], aux_[0]}}, stack_);
    copy(ovl0.data(), ovl0.data()+ovl0.size_block(), ovl0_);
  }
  {
    OverlapBatch ovl1(array<shared_ptr<const Shell>,2>{{aux_[1], aux_[1]}}, stack_);
    copy(ovl1.data(), ovl1.data()+ovl1.size_block(), ovl1_);
  }
  data_ = stack_->get(size_block_*4);
}


SmallNAIBatch::~SmallNAIBatch() {
  stack_->release(size_block_*4, data_);
  stack_->release(aux_[1]->nbasis()*aux_[1]->nbasis(), ovl1_);
  stack_->release(aux_[0]->nbasis()*aux_[0]->nbasis(), ovl0_);
  resources__->release(stack_);
}

void SmallNAIBatch::compute() {
  // first compute uncontracted NAI with auxiliary basis (cartesian)
//OverlapBatch nai(aux_, stack_);
  NAIBatch nai(aux_, geom_, stack_);
  nai.compute();

  // then we need to have momentum integrals
  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int a0size = aux_[0]->nbasis();
  const int a1size = aux_[1]->nbasis();
  double* ints = stack_->get(3 * s0size * a1size);
cout << "ca" << s0size << " " << s1size << " " << a0size << " " << a1size << endl;
  {
    // first half transformation
    // momentum integrals (x,y,z)
    MomentBatch coeff0(array<shared_ptr<const Shell>,2>{{shells_[0]->cartesian_shell(), aux_[0]}}, stack_); 

    // shell[0] runs faster 
    {
    const double* carea = coeff0.data();
    for (int i = 0; i != 3; ++i, carea += coeff0.size_block()) {
      double* const tmparea = stack_->get(s0size * a0size);
      double* const tmparea2 = stack_->get(s0size * a0size);
      if (shells_[0]->spherical()) {
        const int carsphindex = shells_[0]->angular_number() * ANG_HRR_END + 0; // only transform shell
        const int nloop = shells_[0]->num_contracted() * a0size; 
        carsphlist.carsphfunc_call(carsphindex, nloop, carea, tmparea);
      } else {
        assert(coeff0.size_block() == a0size*s0size);
        copy(carea, carea+coeff0.size_block(), tmparea);
      }

      mytranspose_(tmparea, &s0size, &a0size, tmparea2);
      {
//      double* const ipiv = stack_->get(a0size);
        double* const tmp = stack_->get(a0size*a0size);
        copy(ovl0_, ovl0_+a0size*a0size, tmp);
//      int* const ipivi = reinterpret_cast<int* const>(ipiv);
        int* ipivi = new int[a0size];
        int info = 0;
//      dgesv_(a0size, s0size, tmp, a0size, ipivi, tmparea2, a0size, info); 
        if (info) throw runtime_error("DGESV failed in SmallNAIBatch::compute");
        stack_->release(a0size*a0size, tmp);
//      stack_->release(a0size, ipiv);
      }
cout << "done " << i << endl;

      dgemm_("T", "N", s0size, a1size, a0size, 1.0, tmparea2, a0size, nai.data(), a0size, 0.0, ints+i*s0size*a1size, s0size);
      stack_->release(s0size*a0size, tmparea2);
      stack_->release(s0size*a0size, tmparea);
    }
    }
  }
cout << "cb" << endl;
  {
    // second half transformation
    MomentBatch coeff1(array<shared_ptr<const Shell>,2>{{shells_[1]->cartesian_shell(), aux_[1]}}, stack_);
    fill(data_, data_+size_block_*4, 0.0);

    array<double* const,3> data = {{data_+size_block_, data_+size_block_*2, data_+size_block_*3}};
    array<int,3> f = {{1,2,0}};
    array<int,3> b = {{2,0,1}};

    double* tmparea = stack_->get(s1size * a1size);
    const double* carea = coeff1.data();
    for (int i = 0; i != 3; ++i, carea += coeff1.size_block()) {
      if (shells_[1]->spherical()) {
        const int carsphindex = shells_[1]->angular_number() * ANG_HRR_END + 0; // only transform shell
        const int nloop = shells_[1]->num_contracted() * a1size; 
        carsphlist.carsphfunc_call(carsphindex, nloop, carea, tmparea);
      } else {
        assert(coeff1.size_block() == a1size*s1size);
        copy(carea, carea + coeff1.size_block(), tmparea); 
      }
      // slot in appropriate blocks
      // four blocks are needed.
      // 0) x^x + y^y + z^z
      dgemm_("N", "T", s0size, s1size, a1size, 1.0, ints+i*s0size*a1size, s0size, tmparea, s1size, 1.0, data_, s0size);
      // 1) x^y - y^x
      // 2) y^z - z^y
      // 3) z^x - x^z
      dgemm_("N", "T", s0size, s1size, a1size, 1.0, ints+b[i]*s0size*a1size, s0size, tmparea, s1size, 1.0, data[b[i]], s0size);
      dgemm_("N", "T", s0size, s1size, a1size, -1.0, ints+f[i]*s0size*a1size, s0size, tmparea, s1size, 1.0, data[i], s0size);
    }
    stack_->release(s1size*a1size, tmparea);
  }

  stack_->release(3*shells_[0]->nbasis()*aux_[1]->nbasis(), ints);

}
