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


#include <iostream>
#include <src/rel/smallnaibatch.h>
#include <src/osint/momentbatch.h>
#include <src/osint/overlapbatch.h>
#include <src/rysint/carsphlist.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


SmallNAIBatch::SmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Geometry> geom)
  : geom_(geom), shells_(info), aux_inc_{{shells_[0]->kinetic_balance_uncont(1), shells_[1]->kinetic_balance_uncont(1)}},
    size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()), stack_(resources__->get()) {

  aux_dec_ = array<shared_ptr<const Shell>,2>{{shells_[0]->kinetic_balance_uncont(-1), shells_[1]->kinetic_balance_uncont(-1)}};
  if (aux_dec_[0]) {
    ovl0_dec_ = stack_->get(aux_dec_[0]->nbasis()*aux_dec_[0]->nbasis());
    OverlapBatch ovl0_dec(array<shared_ptr<const Shell>,2>{{aux_dec_[0], aux_dec_[0]}}, stack_);
    ovl0_dec.compute();
    copy(ovl0_dec.data(), ovl0_dec.data()+ovl0_dec.size_block(), ovl0_dec_);
  }
  if (aux_dec_[1]) {
    ovl1_dec_ = stack_->get(aux_dec_[1]->nbasis()*aux_dec_[1]->nbasis());
    OverlapBatch ovl1_dec(array<shared_ptr<const Shell>,2>{{aux_dec_[1], aux_dec_[1]}}, stack_);
    ovl1_dec.compute();
    copy(ovl1_dec.data(), ovl1_dec.data()+ovl1_dec.size_block(), ovl1_dec_);
  }

  // TODO this can be done outside
  ovl0_inc_ = stack_->get(aux_inc_[0]->nbasis()*aux_inc_[0]->nbasis());
  ovl1_inc_ = stack_->get(aux_inc_[1]->nbasis()*aux_inc_[1]->nbasis());
  {
    OverlapBatch ovl0_inc(array<shared_ptr<const Shell>,2>{{aux_inc_[0], aux_inc_[0]}}, stack_);
    ovl0_inc.compute();
    copy(ovl0_inc.data(), ovl0_inc.data()+ovl0_inc.size_block(), ovl0_inc_);
  }
  {
    OverlapBatch ovl1_inc(array<shared_ptr<const Shell>,2>{{aux_inc_[1], aux_inc_[1]}}, stack_);
    ovl1_inc.compute();
    copy(ovl1_inc.data(), ovl1_inc.data()+ovl1_inc.size_block(), ovl1_inc_);
  }

  data_ = stack_->get(size_block_*4);
  fill(data_, data_+size_block_*4, 0.0);
}


SmallNAIBatch::~SmallNAIBatch() {
  stack_->release(size_block_*4, data_);
  stack_->release(aux_inc_[1]->nbasis()*aux_inc_[1]->nbasis(), ovl1_inc_);
  stack_->release(aux_inc_[0]->nbasis()*aux_inc_[0]->nbasis(), ovl0_inc_);
  if (aux_dec_[1]) stack_->release(aux_dec_[1]->nbasis()*aux_dec_[1]->nbasis(), ovl1_dec_);
  if (aux_dec_[0]) stack_->release(aux_dec_[0]->nbasis()*aux_dec_[0]->nbasis(), ovl0_dec_);
  resources__->release(stack_);
}

void SmallNAIBatch::compute() {
  // first compute uncontracted NAI with auxiliary basis (cartesian)
#define LOCAL_DEBUG
#ifdef LOCAL_DEBUG
  shared_ptr<OverlapBatch> nai_inc_inc, nai_dec_dec, nai_dec_inc, nai_inc_dec;
  nai_inc_inc = shared_ptr<OverlapBatch>(new OverlapBatch(aux_inc_, stack_));
  nai_inc_inc->compute();
  if (aux_dec_[0] && aux_dec_[1]) {
    nai_dec_dec = shared_ptr<OverlapBatch>(new OverlapBatch(aux_dec_, stack_));
    nai_dec_dec->compute();
  }
  if (aux_dec_[0]) {
    nai_dec_inc = shared_ptr<OverlapBatch>(new OverlapBatch(array<shared_ptr<const Shell>,2>{{aux_dec_[0],aux_inc_[1]}}, stack_));
    nai_dec_inc->compute();
  }
  if (aux_dec_[1]) {
    nai_inc_dec = shared_ptr<OverlapBatch>(new OverlapBatch(array<shared_ptr<const Shell>,2>{{aux_inc_[0],aux_dec_[1]}}, stack_));
    nai_inc_dec->compute();
  }
#else
  assert(false);
//NAIBatch nai(aux_, geom_, stack_);
#endif

  // then we need to have momentum integrals
  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int a0size = aux_inc_[0]->nbasis();
  const int a1size = aux_inc_[1]->nbasis();
  const int a0size_dec = aux_dec_[0] ? aux_dec_[0]->nbasis() : 0;
  const int a1size_dec = aux_dec_[1] ? aux_dec_[1]->nbasis() : 0;

  double* const ints_dec = stack_->get(3 * s0size * a1size_dec);
  double* const ints_inc = stack_->get(3 * s0size * a1size);
  fill(ints_inc, ints_inc+(3 * s0size * a1size), 0.0);
  fill(ints_dec, ints_dec+(3 * s0size * a1size_dec), 0.0);

  {
    // first half transformation
    // momentum integrals (x,y,z)
    MomentBatch coeff0(array<shared_ptr<const Shell>,2>{{shells_[0]->cartesian_shell(), aux_inc_[0]}}, stack_); 
    coeff0.compute();

    // shell[0] runs faster 
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
        double* const ipiv = stack_->get(a0size);
        double* const tmp = stack_->get(a0size*a0size);
        copy(ovl0_inc_, ovl0_inc_+a0size*a0size, tmp);
        int info = 0;
        dgesv_(a0size, s0size, tmp, a0size, (int*)ipiv, tmparea2, a0size, info); 
        if (info) throw runtime_error("DGESV failed in SmallNAIBatch::compute");
        stack_->release(a0size*a0size, tmp);
        stack_->release(a0size, ipiv);
      }

      if (nai_inc_dec) {
        dgemm_("T", "N", s0size, a1size_dec, a0size, 1.0, tmparea2, a0size, nai_inc_dec->data(), a0size, 1.0, ints_dec+i*s0size*a1size_dec, s0size);
      }
      dgemm_("T", "N", s0size, a1size, a0size, 1.0, tmparea2, a0size, nai_inc_inc->data(), a0size, 1.0, ints_inc+i*s0size*a1size, s0size);
      stack_->release(s0size*a0size, tmparea2);
      stack_->release(s0size*a0size, tmparea);
    }
  }

//NEW STUFF!!!

  if (aux_dec_[0]) {
//if (nai_dec_inc) {
    {
      // first half transformation
      // momentum integrals (x,y,z)
      MomentBatch coeff0(array<shared_ptr<const Shell>,2>{{shells_[0]->cartesian_shell(), aux_dec_[0]}}, stack_); 
      coeff0.compute();
  
      // shell[0] runs faster 
      const double* carea = coeff0.data();
      for (int i = 0; i != 3; ++i, carea += coeff0.size_block()) {
        double* const tmparea = stack_->get(s0size * a0size_dec);
        double* const tmparea2 = stack_->get(s0size * a0size_dec);
        if (shells_[0]->spherical()) {
          const int carsphindex = shells_[0]->angular_number() * ANG_HRR_END + 0; // only transform shell
          const int nloop = shells_[0]->num_contracted() * a0size_dec; 
          carsphlist.carsphfunc_call(carsphindex, nloop, carea, tmparea);
        } else {
          assert(coeff0.size_block() == a0size_dec*s0size);
          copy(carea, carea+coeff0.size_block(), tmparea);
        }
  
        mytranspose_(tmparea, &s0size, &a0size_dec, tmparea2);
        {
          double* const ipiv = stack_->get(a0size_dec);
          double* const tmp = stack_->get(a0size_dec*a0size_dec);
          copy(ovl0_dec_, ovl0_dec_+a0size_dec*a0size_dec, tmp);
          int info = 0;
          dgesv_(a0size_dec, s0size, tmp, a0size_dec, (int*)ipiv, tmparea2, a0size_dec, info); 
          if (info) throw runtime_error("DGESV failed in SmallNAIBatch::compute");
          stack_->release(a0size_dec*a0size_dec, tmp);
          stack_->release(a0size_dec, ipiv);
        }
  
        if (nai_dec_dec) {
          dgemm_("T", "N", s0size, a1size_dec, a0size_dec, 1.0, tmparea2, a0size_dec, nai_dec_dec->data(), a0size_dec, 1.0, ints_dec+i*s0size*a1size_dec, s0size);
        } 

        dgemm_("T", "N", s0size, a1size, a0size_dec, 1.0, tmparea2, a0size_dec, nai_dec_inc->data(), a0size_dec, 1.0, ints_inc+i*s0size*a1size, s0size);
      
        stack_->release(s0size*a0size_dec, tmparea2);
        stack_->release(s0size*a0size_dec, tmparea);
      }
    }
  }

  array<double* const,3> data = {{data_+size_block_, data_+size_block_*2, data_+size_block_*3}};
  array<int,3> f = {{1,2,0}};
  array<int,3> b = {{2,0,1}};

  {
    // second half transformation
    MomentBatch coeff1(array<shared_ptr<const Shell>,2>{{shells_[1]->cartesian_shell(), aux_inc_[1]}}, stack_);
    coeff1.compute();

    double* tmparea = stack_->get(s1size * a1size);
    double* const tmparea2 = stack_->get(s1size * a1size);

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

      mytranspose_(tmparea, &s1size, &a1size, tmparea2);
      {
        double* const ipiv = stack_->get(a1size);
        double* const tmp = stack_->get(a1size*a1size);
        copy(ovl1_inc_, ovl1_inc_+a1size*a1size, tmp);
        int info = 0;
        dgesv_(a1size, s1size, tmp, a1size, (int*)ipiv, tmparea2, a1size, info); 
        if (info) throw runtime_error("DGESV failed in SmallNAIBatch::compute");
        stack_->release(a1size*a1size, tmp);
        stack_->release(a1size, ipiv);
      }
      
#if 0
      if (nai_dec_inc) {
        dgemm_("N", "N", s0size, s1size, a1size, 1.0, ints_inc+i*s0size*a1size, s0size, tmparea2, a1size, 1.0, data_, s0size);
        dgemm_("N", "N", s0size, s1size, a1size, 1.0, ints_inc+b[i]*s0size*a1size, s0size, tmparea2, a1size, 1.0, data[b[i]], s0size);
        dgemm_("N", "N", s0size, s1size, a1size, -1.0, ints_inc+f[i]*s0size*a1size, s0size, tmparea2, a1size, 1.0, data[i], s0size);
      }
#endif
      // slot in appropriate blocks
      // four blocks are needed.
      // 0) x^x + y^y + z^z
      dgemm_("N", "N", s0size, s1size, a1size, 1.0, ints_inc+i*s0size*a1size, s0size, tmparea2, a1size, 1.0, data_, s0size);
      // 1) x^y - y^x
      // 2) y^z - z^y
      // 3) z^x - x^z
      dgemm_("N", "N", s0size, s1size, a1size, 1.0, ints_inc+b[i]*s0size*a1size, s0size, tmparea2, a1size, 1.0, data[b[i]], s0size);
      dgemm_("N", "N", s0size, s1size, a1size, -1.0, ints_inc+f[i]*s0size*a1size, s0size, tmparea2, a1size, 1.0, data[i], s0size);
    }
    stack_->release(s1size*a1size, tmparea2);
    stack_->release(s1size*a1size, tmparea);
  }

  if (aux_dec_[1]) {
    {
      // second half transformation
      MomentBatch coeff1(array<shared_ptr<const Shell>,2>{{shells_[1]->cartesian_shell(), aux_dec_[1]}}, stack_);
      coeff1.compute();
  
      double* tmparea = stack_->get(s1size * a1size_dec);
      double* const tmparea2 = stack_->get(s1size * a1size_dec);
  
      const double* carea = coeff1.data();
      for (int i = 0; i != 3; ++i, carea += coeff1.size_block()) {
        if (shells_[1]->spherical()) {
          const int carsphindex = shells_[1]->angular_number() * ANG_HRR_END + 0; // only transform shell
          const int nloop = shells_[1]->num_contracted() * a1size_dec; 
          carsphlist.carsphfunc_call(carsphindex, nloop, carea, tmparea);
        } else {
          assert(coeff1.size_block() == a1size_dec*s1size);
          copy(carea, carea + coeff1.size_block(), tmparea); 
        }
  
        mytranspose_(tmparea, &s1size, &a1size_dec, tmparea2);
        {
          double* const ipiv = stack_->get(a1size_dec);
          double* const tmp = stack_->get(a1size_dec*a1size_dec);
          copy(ovl1_dec_, ovl1_dec_+a1size_dec*a1size_dec, tmp);
          int info = 0;
          dgesv_(a1size_dec, s1size, tmp, a1size_dec, (int*)ipiv, tmparea2, a1size_dec, info); 
          if (info) throw runtime_error("DGESV failed in SmallNAIBatch::compute");
          stack_->release(a1size_dec*a1size_dec, tmp);
          stack_->release(a1size_dec, ipiv);
        }
        
#if 0 
        if (nai_dec_dec) {
          dgemm_("N", "N", s0size, s1size, a1size_dec, 1.0, ints_dec+i*s0size*a1size_dec, s0size, tmparea2, a1size_dec, 1.0, data_, s0size);
          dgemm_("N", "N", s0size, s1size, a1size_dec, 1.0, ints_dec+b[i]*s0size*a1size_dec, s0size, tmparea2, a1size_dec, 1.0, data[b[i]], s0size);
          dgemm_("N", "N", s0size, s1size, a1size_dec, -1.0, ints_dec+f[i]*s0size*a1size_dec, s0size, tmparea2, a1size_dec, 1.0, data[i], s0size);
        }
#endif
        // slot in appropriate blocks
        // four blocks are needed.
        // 0) x^x + y^y + z^z
        dgemm_("N", "N", s0size, s1size, a1size_dec, 1.0, ints_dec+i*s0size*a1size_dec, s0size, tmparea2, a1size_dec, 1.0, data_, s0size);
        // 1) x^y - y^x
        // 2) y^z - z^y
        // 3) z^x - x^z
        dgemm_("N", "N", s0size, s1size, a1size_dec, 1.0, ints_dec+b[i]*s0size*a1size_dec, s0size, tmparea2, a1size_dec, 1.0, data[b[i]], s0size);
        dgemm_("N", "N", s0size, s1size, a1size_dec, -1.0, ints_dec+f[i]*s0size*a1size_dec, s0size, tmparea2, a1size_dec, 1.0, data[i], s0size);
      }
      stack_->release(s1size*a1size_dec, tmparea2);
      stack_->release(s1size*a1size_dec, tmparea);
    }
  }
  stack_->release(3 * s0size * a1size, ints_inc);
  stack_->release(3 * s0size * a1size_dec, ints_dec);

}
