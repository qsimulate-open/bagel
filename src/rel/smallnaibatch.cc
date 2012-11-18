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
#include <iomanip>
#include <src/rel/smallnaibatch.h>
#include <src/osint/momentbatch.h>
#include <src/osint/overlapbatch.h>
#include <src/rysint/carsphlist.h>
#include <src/util/matrix.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


SmallNAIBatch::SmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Geometry> geom)
  : geom_(geom), shells_(info), aux_inc_{{shells_[0]->kinetic_balance_uncont(1), shells_[1]->kinetic_balance_uncont(1)}},
    size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()), stack_(resources__->get()) {

  aux_dec_ = array<shared_ptr<const Shell>,2>{{shells_[0]->kinetic_balance_uncont(-1), shells_[1]->kinetic_balance_uncont(-1)}};
  const std::array<const int,2> a{{aux_inc_[0]->nbasis() + (aux_dec_[0] ? aux_dec_[0]->nbasis() : 0), aux_inc_[1]->nbasis() + (aux_dec_[1] ? aux_dec_[1]->nbasis() : 0)}};
  ovl_[0] = stack_->get(a[0]*a[0]);
  ovl_[1] = stack_->get(a[1]*a[1]);

  // TODO this can be done outside
  ovl_compute(a);

//make matrix
  data_ = stack_->get(size_block_*4);
  fill(data_, data_+size_block_*4, 0.0);
}


SmallNAIBatch::~SmallNAIBatch() {
  const std::array<const int,2> a{{aux_inc_[0]->nbasis() + (aux_dec_[0] ? aux_dec_[0]->nbasis() : 0), aux_inc_[1]->nbasis() + (aux_dec_[1] ? aux_dec_[1]->nbasis() : 0)}};
  stack_->release(size_block_*4, data_);
  stack_->release(a[1]*a[1], ovl_[1]);
  stack_->release(a[0]*a[0], ovl_[0]);
  resources__->release(stack_);
}

void SmallNAIBatch::compute() {
  // then we need to have momentum integrals
  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int a0size_inc = aux_inc_[0]->nbasis();
  const int a1size_inc = aux_inc_[1]->nbasis();
  const int a0size_dec = aux_dec_[0] ? aux_dec_[0]->nbasis() : 0;
  const int a1size_dec = aux_dec_[1] ? aux_dec_[1]->nbasis() : 0;
  const int a0 = a0size_inc + a0size_dec;
  const int a1 = a1size_inc + a1size_dec;


  //Matrix nai
  const std::shared_ptr<Matrix> nai(new Matrix(a0,a1));
  nai->zero();

  // current nai array
  //double* const nai = stack_->get(a0*a1);

  // first compute uncontracted NAI with auxiliary basis (cartesian)
#define LOCAL_DEBUG
#ifdef LOCAL_DEBUG
  //make nai_compute function
  nai_compute(nai);

#else
  assert(false);
//NAIBatch nai(aux_, geom_, stack_);
#endif

#if 0
  double* const ints = stack_->get(3 * s0size * a1);
  fill(ints, ints+(3 * s0size * a1), 0.0);

  {
    // first half transformation
    // momentum integrals (x,y,z)
    shared_ptr<MomentBatch> coeff0(new MomentBatch(array<shared_ptr<const Shell>,2>{{shells_[0]->cartesian_shell(), aux_inc_[0]}}, stack_));
    coeff0->compute();

    shared_ptr<MomentBatch> coeff1;
    if (aux_dec_[0]) {
      coeff1 = shared_ptr<MomentBatch>(new MomentBatch(array<shared_ptr<const Shell>,2>{{shells_[0]->cartesian_shell(), aux_dec_[0]}}, stack_));
      coeff1->compute();
    } else {
      // TODO just to run
      coeff1 = coeff0;
    }

    // shell[0] runs faster 
    const double* carea0 = coeff0->data();
    const double* carea1 = coeff1->data();
    for (int i = 0; i != 3; ++i, carea0 += coeff0->size_block(), carea1 += coeff1->size_block()) {
      double* const tmparea = stack_->get(s0size * a0);
      double* const tmparea2 = stack_->get(s0size * a0);
      if (shells_[0]->spherical()) {
        const int carsphindex = shells_[0]->angular_number() * ANG_HRR_END + 0; // only transform shell
        const int nloop = shells_[0]->num_contracted() * a0size_inc; 
        carsphlist.carsphfunc_call(carsphindex, nloop, carea0, tmparea);
      } else {
        assert(coeff0->size_block() == a0size_inc*s0size);
        copy(carea0, carea0+coeff0->size_block(), tmparea);
      }
      if (aux_dec_[0]) {
        if (shells_[0]->spherical()) {
          const int carsphindex = shells_[0]->angular_number() * ANG_HRR_END + 0; // only transform shell
          const int nloop = shells_[0]->num_contracted() * a0size_dec; 
          carsphlist.carsphfunc_call(carsphindex, nloop, carea1, tmparea+a0size_inc*s0size);
        } else {
          assert(coeff1->size_block() == a0size_dec*s0size);
          copy(carea1, carea1+coeff1->size_block(), tmparea+a0size_inc*s0size);
        }
      }

 //util/matrix.h inverse() 
      mytranspose_(tmparea, &s0size, &a0, tmparea2);
      {
        double* const ipiv = stack_->get(a0);
        double* const tmp = stack_->get(a0*a0);
        copy(ovl_[0], ovl_[0]+a0*a0, tmp);
        int info = 0;
        dgesv_(a0, s0size, tmp, a0, (int*)ipiv, tmparea2, a0, info); 
        if (info) throw runtime_error("DGESV failed in SmallNAIBatch::compute");
        stack_->release(a0*a0, tmp);
        stack_->release(a0, ipiv);
      }

      // -1 because <m|p|n>^dagger = -<n|p|m>  (can be proven by integration by part)
      dgemm_("T", "N", s0size, a1, a0, -1.0, tmparea2, a0, nai, a0, 1.0, ints+i*s0size*a1, s0size);
      stack_->release(s0size*a0, tmparea2);
      stack_->release(s0size*a0, tmparea);
    }
  }

  array<double* const,3> data = {{data_+size_block_, data_+size_block_*2, data_+size_block_*3}};
  array<int,3> f = {{1,2,0}};
  array<int,3> b = {{2,0,1}};

  {
    // second half transformation
    shared_ptr<MomentBatch> coeff0(new MomentBatch(array<shared_ptr<const Shell>,2>{{shells_[1]->cartesian_shell(), aux_inc_[1]}}, stack_));
    coeff0->compute();

    shared_ptr<MomentBatch> coeff1;
    if (aux_dec_[1]) {
      coeff1 = shared_ptr<MomentBatch>(new MomentBatch(array<shared_ptr<const Shell>,2>{{shells_[1]->cartesian_shell(), aux_dec_[1]}}, stack_));
      coeff1->compute();
    } else {
      // TODO just to run
      coeff1 = coeff0;
    }

    double* tmparea = stack_->get(s1size * a1);
    double* const tmparea2 = stack_->get(s1size * a1);

    const double* carea0 = coeff0->data();
    const double* carea1 = coeff1->data();
    for (int i = 0; i != 3; ++i, carea0 += coeff0->size_block(), carea1 += coeff1->size_block()) {
      if (shells_[1]->spherical()) {
        const int carsphindex = shells_[1]->angular_number() * ANG_HRR_END + 0; // only transform shell
        const int nloop = shells_[1]->num_contracted() * a1size_inc;
        carsphlist.carsphfunc_call(carsphindex, nloop, carea0, tmparea);
      } else {
        assert(coeff0->size_block() == a1size_inc*s1size);
        copy(carea0, carea0 + coeff0->size_block(), tmparea); 
      }
      if (aux_dec_[1]) {
        if (shells_[1]->spherical()) {
          const int carsphindex = shells_[1]->angular_number() * ANG_HRR_END + 0; // only transform shell
          const int nloop = shells_[1]->num_contracted() * a1size_dec; 
          carsphlist.carsphfunc_call(carsphindex, nloop, carea1, tmparea+a1size_inc*s1size);
        } else {
          assert(coeff1->size_block() == a1size_dec*s1size);
          copy(carea1, carea1+coeff1->size_block(), tmparea+a1size_inc*s1size);
        }
      }

      mytranspose_(tmparea, &s1size, &a1, tmparea2);
      {
        double* const ipiv = stack_->get(a1);
        double* const tmp = stack_->get(a1*a1);
        copy(ovl_[1], ovl_[1]+a1*a1, tmp);
        int info = 0;
        dgesv_(a1, s1size, tmp, a1, (int*)ipiv, tmparea2, a1, info); 
        if (info) throw runtime_error("DGESV failed in SmallNAIBatch::compute");
        stack_->release(a1*a1, tmp);
        stack_->release(a1, ipiv);
      }
      
      // slot in appropriate blocks
      // four blocks are needed.
      // 0) x^x + y^y + z^z
      dgemm_("N", "N", s0size, s1size, a1, 1.0, ints+i*s0size*a1, s0size, tmparea2, a1, 1.0, data_, s0size);
      // 1) x^y - y^x
      // 2) y^z - z^y
      // 3) z^x - x^z
      dgemm_("N", "N", s0size, s1size, a1, 1.0, ints+b[i]*s0size*a1, s0size, tmparea2, a1, 1.0, data[b[i]], s0size);
      dgemm_("N", "N", s0size, s1size, a1, -1.0, ints+f[i]*s0size*a1, s0size, tmparea2, a1, 1.0, data[i], s0size);
    }
    stack_->release(s1size*a1, tmparea2);
    stack_->release(s1size*a1, tmparea);
  }

  stack_->release(3*s0size*a1, ints);
  stack_->release(a0*a1, nai);
#endif
}


void SmallNAIBatch::ovl_compute(const std::array<const int,2> a) {

  for (int i = 0; i < 2; i++) {
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_inc_[i], aux_inc_[i]}}, stack_);
      ovl.compute();
      for (int k = 0; k != aux_inc_[i]->nbasis(); ++k)
        copy_n(ovl.data() + k*(aux_inc_[i]->nbasis()), aux_inc_[i]->nbasis(), ovl_[i] + k*a[i]);
    }
    if (aux_dec_[i]) {
      {
        OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_dec_[i], aux_dec_[i]}}, stack_);
        ovl.compute();
        for (int k = aux_inc_[i]->nbasis(), j = 0; k != a[i]; ++k, ++j) 
          copy_n(ovl.data() + j*(aux_dec_[i]->nbasis()), aux_dec_[i]->nbasis(), ovl_[i] + k*a[i] + aux_inc_[i]->nbasis());
      }
      {
        OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_inc_[i], aux_dec_[i]}}, stack_);
        ovl.compute();
        for (int k = aux_inc_[i]->nbasis(); k != a[i]; ++k)
          for (int j = 0; j != aux_inc_[i]->nbasis(); ++j)
            *(ovl_[i]+j + k*a[i]) = *(ovl_[i]+k + j*a[i]) = *(ovl.data()+j+aux_inc_[i]->nbasis()*(k-aux_inc_[i]->nbasis()));
      }
    }
  }
}

void SmallNAIBatch::nai_compute(const std::shared_ptr<Matrix>& nai) {

  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int a0size_inc = aux_inc_[0]->nbasis();
  const int a1size_inc = aux_inc_[1]->nbasis();
  const int a0size_dec = aux_dec_[0] ? aux_dec_[0]->nbasis() : 0;
  const int a1size_dec = aux_dec_[1] ? aux_dec_[1]->nbasis() : 0;
  const int a0 = a0size_inc + a0size_dec;
  const int a1 = a1size_inc + a1size_dec;

  {
    shared_ptr<OverlapBatch> naic(new OverlapBatch(aux_inc_, stack_));
    naic->compute();

    nai->copy_block(0, 0, a0size_inc, a1size_inc, naic->data());

    //for (int i = 0; i != a1size_inc; ++i)
    //  copy_n(naic->data() + i*a0size_inc, a0size_inc, nai + i*a0); 
  }
  if (aux_dec_[0] && aux_dec_[1]) {
    shared_ptr<OverlapBatch> naic(new OverlapBatch(aux_dec_, stack_));
    naic->compute();

    nai->copy_block(a0size_inc, a1size_inc, a0size_dec, a1size_dec, naic->data());
    //for (int i = a1size_inc, j = 0; i != a1; ++i, ++j)
    //  copy_n(naic->data() + j*a0size_dec, a0size_dec, nai + i*a0 + a0size_inc); 
  }
  if (aux_dec_[0]) {
    shared_ptr<OverlapBatch> naic(new OverlapBatch(array<shared_ptr<const Shell>,2>{{aux_dec_[0],aux_inc_[1]}}, stack_));
    naic->compute();

    nai->copy_block(a0size_inc, 0, a0size_dec, a1size_inc, naic->data());

    //for (int i = 0; i != a1size_inc; ++i)
    //  copy_n(naic->data() + i*a0size_dec, a0size_dec, nai + i*a0 + a0size_inc); 
  }
  if (aux_dec_[1]) {
    shared_ptr<OverlapBatch> naic(new OverlapBatch(array<shared_ptr<const Shell>,2>{{aux_inc_[0],aux_dec_[1]}}, stack_));
    naic->compute();

    nai->copy_block(0, a1size_inc, a0size_inc, a1size_dec, naic->data());

    //for (int i = a1size_inc, j = 0; i != a1; ++i, ++j)
    //  copy_n(naic->data() + j*a0size_inc, a0size_inc, nai + i*a0); 
  }
}

