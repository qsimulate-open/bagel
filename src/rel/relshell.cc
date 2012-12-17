//
// BAGEL - Parallel electron correlation program.
// Filename: relshell.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Matthew Kelley and Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <array>
#include <src/rel/relshell.h>
#include <src/osint/momentbatch.h>
#include <src/rysint/carsphlist.h>
#include <src/osint/overlapbatch.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


RelShell::RelShell(const shared_ptr<const Shell> o) : Shell(*o), aux_inc_(kinetic_balance_uncont(1)), aux_dec_(kinetic_balance_uncont(-1)) {

  // overlap = S^-1 between auxiliary functions
  shared_ptr<const Matrix> overlap = overlap_compute_();

  // small is a transformation matrix (x,y,z components)
  small_ = moment_compute_(overlap);
}


shared_ptr<const Matrix> RelShell::overlap_compute_() const {

  const int asize_inc = aux_inc_->nbasis();
  const int asize_dec = aux_dec_ ? aux_dec_->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  shared_ptr<Matrix> overlap(new Matrix(a,a));

  {
    OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_inc_, aux_inc_}});
    ovl.compute();
    for (int i = 0; i != aux_inc_->nbasis(); ++i)
      copy_n(ovl.data() + i*(aux_inc_->nbasis()), aux_inc_->nbasis(), overlap->data() + i*a);
  }
  if (aux_dec_) {
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_dec_, aux_dec_}});
      ovl.compute();
      for (int i = aux_inc_->nbasis(), j = 0; i != a; ++i, ++j) 
        copy_n(ovl.data() + j*(aux_dec_->nbasis()), aux_dec_->nbasis(), overlap->data() + i*a + aux_inc_->nbasis());
    }
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_inc_, aux_dec_}});
      ovl.compute();
      for (int i = aux_inc_->nbasis(); i != a; ++i)
        for (int j = 0; j != aux_inc_->nbasis(); ++j)
          overlap->element(j,i) = overlap->element(i,j) = *(ovl.data()+j+aux_inc_->nbasis()*(i-aux_inc_->nbasis()));
    }
  }

  // Make an inverse
  overlap->inverse();
  return overlap;
}


array<shared_ptr<const Matrix>,3> RelShell::moment_compute_(const shared_ptr<const Matrix> overlap) const {
  const int ssize = nbasis();
  const int asize_inc = aux_inc_->nbasis();
  const int asize_dec = aux_dec_ ? aux_dec_->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  shared_ptr<MomentBatch> coeff0(new MomentBatch(array<shared_ptr<const Shell>,2>{{cartesian_shell(), aux_inc_}}));
  coeff0->compute();

  shared_ptr<MomentBatch> coeff1;
  if (aux_dec_) {
    coeff1 = shared_ptr<MomentBatch>(new MomentBatch(array<shared_ptr<const Shell>,2>{{cartesian_shell(), aux_dec_}}));
    coeff1->compute();
  } else {
    // TODO just to run
    coeff1 = coeff0;
  }

  const double* carea0 = coeff0->data();
  const double* carea1 = coeff1->data();

  shared_ptr<Matrix> tmparea(new Matrix(ssize,a));
  array<shared_ptr<const Matrix>,3> out;

  for (int i = 0; i != 3; ++i, carea0 += coeff0->size_block(), carea1 += coeff1->size_block()) {
    if (spherical_) {
      const int carsphindex = angular_number_ * ANG_HRR_END + 0; // only transform shell
      const int nloop = num_contracted() * asize_inc; 
      carsphlist.carsphfunc_call(carsphindex, nloop, carea0, tmparea->data());
    } else {
      assert(coeff0->size_block() == asize_inc*ssize);
      copy(carea0, carea0+coeff0->size_block(), tmparea->data());
    }
    if (aux_dec_) {
      if (spherical_) {
        const int carsphindex = angular_number_ * ANG_HRR_END + 0; // only transform shell
        const int nloop = num_contracted() * asize_dec; 
        carsphlist.carsphfunc_call(carsphindex, nloop, carea1, tmparea->data()+asize_inc*ssize);
      } else {
        assert(coeff1->size_block() == asize_dec*ssize);
        copy(carea1, carea1+coeff1->size_block(), tmparea->data()+asize_inc*ssize);
        }
    }

    out[i] = shared_ptr<const Matrix>(new Matrix(*overlap ^ *tmparea));
  }
  return out;
}
