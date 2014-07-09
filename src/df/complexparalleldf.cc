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


// TODO Might want to generalize this for complex<double> a
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


// TODO Might want to generalize this for complex<double> a
void ComplexParallelDF::scale(const double a) {
  throw runtime_error("ComplexParallelDF::scale should be right for real a, but has not been verified");
  for (auto& i : dfdata_[0]->block_)
    i->scale(a);
  for (auto& i : dfdata_[1]->block_)
    i->scale(a);
}


// TODO Using 4-multiplication
shared_ptr<ZMatrix> ComplexParallelDF::form_2index(shared_ptr<const ComplexParallelDF> o, const double a, const bool swap) const {
  if (swap) throw runtime_error("Still need to set up ComplexParallelDF::form_2index with index swapping.  Should need to take a conjugate somewhere.");
  shared_ptr<Matrix> outr = dfdata_[0]->form_2index(o->get_real(), a, swap);
  shared_ptr<Matrix> outi = dfdata_[0]->form_2index(o->get_imag(), a, swap);
  *outr += *dfdata_[1]->form_2index(o->get_imag(), a, swap);
  *outi += *dfdata_[1]->form_2index(o->get_real(), -a, swap);
  return make_shared<ZMatrix> (*outr, *outi);
}


shared_ptr<ZMatrix> ComplexParallelDF::compute_Jop(const shared_ptr<const ZMatrix> den) const {
  return compute_Jop(this->shared_from_this(), den);
}


shared_ptr<ZMatrix> ComplexParallelDF::compute_Jop(const shared_ptr<const ComplexParallelDF> o, const shared_ptr<const ZMatrix> den, const bool onlyonce) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  shared_ptr<const ZVectorB> tmp0 = o->compute_cd(den, data2_real(), onlyonce);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return compute_Jop_from_cd(tmp0);
}


// TODO Using 4-multiplication
shared_ptr<ZVectorB> ComplexParallelDF::compute_cd(const shared_ptr<const ZMatrix> den, shared_ptr<const Matrix> dat2, const bool onlyonce) const {
  const shared_ptr<Matrix> dr = den->get_real_part();
  const shared_ptr<Matrix> di = den->get_imag_part();
  shared_ptr<VectorB> outr = dfdata_[0]->compute_cd(dr, dat2, onlyonce);
  shared_ptr<VectorB> outi = dfdata_[0]->compute_cd(di, dat2, onlyonce);
  *outr -= *dfdata_[1]->compute_cd(di, dat2, onlyonce);
  *outi += *dfdata_[1]->compute_cd(dr, dat2, onlyonce);
  return make_shared<ZVectorB> (*outr, *outi);
}


// TODO Using 4-multiplication
shared_ptr<ZMatrix> ComplexParallelDF::compute_Jop_from_cd(shared_ptr<const ZVectorB> tmp0) const {
  if (nblock_ != 1) throw logic_error("compute_Jop so far assumes block_.size() == 1");
  const shared_ptr<VectorB> tmpr = tmp0->get_real_part();
  const shared_ptr<VectorB> tmpi = tmp0->get_imag_part();
  shared_ptr<Matrix> outr = dfdata_[0]->compute_Jop_from_cd(tmpr);
  shared_ptr<Matrix> outi = dfdata_[0]->compute_Jop_from_cd(tmpi);
  *outr -= *dfdata_[1]->compute_Jop_from_cd(tmpi);
  *outi += *dfdata_[1]->compute_Jop_from_cd(tmpr);
  return make_shared<ZMatrix> (*outr, *outi);
}

