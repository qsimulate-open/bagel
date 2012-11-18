//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock.cc
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

#include <numeric>
#include <src/df/dfinttask.h>
#include <src/util/taskqueue.h>
#include <src/df/dfblock.h>
#include <src/rysint/libint.h>
#include <src/rysint/eribatch.h>
#include <src/util/f77.h>

using namespace bagel;
using namespace std;


DFBlock::DFBlock(vector<shared_ptr<const Shell> > a, vector<shared_ptr<const Shell> > b1, vector<shared_ptr<const Shell> > b2,
            const int as, const int b1s, const int b2s)
: aux_(a), b1_(b1), b2_(b2),
  asize_(accumulate(a.begin(), a.end(),    0, [](const int& i, const shared_ptr<const Shell>& o) { return i+o->nbasis(); })), 
  b1size_(accumulate(b1.begin(), b1.end(), 0, [](const int& i, const shared_ptr<const Shell>& o) { return i+o->nbasis(); })), 
  b2size_(accumulate(b2.begin(), b2.end(), 0, [](const int& i, const shared_ptr<const Shell>& o) { return i+o->nbasis(); })), 
  astart_(as), b1start_(b1s), b2start_(b2s) {

  // ugly code..
  int tmpa = 0, tmpb1 = 0, tmpb2 = 0;
  for (auto& i : aux_) { aoff_.push_back(tmpa); tmpa += i->nbasis(); }
  for (auto& i : b1_)  { b1off_.push_back(tmpb1); tmpb1 += i->nbasis(); }
  for (auto& i : b2_)  { b2off_.push_back(tmpb2); tmpb2 += i->nbasis(); }
  assert(tmpa == asize_ && tmpb1 == b1size_ && tmpb2 == b2size_);


  ao_init();
}


// protected functions
void DFBlock::ao_init() {
  // allocation of the data area
  data_ = unique_ptr<double[]>(new double[asize_*b1size_*b2size_]);

  const shared_ptr<const Shell> i3(new Shell(b1_.front()->spherical()));

  // making a task list
  vector<DFIntTask> tasks;
  tasks.reserve(b1_.size()*b2_.size()*aux_.size());

  // TODO this is not general, but for the time being I plan to have full basis functions (and limited aux basis functions); 
  assert(b1_ == b2_);
  
  auto j2 = b2off_.begin();
  for (auto& i2 : b2_) { 
    auto j1 = b1off_.begin();
    for (auto& i1 : b1_) { 
      // TODO using symmetry. This assumes that swap(i1, i2) integrals are also located in this block, which might not be the case in general.
      if (*j1 <= *j2) {
        auto j0 = aoff_.begin();
        for (auto& i0 : aux_) { 
          tasks.push_back(DFIntTask(array<shared_ptr<const Shell>,4>{{i3, i0, i1, i2}}, vector<int>{*j2, *j1, *j0}, this));
          ++j0;
        }
      }
      ++j1;
    }
    ++j2;
  }

  TaskQueue<DFIntTask> tq(tasks);
  tq.compute(resources__->max_num_threads());
}


pair<const double*, shared_ptr<RysInt> > DFBlock::compute_batch(array<shared_ptr<const Shell>,4>& input) {
#ifdef LIBINT_INTERFACE
  shared_ptr<Libint> eribatch(new Libint(input));
#else
  shared_ptr<ERIBatch> eribatch(new ERIBatch(input, 2.0));
#endif
  eribatch->compute();
  return make_pair(eribatch->data(), eribatch);
}


shared_ptr<DFBlock> DFBlock::transform_second(const double* const c, const int nocc) const {
  // so far I only consider the following case
  assert(b1start_ == 0);
  unique_ptr<double[]> tmp(new double[asize_*nocc*b2size_]);

  for (size_t i = 0; i != b2size_; ++i)
    dgemm_("N", "N", asize_, nocc, b1size_, 1.0, data_.get()+i*asize_*b1size_, asize_, c, b1size_, 0.0, tmp.get()+i*asize_*nocc, asize_);

  return shared_ptr<DFBlock>(new DFBlock(tmp, asize_, nocc, b2size_, astart_, 0, b2start_));
}


shared_ptr<DFBlock> DFBlock::transform_third(const double* const c, const int nocc) const {
  // so far I only consider the following case
  assert(b2start_ == 0);
  unique_ptr<double[]> tmp(new double[asize_*b1size_*nocc]);

  dgemm_("N", "N", asize_*b1size_, nocc, b2size_, 1.0, data_.get(), asize_*b1size_, c, b2size_, 0.0, tmp.get(), asize_*b1size_);

  return shared_ptr<DFBlock>(new DFBlock(tmp, asize_, b1size_, nocc, astart_, b1start_, 0));
}


shared_ptr<DFBlock> DFBlock::clone() const {
  unique_ptr<double[]> tmp(new double[asize_*b1size_*b2size_]);
  return shared_ptr<DFBlock>(new DFBlock(tmp, asize_, b1size_, b2size_, astart_, b1start_, b2start_));
}


shared_ptr<DFBlock> DFBlock::copy() const {
  unique_ptr<double[]> tmp(new double[asize_*b1size_*b2size_]);
  copy_n(data_.get(), asize_*b1size_*b2size_, tmp.get());
  return shared_ptr<DFBlock>(new DFBlock(tmp, asize_, b1size_, b2size_, astart_, b1start_, b2start_));
}


DFBlock& DFBlock::operator+=(const DFBlock& o) { daxpy( 1.0, o); return *this; }
DFBlock& DFBlock::operator-=(const DFBlock& o) { daxpy(-1.0, o); return *this; }

void DFBlock::daxpy(const double a, const DFBlock& o) {
  if (size() != o.size()) throw logic_error("DFBlock::daxpy called illegally");
  daxpy_(size(), a, o.data_.get(), 1, data_.get(), 1);
}


void DFBlock::scale(const double a) {
  dscal_(size(), a, data_, 1);
}


// TODO not efficient
void DFBlock::symmetrize() {
  if (b1size_ != b2size_) throw logic_error("illegal call of DFBlock::symmetrize()");
  const int n = b1size_;
  for (int i = 0; i != n; ++i)
    for (int j = i; j != n; ++j)
      for (int k = 0; k != asize_; ++k)
        data_[k+asize_*(j+n*i)] = data_[k+asize_*(i+n*j)] = (data_[k+asize_*(j+n*i)] + data_[k+asize_*(i+n*j)]);
}
