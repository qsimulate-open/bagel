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
#include <src/util/constants.h>
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


shared_ptr<DFBlock> DFBlock::transform_second(const double* const c, const int nocc, const bool trans) const {
  // so far I only consider the following case
  assert(b1start_ == 0);
  unique_ptr<double[]> tmp(new double[asize_*nocc*b2size_]);

  for (size_t i = 0; i != b2size_; ++i) {
    if (!trans)
      dgemm_("N", "N", asize_, nocc, b1size_, 1.0, data_.get()+i*asize_*b1size_, asize_, c, b1size_, 0.0, tmp.get()+i*asize_*nocc, asize_);
    else
      dgemm_("N", "T", asize_, nocc, b1size_, 1.0, data_.get()+i*asize_*b1size_, asize_, c, nocc, 0.0, tmp.get()+i*asize_*nocc, asize_);
  }

  return shared_ptr<DFBlock>(new DFBlock(tmp, asize_, nocc, b2size_, astart_, 0, b2start_));
}


shared_ptr<DFBlock> DFBlock::transform_third(const double* const c, const int nocc, const bool trans) const {
  // so far I only consider the following case
  assert(b2start_ == 0);
  unique_ptr<double[]> tmp(new double[asize_*b1size_*nocc]);

  if (!trans)
    dgemm_("N", "N", asize_*b1size_, nocc, b2size_, 1.0, data_.get(), asize_*b1size_, c, b2size_, 0.0, tmp.get(), asize_*b1size_);
  else  // trans -> back transform
    dgemm_("N", "T", asize_*b1size_, nocc, b2size_, 1.0, data_.get(), asize_*b1size_, c, nocc, 0.0, tmp.get(), asize_*b1size_);

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


shared_ptr<DFBlock> DFBlock::apply_rhf_2RDM() const {
  assert(b1size_ == b2size_);
  const int nocc = b1size_;
  shared_ptr<DFBlock> out = clone();
  out->zero();
  // exchange contributions
  out->daxpy(-2.0, *this); 
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[asize_]);
  fill_n(diagsum.get(), asize_, 0.0);
  for (int i = 0; i != nocc; ++i)
    daxpy_(asize_, 1.0, data_.get()+asize_*(i+nocc*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nocc; ++i)
    daxpy_(asize_, 4.0, diagsum.get(), 1, out->get()+asize_*(i+nocc*i), 1);
  return out;
}


// Caution
//   o strictly assuming that we are using natural orbitals.
//
shared_ptr<DFBlock> DFBlock::apply_uhf_2RDM(const double* amat, const double* bmat) const {
  assert(b1size_ == b2size_);
  const int nocc = b1size_;
  shared_ptr<DFBlock> out = clone();
  {
    unique_ptr<double[]> d2(new double[size()]);
    // exchange contributions
    dgemm_("N", "N", asize_*nocc, nocc, nocc, 1.0, data_.get(), asize_*nocc, amat, nocc, 0.0, d2.get(), asize_*nocc);
    for (int i = 0; i != nocc; ++i)
      dgemm_("N", "N", asize_, nocc, nocc, -1.0, d2.get()+asize_*nocc*i, asize_, amat, nocc, 0.0, out->get()+asize_*nocc*i, asize_);
    dgemm_("N", "N", asize_*nocc, nocc, nocc, 1.0, data_.get(), asize_*nocc, bmat, nocc, 0.0, d2.get(), asize_*nocc);
    for (int i = 0; i != nocc; ++i)
      dgemm_("N", "N", asize_, nocc, nocc, -1.0, d2.get()+asize_*nocc*i, asize_, bmat, nocc, 1.0, out->get()+asize_*nocc*i, asize_);
  }

  unique_ptr<double[]> sum(new double[nocc]);
  for (int i = 0; i != nocc; ++i) sum[i] = amat[i+i*nocc] + bmat[i+i*nocc];
  // coulomb contributions (diagonal to diagonal)
  unique_ptr<double[]> diagsum(new double[asize_]);
  fill_n(diagsum.get(), asize_, 0.0);
  for (int i = 0; i != nocc; ++i)
    daxpy_(asize_, sum[i], data_.get()+asize_*(i+nocc*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nocc; ++i)
    daxpy_(asize_, sum[i], diagsum.get(), 1, out->get()+asize_*(i+nocc*i), 1);
  return out;
}



shared_ptr<DFBlock> DFBlock::apply_2RDM(const double* rdm, const double* rdm1, const int nclosed, const int nact) const {
  assert(nclosed+nact == b1size_ && b1size_ == b2size_);
  // checking if natural orbitals...
  {
    const double a = ddot_(nact*nact, rdm1, 1, rdm1, 1);
    double sum = 0.0;
    for (int i = 0; i != nact; ++i) sum += rdm1[i+nact*i]*rdm1[i+nact*i];
    if (fabs(a-sum) > numerical_zero__) throw logic_error("DF_Full::apply_2rdm should be called with natural orbitals");
  }
  shared_ptr<DFBlock> out = clone();
  out->zero();
  // closed-closed part
  // exchange contribution
  for (int i = 0; i != nclosed; ++i)
    for (int j = 0; j != nclosed; ++j)
      daxpy_(asize_, -2.0, data_.get()+asize_*(j+b1size_*i), 1, out->get()+asize_*(j+b1size_*i), 1);
  // coulomb contribution
  unique_ptr<double[]> diagsum(new double[asize_]);
  fill_n(diagsum.get(), asize_, 0.0);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize_, 1.0, data_.get()+asize_*(i+b1size_*i), 1, diagsum.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize_, 4.0, diagsum.get(), 1, out->get()+asize_*(i+b1size_*i), 1);

  // act-act part
  // compress
  unique_ptr<double[]> buf(new double[nact*nact*asize_]);
  unique_ptr<double[]> buf2(new double[nact*nact*asize_]);
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      dcopy_(asize_, data_.get()+asize_*(j+nclosed+b1size_*(i+nclosed)), 1, buf.get()+asize_*(j+nact*i),1);
  // multiply
  dgemm_("N", "N", asize_, nact*nact, nact*nact, 1.0, buf.get(), asize_, rdm, nact*nact, 0.0, buf2.get(), asize_);
  // slot in
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      dcopy_(asize_, buf2.get()+asize_*(j+nact*i),1, out->get()+asize_*(j+nclosed+b1size_*(i+nclosed)), 1);

  // closed-act part
  // coulomb contribution G^ia_ia = 2*gamma_ab
  // ASSUMING natural orbitals
  for (int i = 0; i != nact; ++i)
    daxpy_(asize_, 2.0*rdm1[i+nact*i], diagsum.get(), 1, out->get()+asize_*(i+nclosed+b1size_*(i+nclosed)), 1);
  unique_ptr<double[]> diagsum2(new double[asize_]);
  dgemv_("N", asize_, nact*nact, 1.0, buf.get(), asize_, rdm1, 1, 0.0, diagsum2.get(), 1);
  for (int i = 0; i != nclosed; ++i)
    daxpy_(asize_, 2.0, diagsum2.get(), 1, out->get()+asize_*(i+b1size_*i), 1);
  // exchange contribution
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nclosed; ++j) {
      daxpy_(asize_, -rdm1[i+nact*i], data_.get()+asize_*(j+b1size_*(i+nclosed)), 1, out->get()+asize_*(j+b1size_*(i+nclosed)), 1);
      daxpy_(asize_, -rdm1[i+nact*i], data_.get()+asize_*(i+nclosed+b1size_*j), 1, out->get()+asize_*(i+nclosed+b1size_*j), 1);
    }
  }
  return out;
}


shared_ptr<DFBlock> DFBlock::apply_2RDM(const double* rdm) const {
  shared_ptr<DFBlock> out = clone();
  dgemm_("N", "T", asize_, b1size_*b2size_, b1size_*b2size_, 1.0, data_.get(), asize_, rdm, b1size_*b2size_, 0.0, out->get(), asize_);
  return out;
}


unique_ptr<double[]> DFBlock::form_2index(const shared_ptr<const DFBlock> o, const double a) const {
  if (asize_ != o->asize_ || (b1size_ != o->b1size_ && b2size_ != o->b2size_)) throw logic_error("illegal call of DFBlock::form_2index");
  unique_ptr<double[]> target(new double[b2size_*o->b2size_]);

  if (b1size_ == o->b1size_) {
    dgemm_("T", "N", b2size_, o->b2size_, asize_*b1size_, a, data_.get(), asize_*b1size_, o->data_.get(), asize_*b1size_, 0.0, target.get(), b2size_);
  } else {
    assert(b2size_ == o->b2size_);
    fill_n(target.get(), b2size_*o->b2size_, 0.0);
    for (int i = 0; i != b2size_; ++i)
      dgemm_("T", "N", b1size_, o->b1size_, asize_, a, data_.get()+i*asize_*b1size_, asize_, o->data_.get()+i*asize_*o->b1size_, asize_, 1.0, target.get(), b1size_);
  }

  return target;
}


unique_ptr<double[]> DFBlock::form_4index(const shared_ptr<const DFBlock> o, const double a) const {
  if (asize_ != o->asize_) throw logic_error("illegal call of DFBlock::form_4index");
  unique_ptr<double[]> target(new double[b2size_*o->b2size_*b1size_*o->b1size_]);
  dgemm_("T", "N", b1size_*b2size_, o->b1size_*o->b2size_, asize_, a, data_.get(), asize_, o->data_.get(), asize_, 0.0, target.get(), b1size_*b2size_);
  return target;
}


shared_ptr<Matrix> DFBlock::form_aux_2index(const shared_ptr<const DFBlock> o, const double a) const {
  if (b1size_ != o->b1size_ || b2size_ != o->b2size_) throw logic_error("illegal call of DFBlock::form_aux_2index");
  shared_ptr<Matrix> target(new Matrix(asize_, o->asize_));
  dgemm_("N", "T", asize_, o->asize_, b1size_*b2size_, 1.0, data_.get(), asize_, o->data_.get(), o->asize_, 0.0, target->data(), asize_);
  return target;
}
