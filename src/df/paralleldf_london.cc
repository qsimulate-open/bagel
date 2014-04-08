//
// BAGEL - Parallel electron correlation program.
// Filename: paralleldf_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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

#include <src/df/paralleldf_london.h>

using namespace std;
using namespace bagel;


ParallelDF_London::ParallelDF_London(const size_t naux, const size_t nb1, const size_t nb2, shared_ptr<const ParallelDF_London> df, shared_ptr<ZMatrix> dat)
 : naux_(naux), nindex1_(nb1), nindex2_(nb2), df_(df), data2_(dat), serial_(df ? df->serial_ : false) { }


shared_ptr<ZMatrix> ParallelDF_London::form_2index(shared_ptr<const ParallelDF_London> o, const double a, const bool swap) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<ZMatrix> out = (!swap) ? block_[0]->form_2index(o->block_[0], a) : o->block_[0]->form_2index(block_[0], a);
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<ZMatrix> ParallelDF_London::form_4index(shared_ptr<const ParallelDF_London> o, const double a, const bool swap) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<ZMatrix> out = (!swap) ? block_[0]->form_4index(o->block_[0], a) : o->block_[0]->form_4index(block_[0], a);

  // all reduce
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<ZMatrix> ParallelDF_London::form_aux_2index(shared_ptr<const ParallelDF_London> o, const double a) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
#ifdef HAVE_MPI_H
  if (!serial_) {
    auto work = make_shared<DFDistT>(this->shared_from_this());
    auto work2 = make_shared<DFDistT>(o);
    return work->form_aux_2index(work2, a).front();
  } else {
#else
  {
#endif
    return block_[0]->form_aux_2index(o->block_[0], a);
  }
}


void ParallelDF_London::ax_plus_y(const double a, const shared_ptr<const ParallelDF_London> o) {
  assert(block_.size() == o->block_.size());
  auto j = o->block_.begin();
  for (auto& i : block_)
    i->ax_plus_y(a, *j++);
}


void ParallelDF_London::scale(const double a) {
  for (auto& i : block_)
    i->scale(a);
}


void ParallelDF_London::symmetrize() {
  for (auto& i : block_)
    i->symmetrize();
}


void ParallelDF_London::add_block(shared_ptr<DFBlock_London> o) {
  block_.push_back(o);
}


shared_ptr<ZMatrix> ParallelDF_London::get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const {
  if (block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  // first thing is to find the node
  tuple<size_t, size_t> info = adist_now()->locate(i);

  // date has to be localised in this node
  if (get<0>(info) == mpi__->rank() && !block_[0]->averaged()) {
    return block_[0]->get_block(i, id, j, jd, k, kd);
  } else {
    throw logic_error("ParallelDF_London::get_block is an intra-node function (or bug?)");
  }
  return nullptr;
}


shared_ptr<ZMatrix> ParallelDF_London::compute_Jop_from_cd(shared_ptr<const ZMatrix> tmp0) const {
  if (block_.size() != 1) throw logic_error("compute_Jop so far assumes block_.size() == 1");
  shared_ptr<ZMatrix> out = block_[0]->form_mat(tmp0->data()+block_[0]->astart());
  // all reduce
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<ZMatrix> ParallelDF_London::compute_cd(const shared_ptr<const ZMatrix> den, shared_ptr<const ZMatrix> dat2, const bool onlyonce) const {
  if (!dat2 && !data2_) throw logic_error("ParallelDF::compute_cd was called without 2-index integrals");
  if (!dat2) dat2 = data2_;

  auto tmp0 = make_shared<ZMatrix>(naux_, 1, true);

  // D = (D|rs)*d_rs
  if (block_.size() != 1) throw logic_error("compute_Jop so far assumes block_.size() == 1");
  unique_ptr<complex<double>[]> tmp = block_[0]->form_vec(den);
  copy_n(tmp.get(), block_[0]->asize(), tmp0->data()+block_[0]->astart());
  // All reduce
  if (!serial_)
    tmp0->allreduce();

  tmp0 = make_shared<ZMatrix>(*dat2 * *tmp0);
  if (!onlyonce)
    tmp0 = make_shared<ZMatrix>(*dat2 * *tmp0);
  return tmp0;
}


shared_ptr<ZMatrix> ParallelDF_London::compute_Jop(const shared_ptr<const ZMatrix> den) const {
  return compute_Jop(this->shared_from_this(), den);
}


shared_ptr<ZMatrix> ParallelDF_London::compute_Jop(const shared_ptr<const ParallelDF_London> o, const shared_ptr<const ZMatrix> den, const bool onlyonce) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  shared_ptr<const ZMatrix> tmp0 = o->compute_cd(den, data2_, onlyonce);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return compute_Jop_from_cd(tmp0);
}

