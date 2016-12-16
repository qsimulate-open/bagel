//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: paralleldf.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/df/paralleldf.h>
#include <src/df/dfdistt.h>

using namespace std;
using namespace bagel;


ParallelDF::ParallelDF(const size_t naux, const size_t nb1, const size_t nb2, shared_ptr<const ParallelDF> df, shared_ptr<Matrix> dat, const bool serial)
 : naux_(naux), nindex1_(nb1), nindex2_(nb2), df_(df), data2_(dat), serial_(df ? df->serial_ : serial) {

}


shared_ptr<Matrix> ParallelDF::form_2index(shared_ptr<const ParallelDF> o, const double a, const bool swap) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<Matrix> out = (!swap) ? block_[0]->form_2index(o->block_[0], a) : o->block_[0]->form_2index(block_[0], a);
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<Matrix> ParallelDF::form_4index(shared_ptr<const ParallelDF> o, const double a, const bool swap) const {
  if (block_.size() != 1 || o->block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  shared_ptr<Matrix> out = (!swap) ? block_[0]->form_4index(o->block_[0], a) : o->block_[0]->form_4index(block_[0], a);

  // all reduce
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<Matrix> ParallelDF::form_aux_2index(shared_ptr<const ParallelDF> o, const double a) const {
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


void ParallelDF::add_direct_product(const vector<shared_ptr<const VectorB>> cd, const vector<shared_ptr<const Matrix>> dd, const double a) {
  if (block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  if (cd.size() != dd.size()) throw logic_error("Illegal call of ParallelDF::DFDist");

  auto d = dd.begin();
  for (auto& c : cd) {
    const VecView aslice = c->slice(block_[0]->astart(), block_[0]->astart()+block_[0]->asize());
    block_[0]->add_direct_product(aslice, **d++, a);
  }
  assert(d == dd.end());
}


void ParallelDF::ax_plus_y(const double a, const shared_ptr<const ParallelDF> o) {
  assert(block_.size() == o->block_.size());
  auto j = o->block_.begin();
  for (auto& i : block_)
    i->ax_plus_y(a, *j++);
}


void ParallelDF::scale(const double a) {
  for (auto& i : block_)
    i->scale(a);
}


void ParallelDF::symmetrize() {
  for (auto& i : block_)
    i->symmetrize();
}


void ParallelDF::add_block(shared_ptr<DFBlock> o) {
  block_.push_back(o);
}


shared_ptr<btas::Tensor3<double>> ParallelDF::get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const {
  if (block_.size() != 1) throw logic_error("so far assumes block_.size() == 1");
  // first thing is to find the node
  tuple<size_t, size_t> info = adist_now()->locate(i);

  // date has to be localised in this node
  if (get<0>(info) == mpi__->rank()) {
    return block_[0]->get_block(i, id, j, jd, k, kd);
  } else {
    throw logic_error("ParallelDF::get_block is an intra-node function (or bug?)");
  }
  return nullptr;
}


shared_ptr<Matrix> ParallelDF::compute_Jop_from_cd(shared_ptr<const VectorB> tmp0) const {
  if (block_.size() != 1) throw logic_error("compute_Jop so far assumes block_.size() == 1");
  shared_ptr<Matrix> out = block_[0]->form_mat(tmp0->slice(block_[0]->astart(), block_[0]->astart()+block_[0]->asize()));
  // all reduce
  if (!serial_)
    out->allreduce();
  return out;
}


shared_ptr<VectorB> ParallelDF::compute_cd(const shared_ptr<const Matrix> den, shared_ptr<const Matrix> dat2, const int number_of_j) const {
  if (!dat2 && !data2_) throw logic_error("ParallelDF::compute_cd was called without 2-index integrals");
  if (!dat2) dat2 = data2_;

  auto tmp0 = make_shared<VectorB>(naux_);

  // D = (D|rs)*d_rs
  if (block_.size() != 1) throw logic_error("compute_Jop so far assumes block_.size() == 1");
  shared_ptr<VectorB> tmp = block_[0]->form_vec(den);
  copy_n(tmp->data(), block_[0]->asize(), tmp0->data()+block_[0]->astart());
  // All reduce
  if (!serial_)
    tmp0->allreduce();

  if (number_of_j == 1)
    *tmp0 = *dat2 * *tmp0;
  else if (number_of_j == 2)
    *tmp0 = *dat2 * (*dat2 * *tmp0);
  else if (number_of_j != 0)
    throw logic_error("wrong number of J in ParallelDF::compute_cd");
  return tmp0;
}


shared_ptr<Matrix> ParallelDF::compute_Jop(const shared_ptr<const Matrix> den) const {
  return compute_Jop(this->shared_from_this(), den);
}


shared_ptr<Matrix> ParallelDF::compute_Jop(const shared_ptr<const ParallelDF> o, const shared_ptr<const Matrix> den, const bool onlyonce) const {
  // first compute |E*) = d_rs (D|rs) J^{-1}_DE
  shared_ptr<const VectorB> tmp0 = o->compute_cd(den, data2_, onlyonce ? 1 : 2);
  // then compute J operator J_{rs} = |E*) (E|rs)
  return compute_Jop_from_cd(tmp0);
}

