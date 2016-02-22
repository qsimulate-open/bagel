//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dfdistt.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/df/dfdistt.h>
#include <src/util/simple.h>

using namespace std;
using namespace bagel;

DFDistT::DFDistT(shared_ptr<const ParallelDF> in, shared_ptr<const StaticDist> dist)
 : naux_(in->naux()), nindex1_(in->nindex1()), nindex2_(in->nindex2()),
   dist_(dist ? dist : make_shared<StaticDist>(nindex1_*nindex2_, mpi__->size())),
   bstart_(dist_->start(mpi__->rank())), bsize_(dist_->size(mpi__->rank())), df_(in->df()) {

  vector<int> srequest;
  const int myrank = mpi__->rank();

  // loop over DFBlocks
  for (auto& iblock : in->block()) {
    // Local matrix (true means local)
    auto dat = make_shared<Matrix>(naux_, bsize_, true);

    // second form a matrix
    shared_ptr<Matrix> buf = dat->clone();

    // source block
    shared_ptr<const DFBlock> source = iblock;

    // information on the data layout
    shared_ptr<const StaticDist> adist = df_->adist_now();

    vector<int> rrequest;
    // first issue all the send and receive requests
    for (int i = 0; i != mpi__->size(); ++i) {
      if (i != myrank) {
        srequest.push_back(mpi__->request_send(source->data()+source->asize()*dist_->start(i), source->asize()*dist_->size(i), i, myrank));
        rrequest.push_back(mpi__->request_recv(buf->data()+adist->start(i)*bsize_, adist->size(i)*bsize_, i, i));
      } else {
        assert(source->asize()*dist_->size(i) == adist->size(i)*bsize_);
        copy_n(source->data()+source->asize()*dist_->start(i), source->asize()*dist_->size(i), buf->data()+adist->start(i)*bsize_);
      }
    }
    for (auto& i : rrequest) mpi__->wait(i);

    TaskQueue<CopyBlockTask<double>> task(mpi__->size());
    for (int i = 0; i != mpi__->size(); ++i)
      task.emplace_back(buf->data()+adist->start(i)*bsize_, adist->size(i), dat->data()+adist->start(i), naux_, adist->size(i), bsize_);
    task.compute();

    data_.push_back(dat);
  }

  for (auto& i : srequest) mpi__->wait(i);

}


DFDistT::DFDistT(const size_t naux, shared_ptr<const StaticDist> dist, const size_t nindex1, const size_t nindex2,
                 const shared_ptr<const ParallelDF> p)
 : naux_(naux), nindex1_(nindex1), nindex2_(nindex2), dist_(dist), bstart_(dist->start(mpi__->rank())), bsize_(dist->size(mpi__->rank())), df_(p) {

  const int nblock = p->block().size();
  for (int i = 0; i != nblock; ++i)
    data_.push_back(make_shared<Matrix>(naux_, bsize_, true));

}


shared_ptr<DFDistT> DFDistT::clone() const {
  return make_shared<DFDistT>(naux_, dist_, nindex1_, nindex2_, df_);
}


shared_ptr<DFDistT> DFDistT::apply_J(shared_ptr<const Matrix> d) const {
  shared_ptr<DFDistT> out = clone();
  auto j = data_.begin();
  for (auto& i : out->data_)
    *i = *d % **j++;
  return out;
}


vector<shared_ptr<Matrix>> DFDistT::form_aux_2index(shared_ptr<const DFDistT> o, const double a) const {
  vector<shared_ptr<Matrix>> out;
  auto i = data_.begin();
  for (auto& j : o->data_) {
    auto tmp = make_shared<Matrix>(**i++ ^ *j);
    tmp->scale(a);
    tmp->allreduce();
    out.push_back(tmp);
  }
  return out;
}


void DFDistT::get_paralleldf(shared_ptr<ParallelDF> out) const {

  vector<int> request;
  const int myrank = mpi__->rank();

  // we need buffer n regions (n is the number of blocks)
  assert(out->block().size() == data_.size());
  vector<shared_ptr<Matrix>> bufv;
  for (auto& i : data_) bufv.push_back(i->clone());

  auto dat = data_.begin();
  auto buf = bufv.begin();

  // loop over target blocks
  for (auto& iblock : out->block()) {
    // first, issue all the receive requests
    for (int i = 0; i != mpi__->size(); ++i)
      if (i != myrank)
        request.push_back(mpi__->request_recv(iblock->data()+iblock->asize()*dist_->start(i), iblock->asize()*dist_->size(i), i, i));

    // information on the data layout
    shared_ptr<const StaticDist> adist = df_->adist_now();

    // copy using threads
    TaskQueue<CopyBlockTask<double>> task(mpi__->size());
    for (int i = 0; i != mpi__->size(); ++i)
      task.emplace_back((*dat)->data()+adist->start(i), naux_, (*buf)->data()+adist->start(i)*bsize_, adist->size(i), adist->size(i), bsize_);
    task.compute(resources__->max_num_threads());

    // last, issue all the send requests
    for (int i = 0; i != mpi__->size(); ++i) {
      if (i != myrank) {
        request.push_back(mpi__->request_send((*buf)->data()+adist->start(i)*bsize_, adist->size(i)*bsize_, i, myrank));
      } else {
        copy_n((*buf)->data()+adist->start(i)*bsize_, iblock->asize()*dist_->size(i), iblock->data()+iblock->asize()*dist_->start(i));
      }
    }
    ++dat;
    ++buf;
  }

  for (auto& i : request) mpi__->wait(i);
}


// function that returns a slice of local data
vector<shared_ptr<Matrix>> DFDistT::get_slice(const int start, const int end) const {
  assert(start >= bstart_ && end <= bstart_+bsize_);
  vector<shared_ptr<Matrix>> out;
  for (auto& i : data_)
    out.push_back(i->slice_copy(start-bstart_, end-bstart_));
  return out;
}


shared_ptr<Matrix> DFDistT::replicate(const int n) const {
  auto out = make_shared<Matrix>(naux_, nindex1_*nindex2_, true);
  vector<int> rq;
  for (int i = 0; i != mpi__->size(); ++i) {
    const int start = dist_->start(i);
    const int size = dist_->size(i);
    if (i == mpi__->rank()) {
      assert(start == bstart_ && size == bsize_ && data_[n]->size() == naux_*size);
      copy_n(data_[n]->data(), data_[n]->size(), out->element_ptr(0, start));
    } else {
      // issue send request to all the nodes
      rq.push_back(mpi__->request_send(data_[n]->data(), data_[n]->size(), i, mpi__->rank()*mpi__->size()+i));
      rq.push_back(mpi__->request_recv(out->element_ptr(0, start), naux_*size, i, i*mpi__->size()+mpi__->rank()));
    }
  }
  for (auto& i : rq)
    mpi__->wait(i);
  return out;
}
