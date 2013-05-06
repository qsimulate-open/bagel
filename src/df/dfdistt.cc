//
// BAGEL - Parallel electron correlation program.
// Filename: dfblockt.cc
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

#include <src/df/dfdistt.h>
#include <src/util/simple.h>

using namespace std;
using namespace bagel;

DFDistT::DFDistT(std::shared_ptr<const ParallelDF> in)
 : naux_(in->naux()), nindex1_(in->nindex1()), nindex2_(in->nindex2()), df_(in->df()) {

#ifndef HAVE_MPI_H
  assert(false); // this class should be used with MPI (it works without MPI, but it just transposition of data...)
#endif

  // determine how to distribute...
  const size_t stride = nindex1_*nindex2_ / mpi__->size();
  for (int i = 0; i != mpi__->size(); ++i) {
    tabstart_.push_back(stride*i);
    tabsize_.push_back((i+1 == mpi__->size()) ? nindex1_*nindex2_-tabstart_[i] : stride);
  }
  assert(tabsize_.back() >= 0);
  start_ = tabstart_[mpi__->rank()];
  size_ = tabsize_[mpi__->rank()];

  vector<int> srequest;

  // loop over DFBlocks
  for (auto& iblock : in->block()) {
    // Local matrix (true means local)
    auto dat = make_shared<Matrix>(naux_, size_, true);

    // second form a matrix
    shared_ptr<Matrix> buf = dat->clone();

    // source block
    shared_ptr<const DFBlock> source = iblock;

    // information on the data layout
    vector<pair<size_t, size_t>> atab = df_->adist_now()->atable();

    vector<int> rrequest;
    // first issue all the send and receive requests
    for (int i = 0; i != mpi__->size(); ++i) {
      if (i != mpi__->rank()) {
        srequest.push_back(mpi__->request_send(source->get()+source->asize()*tabstart_[i], source->asize()*tabsize_[i], i, mpi__->rank()));
        rrequest.push_back(mpi__->request_recv(buf->data()+atab[i].first*size_, atab[i].second*size_, i, i));
      } else {
        assert(source->asize()*tabsize_[i] == atab[i].second*size_);
        copy_n(source->get()+source->asize()*tabstart_[i], source->asize()*tabsize_[i], buf->data()+atab[i].first*size_); 
      }
    }
    for (auto& i : rrequest) mpi__->wait(i);

    vector<CopyBlockTask> task;
    task.reserve(mpi__->size());
    for (int i = 0; i != mpi__->size(); ++i)
      task.push_back(CopyBlockTask(buf->data()+atab[i].first*size_, atab[i].second, dat->data()+atab[i].first, naux_, atab[i].second, size_));
    TaskQueue<CopyBlockTask> tq(task);
    tq.compute(resources__->max_num_threads());

    data_.push_back(dat);
  }

  for (auto& i : srequest) mpi__->wait(i);

}


DFDistT::DFDistT(const size_t naux, const vector<size_t> start, const vector<size_t> size, const size_t nindex1, const size_t nindex2,
                 const shared_ptr<const ParallelDF> p)
 : naux_(naux), nindex1_(nindex1), nindex2_(nindex2), start_(start[mpi__->rank()]), size_(size[mpi__->rank()]), tabstart_(start), tabsize_(size), df_(p) {

  const int nblock = p->block().size();
  for (int i = 0; i != nblock; ++i)
    data_.push_back(make_shared<Matrix>(naux_, size_, true));

}


shared_ptr<DFDistT> DFDistT::clone() const {
  return make_shared<DFDistT>(naux_, tabstart_, tabsize_, nindex1_, nindex2_, df_); 
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
    tmp->allreduce();
    out.push_back(tmp);
  }
  return out;
}


void DFDistT::get_paralleldf(std::shared_ptr<ParallelDF> out) const {

  vector<int> request;

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
      if (i != mpi__->rank())
        request.push_back(mpi__->request_recv(iblock->get()+iblock->asize()*tabstart_[i], iblock->asize()*tabsize_[i], i, i));

    // information on the data layout
    vector<pair<size_t, size_t>> atab = df_->adist_now()->atable();

    // copy using threads
    vector<CopyBlockTask> task;
    task.reserve(mpi__->size());
    for (int i = 0; i != mpi__->size(); ++i)
      task.push_back(CopyBlockTask((*dat)->data()+atab[i].first, naux_, (*buf)->data()+atab[i].first*size_, atab[i].second, atab[i].second, size_));
    TaskQueue<CopyBlockTask> tq(task);
    tq.compute(resources__->max_num_threads());

    // last, issue all the send requests
    for (int i = 0; i != mpi__->size(); ++i) {
      if (i != mpi__->rank()) {
        request.push_back(mpi__->request_send((*buf)->data()+atab[i].first*size_, atab[i].second*size_, i, mpi__->rank()));
      } else {
        copy_n((*buf)->data()+atab[i].first*size_, out->block(0)->asize()*tabsize_[i], out->block(0)->get()+out->block(0)->asize()*tabstart_[i]);
      }
    }
    ++dat;
    ++buf;
  }

  for (auto& i : request) mpi__->wait(i);
}
