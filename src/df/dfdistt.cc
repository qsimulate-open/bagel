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

#include <src/util/f77.h>
#include <src/parallel/mpi_interface.h>
#include <src/df/dfdistt.h>

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

  data_ = unique_ptr<double[]>(new double[naux_*size_]);

  // second form a matrix
  unique_ptr<double[]> buf(new double[naux_*size_]);

  // source block
  shared_ptr<const DFBlock> source = in->block();

  // information on the data layout
  vector<pair<int, int> > atab = df_->atable();

  vector<int> srequest, rrequest;
  // first issue all the send and receive requests
  for (int i = 0; i != mpi__->size(); ++i) {
    if (i != mpi__->rank()) {
      srequest.push_back(mpi__->request_send(source->get()+source->asize()*tabstart_[i], source->asize()*tabsize_[i], i));
      rrequest.push_back(mpi__->request_recv(buf.get()+atab[i].first*size_, atab[i].second*size_, i));
    } else {
      assert(source->asize()*tabsize_[i] == atab[i].second*size_);
      copy_n(source->get()+source->asize()*tabstart_[i], source->asize()*tabsize_[i], buf.get()+atab[i].first*size_); 
    }
  }
  for (auto& i : rrequest) mpi__->wait(i);

  // second transpose each block
  for (int i = 0; i != mpi__->size(); ++i) {
    const int o = atab[i].first*size_;
    mytranspose_(buf.get()+o, atab[i].second, size_, data_.get()+o);
  }

  for (auto& i : srequest) mpi__->wait(i);

}


DFDistT::DFDistT(const size_t naux, const vector<int> start, const vector<int> size, const size_t nindex1, const size_t nindex2,
                 const shared_ptr<const ParallelDF> p)
 : data_(new double[naux*size[mpi__->rank()]]), naux_(naux), nindex1_(nindex1), nindex2_(nindex2), start_(start[mpi__->rank()]), size_(size[mpi__->rank()]),
   tabstart_(start), tabsize_(size), df_(p) {

}


shared_ptr<DFDistT> DFDistT::clone() const {
  shared_ptr<DFDistT> out(new DFDistT(naux_, tabstart_, tabsize_, nindex1_, nindex2_, df_)); 
  return out;
}


shared_ptr<DFDistT> DFDistT::apply_J(shared_ptr<const Matrix> d) const {
  shared_ptr<DFDistT> out = clone();
  dgemm_("N", "N", size_, naux_, naux_, 1.0, data_.get(), size_, d->data(), naux_, 0.0, out->data_.get(), size_); 
  return out;
}


void DFDistT::get_paralleldf(std::shared_ptr<ParallelDF> out) const {
  // first, issue all the receive requests
  vector<int> request;
  for (int i = 0; i != mpi__->size(); ++i)
    if (i != mpi__->rank())
      request.push_back(mpi__->request_recv(out->block()->get()+out->block()->asize()*tabstart_[i], out->block()->asize()*tabsize_[i], i));

  // second form a matrix
  unique_ptr<double[]> buf(new double[naux_*size_]);

  // information on the data layout
  vector<pair<int, int> > atab = df_->atable();

  // transpose each block back
  for (int i = 0; i != mpi__->size(); ++i) {
    const int o = atab[i].first*size_;
    mytranspose_(data_.get()+o, size_, atab[i].second, buf.get()+o);
  }

  // last, issue all the send requests
  for (int i = 0; i != mpi__->size(); ++i) {
    if (i != mpi__->rank()) {
      request.push_back(mpi__->request_send(buf.get()+atab[i].first*size_, atab[i].second*size_, i));
    } else {
      copy_n(buf.get()+atab[i].first*size_, out->block()->asize()*tabsize_[i], out->block()->get()+out->block()->asize()*tabstart_[i]);
    }
  }
  for (auto& i : request) mpi__->wait(i);
}
