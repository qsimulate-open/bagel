//
// BAGEL - Parallel electron correlation program.
// Filename: recvrequest.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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

#include <src/parallel/recvrequest.h>
#include <src/util/f77.h>
#include <src/util/constants.h>

using namespace std;
using namespace bagel;


// PutRequest receives probe and returns data
PutRequest::PutRequest(const double* d) : data_(d) {
  for (size_t i = 0; i != pool_size__; ++i)
    init();
  mpi__->barrier();
  turn_on();
}


void PutRequest::init() {
  // we should compute the number of messnages to receive?
  // receives
  unique_ptr<size_t[]> buf(new size_t[4]);
  // receives size,tag,rank
  const int rq = mpi__->request_recv(buf.get(), 4, -1, probe_key2__);
  {
    lock_guard<mutex> lock(block_);
    auto m = calls_.insert(make_pair(rq, move(buf)));
    assert(m.second);
  }
}


PutRequest::~PutRequest() {
  mpi__->soft_barrier();
  turn_off();
  for (auto& i : calls_)
    mpi__->cancel(i.first);
}


void PutRequest::flush() {
  size_t cnt = 0;
  {
    lock_guard<mutex> lock(block_);
    for (auto i = calls_.begin(); i != calls_.end(); ) {
      // if this has already arrived, send data 
      if (mpi__->test(i->first)) {
        const int size = i->second[0];
        const int tag  = i->second[1];
        const int rank = i->second[2];
        const int off  = i->second[3];
        mpi__->request_send(data_+off, size, rank, tag);
        i = calls_.erase(i);
        ++cnt;
      } else {
        ++i;
      }
    }
  }
  for (int i = 0; i != cnt; ++i)
    init();
}


/////////////////////////


RecvRequest::RecvRequest() : counter_(probe_key2__ + mpi__->rank() + 1) {

}


int RecvRequest::request_recv(double* buf, const size_t size, const int dest, const size_t off) {
  lock_guard<mutex> lock(block_);
  // sending size
  shared_ptr<Probe> p(new Probe(size, counter_, mpi__->rank(), dest, off, buf));
  counter_ += mpi__->size();
  const int srq = mpi__->request_send(p->size, 4, dest, probe_key2__);
  probe_.push_back(srq);
  const int rrq = mpi__->request_recv(p->buf, p->size[0], dest, p->tag);
  auto m = request_.insert(make_pair(rrq, p));
  assert(m.second);
  return rrq;
}


bool RecvRequest::test() {
  lock_guard<mutex> lock(block_);
  bool done = true;
  for (auto i = request_.begin(); i != request_.end(); ) {
    if (mpi__->test(i->first)) {
      i = request_.erase(i);
    } else {
      done = false;
      ++i;
    }
  }
  return done;
}

