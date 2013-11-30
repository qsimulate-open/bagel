//
// BAGEL - Parallel electron correlation program.
// Filename: accrequest.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/parallel/accrequest.h>
#include <src/util/f77.h>
#include <src/util/constants.h>

using namespace std;
using namespace bagel;

// SendRequest sends buffer using MPI_send. When completed, releases the buffer. Note that
SendRequest::SendRequest() : counter_(probe_key__+mpi__->rank()+1) {
  turn_on();
}


SendRequest::~SendRequest() {
  turn_off();
}


void SendRequest::request_send(unique_ptr<double[]> buf, const size_t size, const int dest, const size_t off) {
  lock_guard<mutex> lock(mutex_);
  // sending size
  auto p = make_shared<Probe>(size, counter_, mpi__->rank(), dest, off, move(buf));
  counter_ += mpi__->size();
  mpi__->request_send(p->size,    4, dest, probe_key__);
  const int rrq = mpi__->request_recv(&p->target, 1, dest, p->tag);
  auto m = inactive_.insert(make_pair(rrq, p));
  assert(m.second);
}


void SendRequest::flush_() {
  lock_guard<mutex> lock(mutex_);

  // if receive buffer at the destination is created, send the message
  for (auto i = inactive_.begin(); i != inactive_.end(); ) {
    if (mpi__->test(i->first)) {
      shared_ptr<Probe> p = i->second;
      assert(round(p->target) == p->tag);
      const int srq = mpi__->request_send(p->buf.get(), p->size[0], p->targetrank, p->tag);
      auto m = requests_.insert(make_pair(srq, p));
      assert(m.second);
      i = inactive_.erase(i);
    } else {
      ++i;
    }
  }

  // if the message is sent, delete the request
  for (auto i = requests_.begin(); i != requests_.end(); ) {
    if (mpi__->test(i->first))
      i = requests_.erase(i);
    else
      ++i;
  }
}


bool SendRequest::test() {
  lock_guard<mutex> lock(mutex_);
  return requests_.empty() && inactive_.empty();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////

AccRequest::AccRequest(double* const d, vector<mutex>* m) : data_(d), datamutex_(m) {
  for (size_t i = 0; i != pool_size__; ++i)
    init();
  turn_on();
}


AccRequest::~AccRequest() {
  turn_off();
  for (auto& i : calls_)
    mpi__->cancel(i.first);
}


void AccRequest::init() {
  // receives
  auto call = make_shared<Call>();
  // receives size,tag,rank
  const int rq = mpi__->request_recv(call->buf.get(), 4, -1, probe_key__);
  {
    lock_guard<mutex> lock(mutex_);
    auto m = calls_.insert(make_pair(rq, call));
    assert(m.second);
  }
}


void AccRequest::flush_() {
  size_t cnt = 0;
  {
    lock_guard<mutex> lock(mutex_);
    for (auto i = calls_.begin(); i != calls_.end(); ) {
      // if this has already arrived, create a buffer, and return the address.
      if (mpi__->test(i->first)) {
        const size_t size = i->second->buf[0];
        const size_t tag  = i->second->buf[1];
        const size_t rank = i->second->buf[2];
        const size_t off  = i->second->buf[3];
        // allocating buffer
        unique_ptr<double[]> buffer(new double[size]);
        buffer[0] = tag;
        // sending back the buffer address
        mpi__->request_send(buffer.get(), 1, rank, tag);
        const int rq = mpi__->request_recv(buffer.get(), size, rank, tag);
        auto m = requests_.insert(make_pair(rq, make_shared<Prep>(size, off, move(buffer))));
        assert(m.second);
        i = calls_.erase(i);
        ++cnt;
      } else {
        ++i;
      }
    }
  }

  for (int i = 0; i != cnt; ++i)
    init();

  {
    lock_guard<mutex> lock(mutex_);
    for (auto i = requests_.begin(); i != requests_.end(); ) {
      if (mpi__->test(i->first)) {
        shared_ptr<Prep> p = i->second;
        // for the time being, we only consider matrices for mutex
        assert(p->off % p->size == 0 && datamutex_->size() > p->off/p->size);
        lock_guard<mutex> lock((*datamutex_)[p->off/p->size]);
        // perform daxpy
        daxpy_(p->size, 1.0, p->buf.get(), 1, data_+p->off, 1);
        i = requests_.erase(i);
      } else {
        ++i;
      }
    }
  }
}


bool AccRequest::test() {
  lock_guard<mutex> lock(mutex_);
  return requests_.empty();
}
