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

#include <src/parallel/accrequest.h>
#include <src/util/f77.h>
#include <boost/thread/thread.hpp>
#include <src/util/constants.h>

using namespace std;
using namespace bagel;


// SendRequest sends buffer using MPI_send. When completed, releases the buffer. Note that
SendRequest::SendRequest() : counter_(probe_key__/mpi__->size()*mpi__->rank()) {

}

void SendRequest::request_send(unique_ptr<double[]> buf, const size_t size, const int dest, const size_t off) {
  // sending size
  shared_ptr<Probe> p(new Probe(size, counter_, mpi__->rank(), dest, off, buf));
  ++counter_;
  const int srq = mpi__->request_send(p->size,    4, dest, probe_key__);
  const int rrq = mpi__->request_recv(&p->target, 1, dest, p->tag);
  auto m = inactive_.insert(make_pair(rrq, p));
  send_.push_back(srq);
  assert(m.second);
}

void SendRequest::flush() {

  // if receive buffer at the destination is created, send the message
  for (auto i = inactive_.begin(); i != inactive_.end(); ) {
    if (mpi__->test(i->first)) {
      shared_ptr<Probe> p = i->second;
      assert(round(p->target) == p->tag);
      const int srq = mpi__->request_send(p->buf.get(), p->size[0], p->targetrank, p->tag);
      auto m = requests_.insert(make_pair(srq, move(p->buf)));
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

// wait for all calls
bool SendRequest::wait1() {
  bool done = true;
  for (auto i = send_.begin(); i != send_.end(); ) {
    if (mpi__->test(*i)) {
      i = send_.erase(i);
    } else {
      ++i;
      done = false;
    }
  }
  flush();
  return done;
}


void SendRequest::wait2() {
  for (auto& i : inactive_)
    mpi__->wait(i.first);
  flush();
  inactive_.clear();
}


void SendRequest::wait3() {
  for (auto& i : requests_)
    mpi__->wait(i.first);
  flush();
  requests_.clear();
}


void AccRequest::init() {
  // we should compute the number of messnages to receive?
  // receives
  unique_ptr<size_t[]> buf(new size_t[4]);
  // receives size,tag,rank
  const int rq = mpi__->request_recv(buf.get(), 4, -1, probe_key__); 
  auto m = calls_.insert(make_pair(rq, move(buf)));
  assert(m.second);
}


AccRequest::AccRequest(double* const d, vector<boost::mutex>* m) : data_(d), mutex_(m) {
}


AccRequest::~AccRequest() {
  for (auto& i : calls_)
    mpi__->cancel(i.first);
}


void AccRequest::init_request() {
  for (size_t i = 0; i != pool_size__; ++i)
    init();
}


void AccRequest::flush() {
  size_t cnt = 0;
  for (auto i = calls_.begin(); i != calls_.end(); ) {
    // if this has already arrived, create a buffer, and return the address.
    if (mpi__->test(i->first)) {
      const int size = i->second[0];
      const int tag  = i->second[1];
      const int rank = i->second[2];
      const int off  = i->second[3];
      // allocating buffer
      unique_ptr<double[]> buffer(new double[size]);
      buffer[0] = tag;
      // sending back the buffer address
      const int sq = mpi__->request_send(buffer.get(), 1, rank, tag);
      const int rq = mpi__->request_recv(buffer.get(), size, rank, tag); 
      send_.push_back(sq);
      auto m = requests_.insert(make_pair(rq, shared_ptr<Prep>(new Prep(size, off, move(buffer)))));
      assert(m.second);
      i = calls_.erase(i);
      ++cnt;
    } else {
      ++i;
    }
  }
  for (int i = 0; i != cnt; ++i)
    init(); 

  for (auto i = requests_.begin(); i != requests_.end(); ) {
    if (mpi__->test(i->first)) {
      shared_ptr<Prep> p = i->second;
      // for the time being, we only consider matrices for mutex
      assert(p->off % p->size == 0 && mutex_->size() > p->off/p->size);
      boost::lock_guard<boost::mutex> lock((*mutex_)[p->off/p->size]);
      // perform daxpy
      daxpy_(p->size, 1.0, p->buf.get(), 1, data_+p->off, 1);
      i = requests_.erase(i); 
    } else {
      ++i;
    }
  }
}


void AccRequest::wait2() {
  for (auto& i : send_)
    mpi__->wait(i);
  flush();
  send_.clear();
}


void AccRequest::wait3() {
  for (auto& i : requests_)
    mpi__->wait(i.first);
  flush();
  requests_.clear();
}
