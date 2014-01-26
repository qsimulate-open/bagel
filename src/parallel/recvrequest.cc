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

#include <src/parallel/recvrequest.h>
#include <src/util/f77.h>
#include <src/util/constants.h>

using namespace std;
using namespace bagel;


// PutRequest receives probe and returns data
PutRequest::PutRequest(const double* d, const size_t probe_offset) : data_(d), probe_offset_(probe_offset) {
  for (size_t i = 0; i != pool_size__; ++i)
    init();
  turn_on();
}


void PutRequest::init() {
  // we should compute the number of messages to receive?
  // receives
  auto call = make_shared<Call>();
  // receives size,tag,rank
  const int rq = mpi__->request_recv(call->buf.get(), 4, -1, probe_key2__ + probe_offset_);
  {
    lock_guard<mutex> lock(block_);
    auto m = calls_.insert(make_pair(rq, call));
    if (!m.second)
      throw logic_error("PutRequest::init error");
  }
}


PutRequest::~PutRequest() {
  turn_off();
  for (auto& i : calls_)
    mpi__->cancel(i.first);
}


void PutRequest::flush_() {
  size_t cnt = 0;
  {
    lock_guard<mutex> lock(block_);
    for (auto i = calls_.begin(); i != calls_.end(); ) {
      // if this has already arrived, send data
      if (mpi__->test(i->first)) {
        const int size = i->second->buf[0];
        const int tag  = i->second->buf[1];
        const int rank = i->second->buf[2];
        const int off  = i->second->buf[3];
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

// BufferPutRequest receives probe and returns information on data requested
BufferPutRequest::BufferPutRequest(const size_t probe_offset) : probe_offset_(probe_offset) {
  for (size_t i = 0; i != pool_size__; ++i)
    init();
  turn_on();
}

void BufferPutRequest::init() {
  // we should compute the number of messages to receive?
  // receives
  auto call = make_shared<Call>();
  // receives size,tag,rank
  const int rq = mpi__->request_recv(call->buf.get(), 4, -1, probe_key2__ + probe_offset_);
  {
    lock_guard<mutex> lock(block_);
    auto m = calls_.emplace(rq, call);
    if (!m.second)
      throw logic_error("BufferPutRequest::init error");
  }
}


BufferPutRequest::~BufferPutRequest() {
  turn_off();
  for (auto& i : calls_)
    mpi__->cancel(i.first);
}


// cleans out finished sends
void BufferPutRequest::flush_() {
  lock_guard<mutex> lock(block_);
  for (auto i = buffs_.begin(); i != buffs_.end(); ) {
    if (mpi__->test(i->first)) {
      i = buffs_.erase(i);
    }
    else {
      ++i;
    }
  }
}

// takes ownership of buf and initiates a send
void BufferPutRequest::request_send(unique_ptr<double[]>&& buf, const size_t size, const size_t dest, const size_t tag) {
  auto buffer = make_shared<Buff>(move(buf));
  lock_guard<mutex> lock(block_);
  const int srq = mpi__->request_send(buffer->buf.get(), size, dest, tag);
  buffs_.emplace(srq, buffer);
}

// returns all of the requests for data
vector<array<size_t, 4>> BufferPutRequest::get_calls() {
  size_t cnt = 0;
  vector<array<size_t, 4>> out;
  {
    lock_guard<mutex> lock(block_);
    for (auto i = calls_.begin(); i != calls_.end(); ) {
      // if this has already arrived, send data
      if (mpi__->test(i->first)) {
        const size_t size = i->second->buf[0];
        const size_t tag  = i->second->buf[1];
        const size_t rank = i->second->buf[2];
        const size_t off  = i->second->buf[3];
        out.push_back( array<size_t, 4>{{size, tag, rank, off}} );
        i = calls_.erase(i);
        ++cnt;
      } else {
        ++i;
      }
    }
  }
  for (int i = 0; i != cnt; ++i)
    init();

  return out;
}


/////////////////////////

RecvRequest::RecvRequest(const size_t nprobes)
  : counter_(probe_key2__ + nprobes*(mpi__->rank() + 1)), nprobes_(nprobes)
{ }


int RecvRequest::request_recv(double* buf, const size_t size, const int dest, const size_t off, const size_t probe_offset) {
  // sending size
  lock_guard<mutex> lock(block_);
  auto p = make_shared<Probe>(size, counter_, mpi__->rank(), dest, off, buf);
  counter_ += mpi__->size()*nprobes_;
  const int srq = mpi__->request_send(p->size, 4, dest, probe_key2__ + probe_offset);
  probe_.push_back(srq);
  const int rrq = mpi__->request_recv(p->buf, p->size[0], dest, p->tag);
  auto m = request_.emplace(rrq, p);
  if (!m.second)
    throw logic_error("RecvRequest::request_recv error");
  return rrq;
}


bool RecvRequest::test() {
  bool done = true;
  lock_guard<mutex> lock(block_);
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
