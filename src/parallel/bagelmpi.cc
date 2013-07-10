//
// BAGEL - Parallel electron correlation program.
// Filename: bagelmpi.cc
// Copyright (C) 2012 Toru Shiozaki
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

// TODO until GCC fixes this bug
#define _GLIBCXX_USE_NANOSLEEP

#include <cassert>
#include <thread>
#include <algorithm>
//#include <src/util/f77.h>
//#include <src/util/constants.h>
#include <src/parallel/bagelmpi.h>

using namespace std;
using namespace bagel;


void BagelMPI::barrier() const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}


// barrier without locking mutex all the time
void BagelMPI::soft_barrier() {
  vector<int> receive; receive.reserve(size_);
  vector<size_t> msg(size_);
  for (int i = 0; i != size_; ++i) {
    if (i == rank_) continue;
    request_send(&msg[rank_], 1, i, (1<<30));
    receive.push_back(request_recv(&msg[i], 1, i, (1<<30)));
  }
  bool done;
  do {
    done = true;
    for (auto& i : receive)
      if (!test(i)) { done = false; break; }
    if (!done) this_thread::sleep_for(sleeptime__);
  } while (!done);
}


void BagelMPI::reduce(double* a, const size_t size, const int root) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  MPI_Reduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI::SUM, root, MPI_COMM_WORLD);
#endif
}


void BagelMPI::allreduce(double* a, const size_t size) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}


void BagelMPI::allreduce(int* a, const size_t size) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
}


void BagelMPI::allreduce(complex<double>* a, const size_t size) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
}


void BagelMPI::soft_allreduce(size_t* a, const size_t size) {
  vector<size_t> receive;
  vector<size_t> msg(size_*size);

  for (int i = 0; i != size_; ++i) {
    if (i == rank_) continue;
    request_send(a, size, i, (1<<30)+1);
    receive.push_back(request_recv(&msg[i*size], size, i, (1<<30)+1));
  }
  bool done;
  do {
    done = true;
    for (auto& i : receive)
      if (!test(i)) { done = false; break; }
    if (!done) this_thread::sleep_for(sleeptime__);
  } while (!done);

  for (int i = 0; i != size_; ++i)
    if (i != rank_) {
      for (int j = 0; j != size; ++j) a[j] += msg[i*size+j];
    }
}


void BagelMPI::broadcast(size_t* a, const size_t size, const int root) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(long long), "size_t is assumed to be the same size as long long");
  MPI_Bcast(static_cast<void*>(a), size, MPI_LONG_LONG, root, MPI_COMM_WORLD);
#endif
}


void BagelMPI::broadcast(double* a, const size_t size, const int root) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  MPI_Bcast(static_cast<void*>(a), size, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
}


void BagelMPI::broadcast(complex<double>* a, const size_t size, const int root) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  MPI_Bcast(static_cast<void*>(a), size, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD);
#endif
}


void BagelMPI::broadcast_force(const double* a, const size_t size, const int root) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  // sometimes we need to broadcast const objects for consistency...
  double* aa = const_cast<double*>(a);
  MPI_Bcast(static_cast<void*>(aa), size, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
}


void BagelMPI::allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_DOUBLE, static_cast<void*>(rec), rsize, MPI_DOUBLE, MPI_COMM_WORLD);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void BagelMPI::allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(long long), "size_t is assumed to be the same size as long long");
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_LONG_LONG, static_cast<void*>(rec), rsize, MPI_LONG_LONG, MPI_COMM_WORLD);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void BagelMPI::allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_INT, static_cast<void*>(rec), rsize, MPI_INT, MPI_COMM_WORLD);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


int BagelMPI::request_send(const double* sbuf, const size_t size, const int dest, const int tag) {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<double*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &c);
    rq.push_back(c);
  }
  request_.insert(make_pair(cnt_, rq));
#endif
  ++cnt_;
  return cnt_-1;
}


int BagelMPI::request_send(const size_t* sbuf, const size_t size, const int dest, const int tag) {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(long long), "size_t is assumed to be the same size as long long");
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<size_t*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_LONG_LONG, dest, tag, MPI_COMM_WORLD, &c);
    rq.push_back(c);
  }
  request_.insert(make_pair(cnt_, rq));
#endif
  ++cnt_;
  return cnt_-1;
}



int BagelMPI::request_recv(double* rbuf, const size_t size, const int origin, const int tag) {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), MPI_COMM_WORLD, &c);
    rq.push_back(c);
  }
  request_.insert(make_pair(cnt_, rq));
#endif
  ++cnt_;
  return cnt_-1;
}


int BagelMPI::request_recv(size_t* rbuf, const size_t size, const int origin, const int tag) {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_LONG_LONG, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), MPI_COMM_WORLD, &c);
    rq.push_back(c);
  }
  request_.insert(make_pair(cnt_, rq));
#endif
  ++cnt_;
  return cnt_-1;
}


void BagelMPI::wait(const int rq) {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  auto i = request_.find(rq);
  assert(i != request_.end());
  for (auto& j : i->second)
    MPI_Wait(&j, MPI_STATUS_IGNORE);
#endif
}


void BagelMPI::cancel(const int rq) {
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  auto i = request_.find(rq);
  assert(i != request_.end());
  for (auto& j : i->second)
    MPI_Cancel(&j);
#endif
}


bool BagelMPI::test(const int rq) {
  lock_guard<mutex> lock(mpimutex_);
  bool out = true;
#ifdef HAVE_MPI_H
  auto i = request_.find(rq);
  assert(i != request_.end());
  for (auto& j : i->second) {
    int b;
    MPI_Test(&j, &b, MPI_STATUS_IGNORE);
    out &= b;
  }
#endif
  return out;
}
