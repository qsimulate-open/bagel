//
// BAGEL - Parallel electron correlation program.
// Filename: mpi_interface.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <iostream>
#include <cassert>
#include <thread>
#include <stdexcept>
#include <algorithm>
#include <src/util/f77.h>
#include <src/util/constants.h>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

MPI_Interface::MPI_Interface()
 : cnt_(0), nprow_(0), npcol_(0), context_(0), myprow_(0), mypcol_(0), mpimutex_() {

#ifdef HAVE_MPI_H
  int provided;
  MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
  if (provided != MPI_THREAD_MULTIPLE)
    throw runtime_error("MPI_THREAD_MULTIPLE not provided");

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);
#ifdef HAVE_SCALAPACK
  tie(nprow_, npcol_) = numgrid(size());
  if (rank() == 0)
    cout << "  * process grid (" << nprow_ << ", " << npcol_ << ") will be used" << endl;
  sl_init_(context_, nprow_, npcol_);
  blacs_gridinfo_(context_, nprow_, npcol_, myprow_, mypcol_);
#endif
  {
    int flag, *get_val;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &get_val, &flag);
    assert(flag && *get_val >= 32767); // this is what the standard says
    tag_ub_ = *get_val;
  }
#else
  rank_ = 0;
  size_ = 1;
#endif
}


MPI_Interface::~MPI_Interface() {
#ifdef HAVE_MPI_H
#ifndef SCALAPACK
  MPI_Finalize();
#else
  blacs_gridexit_(context_);
  blacs_exit_(0);
#endif
#endif
}


void MPI_Interface::barrier() const {
#ifdef HAVE_MPI_H
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}


// barrier without locking mutex all the time
void MPI_Interface::soft_barrier() {
  vector<int> receive; receive.reserve(size_);
  vector<size_t> msg(size_);
  for (int i = 0; i != size_; ++i) {
    if (i == rank_) continue;
    // using the biggest tag value
    request_send(&msg[rank_], 1, i, tag_ub_);
    receive.push_back(request_recv(&msg[i], 1, i, tag_ub_));
  }
  bool done;
  do {
    done = true;
    for (auto& i : receive)
      if (!test(i)) { done = false; break; }
    if (!done) this_thread::sleep_for(sleeptime__);
  } while (!done);
}


void MPI_Interface::reduce(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI_Reduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI::SUM, root, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::reduce_scatter(double* sendbuf, double* recvbuf, int* recvcnts) const {
#ifdef HAVE_MPI_H
  MPI_Reduce_scatter(sendbuf, recvbuf, recvcnts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}

int MPI_Interface::ireduce_scatter(double* sendbuf, double* recvbuf, int* recvcnts) {
#ifdef HAVE_MPI_H
  vector<MPI_Request> rq;
  MPI_Request c;
  // I hate const_cast. Blame the MPI C binding
  MPI_Ireduce_scatter(sendbuf, recvbuf, recvcnts, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &c);
  rq.push_back(c);
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


void MPI_Interface::allreduce(double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::allreduce(int* a, const size_t size) const {
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::allreduce(complex<double>* a, const size_t size) const {
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
}


int MPI_Interface::iallreduce(size_t* a, const size_t size) {
  static_assert(sizeof(size_t) == sizeof(long long), "size_t is assumed to be the same size as long long");
#ifdef HAVE_MPI_H
  vector<MPI_Request> rq;
  MPI_Request c;
  MPI_Iallreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD, &c);
  rq.push_back(c);
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


void MPI_Interface::soft_allreduce(size_t* a, const size_t size) {
  vector<size_t> receive;
  vector<size_t> msg(size_*size);
  for (int i = 0; i != size_; ++i) {
    if (i == rank_) continue;
    request_send(a, size, i, tag_ub_-1);
    receive.push_back(request_recv(&msg[i*size], size, i, tag_ub_-1));
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


void MPI_Interface::broadcast(size_t* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(long long), "size_t is assumed to be the same size as long long");
  MPI_Bcast(static_cast<void*>(a), size, MPI_LONG_LONG, root, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::broadcast(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI_Bcast(static_cast<void*>(a), size, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::broadcast(complex<double>* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI_Bcast(static_cast<void*>(a), size, MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD);
#endif
}


int MPI_Interface::ibroadcast(double* a, const size_t size, const int root) {
#ifdef HAVE_MPI_H
  vector<MPI_Request> rq;
  MPI_Request c;
  MPI_Ibcast(static_cast<void*>(a), size, MPI_DOUBLE, root, MPI_COMM_WORLD, &c);
  rq.push_back(c);
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


void MPI_Interface::broadcast_force(const double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  // sometimes we need to broadcast const objects for consistency...
  double* aa = const_cast<double*>(a);
  MPI_Bcast(static_cast<void*>(aa), size, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_DOUBLE, static_cast<void*>(rec), rsize, MPI_DOUBLE, MPI_COMM_WORLD);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void MPI_Interface::allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(long long), "size_t is assumed to be the same size as long long");
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_LONG_LONG, static_cast<void*>(rec), rsize, MPI_LONG_LONG, MPI_COMM_WORLD);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void MPI_Interface::allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_INT, static_cast<void*>(rec), rsize, MPI_INT, MPI_COMM_WORLD);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


int MPI_Interface::request_send(const double* sbuf, const size_t size, const int dest, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<double*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


int MPI_Interface::request_send(const size_t* sbuf, const size_t size, const int dest, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  static_assert(sizeof(size_t) == sizeof(long long), "size_t is assumed to be the same size as long long");
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<size_t*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_LONG_LONG, dest, tag, MPI_COMM_WORLD, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}



int MPI_Interface::request_recv(double* rbuf, const size_t size, const int origin, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), MPI_COMM_WORLD, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


int MPI_Interface::request_recv(size_t* rbuf, const size_t size, const int origin, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_LONG_LONG, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), MPI_COMM_WORLD, &c);
    rq.push_back(c);
  }
#endif
  lock_guard<mutex> lock(mpimutex_);
#ifdef HAVE_MPI_H
  request_.emplace(cnt_, rq);
#endif
  ++cnt_;
  return cnt_-1;
}


void MPI_Interface::wait(const int rq) {
#ifdef HAVE_MPI_H
  lock_guard<mutex> lock(mpimutex_);
  auto i = request_.find(rq);
  assert(i != request_.end());
  for (auto& j : i->second)
    MPI_Wait(&j, MPI_STATUS_IGNORE);
#endif
}


void MPI_Interface::cancel(const int rq) {
#ifdef HAVE_MPI_H
  lock_guard<mutex> lock(mpimutex_);
  auto i = request_.find(rq);
  assert(i != request_.end());
  for (auto& j : i->second)
    MPI_Cancel(&j);
#endif
}


bool MPI_Interface::test(const int rq) {
  bool out = true;
#ifdef HAVE_MPI_H
  lock_guard<mutex> lock(mpimutex_);
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


// ScaLapack interfaces

pair<int,int> MPI_Interface::numroc(const int ndim, const int ncol) const {
#ifdef HAVE_SCALAPACK
  return make_pair(numroc_(ndim, blocksize__, myprow_, 0, nprow_),
                   numroc_(ncol, blocksize__, mypcol_, 0, npcol_));
#else
  return make_pair(0,0);
#endif
}

unique_ptr<int[]> MPI_Interface::descinit(const int ndim, const int ncol) const {
  unique_ptr<int[]> desc(new int[9]);
#ifdef HAVE_SCALAPACK
  const int localrow = numroc_(ndim, blocksize__, myprow_, 0, nprow_);
  int info;
  descinit_(desc.get(), ndim, ncol, blocksize__, blocksize__, 0, 0, context_, max(1,localrow), info);
#endif
  return desc;
}

int MPI_Interface::pnum(const int prow, const int pcol) const {
#ifdef HAVE_SCALAPACK
  return blacs_pnum_(context_, prow, pcol);
#else
  return 0;
#endif
}
