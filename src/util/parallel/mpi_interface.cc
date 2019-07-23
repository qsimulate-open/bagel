//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mpi_interface.cc
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

#include <iostream>
#include <iomanip>
#include <cassert>
#include <thread>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <src/util/f77.h>
#include <src/util/constants.h>
#include <src/util/parallel/scalapack.h>
#include <src/util/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

MPI_Interface::MPI_Interface()
 : depth_(0), cnt_(0), nprow_(0), npcol_(0), context_(0), myprow_(0), mypcol_(0), mpimutex_() {

#ifdef HAVE_MPI_H
  int provided;
  MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &provided);
  if (provided != MPI_THREAD_MULTIPLE)
    throw runtime_error("MPI_THREAD_MULTIPLE not provided");

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
  rank_ = world_rank_;
  size_ = world_size_;
#ifdef HAVE_SCALAPACK
  tie(nprow_, npcol_) = numgrid(size());
  if (rank() == 0)
    cout << "  * process grid (" << nprow_ << ", " << npcol_ << ") will be used" << endl;
  sl_init_(context_, nprow_, npcol_);
  blacs_gridinfo_(context_, nprow_, npcol_, myprow_, mypcol_);
#endif

  // print out the node name
  {
    constexpr const size_t maxlen = MPI_MAX_PROCESSOR_NAME;
    int len;
    char name[maxlen];
    MPI_Get_processor_name(name, &len);

    unique_ptr<char[]> buf(new char[maxlen*size_]);
    unique_ptr<int[]> lens(new int[size_]);
    MPI_Gather(static_cast<void*>(name), maxlen, MPI_CHAR, static_cast<void*>(buf.get()), maxlen, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Gather(static_cast<void*>(&len),      1, MPI_INT,  static_cast<void*>(lens.get()),     1, MPI_INT,  0, MPI_COMM_WORLD);
    if (rank() == 0) {
      for (int i = 0; i != size_; ++i)
        cout << left << "    " << setw(32) << string(&buf[i*maxlen], &buf[i*maxlen+lens[i]]) << right << endl;
      cout << endl;
    }
  }

  // obtain the upper bound of tags
  {
    int flag, *get_val;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &get_val, &flag);
    assert(flag && *get_val >= 32767); // this is what the standard says
    tag_ub_ = *get_val;
  }

  // set MPI_COMM_WORLD to mpi_comm_
  mpi_comm_ = MPI_COMM_WORLD;
#else
  world_rank_ = 0;
  world_size_ = 1;
  rank_ = world_rank_;
  size_ = world_size_;
#endif
}


MPI_Interface::~MPI_Interface() {
#ifdef HAVE_MPI_H
#ifndef HAVE_SCALAPACK
  MPI_Finalize();
#else
  blacs_gridexit_(context_);
  blacs_exit_(0);
#endif
#endif
}


void MPI_Interface::barrier() const {
#ifdef HAVE_MPI_H
  MPI_Barrier(mpi_comm_);
#endif
}


void MPI_Interface::allreduce(double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, MPI_SUM, mpi_comm_);
#endif
}


void MPI_Interface::allreduce(int* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_INT, MPI_SUM, mpi_comm_);
#endif
}


void MPI_Interface::allreduce(complex<double>* a, const size_t size) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_);
#endif
}


void MPI_Interface::broadcast(size_t* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(unsigned long long), "size_t is assumed to be the same size as unsigned long long");
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_UNSIGNED_LONG_LONG, root, mpi_comm_);
#endif
}


void MPI_Interface::broadcast(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, root, mpi_comm_);
#endif
}


void MPI_Interface::broadcast(complex<double>* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  assert(size != 0);
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i)
    MPI_Bcast(static_cast<void*>(a+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, root, mpi_comm_);
#endif
}


void MPI_Interface::allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_DOUBLE, static_cast<void*>(rec), rsize, MPI_DOUBLE, mpi_comm_);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void MPI_Interface::allgather(const complex<double>* send, const size_t ssize, complex<double>* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_CXX_DOUBLE_COMPLEX, static_cast<void*>(rec), rsize, MPI_CXX_DOUBLE_COMPLEX, mpi_comm_);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void MPI_Interface::allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  static_assert(sizeof(size_t) == sizeof(unsigned long long), "size_t is assumed to be the same size as unsigned long long");
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_UNSIGNED_LONG_LONG, static_cast<void*>(rec), rsize, MPI_UNSIGNED_LONG_LONG, mpi_comm_);
#else
  assert(ssize == rsize);
  copy_n(send, ssize, rec);
#endif
}


void MPI_Interface::allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_INT, static_cast<void*>(rec), rsize, MPI_INT, mpi_comm_);
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
    MPI_Isend(const_cast<double*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, dest, tag, mpi_comm_, &c);
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


int MPI_Interface::request_send(const complex<double>* sbuf, const size_t size, const int dest, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<complex<double>*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, dest, tag, mpi_comm_, &c);
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
  static_assert(sizeof(size_t) == sizeof(unsigned long long), "size_t is assumed to be the same size as unsigned long long");
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    // I hate const_cast. Blame the MPI C binding
    MPI_Isend(const_cast<size_t*>(sbuf+i*bsize), (i+1 == nbatch ? size-i*bsize : bsize), MPI_UNSIGNED_LONG_LONG, dest, tag, mpi_comm_, &c);
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
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_DOUBLE, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), mpi_comm_, &c);
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


int MPI_Interface::request_recv(complex<double>* rbuf, const size_t size, const int origin, const int tag) {
#ifdef HAVE_MPI_H
  assert(tag <= tag_ub_);
  vector<MPI_Request> rq;
  const int nbatch = (size-1)/bsize  + 1;
  for (int i = 0; i != nbatch; ++i) {
    MPI_Request c;
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_CXX_DOUBLE_COMPLEX, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), mpi_comm_, &c);
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
    MPI_Irecv(rbuf+i*bsize, (i+1 == nbatch ? size-i*bsize : bsize), MPI_UNSIGNED_LONG_LONG, (origin == -1 ? MPI_ANY_SOURCE : origin), (tag==-1 ? MPI_ANY_TAG : tag), mpi_comm_, &c);
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


// MPI Communicators
void MPI_Interface::split(const int n) {
#ifdef HAVE_MPI_H
#ifdef HAVE_SCALAPACK
  // make a map between MPI ranks and process numbers
  if (depth_ == 0) {
    pmap_ = vector<int>(world_size_, 0);
    pmap_[world_rank_] = pnum(myprow_, mypcol_);
    allreduce(pmap_.data(), size_);
  }
#endif
  MPI_Comm new_comm;
  const int icomm = rank_ % n;
  mpi_comm_old_.push_back(pair<MPI_Comm,array<int,5>>(mpi_comm_, {context_, nprow_, npcol_, myprow_, mypcol_}));

  ++depth_;

  MPI_Comm_split(mpi_comm_, icomm, world_rank_, &new_comm);
  mpi_comm_ = new_comm;
  MPI_Comm_rank(mpi_comm_, &rank_);
  MPI_Comm_size(mpi_comm_, &size_);
#ifdef HAVE_SCALAPACK
  tie(nprow_, npcol_) = numgrid(size_);
  vector<int> imap(size_, 0);
  imap[rank_] = pmap_[world_rank_];
  allreduce(imap.data(), size_);
  blacs_get_(0, 0, context_);
  blacs_gridmap_(context_, imap.data(), nprow_, nprow_, npcol_);
  blacs_gridinfo_(context_, nprow_, npcol_, myprow_, mypcol_);
#endif
#endif
}


void MPI_Interface::merge() {
#ifdef HAVE_MPI_H
  MPI_Comm_free(&mpi_comm_);

  --depth_;

  mpi_comm_ = get<0>(mpi_comm_old_[depth_]);
#ifdef HAVE_SCALAPACK
  const array<int,5> scalapack_info = get<1>(mpi_comm_old_[depth_]);
#endif
  MPI_Comm_rank(mpi_comm_, &rank_);
  MPI_Comm_size(mpi_comm_, &size_);

  mpi_comm_old_.pop_back();
#ifdef HAVE_SCALAPACK
  blacs_gridexit_(context_);
  context_ = scalapack_info[0];
  nprow_ = scalapack_info[1];
  npcol_ = scalapack_info[2];
  myprow_ = scalapack_info[3];
  mypcol_ = scalapack_info[4];
#endif
#endif
}


// ScaLapack interfaces

pair<int,int> MPI_Interface::numroc(const int ndim, const int ncol) const {
#ifdef HAVE_SCALAPACK
  return {numroc_(ndim, blocksize__, myprow_, 0, nprow_), numroc_(ncol, blocksize__, mypcol_, 0, npcol_)};
#else
  return {0,0};
#endif
}

vector<int> MPI_Interface::descinit(const int ndim, const int ncol) const {
  vector<int> desc(9);
#ifdef HAVE_SCALAPACK
  const int localrow = numroc_(ndim, blocksize__, myprow_, 0, nprow_);
  int info;
  descinit_(desc.data(), ndim, ncol, blocksize__, blocksize__, 0, 0, context_, max(1,localrow), info);
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
