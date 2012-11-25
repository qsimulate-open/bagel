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

#include <cassert>
#include <stdexcept>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

MPI_Interface::MPI_Interface(int argc, char** argv) : cnt_(0) {
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
#ifdef HAVE_SCALAPACK
  tie(nprow_, npcol_) = numgrid(mpi__->size());
  sl_init_(context_, nprow_, npcol_);
#endif
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


int MPI_Interface::rank() const {
#ifdef HAVE_MPI_H
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  return rank;
#else
  return 0;
#endif
}


int MPI_Interface::size() const {
#ifdef HAVE_MPI_H
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
#else
  return 1;
#endif
}


void MPI_Interface::barrier() const {
#ifdef HAVE_MPI_H
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void MPI_Interface::reduce(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI_Reduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI::SUM, root, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::allreduce(double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::broadcast(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI_Bcast(static_cast<void*>(a), size, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
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
#endif
} 


void MPI_Interface::allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  // I hate const_cast. Blame the MPI C binding
  MPI_Allgather(const_cast<void*>(static_cast<const void*>(send)), ssize, MPI_INT, static_cast<void*>(rec), rsize, MPI_INT, MPI_COMM_WORLD);
#endif
} 


int MPI_Interface::request_send(const double* sbuf, const size_t size, const int dest) {
#ifdef HAVE_MPI_H
  MPI_Request rq;
  // I hate const_cast. Blame the MPI C binding
  MPI_Isend(const_cast<double*>(sbuf), size, MPI_DOUBLE, dest, cnt_, MPI_COMM_WORLD, &rq); 
  request_.insert(make_pair(cnt_, rq));
#endif
  ++cnt_;
  return cnt_-1;
}


int MPI_Interface::request_recv(double* rbuf, const size_t size, const int origin) {
#ifdef HAVE_MPI_H
  MPI_Request rq;
  MPI_Irecv(rbuf, size, MPI_DOUBLE, origin, MPI_ANY_TAG, MPI_COMM_WORLD, &rq); 
  request_.insert(make_pair(cnt_, rq));
#endif
  ++cnt_;
  return cnt_-1;
}


void MPI_Interface::wait(const int rq) {
#ifdef HAVE_MPI_H
  auto i = request_.find(rq);
  assert(i != request_.end());
  MPI_Wait(&i->second, MPI_STATUS_IGNORE);
#endif
}
