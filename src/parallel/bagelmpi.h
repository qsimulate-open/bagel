//
// BAGEL - Parallel electron correlation program.
// Filename: mpi_interface.h
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

#ifndef __SRC_UTIL_BAGELMPI_H
#define __SRC_UTIL_BAGELMPI_H

#include <stddef.h>
#include <config.h>
#include <memory>
#include <complex>
#ifdef HAVE_MPI_H
 #include <mpi.h>
#endif

namespace bagel {

// class for static members only
struct BagelMPI {
  public:
    // collective functions
    // barrier
    static void barrier() const;
    // barrier but with less mutex lock
    static void soft_barrier();
    // sum reduce to the root process
    static void reduce(double*, const size_t size, const int root) const;
    // sum reduce and broadcast to each process
    static void allreduce(int*, const size_t size) const;
    static void allreduce(double*, const size_t size) const;
    static void allreduce(std::complex<double>*, const size_t size) const;
    // all reduce but with less mutex lock
    static void soft_allreduce(size_t*, const size_t size);
    // broadcast
    static void broadcast(size_t*, const size_t size, const int root) const;
    static void broadcast(double*, const size_t size, const int root) const;
    static void broadcast(std::complex<double>*, const size_t size, const int root) const;
    // broadcast of const objects. Use with caution...
    static void broadcast_force(const double*, const size_t size, const int root) const;
    static void allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const;
    static void allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const;
    static void allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const;

    // non-blocking communication with Isend, Irecv
    static int request_send(const double* sbuf, const size_t size, const int dest, const int tag);
    static int request_send(const size_t* sbuf, const size_t size, const int dest, const int tag);
    static int request_recv(double* rbuf, const size_t size, const int source = -1, const int tag = -1);
    static int request_recv(size_t* rbuf, const size_t size, const int source = -1, const int tag = -1);
    static void wait(const int rq);
    static void cancel(const int rq);
    static bool test(const int rq);
};

}

#endif

