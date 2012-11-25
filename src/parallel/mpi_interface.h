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

#ifndef __SRC_PARALLEL_MPI_INTERFACE_H
#define __SRC_PARALLEL_MPI_INTERFACE_H

#include <stddef.h>
#include <config.h>
#include <memory>
#ifdef HAVE_MPI_H
 #include <mpi.h>
#endif

namespace bagel {


class MPI_Interface {
  protected:
    int cnt_;
#ifdef HAVE_MPI_H
    // request handles
    std::map<int, MPI_Request> request_; 
#endif
    int nprow_; 
    int npcol_;
    int context_;
    int myprow_;
    int mypcol_;

  public:
    MPI_Interface(int argc, char** argv);
    ~MPI_Interface();

    int rank() const;
    int size() const;
    bool last() const { return rank() == size()-1; }

    // collective functions
    // barrier
    void barrier() const;
    // sum reduce to the root process
    void reduce(double*, const size_t size, const int root) const;
    // sum reduce and broadcast to each process
    void allreduce(double*, const size_t size) const;
    // broadcast
    void broadcast(double*, const size_t size, const int root) const;
    // broadcast of const objects. Use with caution...
    void broadcast_force(const double*, const size_t size, const int root) const;
    void allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const; 
    void allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const; 

    int request_send(const double* sbuf, const size_t size, const int dest);
    int request_recv(double* rbuf, const size_t size, const int source);
    void wait(const int rq);

    // scalapack
    int nprow() const { return nprow_; }
    int npcol() const { return npcol_; }
    int context() const { return context_; }
    int myprow() const { return myprow_; }
    int mypcol() const { return mypcol_; }

    std::pair<int,int> numroc(const int, const int) const;
    std::unique_ptr<int[]> descinit(const int, const int) const;
};

extern MPI_Interface* mpi__; 

}


#endif

