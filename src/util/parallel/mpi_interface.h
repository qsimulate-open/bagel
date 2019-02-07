//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mpi_interface.h
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

#ifndef __SRC_PARALLEL_MPI_INTERFACE_H
#define __SRC_PARALLEL_MPI_INTERFACE_H

#include <stddef.h>
#include <bagel_config.h>
#include <memory>
#include <complex>
#include <mutex>
#include <vector>
#include <map>
#ifdef HAVE_MPI_H
 #include <mpi.h>
#endif

namespace bagel {

class MPI_Interface {
  protected:
    int world_rank_;
    int world_size_;
    int rank_;
    int size_;
    int depth_;

    int cnt_;
    // request handles
#ifdef HAVE_MPI_H
    MPI_Comm mpi_comm_;
    std::map<int, std::vector<MPI_Request>> request_;
    std::vector<std::pair<MPI_Comm,std::array<int,5>>> mpi_comm_old_;
#endif
#ifdef HAVE_SCALAPACK
    std::vector<int> pmap_;
#endif
    int nprow_;
    int npcol_;
    int context_;
    int myprow_;
    int mypcol_;

    // maximum size of the MPI buffer
    static constexpr size_t bsize = 100000000LU;

    // mutex for isend and irecv
    mutable std::mutex mpimutex_;

    // MPI's internal variables
    int tag_ub_;

  public:
    MPI_Interface();
    ~MPI_Interface();

    int world_rank() const { return world_rank_; }
    int world_size() const { return world_size_; }
    int rank() const { return rank_; }
    int size() const { return size_; }
    int depth() const { return depth_; }
    bool last() const { return rank() == size()-1; }

    // collective functions
    // barrier
    void barrier() const;
    // sum reduce and broadcast to each process
    void allreduce(int*, const size_t size) const;
    void allreduce(double*, const size_t size) const;
    void allreduce(std::complex<double>*, const size_t size) const;
    // broadcast
    void broadcast(size_t*, const size_t size, const int root) const;
    void broadcast(double*, const size_t size, const int root) const;
    void broadcast(std::complex<double>*, const size_t size, const int root) const;
    void allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const;
    void allgather(const std::complex<double>* send, const size_t ssize, std::complex<double>* rec, const size_t rsize) const;
    void allgather(const size_t* send, const size_t ssize, size_t* rec, const size_t rsize) const;
    void allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const;

    // one-sided communication with Isend, Irecv
    int request_send(const double* sbuf, const size_t size, const int dest, const int tag);
    int request_send(const std::complex<double>* sbuf, const size_t size, const int dest, const int tag);
    int request_send(const size_t* sbuf, const size_t size, const int dest, const int tag);
    int request_recv(double* rbuf, const size_t size, const int source = -1, const int tag = -1);
    int request_recv(std::complex<double>* rbuf, const size_t size, const int source = -1, const int tag = -1);
    int request_recv(size_t* rbuf, const size_t size, const int source = -1, const int tag = -1);
    void wait(const int rq);
    void cancel(const int rq);
    bool test(const int rq);

    // scalapack
    int nprow() const { return nprow_; }
    int npcol() const { return npcol_; }
    int context() const { return context_; }
    int myprow() const { return myprow_; }
    int mypcol() const { return mypcol_; }

#ifdef HAVE_MPI_H
    // communicators. n is the number of processes per communicator.
    const MPI_Comm& mpi_comm() const { return mpi_comm_; }
#endif
    void split(const int n);
    void merge();

    int pnum(const int prow, const int pcol) const;
    std::pair<int,int> numroc(const int, const int) const;
    std::vector<int> descinit(const int, const int) const;
};

extern MPI_Interface* mpi__;

}


#endif

