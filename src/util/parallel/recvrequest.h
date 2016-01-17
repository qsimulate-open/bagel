//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: recvrequest.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_PARALLEL_RECVREQUEST_H
#define __SRC_PARALLEL_RECVREQUEST_H

#include <tuple>
#include <cassert>
#include <src/util/serverflush.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/parallel/resources.h>

namespace bagel {

// receives using MPI_irecv and accumulate to the local destination
class PutRequest : public ServerFlush {
  protected:
    struct Call {
      std::unique_ptr<size_t[]> buf;
      Call() : buf(new size_t[4]) { }
    };
    std::map<int, std::shared_ptr<Call>> calls_;

    void init();
    const double* const data_;

    const size_t probe_offset_;

    // this mutex is for MPI calls
    std::mutex block_;

    void flush_() override;

  public:
    PutRequest(const double* d, const size_t probe_offset = 0);
    ~PutRequest();
};

// receives using MPI_irecv and sends while taking ownership of the sent data
class BufferPutRequest : public ServerFlush {
  protected:
    struct Call {
      std::unique_ptr<size_t[]> buf;
      Call() : buf(new size_t[4]) { }
    };
    std::map<int, std::shared_ptr<Call>> calls_;

    struct Buff {
      std::unique_ptr<double[]> buf;
      Buff(std::unique_ptr<double[]>&& b) : buf(std::move(b)) {}
    };
    std::map<int, std::shared_ptr<Buff>> buffs_;

    void init();

    const size_t probe_offset_;

    // this mutex is for MPI calls
    std::mutex block_;

    // not to be used
    void flush_() override;

  public:
    BufferPutRequest(const size_t probe_offset = 0);
    ~BufferPutRequest();

    void request_send(std::unique_ptr<double[]>&& buf, const size_t size, const size_t dest, const size_t tag);
    std::vector<std::array<size_t, 4>> get_calls();
};

class RecvRequest {
  protected:
    struct Probe {
      const size_t size[4];
      const size_t tag;
      double target;
      const size_t myrank;
      const size_t targetrank;
      const size_t off;
      // buf is the target area (which is local)
      double* buf;
      Probe(const size_t s, const size_t c, const size_t r, const size_t t, const size_t o, double* b)
        : size{s,c,r,o}, tag(c), myrank(r), targetrank(t), off(o), buf(b) { }
    };

    size_t counter_;
    const size_t nprobes_;

    // tuple contains: size, if ready, target rank, and buffer
    std::map<int, std::shared_ptr<Probe>> request_;
    std::vector<int> probe_;

    std::mutex block_;

  public:
    RecvRequest(const size_t nprobes = 1);
    // return mpi tag
    int request_recv(double* target, const size_t size, const int dest, const size_t off, const size_t probe_offset = 0);
    bool test();

};

}

#endif
