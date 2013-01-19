//
// BAGEL - Parallel electron correlation program.
// Filename: recvrequest.h
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

#ifndef __SRC_PARALLEL_RECVREQUEST_H
#define __SRC_PARALLEL_RECVREQUEST_H

#include <map>
#include <memory>
#include <mutex>
#include <tuple>
#include <cassert>
#include <src/util/serverflush.h>
#include <src/parallel/mpi_interface.h>
#include <src/parallel/resources.h>

namespace bagel {

// receives using MPI_irecv and accumulate to the local destination
class PutRequest : public ServerFlush {
  protected:
    std::map<int, std::unique_ptr<size_t[]> > calls_;

    void init();
    const double* const data_;

    // this mutex is for MPI calls
    std::mutex block_;

    void flush() override;

  public:
    PutRequest(const double* d);
    ~PutRequest();
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

    // tuple contains: size, if ready, target rank, and buffer 
    std::map<int, std::shared_ptr<Probe> > request_;
    std::vector<int> probe_;

    std::mutex block_;

  public:
    RecvRequest();
    // return mpi tag
    int request_recv(double* target, const size_t size, const int dest, const size_t off);
    bool test();

};

}

#endif

