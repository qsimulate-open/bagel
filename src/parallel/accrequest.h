//
// BAGEL - Parallel electron correlation program.
// Filename: accrequest.h
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

#ifndef __SRC_PARALLEL_ACCREQUEST_H
#define __SRC_PARALLEL_ACCREQUEST_H

#include <tuple>
#include <cassert>
#include <src/parallel/mpi_interface.h>
#include <src/parallel/resources.h>
#include <src/util/serverflush.h>

namespace bagel {


// SendRequest sends buffer using MPI_send. When completed, releases the buffer. Note that
// one needs to periodically call "void flush()"
class SendRequest : public ServerFlush {
  protected:
    struct Probe {
      const size_t size[4];
      const size_t tag;
      double target;
      const size_t myrank;
      const size_t targetrank;
      const size_t off;
      std::unique_ptr<double[]> buf;
      Probe(const size_t s, const size_t c, const size_t r, const size_t t, const size_t o, std::unique_ptr<double[]>&& b)
        : size{s,c,r,o}, tag(c), myrank(r), targetrank(t), off(o), buf(std::move(b)) { }
    };

    size_t counter_;

    // mutex for protected members of this class
    std::mutex mutex_;

    // tuple contains: size, if ready, target rank, and buffer
    std::map<int, std::shared_ptr<Probe>> inactive_;
    std::map<int, std::shared_ptr<Probe>> requests_;

    void flush_() override;

  public:
    SendRequest();
    ~SendRequest();

    void request_send(std::unique_ptr<double[]> buf, const size_t size, const int dest, const size_t off);

    bool test();

};


// receives using MPI_irecv and accumulate to the local destination
class AccRequest : public ServerFlush {
  protected:

    double* const data_;
    // mutex for the target data_ (for each alpha string)
    std::vector<std::mutex>* datamutex_;

    struct Call {
      std::unique_ptr<size_t[]> buf;
      Call() : buf(new size_t[4]) { }
    };

    struct Prep {
      const size_t size;
      const size_t off;
      std::unique_ptr<double[]> buf;
      Prep(const size_t s, const size_t o, std::unique_ptr<double[]> b) : size(s), off(o), buf(std::move(b)) {}
    };

    // mutex for protected members of this class
    std::mutex mutex_;

    // speculative calls to receive probes
    std::map<int, std::shared_ptr<Call>> calls_;
    // requests to receive data from send_
    std::map<int, std::shared_ptr<Prep>> requests_;

    // one speculative call
    void init();
    void flush_();

  public:
    AccRequest(double* const d, std::vector<std::mutex>* m);
    ~AccRequest();

    bool test();

};

}

#endif

