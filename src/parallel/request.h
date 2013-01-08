//
// BAGEL - Parallel electron correlation program.
// Filename: request.h
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

#ifndef __SRC_PARALLEL_REQUEST_H
#define __SRC_PARALLEL_REQUEST_H

#include <map>
#include <memory>
#include <tuple>
#include <cassert>
#include <src/parallel/mpi_interface.h>

namespace bagel {

// SendRequest sends buffer using MPI_send. When completed, releases the buffer. Note that
// one needs to periodically call "void request_test()" 
class SendRequest {
  protected:
    struct Probe {
      const size_t size;
      const size_t tag;
      size_t ready; // 1 ready, 0 not
      const int rank;
      std::unique_ptr<double[]> buf;
      Probe(const size_t s, const size_t c, const size_t r, const int ra, std::unique_ptr<double[]>& b) : size(s), tag(c), ready(r), rank(ra), buf(std::move(b)) { }
    };

    size_t counter_;

    // tuple contains: size, if ready, target rank, and buffer 
    std::map<int, std::shared_ptr<Probe> > inactive_;
    std::map<int, std::unique_ptr<double[]> > requests_;

  public:
    SendRequest() : counter_(0) {}

    void request_send(std::unique_ptr<double[]> buf, const size_t size, const int dest) {
      // sending size
      std::shared_ptr<Probe> p(new Probe(size, counter_, 0, dest, buf));
      ++counter_;
      const int srq = mpi__->request_send(&p->size,  1, dest, p->tag);
      const int rrq = mpi__->request_recv(&p->ready, 1, dest, p->tag);
      auto m = inactive_.insert(std::make_pair(rrq, p));
      assert(m.second);
    }

    void flush() {

      // if receive buffer at the destination is created, send the message
      for (auto i = inactive_.begin(); i != inactive_.end(); ) {
        if (mpi__->test(i->first)) {
          int j = mpi__->request_send(i->second->buf.get(), i->second->size, i->second->rank);
          auto m = requests_.insert(std::make_pair(j, std::move(i->second->buf)));
          assert(m.second);
          i = inactive_.erase(i); 
        } else {
          ++i;
        }
      }

      // if the message is sent, delete the request
      for (auto i = requests_.begin(); i != requests_.end(); ) {
        if (mpi__->test(i->first))
          i = requests_.erase(i);
        else
          ++i;
      }
    }

};


// receives using MPI_irecv and accumulate to the local destination
class AccRequest {
  protected:


  public:
    AccRequest() {}

    void flush() {
      // receives 
    }

};

}

#endif

