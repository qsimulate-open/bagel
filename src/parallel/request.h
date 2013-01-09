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
#include <src/parallel/resources.h>

namespace bagel {


const static size_t probe_key = (1 << 30);

// SendRequest sends buffer using MPI_send. When completed, releases the buffer. Note that
// one needs to periodically call "void request_test()" 
class SendRequest {
  protected:
    struct Probe {
      const size_t size[4];
      const size_t tag;
      double target;
      const size_t myrank;
      const size_t targetrank;
      const size_t off;
      std::unique_ptr<double[]> buf;
      Probe(const size_t s, const size_t c, const size_t r, const size_t t, const size_t o, std::unique_ptr<double[]>& b)
        : size{s,c,r,o}, tag(c), myrank(r), targetrank(t), off(o), buf(std::move(b)) { }
    };

    size_t counter_;

    // tuple contains: size, if ready, target rank, and buffer 
    std::map<int, std::shared_ptr<Probe> > inactive_;
    std::map<int, std::unique_ptr<double[]> > requests_;
    std::vector<int> send_;

  public:
    SendRequest() : counter_(probe_key/mpi__->size()*mpi__->rank()) {}

    void request_send(std::unique_ptr<double[]> buf, const size_t size, const int dest, const size_t off) {
      // sending size
      std::shared_ptr<Probe> p(new Probe(size, counter_, mpi__->rank(), dest, off, buf));
      ++counter_;
      const int srq = mpi__->request_send(p->size,    4, dest, probe_key);
      const int rrq = mpi__->request_recv(&p->target, 1, dest, p->tag);
      auto m = inactive_.insert(std::make_pair(rrq, p));
      send_.push_back(srq);
      assert(m.second);
    }

    void flush() {

      // if receive buffer at the destination is created, send the message
      for (auto i = inactive_.begin(); i != inactive_.end(); ) {
        if (mpi__->test(i->first)) {
          std::shared_ptr<Probe> p = i->second;
          assert(std::round(p->target) == p->tag);
          const int srq = mpi__->request_send(p->buf.get(), p->size[0], p->targetrank, p->tag);
          auto m = requests_.insert(std::make_pair(srq, std::move(p->buf)));
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

    // wait for all calls
    void wait1() {
      for (auto& i : send_)     mpi__->wait(i);
      flush();
    }
    void wait2() {
      for (auto& i : inactive_) mpi__->wait(i.first);
      flush();
    }
    void wait3() {
      for (auto& i : requests_) mpi__->wait(i.first);
      flush();
    }

};


// receives using MPI_irecv and accumulate to the local destination
class AccRequest {
  protected:

    size_t total_;
    double* const data_;

    struct Prep {
      const size_t size;
      const size_t off;
      std::unique_ptr<double[]> buf;
      Prep(const size_t s, const size_t o, std::unique_ptr<double[]> b) : size(s), off(o), buf(std::move(b)) {}
    };

    std::map<int, std::unique_ptr<size_t[]> > calls_;
    // buffer pointer, size, offset at the target, and buffer
    std::map<int, std::shared_ptr<Prep> > requests_;

    std::vector<int> send_;

    void init() {
      // we should compute the number of messnages to receive?
      // receives
      const int dest = 0;
      std::unique_ptr<size_t[]> buf(new size_t[4]);
      // receives size,tag,rank
      const int rq = mpi__->request_recv(buf.get(), 4, -1, probe_key); 
      auto m = calls_.insert(std::make_pair(rq, std::move(buf)));
      assert(m.second);
    }

  public:
    AccRequest(const int total, double* const d) : total_(total), data_(d) {
      for (size_t i = 0; i != total_; ++i) init();
    }

    void flush() {
      for (auto i = calls_.begin(); i != calls_.end(); ) {
        // if this has already arrived, create a buffer, and return the address.
        if (mpi__->test(i->first)) {
          const int size = i->second[0];
          const int tag  = i->second[1];
          const int rank = i->second[2];
          const int off  = i->second[3];
          // allocating buffer
          std::unique_ptr<double[]> buffer(new double[size]);
          buffer[0] = tag;
          // sending back the buffer address
          const int sq = mpi__->request_send(buffer.get(), 1, rank, tag);
          const int rq = mpi__->request_recv(buffer.get(), size, rank, tag); 
          send_.push_back(sq);
          auto m = requests_.insert(std::make_pair(rq, std::shared_ptr<Prep>(new Prep(size, off, std::move(buffer)))));
          assert(m.second);
          i = calls_.erase(i);
        } else {
          ++i;
        }
      }

      for (auto i = requests_.begin(); i != requests_.end(); ) {
        if (mpi__->test(i->first)) {
          std::shared_ptr<Prep> p = i->second;
          daxpy_(p->size, 1.0, p->buf.get(), 1, data_+p->off*p->size, 1);
          i = requests_.erase(i); 
        } else {
          ++i;
        }
      }
    }

    // wait for all calls
    void wait1() {
      for (auto& i : calls_) mpi__->wait(i.first);
      flush();
    }
    void wait2() {
      for (auto& i : send_)  mpi__->wait(i);
      flush();
    }
    void wait3() {
      for (auto& i : requests_) mpi__->wait(i.first);
      flush();
    }


};

}

#endif

