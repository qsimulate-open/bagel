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
      const size_t size[4];
      const size_t tag;
      size_t target; // 1 ready, 0 not
      const size_t rank;
      const size_t off;
      std::unique_ptr<double[]> buf;
      Probe(const size_t s, const size_t c, const size_t t, const size_t ra, const size_t o, std::unique_ptr<double[]>& b)
        : size{s,c,ra,o}, tag(c), target(t), rank(ra), off(o), buf(std::move(b)) { }
    };

    size_t counter_;

    // tuple contains: size, if ready, target rank, and buffer 
    std::map<int, std::shared_ptr<Probe> > inactive_;
    std::map<int, std::unique_ptr<double[]> > requests_;

  public:
    SendRequest() : counter_(0) {}

    void request_send(std::unique_ptr<double[]> buf, const size_t size, const int dest, const size_t off) {
      // sending size
      std::shared_ptr<Probe> p(new Probe(size, counter_, 0, dest, off, buf));
      ++counter_;
      const int srq = mpi__->request_send(p->size,    4, dest, p->tag);
      const int rrq = mpi__->request_recv(&p->target, 1, dest, p->tag);
      auto m = inactive_.insert(std::make_pair(rrq, p));
      assert(m.second);
    }

    void flush() {

      // if receive buffer at the destination is created, send the message
      for (auto i = inactive_.begin(); i != inactive_.end(); ) {
        if (mpi__->test(i->first)) {
          int j = mpi__->request_send(i->second->buf.get(), i->second->size[0], i->second->rank);
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
    double* const data_;

    struct Double_ptr {
      double* ptr;
      Double_ptr(double* p) : ptr(p) { }
    };

    struct Prep {
      const std::shared_ptr<Double_ptr> ptr;
      const size_t size;
      const size_t off;
      std::unique_ptr<double[]> buf;
      Prep(std::shared_ptr<Double_ptr> p, const size_t s, const size_t o, std::unique_ptr<double[]> b) : ptr(p), size(s), off(o), buf(std::move(b)) {}
    };

    std::map<int, std::unique_ptr<size_t[]> > calls_;
    // buffer pointer, size, offset at the target, and buffer
    std::map<int, std::shared_ptr<Prep> > requests_;

    void init() {
      // we should compute the number of messnages to receive?
      // receives
      const int dest = 0;
      std::unique_ptr<size_t[]> buf(new size_t[4]);
      // receives size,tag,rank
      const int rq = mpi__->request_recv(buf.get(), 4); 
      auto m = calls_.insert(std::make_pair(rq, std::move(buf)));
      assert(m.second);
    }

  public:
    AccRequest(double* const d = nullptr) : data_(d) {}

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
          std::shared_ptr<Double_ptr> ptr(new Double_ptr(buffer.get()));
          // sending back the buffer address
          const int sq = mpi__->request_send<Double_ptr>(&*ptr, sizeof(struct Double_ptr), rank, tag);
          const int rq = mpi__->request_recv(buffer.get(), size, rank, tag); 
          auto m = requests_.insert(std::make_pair(rq, std::shared_ptr<Prep>(new Prep(ptr, size, off, std::move(buffer)))));
          assert(m.second);
          i = calls_.erase(i);
        } else {
          ++i;
        }
      }

      for (auto i = requests_.begin(); i != requests_.end(); ++i) {
        if (mpi__->test(i->first)) {
          daxpy_(i->second->size, 1.0, i->second->buf.get(), 1, data_+i->second->off, 1);
          i = requests_.erase(i); 
        } else {
          ++i;
        }
      }
    }

};

}

#endif

