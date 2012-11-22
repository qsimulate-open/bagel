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

// Read only window for one-sided communication
class Window {
  protected:
    const double* const data_;
    const size_t size_;
#ifdef HAVE_MPI_H
    MPI::Win window_;
#endif

  public:
#ifdef HAVE_MPI_H
    Window(const double* a, const size_t size, MPI::Win w) : data_(a), size_(size), window_(w) { }
    ~Window() { window_.Free(); }
#else
    Window(const double* a, const size_t size) : data_(a), size_(size) { }
#endif
    size_t size() { return size_; }
    const double* data() { return data_; }
};


class MPI_Interface {
  protected:

  public:
    MPI_Interface(int argc, char** argv);
    ~MPI_Interface();

    int rank() const;
    int size() const;

    // collective functions
    void barrier() const;
    void reduce(double*, const size_t size, const int root) const;
    void allreduce(double*, const size_t size) const;
    void broadcast(double*, const size_t size, const int root) const;
    void allgather(double* send, const size_t ssize, double* rec, const size_t rsize) const; 
    void allgather(int* send, const size_t ssize, int* rec, const size_t rsize) const; 

    // one sided communication
    std::shared_ptr<Window> create_window(double*, const size_t size) const;
    void get(double*, double*, const size_t, std::shared_ptr<Window>) const; 
};

extern MPI_Interface* mpi__; 

}


#endif

