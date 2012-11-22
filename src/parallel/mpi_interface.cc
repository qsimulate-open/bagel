//
// BAGEL - Parallel electron correlation program.
// Filename: mpi_interface.cc
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

#include <stdexcept>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

MPI_Interface::MPI_Interface(int argc, char** argv) {
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
#endif
}


MPI_Interface::~MPI_Interface() {
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif
}


int MPI_Interface::rank() const {
#ifdef HAVE_MPI_H
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  return rank;
#else
  return 0;
#endif
}


int MPI_Interface::size() const {
#ifdef HAVE_MPI_H
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
#else
  return 1;
#endif
}


void MPI_Interface::barrier() const {
#ifdef HAVE_MPI_H
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void MPI_Interface::reduce(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI_Reduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI::SUM, root, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::allreduce(double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::broadcast(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI_Bcast(static_cast<void*>(a), size, MPI_DOUBLE, root, MPI_COMM_WORLD);
#endif
}


void MPI_Interface::allgather(double* send, const size_t ssize, double* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  MPI_Allgather(static_cast<void*>(send), ssize, MPI_DOUBLE, static_cast<void*>(rec), rsize, MPI_DOUBLE, MPI_COMM_WORLD);
#endif
} 


void MPI_Interface::allgather(int* send, const size_t ssize, int* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  MPI_Allgather(static_cast<void*>(send), ssize, MPI_INT, static_cast<void*>(rec), rsize, MPI_INT, MPI_COMM_WORLD);
#endif
} 


shared_ptr<Window> MPI_Interface::create_window(double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  MPI_Win w;
  MPI_Win_create(static_cast<void*>(a), size*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &w);
  shared_ptr<Window> win(new Window(a, size, w));
#else
  shared_ptr<Window> win(new Window(a, size));
#endif
  return win;
}


void MPI_Interface::get(double* a, double* b, const size_t size, shared_ptr<Window> win) const {
#if 0
#ifdef HAVE_MPI_H
  win->Get(static_cast<void*>(b), size, MPI_DOUBLE,  
#else
  throw logic_error("MPI_Interface::get should not be called without MPI");
#endif
#endif
}
