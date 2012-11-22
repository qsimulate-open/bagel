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

#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

MPI_Interface::MPI_Interface() {
#ifdef HAVE_MPI_H
  MPI::Init();
#endif
}


MPI_Interface::~MPI_Interface() {
#ifdef HAVE_MPI_H
  MPI::Finalize();
#endif
}


int MPI_Interface::rank() const {
#ifdef HAVE_MPI_H
  return MPI::COMM_WORLD.Get_rank(); 
#else
  return 0;
#endif
}


int MPI_Interface::size() const {
#ifdef HAVE_MPI_H
  return MPI::COMM_WORLD.Get_size(); 
#else
  return 1;
#endif
}


void MPI_Interface::barrier() const {
#ifdef HAVE_MPI_H
  MPI::COMM_WORLD.Barrier();
#endif
}

void MPI_Interface::reduce(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI::COMM_WORLD.Reduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI::SUM, root);
#endif
}


void MPI_Interface::allreduce(double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, static_cast<void*>(a), size, MPI_DOUBLE, MPI::SUM);
#endif
}


void MPI_Interface::broadcast(double* a, const size_t size, const int root) const {
#ifdef HAVE_MPI_H
  MPI::COMM_WORLD.Bcast(static_cast<void*>(a), size, MPI_DOUBLE, root);
#endif
}


void MPI_Interface::allgather(const double* send, const size_t ssize, double* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  MPI::COMM_WORLD.Allgather(static_cast<const void*>(send), ssize, MPI_DOUBLE, static_cast<void*>(rec), rsize, MPI_DOUBLE);
#endif
} 


void MPI_Interface::allgather(const int* send, const size_t ssize, int* rec, const size_t rsize) const {
#ifdef HAVE_MPI_H
  MPI::COMM_WORLD.Allgather(static_cast<const void*>(send), ssize, MPI_INT, static_cast<void*>(rec), rsize, MPI_INT);
#endif
} 


shared_ptr<Window> MPI_Interface::create_window(const double* a, const size_t size) const {
#ifdef HAVE_MPI_H
  shared_ptr<Window> win(new Window(a, size, MPI::Win::Create(static_cast<const void*>(a), size*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI::COMM_WORLD)));
#else
  shared_ptr<Window> win(new Window(a, size));
#endif
  return win;
}

