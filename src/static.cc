//
// BAGEL - Parallel electron correlation program.
// Filename: static.cc
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

#ifdef _OPENMP
  #include <omp.h>
#endif
#include <thread>
#include <src/global.h>
#include <src/util/string.h>
#include <src/parallel/mpi_interface.h>
#include <src/parallel/resources.h>

// They are used from other files
namespace bagel{
  Resources* resources__;
  MPI_Interface* mpi__;

  std::unique_ptr<Resources> resources;
  std::unique_ptr<MPI_Interface> mpi;
}

using namespace bagel;
using namespace std;

namespace bagel {

void static_variables() {
  // setup MPI interface. It does nothing for serial runs
  mpi = unique_ptr<MPI_Interface>(new MPI_Interface());
  mpi__ = mpi.get();
  {
    string snum_threads = getenv_multiple("BAGEL_NUM_THREADS", "OMP_NUM_THREADS");
    const int num_threads = snum_threads.empty() ? thread::hardware_concurrency() : lexical_cast<int>(snum_threads);
    if (num_threads < 1)
      throw runtime_error("Set BAGEL_NUM_THREADS for the number of threads used");
    if (mpi__->rank() == 0)
      cout << "  * using " << num_threads << " threads per process" << endl;
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
    resources = unique_ptr<Resources>(new Resources(num_threads));
    resources__ = resources.get();
  }
}

}
