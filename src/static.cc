//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: static.cc
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

#ifdef _OPENMP
  #include <omp.h>
#endif
#include <cfenv>
#include <thread>
#include <src/global.h>
#include <src/util/string_util.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/parallel/resources.h>

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

  // rounding mode in std::rint, std::lrint, and std::llrint
  fesetround(FE_TONEAREST);
}

}
