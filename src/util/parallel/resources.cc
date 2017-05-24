//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: resources.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <src/util/parallel/resources.h>

using namespace std;
using namespace bagel;

StackMem::StackMem() : pointer_(0LU), total_(20000000LU) { // TODO this should not be hardwired
  stack_area_ = unique_ptr<double[]>(new double[total_]);

  // in case we use Libint for ERI
#ifdef LIBINT_INTERFACE
  // TODO 20LU should not be hardwired
  libint_t_ = unique_ptr<Libint_t[]>(new Libint_t[20LU*20LU*20LU*20LU]);
  if (libint2_need_memory_3eri1(LIBINT2_MAX_AM_3eri1) < libint2_need_memory_eri(LIBINT2_MAX_AM_eri))
    LIBINT2_PREFIXED_NAME(libint2_init_eri)(&libint_t_[0], LIBINT2_MAX_AM_eri, 0);
  else
    LIBINT2_PREFIXED_NAME(libint2_init_3eri1)(&libint_t_[0], LIBINT2_MAX_AM_3eri1, 0);
#endif
}


Resources::Resources(const int max) : proc_(make_shared<Process>()), max_num_threads_(max) {
#ifdef LIBINT_INTERFACE
  LIBINT2_PREFIXED_NAME(libint2_static_init)();
#endif
  for (int i = 0; i != max; ++i)
    stackmem_[make_shared<StackMem>()].clear();
}


shared_ptr<StackMem> Resources::get() {
  for (auto& i : stackmem_) {
    if (!i.second.test_and_set())
      return i.first;
  }
  // This error most often occurs if we forget to destruct one integral object before constructing another
  throw runtime_error("Stack Memory exhausted");
  return nullptr;
}


void Resources::release(shared_ptr<StackMem> o) {
  o->clear();
  auto iter = stackmem_.find(o);
  assert(iter != stackmem_.end());
  iter->second.clear();
}

