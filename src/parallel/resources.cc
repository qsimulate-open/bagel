//
// BAGEL - Parallel electron correlation program.
// Filename: stackmem.cc
// Copyright (C) 2011 Toru Shiozaki
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


#include <cassert>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <src/rysint/macros.h>
#include <src/parallel/resources.h>

using namespace std;
using namespace bagel;

StackMem::StackMem() : pointer_(0LU), total_(10000000LU) { // TODO 80MByte
  stack_area_ = unique_ptr<double[]>(new double[total_]);

  // in case we use Libint for ERI
#ifdef LIBINT_INTERFACE
  // TODO 20LU should not be hardwired
  libint_t_ = unique_ptr<Libint_t[]>(new Libint_t[20LU*20LU*20LU*20LU]);
  if (libint2_need_memory_3eri1(LIBINT_MAX_AM) < libint2_need_memory_eri(LIBINT_MAX_AM)) {
    LIBINT2_PREFIXED_NAME(libint2_init_eri)(&libint_t_[0], LIBINT_MAX_AM, 0);
  } else {
    LIBINT2_PREFIXED_NAME(libint2_init_3eri1)(&libint_t_[0], LIBINT_MAX_AM, 0);
  }
#endif
}


double* StackMem::get(const size_t size) {
  assert(pointer_+size < total_);
  double* out = stack_area_.get() + pointer_;
  pointer_ += size;
  return out;
}


void StackMem::release(const size_t size, double* p) {
  pointer_ -= size;
  assert(p == stack_area_.get()+pointer_ || size == 0);
}



Resources::Resources(const int n) : max_num_threads_(n) {
  for (int i = 0; i != max_num_threads_; ++i) {
    flag_.push_back(shared_ptr<atomic_flag>(new atomic_flag(ATOMIC_FLAG_INIT)));
    stackmem_.push_back(shared_ptr<StackMem>(new StackMem()));
  }

#ifdef LIBINT_INTERFACE
  LIBINT2_PREFIXED_NAME(libint2_static_init)();
#endif
}


shared_ptr<StackMem> Resources::get() {
  for (int i = 0; i != max_num_threads_; ++i) {
    bool used = flag_[i]->test_and_set();
    if (!used) {
      return stackmem_[i];
    }
  }
  throw runtime_error("Stack Memory exhausted");
  return stackmem_.front();
}


void Resources::release(shared_ptr<StackMem> o) {
  bool found = false;
  for (int i = 0; i != max_num_threads_; ++i) {
    if (stackmem_[i] == o) {
      o->clear();
      flag_[i]->clear();
      found = true;
      break;
    }
  }
  if (!found) throw logic_error("This should not happen: resources.h");
}
