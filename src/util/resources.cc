//
// Newint - Parallel electron correlation program.
// Filename: stackmem.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <cassert>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <src/util/resources.h>

using namespace std;

StackMem::StackMem() : pointer_(0LU), total_(10000000LU) { // 80MByte
  stack_area_ = unique_ptr<double[]>(new double[total_]);
}


double* StackMem::get(const size_t size) {
  assert(pointer_+size < total_);
  double* out = stack_area_.get() + pointer_;
  pointer_ += size;
  return out;
}


void StackMem::release(const size_t size, double* p) {
  pointer_ -= size; 
  assert(p == stack_area_.get()+pointer_);
}



Resources::Resources(const int n) : n_(n) {
  for (int i = 0; i != n_; ++i) {
    flag_.push_back(shared_ptr<atomic_flag>(new atomic_flag(ATOMIC_FLAG_INIT)));
    stackmem_.push_back(shared_ptr<StackMem>(new StackMem()));
  }
} 


shared_ptr<StackMem> Resources::get() {
  for (int i = 0; i != n_; ++i) {
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
  for (int i = 0; i != n_; ++i) {
    if (stackmem_[i] == o) {
      o->clear();
      flag_[i]->clear(); 
      found = true;
      break;
    }
  }
  if (!found) throw logic_error("This should not happen: resources.h");
}
