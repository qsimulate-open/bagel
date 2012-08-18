//
// Newint - Parallel electron correlation program.
// Filename: stackmem.h
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


#ifndef __SRC_UTIL_RESOURCES_H
#define __SRC_UTIL_RESOURCES_H

// CAUTION last-in-first-out stack to avoid the overhead of new'ing every time

#include <memory>
#include <atomic>
#include <vector>

class StackMem {
  protected:
    std::unique_ptr<double[]> stack_area_;
    size_t pointer_;
    const size_t total_;

  public:
    StackMem();

    double* get(const size_t size);
    void release(const size_t size, double* p);

    void clear() { pointer_ = 0LU; };
};


class Resources {
  private:
    std::vector<std::shared_ptr<StackMem> > stackmem_;
    std::vector<std::shared_ptr<std::atomic_flag> > flag_;
    size_t n_;

  public:
    Resources(const int n);
    ~Resources() {};

    std::shared_ptr<StackMem> get();
    void release(std::shared_ptr<StackMem> o);
};

extern Resources* resources__;

#endif
