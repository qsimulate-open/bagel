//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: resources.h
// Copyright (C) 2011 Quantum Simulation Technologies, Inc.
//
// Author: Toru Shiozaki <shiozaki@qsimulate.com>
// Maintainer: QSimulate
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

#ifndef __SRC_PARALLEL_RESOURCES_H
#define __SRC_PARALLEL_RESOURCES_H

// CAUTION last-in-first-out stack to avoid the overhead of new'ing every time

#include <stddef.h>
#include <deque>
#include <memory>
#include <atomic>
#include <cassert>
#include <stdexcept>
#include <complex>
#include <map>
#ifdef LIBINT_INTERFACE
  #include <libint2.h>
#endif
#include <src/util/constants.h>
#include <src/util/parallel/process.h>

namespace bagel {

class StackMemLine {
  protected:
    std::unique_ptr<double[]> stack_area_;
    size_t pointer_;
    const size_t total_;

  public:
    StackMemLine(const size_t total=20000000LU) : pointer_(0LU), total_(total) {
      stack_area_ = std::unique_ptr<double[]>(new double[total_]);
    }

    template <typename DataType = double>
    DataType* get(const size_t size) {
      assert(size * sizeof(DataType) % sizeof(double) == 0);
      size_t new_pointer = pointer_ + size * sizeof(DataType) / sizeof(double);
      if (new_pointer > total_) { return NULL; }

      DataType* out = reinterpret_cast<DataType*> (stack_area_.get() + pointer_);
      pointer_ = new_pointer;
      return out;
    }

    template <typename DataType = double>
    void release(const size_t size, DataType* p) {
      assert(size * sizeof(DataType) % sizeof(double) == 0);
      pointer_ -= (size * sizeof(DataType) / sizeof(double));
      assert(p == reinterpret_cast<DataType*> (stack_area_.get()+pointer_) || size == 0);
    }

    size_t size() const { return pointer_; }

    bool is_empty() const { return size() == 0LU; }

    size_t max_size() const { return total_; }

    void clear() { pointer_ = 0LU; }
};


class StackMem {
  protected:
    std::deque<StackMemLine> stack_lines_;
    size_t current_size_;
    size_t water_line_;

#ifdef LIBINT_INTERFACE
    std::unique_ptr<Libint_t[]> libint_t_;
#endif

  public:
    StackMem(size_t init_line_size=20000000LU);

    template <typename DataType = double>
    DataType* get(const size_t size) {
      assert((size * sizeof(DataType)) % sizeof(double) == 0);
      const size_t word_size = size * sizeof(DataType) / sizeof(double);
      current_size_ += word_size;
      water_line_ = std::max(water_line_, current_size_);
      DataType* out = stack_lines_.back().template get<DataType>(size);
      if (out == NULL) {
        stack_lines_.push_back(StackMemLine(word_size));
        out = stack_lines_.back().template get<DataType>(size);
        assert(stack_lines_.back().size() == stack_lines_.back().max_size());
      }
      assert(out != NULL);
      return out;
    }

    template <typename DataType = double>
    void release(const size_t size, DataType* p) {
      assert((size * sizeof(DataType)) % sizeof(double) == 0);
      const size_t word_size = size * sizeof(DataType) / sizeof(double);
      current_size_ -= word_size;
      assert(stack_lines_.size() > 0);

      stack_lines_.back().template release<DataType>(size, p);
      if (stack_lines_.size() > 1) {
        assert(stack_lines_.back().size() == 0);
        stack_lines_.pop_back();
      }
      if (is_empty())  { clear(); }
    }

    bool is_empty() const { return current_size_ == 0LU; }

    void clear();

#ifdef LIBINT_INTERFACE
    Libint_t* libint_t_ptr(const int i) { return &libint_t_[i]; }
#endif
};


class Resources {
  private:
    std::shared_ptr<Process> proc_;
    std::map<std::shared_ptr<StackMem>, std::atomic_flag> stackmem_;
    size_t max_num_threads_;

  public:
    Resources(const int max);

    std::shared_ptr<StackMem> get();
    void release(std::shared_ptr<StackMem> o);

    size_t max_num_threads() const { return max_num_threads_; }
    std::shared_ptr<Process> proc() { return proc_; }
};

extern Resources* resources__;
extern bool nocompute__;

}

#endif
