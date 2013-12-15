//
// BAGEL - Parallel electron correlation program.
// Filename: taskqueue.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_UTIL_TASKQUEUE_H
#define __SRC_UTIL_TASKQUEUE_H

#include <stddef.h>
#include <list>
#include <atomic>
#include <memory>
#include <stdexcept>
#include <utility>
#include <thread>

#include <vector>
#include <bagel_config.h>
#ifdef HAVE_MKL_H
  #include "mkl_service.h"
#endif
#include <src/parallel/resources.h>

namespace bagel {

namespace {
  template<typename T> void call_compute(T& task) { task.compute(); }
  template<typename T> void call_compute(std::shared_ptr<T>& task) { task->compute(); }
}

template<typename T>
class TaskQueue {
  protected:
    std::vector<T> task_;
    std::list<std::shared_ptr<std::atomic_flag>> flag_;
    static const int chunck_ = 12;

  public:
    TaskQueue(size_t expected = 0) { task_.reserve(expected); }
    TaskQueue(std::vector<T>&& t) : task_(std::move(t)) { }

    template<typename ...args>
    void emplace_back(args&&... a) { task_.emplace_back(std::forward<args>(a)...); }

    void compute(const int num_threads = resources__->max_num_threads()) {
      if (task_.empty()) return;
#ifdef HAVE_MKL_H
      const int mkl_num = mkl_get_max_threads();
      mkl_set_num_threads(1);
#endif
#ifndef _OPENMP
      for (int i = 0; i != (task_.size()-1)/chunck_+1; ++i)
        flag_.push_back(std::shared_ptr<std::atomic_flag>(new std::atomic_flag(ATOMIC_FLAG_INIT)));
      std::list<std::shared_ptr<std::thread>> threads;
      for (int i = 0; i != num_threads; ++i)
        threads.push_back(std::make_shared<std::thread>(&TaskQueue<T>::compute_one_thread, this));
      for (auto& i : threads) {
        i->join();
      }
#else
      const size_t n = task_.size();
      #pragma omp parallel for schedule(dynamic,chunck_)
      for (size_t i = 0; i < n; ++i)
        call_compute(task_[i]);
#endif
#ifdef HAVE_MKL_H
      mkl_set_num_threads(mkl_num);
#endif
    }

    void compute_one_thread() {
      int j = 0;
      for (auto i = flag_.begin(); i != flag_.end(); ++i, j += chunck_)
        if (!(*i)->test_and_set()) {
          call_compute(task_[j]);
          for (int k = 1; k < chunck_; ++k)
            if (j+k < task_.size()) call_compute(task_[j+k]);
        }
    }
};

}

#endif
