//
// Newint - Parallel electron correlation program.
// Filename: taskqueue.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_UTIL_TASKQUEUE_H
#define __SRC_UTIL_TASKQUEUE_H

#include <list>
#include <atomic>
#include <memory>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>


template<typename T>
class TaskQueue {
  protected:
    std::list<T> task_;
    std::list<std::shared_ptr<std::atomic_flag> > flag_;

  public:
    TaskQueue(std::list<T>& t) : task_(t) {
      for (int i = 0; i != task_.size(); ++i) {
        flag_.push_back(std::shared_ptr<std::atomic_flag>(new std::atomic_flag(ATOMIC_FLAG_INIT)));
      }
    }

    void compute(const int num_threads = 1) {
#if 0
      std::list<std::shared_ptr<boost::thread> > threads;
      for (int i = 0; i != 2; ++i)
//    for (int i = 0; i != num_threads; ++i)
        threads.push_back(std::shared_ptr<boost::thread>(new boost::thread(boost::bind(&TaskQueue<T>::compute_one_thread, this))));

      for (auto i = threads.begin(); i != threads.end(); ++i) (*i)->join();
#else
      for (int i = 0; i < task_.size(); ++i) {
        auto a = task_.begin();
        for (int j = 0; j != i; ++j) ++a;
        a->compute();
      }
#endif
    } 

    void compute_one_thread() {
      auto j = task_.begin();
      for (auto i = flag_.begin(); i != flag_.end(); ++i, ++j) {
        const bool used = (*i)->test_and_set();
        if (!used) {
          j->compute();
        }
      }
    }
};

#endif
