//
// BAGEL - Parallel electron correlation program.
// Filename: distqueue.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_PARALLEL_DISTQUEUE_H
#define __SRC_PARALLEL_DISTQUEUE_H

#include <stddef.h>
#include <list>
#include <memory>
#include <utility>
#include <atomic>
#include <thread>
#include <vector>

#include <src/util/constants.h>
#include <src/parallel/plist.h>
#include <src/parallel/mpi_interface.h>

#include <bagel_config.h>
#ifdef HAVE_MKL_H
  #include "mkl_service.h"
#endif
#include <src/parallel/resources.h>

namespace bagel {

template <typename TaskType, typename... FlushTypes>
class DistQueue : public ParallelList<TaskType> {
  protected:
    std::tuple<FlushTypes...> flushed_;

    std::atomic<bool> done_;

    public:
      DistQueue(FlushTypes... f) : flushed_(std::make_tuple(f...)), done_(false) {}
      ~DistQueue() { this->remove_if( [] (const TaskType& t) { return true; } ); }

      template <typename... Args>
      void emplace_and_compute(Args&&... args) {
        auto t = std::make_shared<TaskType>(std::forward<Args>(args)...);

        if (t->test()) {
          t->compute();
        }
        else {
          this->push_front(t);
        }

        flush();

        while (std::shared_ptr<TaskType> task = pop_first_ready()) {
          task->compute();
        }
      }

      template <typename... Args>
      void emplace(Args&&... args) {
        this->emplace_front(std::forward<Args>(args)...);
      }

      // Only one thread should call this version
      void finish_master() {
        bool done;
        do {
          while (std::shared_ptr<TaskType> t = pop_first_ready())
            t->compute();

          done = this->empty();

          size_t d = done ? 0 : 1;
          mpi__->soft_allreduce(&d, 1);
          done = (d == 0);

          flush();

          if (!done) std::this_thread::sleep_for(sleeptime__);
        } while (!done);
        done_ = true;
      }

      // Similar to above, except it doesn't update done_ and doesn't flush
      // TODO unclear whether this is working
      void finish_worker() {
        do {
          while (std::shared_ptr<TaskType> t = pop_first_ready())
            t->compute();

          if (!done_) std::this_thread::sleep_for(sleeptime__);
        } while (!done_);
      }

      // for convenience
      void finish() { finish_master(); }

      std::shared_ptr<TaskType> pop_first_ready() {
        return this->pop_first_if([] (TaskType& t) { return t.test(); });
      }


    private:
      // Calls flush() on all of the flushable objects
      void flush() { _flush<0, FlushTypes...>(); }

      template <int iter, typename Head, typename... Tail> void _flush() {
        std::get<iter>(flushed_)->flush();
        _flush<iter+1, Tail...>();
      }
      template <int iter> void _flush() {}
};

}

#endif
