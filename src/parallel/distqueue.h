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

#include <src/util/constants.h>
#include <src/parallel/mpi_interface.h>

// TODO until GCC fixes this bug
#ifdef __GNUC__
#if __GNUC__ == 4 && __GNUC_MINOR__ <= 7
#define _GLIBCXX_USE_NANOSLEEP
#endif
#endif
#include <thread>

#include <vector>
#include <bagel_config.h>
#ifdef HAVE_MKL_H
  #include "mkl_service.h"
#endif
#include <src/parallel/resources.h>

namespace bagel {

template <typename TaskType, typename... FlushTypes>
class DistQueue {
  protected:
    std::tuple<FlushTypes...> flushed_;

    std::list<std::shared_ptr<TaskType>> tasks_;
    std::mutex listmut_;

    std::atomic<bool> done_;

    public:
      DistQueue(FlushTypes... f) : flushed_(std::make_tuple(f...)), done_(false) {}

      template <typename... Args>
      void emplace_and_compute(Args&&... args) {
        {
          auto t = std::make_shared<TaskType>(std::forward<Args>(args)...);
          if (t->test()) {
            t->compute();
          }
          else {
            std::lock_guard<std::mutex> lock(listmut_);
            tasks_.push_back(t);
          }
        }

        bool full_pass = false;
        do {
          std::shared_ptr<TaskType> task;
          {
            std::lock_guard<std::mutex> lock(listmut_);
            for (auto i = tasks_.begin(); i != tasks_.end(); ) {
              if ((*i)->test()) {
                task = *i;
                i = tasks_.erase(i);
                break;
              }
              else {
                ++i;
              }
            }
          }
          if (task) task->compute();
          else full_pass = true;
        } while (!full_pass);

        flush();
      }

      template <typename... Args>
      void emplace(Args&&... args) {
        auto t = std::make_shared<TaskType>(std::forward<Args>(args)...);
        if (t->test()) {
          t->compute();
        }
        else {
          std::lock_guard<std::mutex> lock(listmut_);
          tasks_.push_back(t);
        }
      }

      // Only one thread should call this version
      void finish_master() {
        bool done;
        do {
          bool full_pass = false;
          do {
            done = true;
            std::shared_ptr<TaskType> task;
            {
              std::lock_guard<std::mutex> lock(listmut_);
              for (auto i = tasks_.begin(); i != tasks_.end(); ) {
                if ((*i)->test()) {
                  task = *i;
                  i = tasks_.erase(i);
                  break;
                }
                else {
                  ++i;
                  done = false;
                }
              }
            }
            if (task) task->compute();
            else full_pass = true;
          } while (!full_pass);
#ifndef USE_SERVER_THREAD
          size_t d = done ? 0 : 1;
          mpi__->soft_allreduce(&d, 1);
          done = (d == 0);
#endif
          if (!done) flush();
          if (!done) std::this_thread::sleep_for(sleeptime__);
        } while (!done);
        done_ = true;
      }

      // Similar to above, except it doesn't update done_ and doesn't flush
      void finish_worker() {
        do {
          bool full_pass = false;
          do {
            std::shared_ptr<TaskType> task;
            {
              std::lock_guard<std::mutex> lock(listmut_);
              for (auto i = tasks_.begin(); i != tasks_.end(); ) {
                if ((*i)->test()) {
                  task = *i;
                  i = tasks_.erase(i);
                  break;
                }
                else {
                  ++i;
                }
              }
            }
            if (task) task->compute();
            else full_pass = true;
          } while (!full_pass);
          if (!done_) std::this_thread::sleep_for(sleeptime__);
        } while (!done_);
      }

      // for convenience
      void finish() { finish_master(); }

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
