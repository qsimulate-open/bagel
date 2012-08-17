//
// Author : Toru Shiozaki
// Date   : August 2012
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
    TaskQueue(std::list<T> t) : task_(t) {
      for (int i = 0; i != task_.size(); ++i) {
        flag_.push_back(std::shared_ptr<std::atomic_flag>(new std::atomic_flag(ATOMIC_FLAG_INIT)));
      }
    }

    std::pair<T, bool> next() {
      auto j = task_.begin();
      for (auto i = flag_.begin(); i != flag_.end(); ++i, ++j)
        if (!(*i)->test_and_set()) return std::make_pair(*j, true);
      return std::make_pair(T(), false);
    }
};

#endif
