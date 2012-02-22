//
// Author : Toru Shiozaki
// Date   : feb 2012
//

// Impliments a Task queue; in the future this should be replaced by
// some standard interface (boost::graph?)

#ifndef __SRC_SMITH_QUEUE_H
#define __SRC_SMITH_QUEUE_H

#include <src/smith/task.h>
#include <cassert>

namespace SMITH {

template <typename T>
class Queue {
  protected:
    std::list<std::shared_ptr<Task<T> > > tasklist_;
    std::list<std::shared_ptr<Task<T> > >::iterator current_;
    int cnt_;

  public:
    Queue(std::list<std::shared_ptr<Task<T> > > d) : tasklist_(d), cnt_(0) { current_ = tasklist_.begin(); };
    ~Queue() {};


    std::shared_ptr<Task<T> > next() {
      auto i = tasklist_.begin();
      for ( ; i != tasklist_.end(); ++i) {
        if ((*i)->ready()) break;
      }
      assert(i != tasklist_.end());
      ++cnt_;
      return *i;
    };

    bool done() const { return tasklist_.size() == cnt_; };
};

}

#endif
