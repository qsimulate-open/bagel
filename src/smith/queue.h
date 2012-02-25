//
// Newint - Parallel electron correlation program.
// Filename: queue.h
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


// Impliments a Task queue; in the future this should be replaced by
// some standard interface (boost::graph?)

#ifndef __SRC_SMITH_QUEUE_H
#define __SRC_SMITH_QUEUE_H

#include <src/smith/task.h>
#include <cassert>
#include <list>
#include <memory>

namespace SMITH {

template <typename T>
class Queue {
  protected:
    std::list<std::shared_ptr<Task<T> > > tasklist_;
    int cnt_;

  public:
    Queue() : cnt_(0) {};
    Queue(std::list<std::shared_ptr<Task<T> > > d) : tasklist_(d), cnt_(0) { };
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

    void add_task(std::shared_ptr<Task<T> > a) { tasklist_.push_back(a); };

    bool done() const { return tasklist_.size() == cnt_; };

    void initialize() {
      for (auto i = tasklist_.begin(); i != tasklist_.end(); ++i)
        (*i)->initialize();
      cnt_ = 0;
    };
};

}

#endif
