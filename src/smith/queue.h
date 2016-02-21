//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: queue.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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


// Impliments a Task queue; in the future this should be replaced by
// some standard interface (boost::graph?)

#ifndef __SRC_SMITH_QUEUE_H
#define __SRC_SMITH_QUEUE_H

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/task.h>
#include <cassert>
#include <list>
#include <memory>
#include <algorithm>

namespace bagel {
namespace SMITH {

class Queue {
  protected:
    std::list<std::shared_ptr<Task>> tasklist_;

  public:
    Queue() {}
    Queue(const std::list<std::shared_ptr<Task>>& d) : tasklist_(d) { }

    std::shared_ptr<Task> next_compute();

    void add_task(std::shared_ptr<Task> a) { tasklist_.push_back(a); }

    void insert(std::shared_ptr<Queue> b) {
      for (auto& i : b->tasklist_)
        tasklist_.push_back(i);
    }

    bool done() const { return tasklist_.empty(); }

    void initialize() {
      for (auto& i : tasklist_) i->initialize();
    }
};

}
}

#endif
#endif
