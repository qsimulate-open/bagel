//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: task.h
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


#ifndef __SRC_SMITH_TASK_H
#define __SRC_SMITH_TASK_H

#include <stddef.h>
#include <list>
#include <memory>
#include <algorithm>
#include <src/smith/storage.h>

namespace bagel {
namespace SMITH {

// base class for Task objects
// assumes that the operation table is static (not adjustable at runtime).
class Task : public std::enable_shared_from_this<Task> {
  protected:
    std::list<std::shared_ptr<Task>> depend_;
    std::list<std::weak_ptr<Task>> target_;
    bool done_;
    virtual void compute_() = 0;

  public:
    Task() : done_(false) {}
    Task(std::list<std::shared_ptr<Task>>& d) : depend_(d), done_(false) {}
    ~Task() { }
    void compute() {
      compute_();
      done_ = true;
    }

    void add_dep(std::shared_ptr<Task> a) {
      depend_.push_back(a);
      a->set_target(this->shared_from_this());
    }

    void delete_dep(std::shared_ptr<Task> a) {
      depend_.remove(a);
      assert(std::find(depend_.begin(), depend_.end(), a) == depend_.end());
    }

    void set_target(std::shared_ptr<Task> b) { target_.push_back(b); }

    void initialize() { done_ = false; }

    bool done() const { return done_; }

    bool ready() const { return !done_ && std::all_of(depend_.begin(), depend_.end(), [](std::shared_ptr<Task> o) { return o->done(); }); }

    virtual double target() const { return 0.0; }
};


class AccTask : public Task {
  protected:
    double target_;
  public:
    AccTask() {}
    ~AccTask() {}
    double target() const { return target_; }

};

}
}

#endif
