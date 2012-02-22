//
// Newint - Parallel electron correlation program.
// Filename: task.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#ifndef __SRC_SMITH_TASK_H
#define __SRC_SMITH_TASK_H

namespace SMITH {

// base class for Task objects
// assumes that the operation table is static (not adjustable at runtime). 
template <typename T>
class Task {
  protected:
    std::list<std::shared_ptr<Task<T> > > depend_;
    std::list<Task<T>*> target_;
    bool done_;
    virtual void compute_() = 0;

  public:
    Task() : done_(false) {};
    Task(std::list<std::shared_ptr<Task<T> > >& d) : depend_(d), done_(false) {};
    ~Task() { };
    void compute() {
      compute_();
      done_ = true;
    };

    void add_dep(std::shared_ptr<Task<T> > a) {
      depend_.push_back(a);
      a->set_target(this);
    };
    void set_target(Task<T>* b) { target_.push_back(b); };

    bool ready() const {
      bool out = true;
      for (auto i = depend_.begin(); i != depend_.end(); ++i) out &= (*i)->ready();
      return out && !done_;
    };
};

}

#endif
