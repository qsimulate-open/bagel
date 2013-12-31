//
// BAGEL - Parallel electron correlation program.
// Filename: task.h
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


#ifndef __SRC_SMITH_TASK_H
#define __SRC_SMITH_TASK_H

#include <stddef.h>
#include <list>
#include <memory>

namespace bagel {
namespace SMITH {

// base class for Task objects
// assumes that the operation table is static (not adjustable at runtime).
template <typename T>
class Task : public std::enable_shared_from_this<Task<T>> {
  protected:
    std::list<std::shared_ptr<Task<T>>> depend_;
    std::list<std::weak_ptr<Task<T>>> target_;
    bool done_;
    virtual void compute_() = 0;

  public:
    Task() : done_(false) {}
    Task(std::list<std::shared_ptr<Task<T>>>& d) : depend_(d), done_(false) {}
    ~Task() { }
    void compute() {
      compute_();
      done_ = true;
    }

    void add_dep(std::shared_ptr<Task<T>> a) {
      depend_.push_back(a);
      a->set_target(this->shared_from_this());
    }

    void delete_dep(std::shared_ptr<Task<T>> a) {
      auto iter = std::find(depend_.begin(), depend_.end(), a);
      if (iter != depend_.end()) depend_.erase(iter);
      // depend_ should not have duplicated records.
      assert(std::find(depend_.begin(), depend_.end(), a) == depend_.end());
    }

    void set_target(std::shared_ptr<Task<T>> b) { target_.push_back(b); }

    void initialize() { done_ = false; }

    bool done() const { return done_; }

    bool ready() const { return !done_ && std::all_of(depend_.begin(), depend_.end(), [](std::shared_ptr<Task<T>> o) { return o->done(); }); }

    virtual double energy() const { return 0.0; }
    virtual double dedci() const { return 0.0; }
    virtual double correction() const { return 0.0; }
};


template <typename T>
class EnergyTask : public Task<T> {
  protected:
    double energy_;
  public:
    EnergyTask() : Task<T>() {}
    ~EnergyTask() {}
    double energy() const { return energy_; }

};

template <typename T>
class DedciTask : public Task<T> {
  protected:
    double dedci_;
  public:
    DedciTask() : Task<T>() {}
    ~DedciTask() {}
    double dedci() const { return dedci_; }

};

template <typename T>
class CorrectionTask : public Task<T> {
  protected:
    double correction_;
  public:
    CorrectionTask() : Task<T>() {}
    ~CorrectionTask() {}
    double correction() const { return correction_; }

};

template <typename T>
class DensityTask : public Task<T> {
  protected:
  public:
    DensityTask() : Task<T>() {}
    ~DensityTask() {}

};

template <typename T>
class Density2Task : public Task<T> {
  protected:
  public:
    Density2Task() : Task<T>() {}
    ~Density2Task() {}

};

}
}

#endif
