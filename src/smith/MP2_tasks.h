//
// Newint - Parallel electron correlation program.
// Filename: MP2_tasks.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and/or modify
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


#ifndef __SRC_SMITH_MP2_TASKS_H 
#define __SRC_SMITH_MP2_TASKS_H 

#include <memory>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <vector>

namespace SMITH {
namespace MP2{

template <typename T>
class Task0 : public Task<T> {
  protected:
    std::shared_ptr<Tensor<T> > r2_;
    IndexRange closed_;
    IndexRange act_;
    IndexRange virt_;

    void compute_() {
      r2_->zero();
    };  

  public:
    Task0(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      if (t.size() != 1) throw std::logic_error("Task0 error");
      r2_ =  t[0];
      closed_ = i[0];
      act_    = i[1];
      virt_   = i[2];
    };  
    ~Task0() {}; 
};

template <typename T>
class Task1 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange act_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I0;

    void compute_() {
    };

  public:
    Task1(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      act_    = i[1];
      virt_   = i[2];
      r = t[0];
      I0 = t[1];
    };
    ~Task1() {};
};

template <typename T>
class Task2 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange act_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I1;

    void compute_() {
    };

  public:
    Task2(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      act_    = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I1 = t[2];
    };
    ~Task2() {};
};

template <typename T>
class Task3 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange act_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I0;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > I3;

    void compute_() {
    };

  public:
    Task3(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      act_    = i[1];
      virt_   = i[2];
      I0 = t[0];
      t2 = t[1];
      I3 = t[2];
    };
    ~Task3() {};
};

template <typename T>
class Task4 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange act_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > r;
    std::shared_ptr<Tensor<T> > I4;

    void compute_() {
    };

  public:
    Task4(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      act_    = i[1];
      virt_   = i[2];
      r = t[0];
      I4 = t[1];
    };
    ~Task4() {};
};

template <typename T>
class Task5 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange act_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I4;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I5;

    void compute_() {
    };

  public:
    Task5(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      act_    = i[1];
      virt_   = i[2];
      I4 = t[0];
      f1 = t[1];
      I5 = t[2];
    };
    ~Task5() {};
};

template <typename T>
class Task6 : public Task<T> {
  protected:
    IndexRange closed_;
    IndexRange act_;
    IndexRange virt_;
    std::shared_ptr<Tensor<T> > I4;
    std::shared_ptr<Tensor<T> > f1;
    std::shared_ptr<Tensor<T> > I9;

    void compute_() {
    };

  public:
    Task6(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      closed_ = i[0];
      act_    = i[1];
      virt_   = i[2];
      I4 = t[0];
      f1 = t[1];
      I9 = t[2];
    };
    ~Task6() {};
};


}
}
#endif

