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
      for (auto a2 = virt_.begin(); a2 != virt_.end(); ++a2) {
        for (auto c3 = closed_.begin(); c3 != closed_.end(); ++c3) {
          for (auto a4 = virt_.begin(); a4 != virt_.end(); ++a4) {
            for (auto c1 = closed_.begin(); c1 != closed_.end(); ++c1) {
              std::vector<size_t> ohash = vec(c1->key(), a4->key(), c3->key(), a2->key());
              std::unique_ptr<double[]> odata = I0->move_block(ohash);

              std::vector<size_t> i0hash = vec(c1->key(), a4->key(), c3->key(), a2->key());
              std::unique_ptr<double[]> i0data = t2->get_block(i0hash);

              std::vector<size_t> i1hash = vec();
              std::unique_ptr<double[]> i1data = I1->get_block(i1hash);


              I0->put_block(ohash, odata);
            }
          }
        }
      }
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
      for (auto a2 = virt_.begin(); a2 != virt_.end(); ++a2) {
        for (auto c3 = closed_.begin(); c3 != closed_.end(); ++c3) {
          for (auto a4 = virt_.begin(); a4 != virt_.end(); ++a4) {
            for (auto c1 = closed_.begin(); c1 != closed_.end(); ++c1) {
              std::vector<size_t> ohash = vec(c1->key(), a4->key(), c3->key(), a2->key());
              std::unique_ptr<double[]> odata = I0->move_block(ohash);

              std::vector<size_t> i0hash = vec(c1->key(), a2->key(), c3->key(), a4->key());
              std::unique_ptr<double[]> i0data = t2->get_block(i0hash);

              std::vector<size_t> i1hash = vec();
              std::unique_ptr<double[]> i1data = I3->get_block(i1hash);


              I0->put_block(ohash, odata);
            }
          }
        }
      }
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
      for (auto c3 = closed_.begin(); c3 != closed_.end(); ++c3) {
        for (auto a2 = virt_.begin(); a2 != virt_.end(); ++a2) {
          for (auto a4 = virt_.begin(); a4 != virt_.end(); ++a4) {
            for (auto c1 = closed_.begin(); c1 != closed_.end(); ++c1) {
              std::vector<size_t> ohash = vec(c1->key(), a4->key(), a2->key(), c3->key());
              std::unique_ptr<double[]> odata = I4->move_block(ohash);

              for (auto c5 = closed_.begin(); c5 != closed_.end(); ++c5) {
                std::vector<size_t> i0hash = vec(c3->key(), c5->key());
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);

                std::vector<size_t> i1hash = vec(c1->key(), a4->key(), c5->key(), a2->key());
                std::unique_ptr<double[]> i1data = I5->get_block(i1hash);

              }

              I4->put_block(ohash, odata);
            }
          }
        }
      }
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
      for (auto c3 = closed_.begin(); c3 != closed_.end(); ++c3) {
        for (auto a2 = virt_.begin(); a2 != virt_.end(); ++a2) {
          for (auto a4 = virt_.begin(); a4 != virt_.end(); ++a4) {
            for (auto c1 = closed_.begin(); c1 != closed_.end(); ++c1) {
              std::vector<size_t> ohash = vec(c1->key(), a4->key(), a2->key(), c3->key());
              std::unique_ptr<double[]> odata = I4->move_block(ohash);

              for (auto a5 = virt_.begin(); a5 != virt_.end(); ++a5) {
                std::vector<size_t> i0hash = vec(a5->key(), a4->key());
                std::unique_ptr<double[]> i0data = f1->get_block(i0hash);

                std::vector<size_t> i1hash = vec(c1->key(), a5->key(), c3->key(), a2->key());
                std::unique_ptr<double[]> i1data = I9->get_block(i1hash);

              }

              I4->put_block(ohash, odata);
            }
          }
        }
      }
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

