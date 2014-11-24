//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks.h
// Copyright (C) 2012 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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


#ifndef __SRC_SMITH_CAS_test_TASKS_H
#define __SRC_SMITH_CAS_test_TASKS_H

#include <memory>
#include <algorithm>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <src/smith/subtask.h>
#include <vector>

namespace bagel {
namespace SMITH {
namespace CAS_test{

template <typename T>
class Task0 : public Task<T> {
  protected:
    std::shared_ptr<Tensor<T>> r_;
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;

    void compute_() {
      r_->zero();
    };

  public:
    Task0(std::vector<std::shared_ptr<Tensor<T>>> t) : Task<T>() {
      r_ =  t[0];
    };
    ~Task0() {};
};

template <typename T>
class Task1 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          // std::shared_ptr<Tensor<T> > Gamma0;
          // std::shared_ptr<Tensor<T> > rdm1;
          // std::shared_ptr<Tensor<T> > f1;

          // scalar
          // tensor label: Gamma0
          std::unique_ptr<double[]> odata = out()->move_block();
          // associated with merged
          std::unique_ptr<double[]> fdata = in(1)->get_block(x1, x0);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            for (int i0 = 0; i0 != x0.size(); ++i0) {
              for (int i1 = 0; i1 != x1.size(); ++i1) {
                odata[0]
                  += (1.0) * i0data[i1+x1.size()*(i0)] * fdata[i1+x1.size()*(i0)];
              }
            }
          }
          out()->put_block(odata);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task1(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task1() {};
};

template <typename T>
class Task2 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          // std::shared_ptr<Tensor<T> > Gamma2;
          // std::shared_ptr<Tensor<T> > rdm1;

          // tensor label: Gamma2
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task2(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task2() {};
};

template <typename T>
class Task3 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x3 = b(1);
          const Index x1 = b(2);
          const Index x2 = b(3);
          // std::shared_ptr<Tensor<T> > Gamma6;
          // std::shared_ptr<Tensor<T> > rdm2;
          // std::shared_ptr<Tensor<T> > f1;

          // tensor label: Gamma6
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
          // associated with merged
          std::unique_ptr<double[]> fdata = in(1)->get_block(x2, x1);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            for (int i1 = 0; i1 != x1.size(); ++i1) {
              for (int i2 = 0; i2 != x2.size(); ++i2) {
                for (int i0 = 0; i0 != x0.size(); ++i0) {
                  for (int i3 = 0; i3 != x3.size(); ++i3) {
                    odata[i3+x3.size()*(i0)]
                      += (1.0) * i0data[i3+x3.size()*(i0+x0.size()*(i2+x2.size()*(i1)))] * fdata[i2+x2.size()*(i1)];
                  }
                }
              }
            }
          }
          out()->put_block(odata, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task3(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x1 : *range[1])
                for (auto& x2 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task3() {};
};

template <typename T>
class Task4 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x2 = b(1);
          const Index x0 = b(2);
          const Index x3 = b(3);
          // std::shared_ptr<Tensor<T> > Gamma14;
          // std::shared_ptr<Tensor<T> > rdm2;

          // tensor label: Gamma14
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task4(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              for (auto& x0 : *range[1])
                for (auto& x3 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task4() {};
};

template <typename T>
class Task5 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x4 = b(1);
          const Index x0 = b(2);
          const Index x5 = b(3);
          const Index x2 = b(4);
          const Index x3 = b(5);
          // std::shared_ptr<Tensor<T> > Gamma16;
          // std::shared_ptr<Tensor<T> > rdm3;
          // std::shared_ptr<Tensor<T> > f1;

          // tensor label: Gamma16
          std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
          // associated with merged
          std::unique_ptr<double[]> fdata = in(1)->get_block(x3, x2);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
            for (int i2 = 0; i2 != x2.size(); ++i2) {
              for (int i3 = 0; i3 != x3.size(); ++i3) {
                for (int i1 = 0; i1 != x1.size(); ++i1) {
                  for (int i4 = 0; i4 != x4.size(); ++i4) {
                    for (int i0 = 0; i0 != x0.size(); ++i0) {
                      for (int i5 = 0; i5 != x5.size(); ++i5) {
                        odata[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1)))]
                          += (1.0) * i0data[i5+x5.size()*(i0+x0.size()*(i4+x4.size()*(i1+x1.size()*(i3+x3.size()*(i2)))))] * fdata[i3+x3.size()*(i2)];
                      }
                    }
                  }
                }
              }
            }
          }
          out()->put_block(odata, x5, x0, x4, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task5(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
          for (auto& x1 : *range[1])
            for (auto& x4 : *range[1])
              for (auto& x0 : *range[1])
                for (auto& x5 : *range[1])
                  for (auto& x2 : *range[1])
                    for (auto& x3 : *range[1])
                      subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x5, x0, x4, x1, x3, x2}}, in, t[0], range)));
    };
    ~Task5() {};
};

template <typename T>
class Task6 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<6,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x2 = b(0);
          const Index x3 = b(1);
          const Index x1 = b(2);
          const Index x4 = b(3);
          const Index x0 = b(4);
          const Index x5 = b(5);
          // std::shared_ptr<Tensor<T> > Gamma67;
          // std::shared_ptr<Tensor<T> > rdm3;

          // tensor label: Gamma67
          std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1, x3, x2);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
            sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
          }
          out()->put_block(odata, x5, x0, x4, x1, x3, x2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task6(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x1 : *range[1])
                for (auto& x4 : *range[1])
                  for (auto& x0 : *range[1])
                    for (auto& x5 : *range[1])
                      subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x5, x0, x4, x1, x3, x2}}, in, t[0], range)));
    };
    ~Task6() {};
};

template <typename T>
class Task7 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<3,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);
          // std::shared_ptr<Tensor<T> > Gamma72;
          // std::shared_ptr<Tensor<T> > rdm1I0;
          // std::shared_ptr<Tensor<T> > f1;

          // tensor label: Gamma72
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          // associated with merged
          std::unique_ptr<double[]> fdata = in(1)->get_block(x1, x0);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
              for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
                for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                  odata[ici0]
                    += (1.0) * i0data[ici0+ci0.size()*(ix1+x1.size()*(ix0))] * fdata[ix1+x1.size()*(ix0)];
                }
              }
            }
          }
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task7(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,4> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
          for (auto& ci0 : *range[3])
            for (auto& x0 : *range[1])
              for (auto& x1 : *range[1])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task7() {};
};

template <typename T>
class Task8 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);
          // std::shared_ptr<Tensor<T> > Gamma74;
          // std::shared_ptr<Tensor<T> > rdm1I0;

          // tensor label: Gamma74
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task8(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,4> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task8() {};
};

template <typename T>
class Task9 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);
          // std::shared_ptr<Tensor<T> > Gamma78;
          // std::shared_ptr<Tensor<T> > rdm2I0;
          // std::shared_ptr<Tensor<T> > f1;

          // tensor label: Gamma78
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
          // associated with merged
          std::unique_ptr<double[]> fdata = in(1)->get_block(x2, x1);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
              for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
                for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
                  for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                    for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                      odata[ici0+ci0.size()*(ix3+x3.size()*(ix0))]
                        += (1.0) * i0data[ici0+ci0.size()*(ix3+x3.size()*(ix0+x0.size()*(ix2+x2.size()*(ix1))))] * fdata[ix2+x2.size()*(ix1)];
                    }
                  }
                }
              }
            }
          }
          out()->put_block(odata, ci0, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task9(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,4> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                for (auto& x1 : *range[1])
                  for (auto& x2 : *range[1])
                    subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task9() {};
};

template <typename T>
class Task10 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);
          // std::shared_ptr<Tensor<T> > Gamma86;
          // std::shared_ptr<Tensor<T> > rdm2I0;

          // tensor label: Gamma86
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task10(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,4> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              for (auto& x0 : *range[1])
                for (auto& x3 : *range[1])
                  for (auto& ci0 : *range[3])
                    subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task10() {};
};

template <typename T>
class Task11 : public Task<T> {  // associated with gamma
  protected:
    class Task_local : public SubTask<7,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x5 = b(1);
          const Index x0 = b(2);
          const Index x4 = b(3);
          const Index x1 = b(4);
          const Index x3 = b(5);
          const Index x2 = b(6);
          // std::shared_ptr<Tensor<T> > Gamma88;
          // std::shared_ptr<Tensor<T> > rdm3I0;
          // std::shared_ptr<Tensor<T> > f1;

          // tensor label: Gamma88
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
          // associated with merged
          std::unique_ptr<double[]> fdata = in(1)->get_block(x3, x2);
          {
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1, x3, x2);
            for (int ix2 = 0; ix2 != x2.size(); ++ix2) {
              for (int ix3 = 0; ix3 != x3.size(); ++ix3) {
                for (int ix1 = 0; ix1 != x1.size(); ++ix1) {
                  for (int ix4 = 0; ix4 != x4.size(); ++ix4) {
                    for (int ix0 = 0; ix0 != x0.size(); ++ix0) {
                      for (int ix5 = 0; ix5 != x5.size(); ++ix5) {
                        for (int ici0 = 0; ici0 != ci0.size(); ++ici0) {
                          odata[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1))))]
                            += (1.0) * i0data[ici0+ci0.size()*(ix5+x5.size()*(ix0+x0.size()*(ix4+x4.size()*(ix1+x1.size()*(ix3+x3.size()*(ix2))))))] * fdata[ix3+x3.size()*(ix2)];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          out()->put_block(odata, ci0, x5, x0, x4, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task11(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,4> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
          for (auto& x1 : *range[1])
            for (auto& x4 : *range[1])
              for (auto& x0 : *range[1])
                for (auto& x5 : *range[1])
                  for (auto& ci0 : *range[3])
                    for (auto& x2 : *range[1])
                      for (auto& x3 : *range[1])
                        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,7>{{ci0, x5, x0, x4, x1, x3, x2}}, in, t[0], range)));
    };
    ~Task11() {};
};

template <typename T>
class Task12 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c3 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: r
          std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
          {
            // tensor label: I0
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
            sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          out()->put_block(odata, c3, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task12(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c3 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task12() {};
};

template <typename T>
class Task13 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        const double e0_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<4,2,T>(block, in, out), range_(ran), e0_(e) { }

        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I0
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
            dscal_(c1.size()*a4.size()*c3.size()*a2.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
            dscal_(c1.size()*a2.size()*c3.size()*a4.size(), e0_, i1data.get(), 1);
            sort_indices<0,3,2,1,1,1,-8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
          }
          {
            // tensor label: v2
            std::unique_ptr<double[]> i2data = in(1)->get_block(c1, a4, c3, a2);
            sort_indices<0,1,2,3,1,1,-4,1>(i2data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          {
            // tensor label: v2
            std::unique_ptr<double[]> i3data = in(1)->get_block(c1, a2, c3, a4);
            sort_indices<0,3,2,1,1,1,8,1>(i3data, odata, c1.size(), a2.size(), c3.size(), a4.size());
          }
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task13(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range, e)));
    };
    ~Task13() {};
};

template <typename T>
class Task14 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I0
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

          // tensor label: I1
          std::unique_ptr<double[]> i1data = in(1)->get_block();
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
          sort_indices<0,1,1,1>(i1data, i1data_sorted);

          dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), 1, 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

          sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task14(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task14() {};
};

template <typename T>
class Task15 : public Task<T> {
  protected:
    class Task_local : public SubTask<1,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<1,1,T>(std::array<const Index, 1>(), in, out), range_(ran) { }

        void compute() override {

          // tensor label: I1
          std::unique_ptr<double[]> odata = out()->move_block();
          {
            // scalar
            // tensor label: Gamma0
            std::unique_ptr<double[]> i0data = in(0)->get_block();
            sort_indices<1,1,-4,1>(i0data, odata);
          }
          out()->put_block(odata);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task15(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(in, t[0], range)));
    };
    ~Task15() {};
};

template <typename T>
class Task16 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I0
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I3
          std::unique_ptr<double[]> i1data = in(1)->get_block();
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
          sort_indices<0,1,1,1>(i1data, i1data_sorted);

          dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), 1, 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

          sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task16(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task16() {};
};

template <typename T>
class Task17 : public Task<T> {
  protected:
    class Task_local : public SubTask<1,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<1,1,T>(std::array<const Index, 1>(), in, out), range_(ran) { }

        void compute() override {

          // tensor label: I3
          std::unique_ptr<double[]> odata = out()->move_block();
          {
            // scalar
            // tensor label: Gamma0
            std::unique_ptr<double[]> i0data = in(0)->get_block();
            sort_indices<1,1,8,1>(i0data, odata);
          }
          out()->put_block(odata);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task17(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(in, t[0], range)));
    };
    ~Task17() {};
};

template <typename T>
class Task18 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c3 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: r
          std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
          {
            // tensor label: I4
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, a2, c3);
            sort_indices<3,1,0,2,1,1,1,1>(i0data, odata, c1.size(), a4.size(), a2.size(), c3.size());
          }
          {
            // tensor label: I4
            std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, a4, c1);
            sort_indices<0,2,3,1,1,1,1,1>(i0data, odata, c3.size(), a2.size(), a4.size(), c1.size());
          }
          out()->put_block(odata, c3, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task18(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c3 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task18() {};
};

template <typename T>
class Task19 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index a2 = b(2);
          const Index c3 = b(3);

          // tensor label: I4
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, a2, c3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, a2, c3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, a2, c3), 0.0);

          for (auto& c5 : *range_[0]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c5)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c5.size());

            // tensor label: I5
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
            sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());

            dgemm_("T", "N", c3.size(), c1.size()*a4.size()*a2.size(), c5.size(),
                   1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
                   1.0, odata_sorted, c3.size());
          }

          sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a4.size(), a2.size());
          out()->put_block(odata, c1, a4, a2, c3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task19(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, a2, c3}}, in, t[0], range)));
    };
    ~Task19() {};
};

template <typename T>
class Task20 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c5 = b(2);
          const Index a2 = b(3);

          // tensor label: I5
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
            sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c5, a4);
            sort_indices<0,3,2,1,1,1,-8,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
          }
          out()->put_block(odata, c1, a4, c5, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task20(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c5 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c5, a2}}, in, t[0], range)));
    };
    ~Task20() {};
};

template <typename T>
class Task21 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index a2 = b(2);
          const Index c3 = b(3);

          // tensor label: I4
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, a2, c3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, a2, c3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, a2, c3), 0.0);

          for (auto& a5 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a5, a4)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a5.size(), a4.size());

            // tensor label: I9
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, a2);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, a2)]);
            sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());

            dgemm_("T", "N", a4.size(), c1.size()*c3.size()*a2.size(), a5.size(),
                   1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
                   1.0, odata_sorted, a4.size());
          }

          sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
          out()->put_block(odata, c1, a4, a2, c3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task21(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, a2, c3}}, in, t[0], range)));
    };
    ~Task21() {};
};

template <typename T>
class Task22 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c1 = b(0);
          const Index a5 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I9
          std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
            sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a5);
            sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
          }
          out()->put_block(odata, c1, a5, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task22(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a5 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a5, c3, a2}}, in, t[0], range)));
    };
    ~Task22() {};
};

template <typename T>
class Task23 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index a2 = b(2);
          const Index c3 = b(3);

          // tensor label: I4
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, a2, c3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, a2, c3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, a2, c3), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size());

            // tensor label: I13
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c1, a2);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c1, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

            dgemm_("T", "N", c3.size(), a4.size()*c1.size()*a2.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, c3.size());
          }

          sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), a2.size());
          out()->put_block(odata, c1, a4, a2, c3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task23(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, a2, c3}}, in, t[0], range)));
    };
    ~Task23() {};
};

template <typename T>
class Task24 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: I13
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

            // tensor label: I14
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c1.size()*a2.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
          out()->put_block(odata, x0, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task24(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task24() {};
};

template <typename T>
class Task25 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I14
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task25(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task25() {};
};

template <typename T>
class Task26 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: I13
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

            // tensor label: I17
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a2.size()*c1.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task26(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task26() {};
};

template <typename T>
class Task27 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I17
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task27(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task27() {};
};

template <typename T>
class Task28 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index c2 = b(0);
          const Index a3 = b(1);
          const Index x0 = b(2);
          const Index a1 = b(3);

          // tensor label: r
          std::unique_ptr<double[]> odata = out()->move_block(c2, a3, x0, a1);
          {
            // tensor label: I18
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c2, a3, a1);
            sort_indices<1,2,0,3,1,1,1,1>(i0data, odata, x0.size(), c2.size(), a3.size(), a1.size());
          }
          out()->put_block(odata, c2, a3, x0, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task28(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c2, a3, x0, a1}}, in, t[0], range)));
    };
    ~Task28() {};
};

template <typename T>
class Task29 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            for (auto& c4 : *range_[0]) {
              // tensor label: f1
              std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
              sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());

              // tensor label: I19
              std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a3, c4, a1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a3, c4, a1)]);
              sort_indices<0,4,1,2,3,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

              dgemm_("T", "N", 1, x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
                     1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
                     1.0, odata_sorted, 1);
            }
          }

          sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task29(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task29() {};
};

template <typename T>
class Task30 : public Task<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index c4 = b(4);
          const Index a1 = b(5);

          // tensor label: I19
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

          // tensor label: I20
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
          sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

          sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), x1.size(), x0.size());
          out()->put_block(odata, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task30(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task30() {};
};

template <typename T>
class Task31 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I20
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-4,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task31(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task31() {};
};

template <typename T>
class Task32 : public Task<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index c4 = b(4);
          const Index a1 = b(5);

          // tensor label: I19
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

          // tensor label: I23
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
          sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

          sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), x1.size(), x0.size());
          out()->put_block(odata, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task32(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task32() {};
};

template <typename T>
class Task33 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I23
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task33(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task33() {};
};

template <typename T>
class Task34 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I25
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task34(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task34() {};
};

template <typename T>
class Task35 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);

          // tensor label: I25
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
          {
            // tensor label: Gamma6
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x3.size(), x0.size());
          }
          out()->put_block(odata, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task35(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x3, x0}}, in, t[0], range)));
    };
    ~Task35() {};
};

template <typename T>
class Task36 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I27
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task36(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task36() {};
};

template <typename T>
class Task37 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);

          // tensor label: I27
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
          {
            // tensor label: Gamma6
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x3.size(), x0.size());
          }
          out()->put_block(odata, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task37(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x3, x0}}, in, t[0], range)));
    };
    ~Task37() {};
};

template <typename T>
class Task38 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& c4 : *range_[0]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());

            // tensor label: I29
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c4, a1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c4, a1)]);
            sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a1.size());

            dgemm_("T", "N", c2.size(), x0.size()*a3.size()*a1.size(), c4.size(),
                   1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
                   1.0, odata_sorted, c2.size());
          }

          sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a3.size(), a1.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task38(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task38() {};
};

template <typename T>
class Task39 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c4 = b(2);
          const Index a1 = b(3);

          // tensor label: I29
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());

            // tensor label: I30
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c4.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c4.size()*a1.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task39(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task39() {};
};

template <typename T>
class Task40 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I30
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task40(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task40() {};
};

template <typename T>
class Task41 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c4 = b(2);
          const Index a1 = b(3);

          // tensor label: I29
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());

            // tensor label: I33
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c4.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c4.size()*a3.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task41(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task41() {};
};

template <typename T>
class Task42 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I33
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task42(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task42() {};
};

template <typename T>
class Task43 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& a4 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());

            // tensor label: I35
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a3)]);
            sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a3.size());

            dgemm_("T", "N", a1.size(), x0.size()*c2.size()*a3.size(), a4.size(),
                   1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                   1.0, odata_sorted, a1.size());
          }

          sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), c2.size(), a3.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task43(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task43() {};
};

template <typename T>
class Task44 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);

          // tensor label: I35
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());

            // tensor label: I36
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a3.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task44(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task44() {};
};

template <typename T>
class Task45 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I36
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task45(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task45() {};
};

template <typename T>
class Task46 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);

          // tensor label: I35
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());

            // tensor label: I39
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task46(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task46() {};
};

template <typename T>
class Task47 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I39
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task47(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task47() {};
};

template <typename T>
class Task48 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& a4 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

            // tensor label: I41
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a1)]);
            sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a1.size());

            dgemm_("T", "N", a3.size(), x0.size()*c2.size()*a1.size(), a4.size(),
                   1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                   1.0, odata_sorted, a3.size());
          }

          sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), c2.size(), a1.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task48(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task48() {};
};

template <typename T>
class Task49 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I41
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());

            // tensor label: I42
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a1.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task49(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task49() {};
};

template <typename T>
class Task50 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I42
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task50(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task50() {};
};

template <typename T>
class Task51 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I41
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());

            // tensor label: I45
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task51(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task51() {};
};

template <typename T>
class Task52 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I45
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task52(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task52() {};
};

template <typename T>
class Task53 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());

            // tensor label: I47
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
            sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

            dgemm_("T", "N", c2.size(), x0.size()*a1.size()*a3.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, c2.size());
          }

          sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a1.size(), a3.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task53(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task53() {};
};

template <typename T>
class Task54 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a3 = b(3);

          // tensor label: I47
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

              // tensor label: I48
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task54(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range)));
    };
    ~Task54() {};
};

template <typename T>
class Task55 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I48
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task55(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task55() {};
};

template <typename T>
class Task56 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I60
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task56(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task56() {};
};

template <typename T>
class Task57 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        const double e0_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<2,1,T>(block, in, out), range_(ran), e0_(e) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I60
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task57(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range, e)));
    };
    ~Task57() {};
};

template <typename T>
class Task58 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I62
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task58(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task58() {};
};

template <typename T>
class Task59 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        const double e0_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<2,1,T>(block, in, out), range_(ran), e0_(e) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I62
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task59(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range, e)));
    };
    ~Task59() {};
};

template <typename T>
class Task60 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I68
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task60(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task60() {};
};

template <typename T>
class Task61 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I68
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task61(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task61() {};
};

template <typename T>
class Task62 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I18
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I70
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task62(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task62() {};
};

template <typename T>
class Task63 : public Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I70
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,4,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task63(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task63() {};
};

template <typename T>
class Task64 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index a2 = b(1);
          const Index x0 = b(2);
          const Index a1 = b(3);

          // tensor label: r
          std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
          {
            // tensor label: I49
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
            sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
          }
          {
            // tensor label: I49
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, a2, a1);
            sort_indices<0,2,1,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), a2.size(), a1.size());
          }
          out()->put_block(odata, x1, a2, x0, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task64(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& a2 : *range[2])
            for (auto& x1 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range)));
    };
    ~Task64() {};
};

template <typename T>
class Task65 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I49
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& c3 : *range_[0]) {
              // tensor label: f1
              std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
              sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());

              // tensor label: I50
              std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1, c3, a2);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1, c3, a2)]);
              sort_indices<1,4,0,2,3,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

              dgemm_("T", "N", 1, x0.size()*x1.size()*a1.size()*a2.size(), x2.size()*c3.size(),
                     1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
                     1.0, odata_sorted, 1);
            }
          }

          sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task65(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task65() {};
};

template <typename T>
class Task66 : public Task<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index x2 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index c3 = b(4);
          const Index a2 = b(5);

          // tensor label: I50
          std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1, c3, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());

            // tensor label: I51
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

            dgemm_("T", "N", a1.size()*c3.size()*a2.size(), x0.size()*x2.size()*x1.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a1.size()*c3.size()*a2.size());
          }

          sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), x0.size(), x2.size(), x1.size());
          out()->put_block(odata, x0, x2, x1, a1, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task66(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a1 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& x2 : *range[1])
                for (auto& x0 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x0, x2, x1, a1, c3, a2}}, in, t[0], range)));
    };
    ~Task66() {};
};

template <typename T>
class Task67 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I51
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task67(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task67() {};
};

template <typename T>
class Task68 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I49
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& a3 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

            // tensor label: I55
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
            sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

            dgemm_("T", "N", a2.size(), x0.size()*x1.size()*a1.size(), a3.size(),
                   1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
                   1.0, odata_sorted, a2.size());
          }

          sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x1.size(), a1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task68(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task68() {};
};

template <typename T>
class Task69 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a3 = b(3);

          // tensor label: I55
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

              // tensor label: I56
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task69(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range)));
    };
    ~Task69() {};
};

template <typename T>
class Task70 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I56
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task70(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task70() {};
};

template <typename T>
class Task71 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x1 = b(0);
          const Index a2 = b(1);
          const Index x0 = b(2);
          const Index a1 = b(3);

          // tensor label: r
          std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
          {
            // tensor label: I52
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
            sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
          }
          out()->put_block(odata, x1, a2, x0, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task71(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& a2 : *range[2])
            for (auto& x1 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range)));
    };
    ~Task71() {};
};

template <typename T>
class Task72 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I52
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x5 : *range_[1]) {
            for (auto& x4 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

              // tensor label: I53
              std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x5.size()*x4.size(),
                     1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task72(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task72() {};
};

template <typename T>
class Task73 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x5 = b(0);
          const Index x0 = b(1);
          const Index x4 = b(2);
          const Index x1 = b(3);

          // tensor label: I53
          std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
          {
            // tensor label: Gamma16
            std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
            sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size());
          }
          out()->put_block(odata, x5, x0, x4, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task73(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x5, x0, x4, x1}}, in, t[0], range)));
    };
    ~Task73() {};
};

template <typename T>
class Task74 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I52
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

              // tensor label: I64
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task74(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task74() {};
};

template <typename T>
class Task75 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        const double e0_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<4,1,T>(block, in, out), range_(ran), e0_(e) { }

        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I64
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            dscal_(x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task75(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range, e)));
    };
    ~Task75() {};
};

template <typename T>
class Task76 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I52
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: v2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

              // tensor label: I72
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task76(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task76() {};
};

template <typename T>
class Task77 : public Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I72
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) i->compute();
    }

  public:
    Task77(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task77() {};
};

template <typename T>
class Task78 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a2 = b(1);
          const Index c3 = b(2);
          const Index a4 = b(3);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I74
          std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c3, a2)]);
          sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

          energy_ += ddot_(c1.size()*a4.size()*c3.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task78(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a2 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range)));
    };
    ~Task78() {};
};

template <typename T>
class Task79 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;
        double e0_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<4,2,T>(block, in, out), range_(ran), e0_(e) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I74
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
            dscal_(c1.size()*a4.size()*c3.size()*a2.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
            dscal_(c1.size()*a2.size()*c3.size()*a4.size(), e0_, i1data.get(), 1);
            sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
          }
          {
            // tensor label: v2
            std::unique_ptr<double[]> i2data = in(1)->get_block(c1, a4, c3, a2);
            sort_indices<0,1,2,3,1,1,-2,1>(i2data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          {
            // tensor label: v2
            std::unique_ptr<double[]> i3data = in(1)->get_block(c1, a2, c3, a4);
            sort_indices<0,3,2,1,1,1,4,1>(i3data, odata, c1.size(), a2.size(), c3.size(), a4.size());
          }
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task79(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range, e)));
    };
    ~Task79() {};
};

template <typename T>
class Task80 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I74
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

          // tensor label: I75
          std::unique_ptr<double[]> i1data = in(1)->get_block();
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
          sort_indices<0,1,1,1>(i1data, i1data_sorted);

          dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), 1, 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

          sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task80(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task80() {};
};

template <typename T>
class Task81 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<1,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<1,1,T>(std::array<const Index, 1>(), in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;

          // tensor label: I75
          std::unique_ptr<double[]> odata = out()->move_block();
          {
            // scalar
            // tensor label: Gamma0
            std::unique_ptr<double[]> i0data = in(0)->get_block();
            sort_indices<1,1,-1,1>(i0data, odata);
          }
          out()->put_block(odata);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task81(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(in, t[0], range)));
    };
    ~Task81() {};
};

template <typename T>
class Task82 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I74
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I78
          std::unique_ptr<double[]> i1data = in(1)->get_block();
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
          sort_indices<0,1,1,1>(i1data, i1data_sorted);

          dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), 1, 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

          sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task82(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task82() {};
};

template <typename T>
class Task83 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<1,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<1,1,T>(std::array<const Index, 1>(), in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;

          // tensor label: I78
          std::unique_ptr<double[]> odata = out()->move_block();
          {
            // scalar
            // tensor label: Gamma0
            std::unique_ptr<double[]> i0data = in(0)->get_block();
            sort_indices<1,1,2,1>(i0data, odata);
          }
          out()->put_block(odata);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task83(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(in, t[0], range)));
    };
    ~Task83() {};
};

template <typename T>
class Task84 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I74
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

          for (auto& c5 : *range_[0]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c5)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c5.size());

            // tensor label: I81
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
            sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());

            dgemm_("T", "N", c3.size(), c1.size()*a4.size()*a2.size(), c5.size(),
                   1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
                   1.0, odata_sorted, c3.size());
          }

          sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a4.size(), a2.size());
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task84(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task84() {};
};

template <typename T>
class Task85 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c5 = b(2);
          const Index a2 = b(3);

          // tensor label: I81
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
            sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c5, a4);
            sort_indices<0,3,2,1,1,1,-4,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
          }
          out()->put_block(odata, c1, a4, c5, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task85(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c5 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c5, a2}}, in, t[0], range)));
    };
    ~Task85() {};
};

template <typename T>
class Task86 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I74
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

          for (auto& a5 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a5, a4)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a5.size(), a4.size());

            // tensor label: I87
            std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, a2);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, a2)]);
            sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());

            dgemm_("T", "N", a4.size(), c1.size()*c3.size()*a2.size(), a5.size(),
                   1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
                   1.0, odata_sorted, a4.size());
          }

          sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task86(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task86() {};
};

template <typename T>
class Task87 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a5 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I87
          std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
            sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a5);
            sort_indices<0,3,2,1,1,1,4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
          }
          out()->put_block(odata, c1, a5, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task87(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a5 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a5, c3, a2}}, in, t[0], range)));
    };
    ~Task87() {};
};

template <typename T>
class Task88 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I74
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size());

            // tensor label: I93
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c1, a2);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c1, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

            dgemm_("T", "N", c3.size(), a4.size()*c1.size()*a2.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, c3.size());
          }

          sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), a2.size());
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task88(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task88() {};
};

template <typename T>
class Task89 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: I93
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

            // tensor label: I94
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c1.size()*a2.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
          out()->put_block(odata, x0, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task89(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task89() {};
};

template <typename T>
class Task90 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I94
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task90(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task90() {};
};

template <typename T>
class Task91 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: I93
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

            // tensor label: I98
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a2.size()*c1.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task91(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task91() {};
};

template <typename T>
class Task92 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I98
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task92(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task92() {};
};

template <typename T>
class Task93 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
          sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

          // tensor label: I100
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c2, a3, a1)]);
          sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c2.size(), a3.size(), a1.size());

          energy_ += ddot_(x0.size()*c2.size()*a3.size()*a1.size(), i0data_sorted, 1, i1data_sorted, 1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task93(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a1 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a1, c2, a3}}, in, t[0], range)));
    };
    ~Task93() {};
};

template <typename T>
class Task94 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            for (auto& c4 : *range_[0]) {
              // tensor label: f1
              std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
              sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());

              // tensor label: I101
              std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a3, c4, a1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a3, c4, a1)]);
              sort_indices<0,4,1,2,3,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

              dgemm_("T", "N", 1, x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
                     1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
                     1.0, odata_sorted, 1);
            }
          }

          sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task94(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task94() {};
};

template <typename T>
class Task95 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index c4 = b(4);
          const Index a1 = b(5);

          // tensor label: I101
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

          // tensor label: I102
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
          sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

          sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), x1.size(), x0.size());
          out()->put_block(odata, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task95(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task95() {};
};

template <typename T>
class Task96 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I102
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task96(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task96() {};
};

template <typename T>
class Task97 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index c4 = b(4);
          const Index a1 = b(5);

          // tensor label: I101
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

          // tensor label: I106
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
          sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

          sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), x1.size(), x0.size());
          out()->put_block(odata, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task97(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task97() {};
};

template <typename T>
class Task98 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I106
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task98(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task98() {};
};

template <typename T>
class Task99 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I109
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task99(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task99() {};
};

template <typename T>
class Task100 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x3 = b(0);
          const Index x0 = b(1);

          // tensor label: I109
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
          {
            // tensor label: Gamma6
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
            sort_indices<0,1,1,1,-1,4>(i0data, odata, x3.size(), x0.size());
          }
          out()->put_block(odata, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task100(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x3, x0}}, in, t[0], range)));
    };
    ~Task100() {};
};

template <typename T>
class Task101 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I112
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task101(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task101() {};
};

template <typename T>
class Task102 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x3 = b(0);
          const Index x0 = b(1);

          // tensor label: I112
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
          {
            // tensor label: Gamma6
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x3.size(), x0.size());
          }
          out()->put_block(odata, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task102(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x3, x0}}, in, t[0], range)));
    };
    ~Task102() {};
};

template <typename T>
class Task103 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& c4 : *range_[0]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());

            // tensor label: I115
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c4, a1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c4, a1)]);
            sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a1.size());

            dgemm_("T", "N", c2.size(), x0.size()*a3.size()*a1.size(), c4.size(),
                   1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
                   1.0, odata_sorted, c2.size());
          }

          sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a3.size(), a1.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task103(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task103() {};
};

template <typename T>
class Task104 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c4 = b(2);
          const Index a1 = b(3);

          // tensor label: I115
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());

            // tensor label: I116
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c4.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c4.size()*a1.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task104(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task104() {};
};

template <typename T>
class Task105 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I116
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,4>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task105(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task105() {};
};

template <typename T>
class Task106 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c4 = b(2);
          const Index a1 = b(3);

          // tensor label: I115
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());

            // tensor label: I120
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c4.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c4.size()*a3.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task106(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task106() {};
};

template <typename T>
class Task107 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I120
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task107(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task107() {};
};

template <typename T>
class Task108 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& a4 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());

            // tensor label: I123
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a3)]);
            sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a3.size());

            dgemm_("T", "N", a1.size(), x0.size()*c2.size()*a3.size(), a4.size(),
                   1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                   1.0, odata_sorted, a1.size());
          }

          sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), c2.size(), a3.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task108(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task108() {};
};

template <typename T>
class Task109 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);

          // tensor label: I123
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());

            // tensor label: I124
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a3.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task109(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task109() {};
};

template <typename T>
class Task110 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I124
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task110(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task110() {};
};

template <typename T>
class Task111 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);

          // tensor label: I123
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());

            // tensor label: I128
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task111(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task111() {};
};

template <typename T>
class Task112 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I128
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task112(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task112() {};
};

template <typename T>
class Task113 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& a4 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

            // tensor label: I131
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a1)]);
            sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a1.size());

            dgemm_("T", "N", a3.size(), x0.size()*c2.size()*a1.size(), a4.size(),
                   1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                   1.0, odata_sorted, a3.size());
          }

          sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), c2.size(), a1.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task113(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task113() {};
};

template <typename T>
class Task114 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I131
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());

            // tensor label: I132
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a1.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task114(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task114() {};
};

template <typename T>
class Task115 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I132
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task115(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task115() {};
};

template <typename T>
class Task116 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I131
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());

            // tensor label: I136
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task116(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task116() {};
};

template <typename T>
class Task117 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I136
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task117(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task117() {};
};

template <typename T>
class Task118 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());

            // tensor label: I139
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
            sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

            dgemm_("T", "N", c2.size(), x0.size()*a1.size()*a3.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, c2.size());
          }

          sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), a1.size(), a3.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task118(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task118() {};
};

template <typename T>
class Task119 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a3 = b(3);

          // tensor label: I139
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

              // tensor label: I140
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task119(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range)));
    };
    ~Task119() {};
};

template <typename T>
class Task120 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I140
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task120(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task120() {};
};

template <typename T>
class Task121 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I158
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task121(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task121() {};
};

template <typename T>
class Task122 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;
        double e0_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<2,1,T>(block, in, out), range_(ran), e0_(e) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I158
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,1,1,1,4>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task122(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range, e)));
    };
    ~Task122() {};
};

template <typename T>
class Task123 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I161
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task123(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task123() {};
};

template <typename T>
class Task124 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;
        double e0_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<2,1,T>(block, in, out), range_(ran), e0_(e) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I161
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            dscal_(x1.size()*x0.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task124(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range, e)));
    };
    ~Task124() {};
};

template <typename T>
class Task125 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I171
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task125(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task125() {};
};

template <typename T>
class Task126 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I171
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task126(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task126() {};
};

template <typename T>
class Task127 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index c2 = b(1);
          const Index a3 = b(2);
          const Index a1 = b(3);

          // tensor label: I100
          std::unique_ptr<double[]> odata = out()->move_block(x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I174
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task127(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task127() {};
};

template <typename T>
class Task128 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I174
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task128(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task128() {};
};

template <typename T>
class Task129 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index a1 = b(1);
          const Index x1 = b(2);
          const Index a2 = b(3);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
          sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

          // tensor label: I142
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a2)]);
          sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a2.size());

          energy_ += ddot_(x0.size()*x1.size()*a1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task129(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& a1 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a1, x1, a2}}, in, t[0], range)));
    };
    ~Task129() {};
};

template <typename T>
class Task130 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I142
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& c3 : *range_[0]) {
              // tensor label: f1
              std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
              sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());

              // tensor label: I143
              std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1, c3, a2);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1, c3, a2)]);
              sort_indices<1,4,0,2,3,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

              dgemm_("T", "N", 1, x0.size()*x1.size()*a1.size()*a2.size(), x2.size()*c3.size(),
                     1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
                     1.0, odata_sorted, 1);
            }
          }

          sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task130(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task130() {};
};

template <typename T>
class Task131 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index x2 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index c3 = b(4);
          const Index a2 = b(5);

          // tensor label: I143
          std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1, c3, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());

            // tensor label: I144
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

            dgemm_("T", "N", a1.size()*c3.size()*a2.size(), x0.size()*x2.size()*x1.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a1.size()*c3.size()*a2.size());
          }

          sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), x0.size(), x2.size(), x1.size());
          out()->put_block(odata, x0, x2, x1, a1, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task131(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a1 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& x2 : *range[1])
                for (auto& x0 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x0, x2, x1, a1, c3, a2}}, in, t[0], range)));
    };
    ~Task131() {};
};

template <typename T>
class Task132 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I144
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task132(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task132() {};
};

template <typename T>
class Task133 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I142
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x5 : *range_[1]) {
            for (auto& x4 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

              // tensor label: I147
              std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x5.size()*x4.size(),
                     1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task133(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task133() {};
};

template <typename T>
class Task134 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x5 = b(0);
          const Index x0 = b(1);
          const Index x4 = b(2);
          const Index x1 = b(3);

          // tensor label: I147
          std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1);
          {
            // tensor label: Gamma16
            std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
            sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size());
          }
          out()->put_block(odata, x5, x0, x4, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task134(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x5, x0, x4, x1}}, in, t[0], range)));
    };
    ~Task134() {};
};

template <typename T>
class Task135 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I142
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& a3 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

            // tensor label: I150
            std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
            sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

            dgemm_("T", "N", a2.size(), x0.size()*x1.size()*a1.size(), a3.size(),
                   1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
                   1.0, odata_sorted, a2.size());
          }

          sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), x1.size(), a1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task135(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task135() {};
};

template <typename T>
class Task136 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a3 = b(3);

          // tensor label: I150
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

              // tensor label: I151
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task136(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range)));
    };
    ~Task136() {};
};

template <typename T>
class Task137 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I151
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task137(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task137() {};
};

template <typename T>
class Task138 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I142
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

              // tensor label: I164
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task138(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task138() {};
};

template <typename T>
class Task139 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;
        double e0_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran, const double e)
          : SubTask<4,1,T>(block, in, out), range_(ran), e0_(e) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I164
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            dscal_(x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task139(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range, e)));
    };
    ~Task139() {};
};

template <typename T>
class Task140 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I142
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: v2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

              // tensor label: I177
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task140(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task140() {};
};

template <typename T>
class Task141 : public EnergyTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double energy_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double energy() const { return energy_; }

        void compute() override {
          energy_ = 0.0;
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I177
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->energy_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->energy_ += i->energy();
      }
    }

  public:
    Task141(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task141() {};
};

template <typename T>
class Task142 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index c1 = b(0);
          const Index a2 = b(1);
          const Index c3 = b(2);
          const Index a4 = b(3);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I179
          std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c3, a2)]);
          sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

          correction_ += ddot_(c1.size()*a4.size()*c3.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task142(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a2 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range)));
    };
    ~Task142() {};
};

template <typename T>
class Task143 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I179
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
            sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
            sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
          }
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task143(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task143() {};
};

template <typename T>
class Task144 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index x0 = b(0);
          const Index a1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
          sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

          // tensor label: I183
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c2, a1)]);
          sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c2.size(), a1.size());

          correction_ += ddot_(x0.size()*a3.size()*c2.size()*a1.size(), i0data_sorted, 1, i1data_sorted, 1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task144(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a1 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a1, c2, a3}}, in, t[0], range)));
    };
    ~Task144() {};
};

template <typename T>
class Task145 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I183
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I184
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task145(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task145() {};
};

template <typename T>
class Task146 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I184
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task146(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task146() {};
};

template <typename T>
class Task147 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I183
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I187
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task147(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task147() {};
};

template <typename T>
class Task148 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I187
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task148(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task148() {};
};

template <typename T>
class Task149 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index x0 = b(0);
          const Index a1 = b(1);
          const Index x1 = b(2);
          const Index a2 = b(3);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
          sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

          // tensor label: I189
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a2)]);
          sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a2.size());

          correction_ += ddot_(x0.size()*x1.size()*a1.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task149(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& a1 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a1, x1, a2}}, in, t[0], range)));
    };
    ~Task149() {};
};

template <typename T>
class Task150 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I189
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
              sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

              // tensor label: I190
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task150(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task150() {};
};

template <typename T>
class Task151 : public CorrectionTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }
        double correction_;

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }

        double correction() const { return correction_; }

        void compute() override {
          correction_ = 0.0;
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I190
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      this->correction_ = 0.0;
      for (auto& i : subtasks_) {
        i->compute();
        this->correction_ += i->correction();
      }
    }

  public:
    Task151(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task151() {};
};

template <typename T>
class Task152 : public DensityTask<T> {
  protected:
    std::shared_ptr<Tensor<T>> d_;
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;

    void compute_() {
      d_->zero();
    };

  public:
    Task152(std::vector<std::shared_ptr<Tensor<T>>> t) : DensityTask<T>() {
      d_ =  t[0];
    };
    ~Task152() {};
};

template <typename T>
class Task153 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
          {
            // tensor label: I191
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x0, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task153(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x0, x1}}, in, t[0], range)));
    };
    ~Task153() {};
};

template <typename T>
class Task154 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I191
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);

          for (auto& c1 : *range_[0]) {
            for (auto& a2 : *range_[2]) {
              for (auto& c3 : *range_[0]) {
                for (auto& a4 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

                  // tensor label: I192
                  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, a4, c3, a2);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, a4, c3, a2)]);
                  sort_indices<2,5,4,3,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), a4.size(), c3.size(), a2.size());

                  dgemm_("T", "N", 1, x1.size()*x0.size(), c1.size()*a4.size()*c3.size()*a2.size(),
                         1.0, i0data_sorted, c1.size()*a4.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*c3.size()*a2.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task154(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task154() {};
};

template <typename T>
class Task155 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);
          const Index c1 = b(2);
          const Index a4 = b(3);
          const Index c3 = b(4);
          const Index a2 = b(5);

          // tensor label: I192
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

          // tensor label: I193
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
          sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

          dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

          sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), x1.size(), x0.size());
          out()->put_block(odata, x1, x0, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task155(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x1, x0, c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task155() {};
};

template <typename T>
class Task156 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I193
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task156(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task156() {};
};

template <typename T>
class Task157 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);
          const Index c1 = b(2);
          const Index a4 = b(3);
          const Index c3 = b(4);
          const Index a2 = b(5);

          // tensor label: I192
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I196
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
          sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

          dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

          sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), x1.size(), x0.size());
          out()->put_block(odata, x1, x0, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task157(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x1, x0, c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task157() {};
};

template <typename T>
class Task158 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I196
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task158(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task158() {};
};

template <typename T>
class Task159 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c5 = b(0);
          const Index c3 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
          {
            // tensor label: I197
            std::unique_ptr<double[]> i0data = in(0)->get_block(c5, c3);
            sort_indices<0,1,1,1,1,1>(i0data, odata, c5.size(), c3.size());
          }
          out()->put_block(odata, c5, c3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task159(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& c5 : *range[0])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{c5, c3}}, in, t[0], range)));
    };
    ~Task159() {};
};

template <typename T>
class Task160 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c5 = b(0);
          const Index c3 = b(1);

          // tensor label: I197
          std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c5, c3), 0.0);

          for (auto& c1 : *range_[0]) {
            for (auto& a2 : *range_[2]) {
              for (auto& a4 : *range_[2]) {
                // tensor label: t2
                std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
                std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
                sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

                // tensor label: I198
                std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
                std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
                sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());

                dgemm_("T", "N", c3.size(), c5.size(), c1.size()*a4.size()*a2.size(),
                       1.0, i0data_sorted, c1.size()*a4.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*a2.size(),
                       1.0, odata_sorted, c3.size());
              }
            }
          }

          sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
          out()->put_block(odata, c5, c3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task160(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& c5 : *range[0])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{c5, c3}}, in, t[0], range)));
    };
    ~Task160() {};
};

template <typename T>
class Task161 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c5 = b(2);
          const Index a2 = b(3);

          // tensor label: I198
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
            sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c5, a4);
            sort_indices<0,3,2,1,1,1,-4,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
          }
          out()->put_block(odata, c1, a4, c5, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task161(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c5 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c5, a2}}, in, t[0], range)));
    };
    ~Task161() {};
};

template <typename T>
class Task162 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index a4 = b(0);
          const Index a5 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(a4, a5);
          {
            // tensor label: I201
            std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
            sort_indices<1,0,1,1,1,1>(i0data, odata, a5.size(), a4.size());
          }
          out()->put_block(odata, a4, a5);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task162(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a5 : *range[2])
        for (auto& a4 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a4, a5}}, in, t[0], range)));
    };
    ~Task162() {};
};

template <typename T>
class Task163 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index a5 = b(0);
          const Index a4 = b(1);

          // tensor label: I201
          std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
          std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);

          for (auto& c1 : *range_[0]) {
            for (auto& a2 : *range_[2]) {
              for (auto& c3 : *range_[0]) {
                // tensor label: t2
                std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
                std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

                // tensor label: I202
                std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, a2);
                std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, a2)]);
                sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());

                dgemm_("T", "N", a4.size(), a5.size(), c1.size()*c3.size()*a2.size(),
                       1.0, i0data_sorted, c1.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*c3.size()*a2.size(),
                       1.0, odata_sorted, a4.size());
              }
            }
          }

          sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a5.size());
          out()->put_block(odata, a5, a4);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task163(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a4 : *range[2])
        for (auto& a5 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a5, a4}}, in, t[0], range)));
    };
    ~Task163() {};
};

template <typename T>
class Task164 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c1 = b(0);
          const Index a5 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I202
          std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
            sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a5);
            sort_indices<0,3,2,1,1,1,4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
          }
          out()->put_block(odata, c1, a5, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task164(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a5 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a5, c3, a2}}, in, t[0], range)));
    };
    ~Task164() {};
};

template <typename T>
class Task165 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index c3 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
          {
            // tensor label: I205
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c3);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x0.size(), c3.size());
          }
          out()->put_block(odata, x0, c3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task165(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x0, c3}}, in, t[0], range)));
    };
    ~Task165() {};
};

template <typename T>
class Task166 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index c3 = b(1);

          // tensor label: I205
          std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, c3), 0.0);

          for (auto& c1 : *range_[0]) {
            for (auto& a2 : *range_[2]) {
              for (auto& a4 : *range_[2]) {
                // tensor label: t2
                std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
                std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
                sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

                // tensor label: I206
                std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c1, a2);
                std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c1, a2)]);
                sort_indices<2,3,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

                dgemm_("T", "N", c3.size(), x0.size(), a4.size()*c1.size()*a2.size(),
                       1.0, i0data_sorted, a4.size()*c1.size()*a2.size(), i1data_sorted, a4.size()*c1.size()*a2.size(),
                       1.0, odata_sorted, c3.size());
              }
            }
          }

          sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size());
          out()->put_block(odata, x0, c3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task166(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x0, c3}}, in, t[0], range)));
    };
    ~Task166() {};
};

template <typename T>
class Task167 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: I206
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

            // tensor label: I207
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c1.size()*a2.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
          out()->put_block(odata, x0, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task167(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task167() {};
};

template <typename T>
class Task168 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I207
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task168(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task168() {};
};

template <typename T>
class Task169 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: I206
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

            // tensor label: I210
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a2.size()*c1.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task169(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task169() {};
};

template <typename T>
class Task170 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I210
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task170(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task170() {};
};

template <typename T>
class Task171 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c4 = b(0);
          const Index x1 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(c4, x1);
          {
            // tensor label: I211
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
            sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), c4.size());
          }
          out()->put_block(odata, c4, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task171(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
      for (auto& x1 : *range[1])
        for (auto& c4 : *range[0])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{c4, x1}}, in, t[0], range)));
    };
    ~Task171() {};
};

template <typename T>
class Task172 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index c4 = b(1);

          // tensor label: I211
          std::unique_ptr<double[]> odata = out()->move_block(x1, c4);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c4)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, c4), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& c2 : *range_[0]) {
                for (auto& a3 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

                  // tensor label: I212
                  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a3, c4, a1);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a3, c4, a1)]);
                  sort_indices<1,5,2,3,0,4,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

                  dgemm_("T", "N", 1, x1.size()*c4.size(), x0.size()*c2.size()*a3.size()*a1.size(),
                         1.0, i0data_sorted, x0.size()*c2.size()*a3.size()*a1.size(), i1data_sorted, x0.size()*c2.size()*a3.size()*a1.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), c4.size());
          out()->put_block(odata, x1, c4);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task172(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
      for (auto& c4 : *range[0])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, c4}}, in, t[0], range)));
    };
    ~Task172() {};
};

template <typename T>
class Task173 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index c4 = b(4);
          const Index a1 = b(5);

          // tensor label: I212
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

          // tensor label: I213
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
          sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

          sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), x1.size(), x0.size());
          out()->put_block(odata, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task173(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task173() {};
};

template <typename T>
class Task174 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I213
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task174(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task174() {};
};

template <typename T>
class Task175 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index c4 = b(4);
          const Index a1 = b(5);

          // tensor label: I212
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

          // tensor label: I216
          std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
          sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

          sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), x1.size(), x0.size());
          out()->put_block(odata, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task175(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task175() {};
};

template <typename T>
class Task176 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I216
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task176(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task176() {};
};

template <typename T>
class Task177 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x2 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(x1, x2);
          {
            // tensor label: I217
            std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1);
            sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), x1.size());
          }
          out()->put_block(odata, x1, x2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task177(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x2}}, in, t[0], range)));
    };
    ~Task177() {};
};

template <typename T>
class Task178 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x2 = b(0);
          const Index x1 = b(1);

          // tensor label: I217
          std::unique_ptr<double[]> odata = out()->move_block(x2, x1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& c2 : *range_[0]) {
                for (auto& a3 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

                  // tensor label: I218
                  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a3, c2, a1);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a3, c2, a1)]);
                  sort_indices<0,5,4,3,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a3.size(), c2.size(), a1.size());

                  dgemm_("T", "N", 1, x2.size()*x1.size(), x0.size()*a3.size()*c2.size()*a1.size(),
                         1.0, i0data_sorted, x0.size()*a3.size()*c2.size()*a1.size(), i1data_sorted, x0.size()*a3.size()*c2.size()*a1.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
          out()->put_block(odata, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task178(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x2, x1}}, in, t[0], range)));
    };
    ~Task178() {};
};

template <typename T>
class Task179 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x2 = b(1);
          const Index x1 = b(2);
          const Index a3 = b(3);
          const Index c2 = b(4);
          const Index a1 = b(5);

          // tensor label: I218
          std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a3, c2, a1), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I219
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size(), x2.size(), x1.size());
          out()->put_block(odata, x0, x2, x1, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task179(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& x2 : *range[1])
                for (auto& x0 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x0, x2, x1, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task179() {};
};

template <typename T>
class Task180 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I219
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task180(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task180() {};
};

template <typename T>
class Task181 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x2 = b(1);
          const Index x1 = b(2);
          const Index a3 = b(3);
          const Index c2 = b(4);
          const Index a1 = b(5);

          // tensor label: I218
          std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a3, c2, a1), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I222
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size()*x2.size()*x1.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size(), x2.size(), x1.size());
          out()->put_block(odata, x0, x2, x1, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task181(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& x2 : *range[1])
                for (auto& x0 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x0, x2, x1, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task181() {};
};

template <typename T>
class Task182 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I222
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task182(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task182() {};
};

template <typename T>
class Task183 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c4 = b(0);
          const Index c2 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(c4, c2);
          {
            // tensor label: I223
            std::unique_ptr<double[]> i0data = in(0)->get_block(c4, c2);
            sort_indices<0,1,1,1,1,1>(i0data, odata, c4.size(), c2.size());
          }
          out()->put_block(odata, c4, c2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task183(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
      for (auto& c2 : *range[0])
        for (auto& c4 : *range[0])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{c4, c2}}, in, t[0], range)));
    };
    ~Task183() {};
};

template <typename T>
class Task184 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c4 = b(0);
          const Index c2 = b(1);

          // tensor label: I223
          std::unique_ptr<double[]> odata = out()->move_block(c4, c2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, c2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(c4, c2), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& a3 : *range_[2]) {
                // tensor label: t2
                std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
                std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
                sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

                // tensor label: I224
                std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3, c4, a1);
                std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3, c4, a1)]);
                sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size(), c4.size(), a1.size());

                dgemm_("T", "N", c2.size(), c4.size(), x0.size()*a3.size()*a1.size(),
                       1.0, i0data_sorted, x0.size()*a3.size()*a1.size(), i1data_sorted, x0.size()*a3.size()*a1.size(),
                       1.0, odata_sorted, c2.size());
              }
            }
          }

          sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), c4.size());
          out()->put_block(odata, c4, c2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task184(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
      for (auto& c2 : *range[0])
        for (auto& c4 : *range[0])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{c4, c2}}, in, t[0], range)));
    };
    ~Task184() {};
};

template <typename T>
class Task185 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c4 = b(2);
          const Index a1 = b(3);

          // tensor label: I224
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());

            // tensor label: I225
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c4.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c4.size()*a1.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task185(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task185() {};
};

template <typename T>
class Task186 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I225
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,4>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task186(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task186() {};
};

template <typename T>
class Task187 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c4 = b(2);
          const Index a1 = b(3);

          // tensor label: I224
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c4, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());

            // tensor label: I228
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c4.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c4.size()*a3.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task187(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task187() {};
};

template <typename T>
class Task188 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I228
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task188(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task188() {};
};

template <typename T>
class Task189 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index a1 = b(0);
          const Index a4 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(a1, a4);
          {
            // tensor label: I229
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
            sort_indices<1,0,1,1,1,1>(i0data, odata, a4.size(), a1.size());
          }
          out()->put_block(odata, a1, a4);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task189(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a4 : *range[2])
        for (auto& a1 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a1, a4}}, in, t[0], range)));
    };
    ~Task189() {};
};

template <typename T>
class Task190 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index a4 = b(0);
          const Index a1 = b(1);

          // tensor label: I229
          std::unique_ptr<double[]> odata = out()->move_block(a4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(a4, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& c2 : *range_[0]) {
              for (auto& a3 : *range_[2]) {
                // tensor label: t2
                std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
                std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
                sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

                // tensor label: I230
                std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a3);
                std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a3)]);
                sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a3.size());

                dgemm_("T", "N", a1.size(), a4.size(), x0.size()*c2.size()*a3.size(),
                       1.0, i0data_sorted, x0.size()*c2.size()*a3.size(), i1data_sorted, x0.size()*c2.size()*a3.size(),
                       1.0, odata_sorted, a1.size());
              }
            }
          }

          sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), a4.size());
          out()->put_block(odata, a4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task190(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a4 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a4, a1}}, in, t[0], range)));
    };
    ~Task190() {};
};

template <typename T>
class Task191 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);

          // tensor label: I230
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());

            // tensor label: I231
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a3.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task191(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task191() {};
};

template <typename T>
class Task192 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I231
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task192(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task192() {};
};

template <typename T>
class Task193 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);

          // tensor label: I230
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());

            // tensor label: I234
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task193(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task193() {};
};

template <typename T>
class Task194 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I234
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task194(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task194() {};
};

template <typename T>
class Task195 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index a3 = b(0);
          const Index a4 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(a3, a4);
          {
            // tensor label: I235
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
            sort_indices<1,0,1,1,1,1>(i0data, odata, a4.size(), a3.size());
          }
          out()->put_block(odata, a3, a4);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task195(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a4 : *range[2])
        for (auto& a3 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a3, a4}}, in, t[0], range)));
    };
    ~Task195() {};
};

template <typename T>
class Task196 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index a4 = b(0);
          const Index a3 = b(1);

          // tensor label: I235
          std::unique_ptr<double[]> odata = out()->move_block(a4, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(a4, a3), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& c2 : *range_[0]) {
                // tensor label: t2
                std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
                std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

                // tensor label: I236
                std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c2, a1);
                std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c2, a1)]);
                sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c2.size(), a1.size());

                dgemm_("T", "N", a3.size(), a4.size(), x0.size()*c2.size()*a1.size(),
                       1.0, i0data_sorted, x0.size()*c2.size()*a1.size(), i1data_sorted, x0.size()*c2.size()*a1.size(),
                       1.0, odata_sorted, a3.size());
              }
            }
          }

          sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), a4.size());
          out()->put_block(odata, a4, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task196(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a4 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a4, a3}}, in, t[0], range)));
    };
    ~Task196() {};
};

template <typename T>
class Task197 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I236
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());

            // tensor label: I237
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a1.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task197(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task197() {};
};

template <typename T>
class Task198 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I237
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,4>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task198(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task198() {};
};

template <typename T>
class Task199 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a4 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I236
          std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());

            // tensor label: I240
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a4.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a4.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), x0.size());
          out()->put_block(odata, x0, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task199(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task199() {};
};

template <typename T>
class Task200 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I240
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task200(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task200() {};
};

template <typename T>
class Task201 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index c2 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(x1, c2);
          {
            // tensor label: I241
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c2);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c2.size());
          }
          out()->put_block(odata, x1, c2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task201(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, c2}}, in, t[0], range)));
    };
    ~Task201() {};
};

template <typename T>
class Task202 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index c2 = b(1);

          // tensor label: I241
          std::unique_ptr<double[]> odata = out()->move_block(x1, c2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, c2), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& a3 : *range_[2]) {
                // tensor label: t2
                std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
                std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
                sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

                // tensor label: I242
                std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
                std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
                sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

                dgemm_("T", "N", c2.size(), x1.size(), x0.size()*a1.size()*a3.size(),
                       1.0, i0data_sorted, x0.size()*a1.size()*a3.size(), i1data_sorted, x0.size()*a1.size()*a3.size(),
                       1.0, odata_sorted, c2.size());
              }
            }
          }

          sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c2.size(), x1.size());
          out()->put_block(odata, x1, c2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task202(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, c2}}, in, t[0], range)));
    };
    ~Task202() {};
};

template <typename T>
class Task203 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a3 = b(3);

          // tensor label: I242
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& x3 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

              // tensor label: I243
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task203(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range)));
    };
    ~Task203() {};
};

template <typename T>
class Task204 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I243
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task204(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task204() {};
};

template <typename T>
class Task205 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c3 = b(0);
          const Index x2 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(c3, x2);
          {
            // tensor label: I244
            std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
            sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), c3.size());
          }
          out()->put_block(odata, c3, x2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task205(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
      for (auto& x2 : *range[1])
        for (auto& c3 : *range[0])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{c3, x2}}, in, t[0], range)));
    };
    ~Task205() {};
};

template <typename T>
class Task206 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x2 = b(0);
          const Index c3 = b(1);

          // tensor label: I244
          std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x2, c3), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& x1 : *range_[1]) {
                for (auto& a2 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

                  // tensor label: I245
                  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a1, c3, a2);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a1, c3, a2)]);
                  sort_indices<0,3,2,5,1,4,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

                  dgemm_("T", "N", 1, x2.size()*c3.size(), x0.size()*x1.size()*a1.size()*a2.size(),
                         1.0, i0data_sorted, x0.size()*x1.size()*a1.size()*a2.size(), i1data_sorted, x0.size()*x1.size()*a1.size()*a2.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), c3.size());
          out()->put_block(odata, x2, c3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task206(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
      for (auto& c3 : *range[0])
        for (auto& x2 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x2, c3}}, in, t[0], range)));
    };
    ~Task206() {};
};

template <typename T>
class Task207 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x2 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index c3 = b(4);
          const Index a2 = b(5);

          // tensor label: I245
          std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1, c3, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());

            // tensor label: I246
            std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

            dgemm_("T", "N", a1.size()*c3.size()*a2.size(), x0.size()*x2.size()*x1.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a1.size()*c3.size()*a2.size());
          }

          sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), x0.size(), x2.size(), x1.size());
          out()->put_block(odata, x0, x2, x1, a1, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task207(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a1 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& x2 : *range[1])
                for (auto& x0 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x0, x2, x1, a1, c3, a2}}, in, t[0], range)));
    };
    ~Task207() {};
};

template <typename T>
class Task208 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I246
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task208(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task208() {};
};

template <typename T>
class Task209 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x2 = b(0);
          const Index x3 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(x2, x3);
          {
            // tensor label: I247
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2);
            sort_indices<1,0,1,1,1,1>(i0data, odata, x3.size(), x2.size());
          }
          out()->put_block(odata, x2, x3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task209(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x2, x3}}, in, t[0], range)));
    };
    ~Task209() {};
};

template <typename T>
class Task210 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x3 = b(0);
          const Index x2 = b(1);

          // tensor label: I247
          std::unique_ptr<double[]> odata = out()->move_block(x3, x2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x3, x2), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& x1 : *range_[1]) {
                for (auto& a2 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

                  // tensor label: I248
                  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x3, x2, a1, a2);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x3, x2, a1, a2)]);
                  sort_indices<0,4,1,5,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x3.size(), x2.size(), a1.size(), a2.size());

                  dgemm_("T", "N", 1, x3.size()*x2.size(), x0.size()*x1.size()*a1.size()*a2.size(),
                         1.0, i0data_sorted, x0.size()*x1.size()*a1.size()*a2.size(), i1data_sorted, x0.size()*x1.size()*a1.size()*a2.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size());
          out()->put_block(odata, x3, x2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task210(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x3, x2}}, in, t[0], range)));
    };
    ~Task210() {};
};

template <typename T>
class Task211 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index x3 = b(2);
          const Index x2 = b(3);
          const Index a1 = b(4);
          const Index a2 = b(5);

          // tensor label: I248
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, x3, x2, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x3, x2, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x3, x2, a1, a2), 0.0);

          for (auto& x4 : *range_[1]) {
            for (auto& x5 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

              // tensor label: I249
              std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x1, x3, x2);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x1, x3, x2)]);
              sort_indices<2,0,1,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size()*x3.size()*x2.size(), x5.size()*x4.size(),
                     1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size(), x3.size(), x2.size());
          out()->put_block(odata, x0, x1, x3, x2, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task211(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x1 : *range[1])
                for (auto& x0 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x0, x1, x3, x2, a1, a2}}, in, t[0], range)));
    };
    ~Task211() {};
};

template <typename T>
class Task212 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<6,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,6>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<6,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x5 = b(0);
          const Index x0 = b(1);
          const Index x4 = b(2);
          const Index x1 = b(3);
          const Index x3 = b(4);
          const Index x2 = b(5);

          // tensor label: I249
          std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x1, x3, x2);
          {
            // tensor label: Gamma67
            std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1, x3, x2);
            sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x0.size(), x4.size(), x1.size(), x3.size(), x2.size());
          }
          out()->put_block(odata, x5, x0, x4, x1, x3, x2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task212(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x4 : *range[1])
              for (auto& x0 : *range[1])
                for (auto& x5 : *range[1])
                  subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,6>{{x5, x0, x4, x1, x3, x2}}, in, t[0], range)));
    };
    ~Task212() {};
};

template <typename T>
class Task213 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index a2 = b(0);
          const Index a3 = b(1);

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(a2, a3);
          {
            // tensor label: I250
            std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
            sort_indices<1,0,1,1,1,1>(i0data, odata, a3.size(), a2.size());
          }
          out()->put_block(odata, a2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task213(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a2 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a2, a3}}, in, t[0], range)));
    };
    ~Task213() {};
};

template <typename T>
class Task214 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<2,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index a3 = b(0);
          const Index a2 = b(1);

          // tensor label: I250
          std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(a3, a2), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& x1 : *range_[1]) {
                // tensor label: t2
                std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
                std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
                sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

                // tensor label: I251
                std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, a3);
                std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, a3)]);
                sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), a3.size());

                dgemm_("T", "N", a2.size(), a3.size(), x0.size()*x1.size()*a1.size(),
                       1.0, i0data_sorted, x0.size()*x1.size()*a1.size(), i1data_sorted, x0.size()*x1.size()*a1.size(),
                       1.0, odata_sorted, a2.size());
              }
            }
          }

          sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), a3.size());
          out()->put_block(odata, a3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task214(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a3 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a3, a2}}, in, t[0], range)));
    };
    ~Task214() {};
};

template <typename T>
class Task215 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a3 = b(3);

          // tensor label: I251
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a3), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& x3 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

              // tensor label: I252
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task215(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range)));
    };
    ~Task215() {};
};

template <typename T>
class Task216 : public DensityTask<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I252
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task216(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task216() {};
};

template <typename T>
class Task217 : public Density2Task<T> {
  protected:
    std::shared_ptr<Tensor<T>> d2_;
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;

    void compute_() {
      d2_->zero();
    };

  public:
    Task217(std::vector<std::shared_ptr<Tensor<T>>> t) : Density2Task<T>() {
      d2_ =  t[0];
    };
    ~Task217() {};
};

template <typename T>
class Task218 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c3 = b(0);
          const Index a4 = b(1);
          const Index c1 = b(2);
          const Index a2 = b(3);

          // tensor label: Den1
          std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
          {
            // tensor label: I253
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
            sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          out()->put_block(odata, c3, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task218(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c3 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task218() {};
};

template <typename T>
class Task219 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c1 = b(0);
          const Index a4 = b(1);
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I253
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
            sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
            sort_indices<0,3,2,1,1,1,4,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
          }
          out()->put_block(odata, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task219(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task219() {};
};

template <typename T>
class Task220 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index c2 = b(0);
          const Index a3 = b(1);
          const Index x0 = b(2);
          const Index a1 = b(3);

          // tensor label: Den1
          std::unique_ptr<double[]> odata = out()->move_block(c2, a3, x0, a1);
          {
            // tensor label: I255
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
            sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, x0.size(), a3.size(), c2.size(), a1.size());
          }
          out()->put_block(odata, c2, a3, x0, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task220(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c2, a3, x0, a1}}, in, t[0], range)));
    };
    ~Task220() {};
};

template <typename T>
class Task221 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I255
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I256
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
          out()->put_block(odata, x0, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task221(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task221() {};
};

template <typename T>
class Task222 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I256
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task222(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task222() {};
};

template <typename T>
class Task223 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index a3 = b(1);
          const Index c2 = b(2);
          const Index a1 = b(3);

          // tensor label: I255
          std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I258
            std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
            sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
          out()->put_block(odata, x0, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task223(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task223() {};
};

template <typename T>
class Task224 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<2,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,2>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<2,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index x0 = b(1);

          // tensor label: I258
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
          {
            // tensor label: Gamma2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
            sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
          }
          out()->put_block(odata, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task224(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task224() {};
};

template <typename T>
class Task225 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x1 = b(0);
          const Index a2 = b(1);
          const Index x0 = b(2);
          const Index a1 = b(3);

          // tensor label: Den1
          std::unique_ptr<double[]> odata = out()->move_block(x1, a2, x0, a1);
          {
            // tensor label: I259
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, a1, a2);
            sort_indices<1,3,0,2,1,1,1,1>(i0data, odata, x0.size(), x1.size(), a1.size(), a2.size());
          }
          out()->put_block(odata, x1, a2, x0, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task225(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& a2 : *range[2])
            for (auto& x1 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range)));
    };
    ~Task225() {};
};

template <typename T>
class Task226 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<4,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x0 = b(0);
          const Index x1 = b(1);
          const Index a1 = b(2);
          const Index a2 = b(3);

          // tensor label: I259
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, a2), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& x3 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

              // tensor label: I260
              std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
              sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), x0.size(), x1.size());
          out()->put_block(odata, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task226(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task226() {};
};

template <typename T>
class Task227 : public Density2Task<T> {
  protected:
    class Task_local : public SubTask<4,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,3> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,4>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,3>& ran)
          : SubTask<4,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index x3 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);

          // tensor label: I260
          std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
          {
            // tensor label: Gamma14
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
            sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task227(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task227() {};
};

template <typename T>
class Task228 : public DedciTask<T> {
  protected:
    std::shared_ptr<Tensor<T>> dec_;
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;
    IndexRange ci_;

    void compute_() {
      dec_->zero();
    };

  public:
    Task228(std::vector<std::shared_ptr<Tensor<T>>> t) : DedciTask<T>() {
      dec_ =  t[0];
    };
    ~Task228() {};
};

template <typename T>
class Task229 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: deci
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          {
            // tensor label: I261
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
            sort_indices<0,1,1,1,1>(i0data, odata, ci0.size());
          }
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task229(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task229() {};
};

template <typename T>
class Task230 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I261
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

          for (auto& c1 : *range_[0]) {
            for (auto& a2 : *range_[2]) {
              for (auto& c3 : *range_[0]) {
                for (auto& a4 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

                  // tensor label: I262
                  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, c1, a4, c3, a2);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, c1, a4, c3, a2)]);
                  sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), c1.size(), a4.size(), c3.size(), a2.size());

                  dgemm_("T", "N", 1, ci0.size(), c1.size()*a4.size()*c3.size()*a2.size(),
                         1.0, i0data_sorted, c1.size()*a4.size()*c3.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*c3.size()*a2.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task230(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task230() {};
};

template <typename T>
class Task231 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index c1 = b(1);
          const Index a4 = b(2);
          const Index c3 = b(3);
          const Index a2 = b(4);

          // tensor label: I262
          std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

          // tensor label: I263
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0)]);
          sort_indices<0,0,1,1,1>(i1data, i1data_sorted, ci0.size());

          dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), ci0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());

          sort_indices<4,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), ci0.size());
          out()->put_block(odata, ci0, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task231(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task231() {};
};

template <typename T>
class Task232 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I263
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          {
            // tensor label: Gamma72
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
            sort_indices<0,1,1,-1,1>(i0data, odata, ci0.size());
          }
          {
            // tensor label: Gamma72
            std::unique_ptr<double[]> i1data = in(0)->get_block(ci0);
            sort_indices<0,1,1,-1,1>(i1data, odata, ci0.size());
          }
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task232(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task232() {};
};

template <typename T>
class Task233 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index c1 = b(1);
          const Index a4 = b(2);
          const Index c3 = b(3);
          const Index a2 = b(4);

          // tensor label: I262
          std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I266
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0)]);
          sort_indices<0,0,1,1,1>(i1data, i1data_sorted, ci0.size());

          dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), ci0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());

          sort_indices<4,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), ci0.size());
          out()->put_block(odata, ci0, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task233(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task233() {};
};

template <typename T>
class Task234 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I266
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          {
            // tensor label: Gamma72
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
            sort_indices<0,1,1,2,1>(i0data, odata, ci0.size());
          }
          {
            // tensor label: Gamma72
            std::unique_ptr<double[]> i1data = in(0)->get_block(ci0);
            sort_indices<0,1,1,2,1>(i1data, odata, ci0.size());
          }
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task234(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task234() {};
};

template <typename T>
class Task235 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index c1 = b(1);
          const Index a4 = b(2);
          const Index c3 = b(3);
          const Index a2 = b(4);

          // tensor label: I262
          std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size());

            // tensor label: I269
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c1, a2);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c1, a2)]);
            sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c1.size(), a2.size());

            dgemm_("T", "N", c3.size(), ci0.size()*a4.size()*c1.size()*a2.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, c3.size());
          }

          sort_indices<1,3,2,0,4,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), a4.size(), c1.size(), a2.size());
          out()->put_block(odata, ci0, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task235(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task235() {};
};

template <typename T>
class Task236 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index a4 = b(2);
          const Index c1 = b(3);
          const Index a2 = b(4);

          // tensor label: I269
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());

            // tensor label: I270
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c1.size()*a2.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c1.size()*a2.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task236(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task236() {};
};

template <typename T>
class Task237 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I270
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task237(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task237() {};
};

template <typename T>
class Task238 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index a4 = b(2);
          const Index c1 = b(3);
          const Index a2 = b(4);

          // tensor label: I269
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());

            // tensor label: I274
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a2.size()*c1.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a2.size()*c1.size()*a4.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task238(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task238() {};
};

template <typename T>
class Task239 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I274
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task239(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task239() {};
};

template <typename T>
class Task240 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index c1 = b(1);
          const Index a4 = b(2);
          const Index c3 = b(3);
          const Index a2 = b(4);

          // tensor label: I262
          std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());

            // tensor label: I336
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c1, a2);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c1, a2)]);
            sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c1.size(), a2.size());

            dgemm_("T", "N", c3.size(), ci0.size()*a4.size()*c1.size()*a2.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, c3.size());
          }

          sort_indices<1,3,2,0,4,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), a4.size(), c1.size(), a2.size());
          out()->put_block(odata, ci0, c1, a4, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task240(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task240() {};
};

template <typename T>
class Task241 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index a4 = b(2);
          const Index c1 = b(3);
          const Index a2 = b(4);

          // tensor label: I336
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c1, a2), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c1, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c1, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c1.size(), a2.size());

            // tensor label: I337
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c1.size()*a2.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a4.size()*c1.size()*a2.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task241(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task241() {};
};

template <typename T>
class Task242 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I337
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task242(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task242() {};
};

template <typename T>
class Task243 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index a4 = b(2);
          const Index c1 = b(3);
          const Index a2 = b(4);

          // tensor label: I336
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c1, a2), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), a4.size());

            // tensor label: I341
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a2.size()*c1.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a2.size()*c1.size()*a4.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, a4, c1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task243(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task243() {};
};

template <typename T>
class Task244 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I341
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task244(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task244() {};
};

template <typename T>
class Task245 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I261
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& c2 : *range_[0]) {
                for (auto& a3 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

                  // tensor label: I276
                  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, c2, a3, a1);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, c2, a3, a1)]);
                  sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), c2.size(), a3.size(), a1.size());

                  dgemm_("T", "N", 1, ci0.size(), x0.size()*c2.size()*a3.size()*a1.size(),
                         1.0, i0data_sorted, x0.size()*c2.size()*a3.size()*a1.size(), i1data_sorted, x0.size()*c2.size()*a3.size()*a1.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task245(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task245() {};
};

template <typename T>
class Task246 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& c4 : *range_[0]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: f1
              std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
              sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());

              // tensor label: I277
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a3, c4, a1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
              sort_indices<5,1,0,2,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

              dgemm_("T", "N", 1, ci0.size()*x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
                     1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
                     1.0, odata_sorted, 1);
            }
          }

          sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), c2.size(), a3.size(), a1.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task246(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task246() {};
};

template <typename T>
class Task247 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<7,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);
          const Index c2 = b(3);
          const Index a3 = b(4);
          const Index c4 = b(5);
          const Index a1 = b(6);

          // tensor label: I277
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

          // tensor label: I278
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
          sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

          sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), ci0.size(), x1.size(), x0.size());
          out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task247(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  for (auto& ci0 : *range[3])
                    subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,7>{{ci0, x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task247() {};
};

template <typename T>
class Task248 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I278
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task248(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task248() {};
};

template <typename T>
class Task249 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<7,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);
          const Index c2 = b(3);
          const Index a3 = b(4);
          const Index c4 = b(5);
          const Index a1 = b(6);

          // tensor label: I277
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

          // tensor label: I282
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
          sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), ci0.size()*x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

          sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), ci0.size(), x1.size(), x0.size());
          out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task249(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  for (auto& ci0 : *range[3])
                    subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,7>{{ci0, x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task249() {};
};

template <typename T>
class Task250 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I282
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task250(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task250() {};
};

template <typename T>
class Task251 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I285
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task251(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task251() {};
};

template <typename T>
class Task252 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);

          // tensor label: I285
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
          {
            // tensor label: Gamma78
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
            sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size());
          }
          out()->put_block(odata, ci0, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task252(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x3, x0}}, in, t[0], range)));
    };
    ~Task252() {};
};

template <typename T>
class Task253 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I288
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task253(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task253() {};
};

template <typename T>
class Task254 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);

          // tensor label: I288
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
          {
            // tensor label: Gamma78
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size());
          }
          out()->put_block(odata, ci0, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task254(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x3, x0}}, in, t[0], range)));
    };
    ~Task254() {};
};

template <typename T>
class Task255 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& c4 : *range_[0]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());

            // tensor label: I291
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a3, c4, a1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a3, c4, a1)]);
            sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a3.size(), c4.size(), a1.size());

            dgemm_("T", "N", c2.size(), ci0.size()*x0.size()*a3.size()*a1.size(), c4.size(),
                   1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
                   1.0, odata_sorted, c2.size());
          }

          sort_indices<1,2,0,3,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x0.size(), a3.size(), a1.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task255(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task255() {};
};

template <typename T>
class Task256 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index a3 = b(2);
          const Index c4 = b(3);
          const Index a1 = b(4);

          // tensor label: I291
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c4, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());

            // tensor label: I292
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c4.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c4.size()*a1.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task256(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task256() {};
};

template <typename T>
class Task257 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I292
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task257(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task257() {};
};

template <typename T>
class Task258 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index a3 = b(2);
          const Index c4 = b(3);
          const Index a1 = b(4);

          // tensor label: I291
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c4, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());

            // tensor label: I296
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c4.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c4.size()*a3.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task258(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task258() {};
};

template <typename T>
class Task259 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I296
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task259(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task259() {};
};

template <typename T>
class Task260 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& a4 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());

            // tensor label: I299
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c2, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c2, a3)]);
            sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c2.size(), a3.size());

            dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*c2.size()*a3.size(), a4.size(),
                   1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                   1.0, odata_sorted, a1.size());
          }

          sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), c2.size(), a3.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task260(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task260() {};
};

template <typename T>
class Task261 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index a4 = b(2);
          const Index c2 = b(3);
          const Index a3 = b(4);

          // tensor label: I299
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());

            // tensor label: I300
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task261(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task261() {};
};

template <typename T>
class Task262 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I300
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task262(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task262() {};
};

template <typename T>
class Task263 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index a4 = b(2);
          const Index c2 = b(3);
          const Index a3 = b(4);

          // tensor label: I299
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());

            // tensor label: I304
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a4.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task263(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task263() {};
};

template <typename T>
class Task264 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I304
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task264(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task264() {};
};

template <typename T>
class Task265 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& a4 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

            // tensor label: I307
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c2, a1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c2, a1)]);
            sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c2.size(), a1.size());

            dgemm_("T", "N", a3.size(), ci0.size()*x0.size()*c2.size()*a1.size(), a4.size(),
                   1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                   1.0, odata_sorted, a3.size());
          }

          sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x0.size(), c2.size(), a1.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task265(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task265() {};
};

template <typename T>
class Task266 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index a4 = b(2);
          const Index c2 = b(3);
          const Index a1 = b(4);

          // tensor label: I307
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());

            // tensor label: I308
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task266(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task266() {};
};

template <typename T>
class Task267 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I308
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task267(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task267() {};
};

template <typename T>
class Task268 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index a4 = b(2);
          const Index c2 = b(3);
          const Index a1 = b(4);

          // tensor label: I307
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());

            // tensor label: I312
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a4.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task268(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task268() {};
};

template <typename T>
class Task269 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I312
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task269(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task269() {};
};

template <typename T>
class Task270 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());

            // tensor label: I315
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a3)]);
            sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a3.size());

            dgemm_("T", "N", c2.size(), ci0.size()*x0.size()*a1.size()*a3.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, c2.size());
          }

          sort_indices<1,2,0,4,3,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x0.size(), a1.size(), a3.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task270(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task270() {};
};

template <typename T>
class Task271 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index a3 = b(4);

          // tensor label: I315
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a3), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& x3 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

              // tensor label: I316
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
              sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x0.size(), x1.size());
          out()->put_block(odata, ci0, x0, x1, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task271(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, x1, a1, a3}}, in, t[0], range)));
    };
    ~Task271() {};
};

template <typename T>
class Task272 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I316
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task272(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task272() {};
};

template <typename T>
class Task273 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I397
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task273(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task273() {};
};

template <typename T>
class Task274 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

        double e0_;

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<3,1,T>(block, in, out), range_(ran), e0_(e) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I397
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task274(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range, e)));
    };
    ~Task274() {};
};

template <typename T>
class Task275 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I400
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task275(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task275() {};
};

template <typename T>
class Task276 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

        double e0_;

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<3,1,T>(block, in, out), range_(ran), e0_(e) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I400
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task276(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range, e)));
    };
    ~Task276() {};
};

template <typename T>
class Task277 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I415
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task277(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task277() {};
};

template <typename T>
class Task278 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I415
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task278(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task278() {};
};

template <typename T>
class Task279 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I276
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, c2, a3, a1), 0.0);

          for (auto& x1 : *range_[1]) {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I418
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
                   1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
          out()->put_block(odata, ci0, x0, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task279(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task279() {};
};

template <typename T>
class Task280 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I418
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task280(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task280() {};
};

template <typename T>
class Task281 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I261
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

          for (auto& x0 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& x1 : *range_[1]) {
                for (auto& a2 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

                  // tensor label: I318
                  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a2);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a2)]);
                  sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a2.size());

                  dgemm_("T", "N", 1, ci0.size(), x0.size()*x1.size()*a1.size()*a2.size(),
                         1.0, i0data_sorted, x0.size()*x1.size()*a1.size()*a2.size(), i1data_sorted, x0.size()*x1.size()*a1.size()*a2.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task281(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task281() {};
};

template <typename T>
class Task282 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I318
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

          for (auto& c3 : *range_[0]) {
            for (auto& x2 : *range_[1]) {
              // tensor label: f1
              std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
              sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());

              // tensor label: I319
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a1, c3, a2);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a1, c3, a2)]);
              sort_indices<5,2,0,1,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

              dgemm_("T", "N", 1, ci0.size()*x0.size()*x1.size()*a1.size()*a2.size(), x2.size()*c3.size(),
                     1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
                     1.0, odata_sorted, 1);
            }
          }

          sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), x1.size(), a1.size(), a2.size());
          out()->put_block(odata, ci0, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task282(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task282() {};
};

template <typename T>
class Task283 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<7,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);
          const Index a1 = b(4);
          const Index c3 = b(5);
          const Index a2 = b(6);

          // tensor label: I319
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1, c3, a2), 0.0);

          for (auto& x3 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());

            // tensor label: I320
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
            sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

            dgemm_("T", "N", a1.size()*c3.size()*a2.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x3.size(),
                   1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
                   1.0, odata_sorted, a1.size()*c3.size()*a2.size());
          }

          sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), ci0.size(), x0.size(), x2.size(), x1.size());
          out()->put_block(odata, ci0, x0, x2, x1, a1, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task283(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a1 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& x2 : *range[1])
                for (auto& x0 : *range[1])
                  for (auto& ci0 : *range[3])
                    subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,7>{{ci0, x0, x2, x1, a1, c3, a2}}, in, t[0], range)));
    };
    ~Task283() {};
};

template <typename T>
class Task284 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I320
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task284(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task284() {};
};

template <typename T>
class Task285 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I318
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

          for (auto& x4 : *range_[1]) {
            for (auto& x5 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

              // tensor label: I323
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x1)]);
              sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x5.size()*x4.size(),
                     1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
          out()->put_block(odata, ci0, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task285(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task285() {};
};

template <typename T>
class Task286 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x5 = b(1);
          const Index x0 = b(2);
          const Index x4 = b(3);
          const Index x1 = b(4);

          // tensor label: I323
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
          {
            // tensor label: Gamma88
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
            sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
          }
          out()->put_block(odata, ci0, x5, x0, x4, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task286(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x5, x0, x4, x1}}, in, t[0], range)));
    };
    ~Task286() {};
};

template <typename T>
class Task287 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I318
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

          for (auto& a3 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

            // tensor label: I326
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a3)]);
            sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a3.size());

            dgemm_("T", "N", a2.size(), ci0.size()*x0.size()*x1.size()*a1.size(), a3.size(),
                   1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
                   1.0, odata_sorted, a2.size());
          }

          sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x0.size(), x1.size(), a1.size());
          out()->put_block(odata, ci0, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task287(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task287() {};
};

template <typename T>
class Task288 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index a3 = b(4);

          // tensor label: I326
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a3), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& x3 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());

              // tensor label: I327
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
              sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x0.size(), x1.size());
          out()->put_block(odata, ci0, x0, x1, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task288(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, x1, a1, a3}}, in, t[0], range)));
    };
    ~Task288() {};
};

template <typename T>
class Task289 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I327
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task289(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task289() {};
};

template <typename T>
class Task290 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I318
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& x3 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

              // tensor label: I403
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
              sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
          out()->put_block(odata, ci0, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task290(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task290() {};
};

template <typename T>
class Task291 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

        double e0_;

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<5,1,T>(block, in, out), range_(ran), e0_(e) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I403
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            dscal_(ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task291(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range, e)));
    };
    ~Task291() {};
};

template <typename T>
class Task292 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x0 = b(1);
          const Index x1 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I318
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);

          for (auto& x2 : *range_[1]) {
            for (auto& x3 : *range_[1]) {
              // tensor label: v2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

              // tensor label: I421
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
              sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
                     1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
          out()->put_block(odata, ci0, x0, x1, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task292(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range)));
    };
    ~Task292() {};
};

template <typename T>
class Task293 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I421
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task293(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task293() {};
};

template <typename T>
class Task294 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I261
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

          for (auto& x1 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& c2 : *range_[0]) {
                for (auto& a3 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());

                  // tensor label: I343
                  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, c2, a3, a1);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, c2, a3, a1)]);
                  sort_indices<1,4,2,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), c2.size(), a3.size(), a1.size());

                  dgemm_("T", "N", 1, ci0.size(), x1.size()*c2.size()*a3.size()*a1.size(),
                         1.0, i0data_sorted, x1.size()*c2.size()*a3.size()*a1.size(), i1data_sorted, x1.size()*c2.size()*a3.size()*a1.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task294(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task294() {};
};

template <typename T>
class Task295 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I343
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

          for (auto& c4 : *range_[0]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: f1
              std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c4);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, c4)]);
              sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), c4.size());

              // tensor label: I344
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a3, c4, a1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
              sort_indices<5,2,0,1,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());

              dgemm_("T", "N", 1, ci0.size()*x1.size()*c2.size()*a3.size()*a1.size(), x0.size()*c4.size(),
                     1.0, i0data_sorted, x0.size()*c4.size(), i1data_sorted, x0.size()*c4.size(),
                     1.0, odata_sorted, 1);
            }
          }

          sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), c2.size(), a3.size(), a1.size());
          out()->put_block(odata, ci0, x1, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task295(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task295() {};
};

template <typename T>
class Task296 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<7,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);
          const Index c2 = b(3);
          const Index a3 = b(4);
          const Index c4 = b(5);
          const Index a1 = b(6);

          // tensor label: I344
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());

          // tensor label: I345
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
          sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());

          sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), ci0.size(), x1.size(), x0.size());
          out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task296(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  for (auto& ci0 : *range[3])
                    subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,7>{{ci0, x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task296() {};
};

template <typename T>
class Task297 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I345
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task297(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task297() {};
};

template <typename T>
class Task298 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<7,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);
          const Index c2 = b(3);
          const Index a3 = b(4);
          const Index c4 = b(5);
          const Index a1 = b(6);

          // tensor label: I344
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());

          // tensor label: I349
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
          sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

          dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), ci0.size()*x1.size()*x0.size(), 1,
                 1.0, i0data_sorted, 1, i1data_sorted, 1,
                 1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());

          sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), ci0.size(), x1.size(), x0.size());
          out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task298(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& c2 : *range[0])
              for (auto& x0 : *range[1])
                for (auto& x1 : *range[1])
                  for (auto& ci0 : *range[3])
                    subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,7>{{ci0, x1, x0, c2, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task298() {};
};

template <typename T>
class Task299 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I349
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task299(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task299() {};
};

template <typename T>
class Task300 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I343
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

          for (auto& c4 : *range_[0]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());

            // tensor label: I358
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a3, c4, a1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a3, c4, a1)]);
            sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a3.size(), c4.size(), a1.size());

            dgemm_("T", "N", c2.size(), ci0.size()*x1.size()*a3.size()*a1.size(), c4.size(),
                   1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
                   1.0, odata_sorted, c2.size());
          }

          sort_indices<1,2,0,3,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x1.size(), a3.size(), a1.size());
          out()->put_block(odata, ci0, x1, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task300(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task300() {};
};

template <typename T>
class Task301 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index a3 = b(2);
          const Index c4 = b(3);
          const Index a1 = b(4);

          // tensor label: I358
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a3, c4, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c4, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c4.size(), a1.size());

            // tensor label: I359
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c4.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a3.size()*c4.size()*a1.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task301(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task301() {};
};

template <typename T>
class Task302 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I359
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task302(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task302() {};
};

template <typename T>
class Task303 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index a3 = b(2);
          const Index c4 = b(3);
          const Index a1 = b(4);

          // tensor label: I358
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a3, c4, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a3, c4, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a3, c4, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c4, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c4.size(), a3.size());

            // tensor label: I363
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c4.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a1.size()*c4.size()*a3.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, a3, c4, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task303(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c4 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, a3, c4, a1}}, in, t[0], range)));
    };
    ~Task303() {};
};

template <typename T>
class Task304 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I363
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task304(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task304() {};
};

template <typename T>
class Task305 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I343
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

          for (auto& a4 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());

            // tensor label: I366
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c2, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c2, a3)]);
            sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c2.size(), a3.size());

            dgemm_("T", "N", a1.size(), ci0.size()*x1.size()*c2.size()*a3.size(), a4.size(),
                   1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                   1.0, odata_sorted, a1.size());
          }

          sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x1.size(), c2.size(), a3.size());
          out()->put_block(odata, ci0, x1, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task305(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task305() {};
};

template <typename T>
class Task306 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index a4 = b(2);
          const Index c2 = b(3);
          const Index a3 = b(4);

          // tensor label: I366
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a3), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c2.size(), a3.size());

            // tensor label: I367
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task306(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task306() {};
};

template <typename T>
class Task307 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I367
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task307(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task307() {};
};

template <typename T>
class Task308 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index a4 = b(2);
          const Index c2 = b(3);
          const Index a3 = b(4);

          // tensor label: I366
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a3), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a4.size());

            // tensor label: I371
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a4.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, a4, c2, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task308(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, a4, c2, a3}}, in, t[0], range)));
    };
    ~Task308() {};
};

template <typename T>
class Task309 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I371
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task309(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task309() {};
};

template <typename T>
class Task310 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I343
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

          for (auto& a4 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());

            // tensor label: I374
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c2, a1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c2, a1)]);
            sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c2.size(), a1.size());

            dgemm_("T", "N", a3.size(), ci0.size()*x1.size()*c2.size()*a1.size(), a4.size(),
                   1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
                   1.0, odata_sorted, a3.size());
          }

          sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x1.size(), c2.size(), a1.size());
          out()->put_block(odata, ci0, x1, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task310(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task310() {};
};

template <typename T>
class Task311 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index a4 = b(2);
          const Index c2 = b(3);
          const Index a1 = b(4);

          // tensor label: I374
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c2.size(), a1.size());

            // tensor label: I375
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a4.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a4.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task311(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task311() {};
};

template <typename T>
class Task312 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I375
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task312(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task312() {};
};

template <typename T>
class Task313 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index a4 = b(2);
          const Index c2 = b(3);
          const Index a1 = b(4);

          // tensor label: I374
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a4);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a4)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a4.size());

            // tensor label: I379
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a4.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, a4, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task313(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, a4, c2, a1}}, in, t[0], range)));
    };
    ~Task313() {};
};

template <typename T>
class Task314 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I379
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task314(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task314() {};
};

template <typename T>
class Task315 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I343
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I406
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task315(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task315() {};
};

template <typename T>
class Task316 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

        double e0_;

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<3,1,T>(block, in, out), range_(ran), e0_(e) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I406
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task316(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range, e)));
    };
    ~Task316() {};
};

template <typename T>
class Task317 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I343
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I409
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task317(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task317() {};
};

template <typename T>
class Task318 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

        double e0_;

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<3,1,T>(block, in, out), range_(ran), e0_(e) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I409
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task318(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range, e)));
    };
    ~Task318() {};
};

template <typename T>
class Task319 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I343
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I424
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task319(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task319() {};
};

template <typename T>
class Task320 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I424
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task320(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task320() {};
};

template <typename T>
class Task321 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index c2 = b(2);
          const Index a3 = b(3);
          const Index a1 = b(4);

          // tensor label: I343
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I427
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x1.size());
          out()->put_block(odata, ci0, x1, c2, a3, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task321(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          for (auto& c2 : *range[0])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range)));
    };
    ~Task321() {};
};

template <typename T>
class Task322 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x1 = b(1);
          const Index x0 = b(2);

          // tensor label: I427
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
          {
            // tensor label: Gamma74
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
            sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
          }
          out()->put_block(odata, ci0, x1, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task322(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task322() {};
};

template <typename T>
class Task323 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I261
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& c2 : *range_[0]) {
                for (auto& a3 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());

                  // tensor label: I351
                  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, a3, c2, a1);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, a3, c2, a1)]);
                  sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), a3.size(), c2.size(), a1.size());

                  dgemm_("T", "N", 1, ci0.size(), x3.size()*a3.size()*c2.size()*a1.size(),
                         1.0, i0data_sorted, x3.size()*a3.size()*c2.size()*a1.size(), i1data_sorted, x3.size()*a3.size()*c2.size()*a1.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task323(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task323() {};
};

template <typename T>
class Task324 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index a3 = b(2);
          const Index c2 = b(3);
          const Index a1 = b(4);

          // tensor label: I351
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());

            // tensor label: I352
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());

            dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x3.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a3.size()*c2.size()*a1.size());
          }

          sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x3.size());
          out()->put_block(odata, ci0, x3, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task324(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task324() {};
};

template <typename T>
class Task325 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);

          // tensor label: I352
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
          {
            // tensor label: Gamma78
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
            sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size());
          }
          out()->put_block(odata, ci0, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task325(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x3, x0}}, in, t[0], range)));
    };
    ~Task325() {};
};

template <typename T>
class Task326 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index a3 = b(2);
          const Index c2 = b(3);
          const Index a1 = b(4);

          // tensor label: I351
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());

            // tensor label: I355
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
            sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());

            dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x3.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a1.size()*c2.size()*a3.size());
          }

          sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x3.size());
          out()->put_block(odata, ci0, x3, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task326(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task326() {};
};

template <typename T>
class Task327 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<3,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,3>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<3,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);

          // tensor label: I355
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
          {
            // tensor label: Gamma78
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
            sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size());
          }
          out()->put_block(odata, ci0, x3, x0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task327(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x3, x0}}, in, t[0], range)));
    };
    ~Task327() {};
};

template <typename T>
class Task328 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index a3 = b(2);
          const Index c2 = b(3);
          const Index a1 = b(4);

          // tensor label: I351
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, a3, c2, a1);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, a3, c2, a1)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, a3, c2, a1), 0.0);

          for (auto& x2 : *range_[1]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
            sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());

            // tensor label: I382
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a3)]);
            sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a3.size());

            dgemm_("T", "N", c2.size(), ci0.size()*x3.size()*a1.size()*a3.size(), x2.size(),
                   1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
                   1.0, odata_sorted, c2.size());
          }

          sort_indices<1,2,4,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x3.size(), a1.size(), a3.size());
          out()->put_block(odata, ci0, x3, a3, c2, a1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task328(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& a3 : *range[2])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, a3, c2, a1}}, in, t[0], range)));
    };
    ~Task328() {};
};

template <typename T>
class Task329 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x2 = b(2);
          const Index a1 = b(3);
          const Index a3 = b(4);

          // tensor label: I382
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a3)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a3.size());

              // tensor label: I383
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
              sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
                     1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x3.size(), x2.size());
          out()->put_block(odata, ci0, x3, x2, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task329(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x2, a1, a3}}, in, t[0], range)));
    };
    ~Task329() {};
};

template <typename T>
class Task330 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I383
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task330(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task330() {};
};

template <typename T>
class Task331 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I261
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

          for (auto& x3 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& x2 : *range_[1]) {
                for (auto& a2 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());

                  // tensor label: I385
                  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a2);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a2)]);
                  sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a2.size());

                  dgemm_("T", "N", 1, ci0.size(), x3.size()*x2.size()*a1.size()*a2.size(),
                         1.0, i0data_sorted, x3.size()*x2.size()*a1.size()*a2.size(), i1data_sorted, x3.size()*x2.size()*a1.size()*a2.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task331(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task331() {};
};

template <typename T>
class Task332 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x2 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I385
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

          for (auto& c3 : *range_[0]) {
            for (auto& x1 : *range_[1]) {
              // tensor label: f1
              std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c3)]);
              sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c3.size());

              // tensor label: I386
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, a1, c3, a2);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, a1, c3, a2)]);
              sort_indices<5,3,0,1,2,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());

              dgemm_("T", "N", 1, ci0.size()*x3.size()*x2.size()*a1.size()*a2.size(), x1.size()*c3.size(),
                     1.0, i0data_sorted, x1.size()*c3.size(), i1data_sorted, x1.size()*c3.size(),
                     1.0, odata_sorted, 1);
            }
          }

          sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x3.size(), x2.size(), a1.size(), a2.size());
          out()->put_block(odata, ci0, x3, x2, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task332(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x2, a1, a2}}, in, t[0], range)));
    };
    ~Task332() {};
};

template <typename T>
class Task333 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<7,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,7>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<7,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x2 = b(2);
          const Index x1 = b(3);
          const Index a1 = b(4);
          const Index c3 = b(5);
          const Index a2 = b(6);

          // tensor label: I386
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a1, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a1, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a1, c3, a2), 0.0);

          for (auto& x0 : *range_[1]) {
            // tensor label: t2
            std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c3, a2)]);
            sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c3.size(), a2.size());

            // tensor label: I387
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
            sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

            dgemm_("T", "N", a1.size()*c3.size()*a2.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
                   1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
                   1.0, odata_sorted, a1.size()*c3.size()*a2.size());
          }

          sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), ci0.size(), x3.size(), x2.size(), x1.size());
          out()->put_block(odata, ci0, x3, x2, x1, a1, c3, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task333(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a1 : *range[2])
            for (auto& x1 : *range[1])
              for (auto& x2 : *range[1])
                for (auto& x3 : *range[1])
                  for (auto& ci0 : *range[3])
                    subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,7>{{ci0, x3, x2, x1, a1, c3, a2}}, in, t[0], range)));
    };
    ~Task333() {};
};

template <typename T>
class Task334 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I387
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task334(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task334() {};
};

template <typename T>
class Task335 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x2 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I385
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

          for (auto& a3 : *range_[2]) {
            // tensor label: f1
            std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
            std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
            sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());

            // tensor label: I393
            std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a3);
            std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a3)]);
            sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a3.size());

            dgemm_("T", "N", a2.size(), ci0.size()*x3.size()*x2.size()*a1.size(), a3.size(),
                   1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
                   1.0, odata_sorted, a2.size());
          }

          sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x3.size(), x2.size(), a1.size());
          out()->put_block(odata, ci0, x3, x2, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task335(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x2, a1, a2}}, in, t[0], range)));
    };
    ~Task335() {};
};

template <typename T>
class Task336 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x2 = b(2);
          const Index a1 = b(3);
          const Index a3 = b(4);

          // tensor label: I393
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a3);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a3)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a3), 0.0);

          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a3);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a3)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a3.size());

              // tensor label: I394
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
              sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
                     1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
                     1.0, odata_sorted, a1.size()*a3.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x3.size(), x2.size());
          out()->put_block(odata, ci0, x3, x2, a1, a3);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task336(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x2, a1, a3}}, in, t[0], range)));
    };
    ~Task336() {};
};

template <typename T>
class Task337 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I394
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task337(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task337() {};
};

template <typename T>
class Task338 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x2 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I385
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

              // tensor label: I412
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
              sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
                     1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x3.size(), x2.size());
          out()->put_block(odata, ci0, x3, x2, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task338(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x2, a1, a2}}, in, t[0], range)));
    };
    ~Task338() {};
};

template <typename T>
class Task339 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

        double e0_;

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran, const double e)
          : SubTask<5,1,T>(block, in, out), range_(ran), e0_(e) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I412
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            dscal_(ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
            sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task339(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range, double e) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range, e)));
    };
    ~Task339() {};
};

template <typename T>
class Task340 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x2 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I385
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: v2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

              // tensor label: I430
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
              sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
                     1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x3.size(), x2.size());
          out()->put_block(odata, ci0, x3, x2, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task340(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x2, a1, a2}}, in, t[0], range)));
    };
    ~Task340() {};
};

template <typename T>
class Task341 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x3 = b(1);
          const Index x0 = b(2);
          const Index x2 = b(3);
          const Index x1 = b(4);

          // tensor label: I430
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
          {
            // tensor label: Gamma86
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
            sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
          }
          out()->put_block(odata, ci0, x3, x0, x2, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task341(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range)));
    };
    ~Task341() {};
};

template <typename T>
class Task342 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<1,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,1>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<1,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);

          // tensor label: I261
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);

          for (auto& x5 : *range_[1]) {
            for (auto& a1 : *range_[2]) {
              for (auto& x4 : *range_[1]) {
                for (auto& a2 : *range_[2]) {
                  // tensor label: t2
                  std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
                  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
                  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());

                  // tensor label: I389
                  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, a1, a2);
                  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, a1, a2)]);
                  sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), a1.size(), a2.size());

                  dgemm_("T", "N", 1, ci0.size(), x5.size()*x4.size()*a1.size()*a2.size(),
                         1.0, i0data_sorted, x5.size()*x4.size()*a1.size()*a2.size(), i1data_sorted, x5.size()*x4.size()*a1.size()*a2.size(),
                         1.0, odata_sorted, 1);
                }
              }
            }
          }

          sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
          out()->put_block(odata, ci0);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task342(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task342() {};
};

template <typename T>
class Task343 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,2,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,2>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,2,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x5 = b(1);
          const Index x4 = b(2);
          const Index a1 = b(3);
          const Index a2 = b(4);

          // tensor label: I389
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, a1, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, a1, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, a1, a2), 0.0);

          for (auto& x1 : *range_[1]) {
            for (auto& x0 : *range_[1]) {
              // tensor label: t2
              std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
              std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
              sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());

              // tensor label: I390
              std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x1);
              std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x1)]);
              sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());

              dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x5.size()*x4.size(), x0.size()*x1.size(),
                     1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
                     1.0, odata_sorted, a1.size()*a2.size());
            }
          }

          sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x5.size(), x4.size());
          out()->put_block(odata, ci0, x5, x4, a1, a2);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task343(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x5, x4, a1, a2}}, in, t[0], range)));
    };
    ~Task343() {};
};

template <typename T>
class Task344 : public DedciTask<T> {
  protected:
    class Task_local : public SubTask<5,1,T> {
      protected:
        const std::array<std::shared_ptr<const IndexRange>,4> range_;

        const Index& b(const size_t& i) const { return this->block(i); }
        const std::shared_ptr<const Tensor<T>>& in(const size_t& i) const { return this->in_tensor(i); }
        const std::shared_ptr<Tensor<T>>& out() const { return this->out_tensor(); }

      public:
        Task_local(const std::array<const Index,5>& block, const std::array<std::shared_ptr<const Tensor<T>>,1>& in, std::shared_ptr<Tensor<T>>& out,
                   std::array<std::shared_ptr<const IndexRange>,4>& ran)
          : SubTask<5,1,T>(block, in, out), range_(ran) { }


        void compute() override {
          const Index ci0 = b(0);
          const Index x5 = b(1);
          const Index x0 = b(2);
          const Index x4 = b(3);
          const Index x1 = b(4);

          // tensor label: I390
          std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
          {
            // tensor label: Gamma88
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
            sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
          }
          out()->put_block(odata, ci0, x5, x0, x4, x1);
        }
    };

    std::vector<std::shared_ptr<Task_local>> subtasks_;

    void compute_() override {
      for (auto& i : subtasks_) {
        i->compute();
      }
    }

  public:
    Task344(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, x5, x0, x4, x1}}, in, t[0], range)));
    };
    ~Task344() {};
};


}
}
}
#endif

