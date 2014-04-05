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
          // std::shared_ptr<Tensor<T> > Gamma4;
          // std::shared_ptr<Tensor<T> > rdm1I0;
          // std::shared_ptr<Tensor<T> > f1;

          // tensor label: Gamma4
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
    Task3(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,4> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
          for (auto& ci0 : *range[3])
            for (auto& x0 : *range[1])
              for (auto& x1 : *range[1])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range)));
    };
    ~Task3() {};
};

template <typename T>
class Task4 : public Task<T> {
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
    Task4(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c3 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task4() {};
};

template <typename T>
class Task5 : public Task<T> {
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
            dscal_(c1.size()*a4.size()*c3.size()*a2.size(), -e0_, i0data.get(), 1);
            sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          {
            // tensor label: t2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
            dscal_(c1.size()*a2.size()*c3.size()*a4.size(), -e0_, i1data.get(), 1);
            sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
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
    Task5(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range, const double e) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range, e)));
    };
    ~Task5() {};
};

template <typename T>
class Task6 : public Task<T> {
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
    Task6(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task6() {};
};

template <typename T>
class Task7 : public Task<T> {
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
    Task7(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(in, t[0], range)));
    };
    ~Task7() {};
};

template <typename T>
class Task8 : public Task<T> {
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
    Task8(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task8() {};
};

template <typename T>
class Task9 : public Task<T> {
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
    Task9(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(in, t[0], range)));
    };
    ~Task9() {};
};

template <typename T>
class Task10 : public Task<T> {
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
    Task10(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c3 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task10() {};
};

template <typename T>
class Task11 : public Task<T> {
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
    Task11(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, a2, c3}}, in, t[0], range)));
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
    Task12(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c5 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c5, a2}}, in, t[0], range)));
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
    Task13(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, a2, c3}}, in, t[0], range)));
    };
    ~Task13() {};
};

template <typename T>
class Task14 : public Task<T> {
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
    Task14(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a5 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a5, c3, a2}}, in, t[0], range)));
    };
    ~Task14() {};
};

template <typename T>
class Task15 : public EnergyTask<T> {
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

          // tensor label: I17
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
    Task15(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a2 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range)));
    };
    ~Task15() {};
};

template <typename T>
class Task16 : public EnergyTask<T> {
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
          const Index c3 = b(2);
          const Index a2 = b(3);

          // tensor label: I17
          std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
          {
            // tensor label: v2
            std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
            sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
          }
          {
            // tensor label: v2
            std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
            sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
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
    Task16(std::vector<std::shared_ptr<Tensor<T>>> t,  std::array<std::shared_ptr<const IndexRange>,3> range) : EnergyTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

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
class Task17 : public CorrectionTask<T> {
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

          // tensor label: I21
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
    Task17(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a2 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range)));
    };
    ~Task17() {};
};

template <typename T>
class Task18 : public CorrectionTask<T> {
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

          // tensor label: I21
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
    Task18(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : CorrectionTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task18() {};
};

template <typename T>
class Task19 : public DensityTask<T> {
  protected:
    std::shared_ptr<Tensor<T>> d_;
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;

    void compute_() {
      d_->zero();
    };

  public:
    Task19(std::vector<std::shared_ptr<Tensor<T>>> t) : DensityTask<T>() {
      d_ =  t[0];
    };
    ~Task19() {};
};

template <typename T>
class Task20 : public DensityTask<T> {
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

          // tensor label: den1
          std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
          {
            // tensor label: I24
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
    Task20(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x0, x1}}, in, t[0], range)));
    };
    ~Task20() {};
};

template <typename T>
class Task21 : public DensityTask<T> {
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

          // tensor label: I24
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

                  // tensor label: I25
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
    Task21(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task21() {};
};

template <typename T>
class Task22 : public DensityTask<T> {
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

          // tensor label: I25
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

          // tensor label: I26
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
    Task22(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
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
    ~Task22() {};
};

template <typename T>
class Task23 : public DensityTask<T> {
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

          // tensor label: I26
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
    Task23(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task23() {};
};

template <typename T>
class Task24 : public DensityTask<T> {
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

          // tensor label: I25
          std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I29
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
    Task24(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
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
    ~Task24() {};
};

template <typename T>
class Task25 : public DensityTask<T> {
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

          // tensor label: I29
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
    Task25(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{x1, x0}}, in, t[0], range)));
    };
    ~Task25() {};
};

template <typename T>
class Task26 : public DensityTask<T> {
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

          // tensor label: den1
          std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
          {
            // tensor label: I30
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
    Task26(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& c5 : *range[0])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{c5, c3}}, in, t[0], range)));
    };
    ~Task26() {};
};

template <typename T>
class Task27 : public DensityTask<T> {
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

          // tensor label: I30
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

                // tensor label: I31
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
    Task27(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
      for (auto& c3 : *range[0])
        for (auto& c5 : *range[0])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{c5, c3}}, in, t[0], range)));
    };
    ~Task27() {};
};

template <typename T>
class Task28 : public DensityTask<T> {
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

          // tensor label: I31
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
    Task28(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c5 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c5, a2}}, in, t[0], range)));
    };
    ~Task28() {};
};

template <typename T>
class Task29 : public DensityTask<T> {
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

          // tensor label: den1
          std::unique_ptr<double[]> odata = out()->move_block(a4, a5);
          {
            // tensor label: I34
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
    Task29(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a5 : *range[2])
        for (auto& a4 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a4, a5}}, in, t[0], range)));
    };
    ~Task29() {};
};

template <typename T>
class Task30 : public DensityTask<T> {
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

          // tensor label: I34
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

                // tensor label: I35
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
    Task30(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
      for (auto& a4 : *range[2])
        for (auto& a5 : *range[2])
          subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,2>{{a5, a4}}, in, t[0], range)));
    };
    ~Task30() {};
};

template <typename T>
class Task31 : public DensityTask<T> {
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

          // tensor label: I35
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
    Task31(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : DensityTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a5 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a5, c3, a2}}, in, t[0], range)));
    };
    ~Task31() {};
};

template <typename T>
class Task32 : public Density2Task<T> {
  protected:
    std::shared_ptr<Tensor<T>> d2_;
    IndexRange closed_;
    IndexRange active_;
    IndexRange virt_;

    void compute_() {
      d2_->zero();
    };

  public:
    Task32(std::vector<std::shared_ptr<Tensor<T>>> t) : Density2Task<T>() {
      d2_ =  t[0];
    };
    ~Task32() {};
};

template <typename T>
class Task33 : public Density2Task<T> {
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

          // tensor label: den2
          std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
          {
            // tensor label: I38
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
    Task33(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c3 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range)));
    };
    ~Task33() {};
};

template <typename T>
class Task34 : public Density2Task<T> {
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

          // tensor label: I38
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
    Task34(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,3> range) : Density2Task<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task34() {};
};

template <typename T>
class Task35 : public DedciTask<T> {
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
    Task35(std::vector<std::shared_ptr<Tensor<T>>> t) : DedciTask<T>() {
      dec_ =  t[0];
    };
    ~Task35() {};
};

template <typename T>
class Task36 : public DedciTask<T> {
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
            // tensor label: I40
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
    Task36(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task36() {};
};

template <typename T>
class Task37 : public DedciTask<T> {
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

          // tensor label: I40
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

                  // tensor label: I41
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
    Task37(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task37() {};
};

template <typename T>
class Task38 : public DedciTask<T> {
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

          // tensor label: I41
          std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());

          // tensor label: I42
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
    Task38(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task38() {};
};

template <typename T>
class Task39 : public DedciTask<T> {
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

          // tensor label: I42
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          {
            // tensor label: Gamma4
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
            sort_indices<0,1,1,-1,2>(i0data, odata, ci0.size());
          }
          {
            // tensor label: Gamma4
            std::unique_ptr<double[]> i1data = in(0)->get_block(ci0);
            sort_indices<0,1,1,-1,2>(i1data, odata, ci0.size());
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
    Task39(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task39() {};
};

template <typename T>
class Task40 : public DedciTask<T> {
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

          // tensor label: I41
          std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
          std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
          std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);

          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());

          // tensor label: I45
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
    Task40(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,2> in = {{t[1], t[2]}};

      subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          for (auto& a4 : *range[2])
            for (auto& c1 : *range[0])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range)));
    };
    ~Task40() {};
};

template <typename T>
class Task41 : public DedciTask<T> {
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

          // tensor label: I45
          std::unique_ptr<double[]> odata = out()->move_block(ci0);
          {
            // tensor label: Gamma4
            std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
            sort_indices<0,1,1,1,1>(i0data, odata, ci0.size());
          }
          {
            // tensor label: Gamma4
            std::unique_ptr<double[]> i1data = in(0)->get_block(ci0);
            sort_indices<0,1,1,1,1>(i1data, odata, ci0.size());
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
    Task41(std::vector<std::shared_ptr<Tensor<T>>> t, std::array<std::shared_ptr<const IndexRange>,4> range) : DedciTask<T>() {
      std::array<std::shared_ptr<const Tensor<T>>,1> in = {{t[1]}};

      subtasks_.reserve(range[3]->nblock());
      for (auto& ci0 : *range[3])
        subtasks_.push_back(std::shared_ptr<Task_local>(new Task_local(std::array<const Index,1>{{ci0}}, in, t[0], range)));
    };
    ~Task41() {};
};


}
}
}
#endif

