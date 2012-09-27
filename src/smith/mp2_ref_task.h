//
// BAGEL - Parallel electron correlation program.
// Filename: mp2_ref_task.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

// This file will be replaced by a generated version

#ifndef __SRC_SMITH_MP2_REF_TASK_H
#define __SRC_SMITH_MP2_REF_TASK_H

#include <memory>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/task.h>
#include <vector>

namespace bagel {
namespace SMITH {

template <typename T>
class Task0 : public Task<T> {
  protected:
    std::shared_ptr<Tensor<T> > r2_;
    std::shared_ptr<Tensor<T> > f1_;
    std::shared_ptr<Tensor<T> > t2_;
    IndexRange closed_;
    IndexRange virt_;

    void compute_() {
      for (auto i3 = virt_.begin(); i3 != virt_.end(); ++i3) {
        for (auto i2 = closed_.begin(); i2 != closed_.end(); ++i2) {
          for (auto i1 = virt_.begin(); i1 != virt_.end(); ++i1) {
            for (auto i0 = closed_.begin(); i0 != closed_.end(); ++i0) {
              std::vector<size_t> ohash = {i0->key(), i1->key(), i2->key(), i3->key()};
              std::unique_ptr<double[]> odata = r2_->move_block(ohash);

              for (auto c0 = closed_.begin(); c0 != closed_.end(); ++c0) {
                std::vector<size_t> ihash0 = {c0->key(), i0->key()};
                const std::unique_ptr<double[]> idata0 = f1_->get_block(ihash0);

                std::vector<size_t> ihash1 = {c0->key(), i1->key(), i2->key(), i3->key()};
                std::vector<size_t> ihash2 = {c0->key(), i3->key(), i2->key(), i1->key()};
                const std::unique_ptr<double[]> idata1 = t2_->get_block(ihash1);
                const std::unique_ptr<double[]> idata2 = t2_->get_block(ihash2);
                std::unique_ptr<double[]> idata3(new double[t2_->get_size(ihash1)]);

                assert(t2_->get_size(ihash1) == t2_->get_size(ihash2));

                sort_indices<0,3,2,1,0,1,-1,1>(idata2, idata3, c0->size(), i3->size(), i2->size(), i1->size());
                daxpy_(t2_->get_size(ihash1), 2.0, idata1, 1, idata3, 1);

                const int common = c0->size();
                const int isize0 = i0->size();
                const int isize1 = i1->size() * i2->size() * i3->size();
                dgemm_("T", "N", isize0, isize1, common, -1.0, idata0, common, idata3, common, 1.0, odata, isize0);
              }
              r2_->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task0(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      if (t.size() != 3 || i.size() != 2) throw std::logic_error("Task0 error");
      r2_ =  t[0];
      f1_ =  t[1];
      t2_ =  t[2];
      closed_ = i[0];
      virt_   = i[1];
    };
    ~Task0() {};

};


template <typename T>
class Task1 : public Task<T> {
  protected:
    std::shared_ptr<Tensor<T> > r2_;
    std::shared_ptr<Tensor<T> > f1_;
    std::shared_ptr<Tensor<T> > t2_;
    IndexRange closed_;
    IndexRange virt_;

    void compute_() {
      for (auto i3 = virt_.begin(); i3 != virt_.end(); ++i3) {
        for (auto i2 = closed_.begin(); i2 != closed_.end(); ++i2) {
          for (auto i1 = virt_.begin(); i1 != virt_.end(); ++i1) {
            for (auto i0 = closed_.begin(); i0 != closed_.end(); ++i0) {
              std::vector<size_t> ohash = {i0->key(), i1->key(), i2->key(), i3->key()};
              std::unique_ptr<double[]> odata = r2_->move_block(ohash);
              std::unique_ptr<double[]> idata4(new double[r2_->get_size(ohash)]);

              for (auto c0 = virt_.begin(); c0 != virt_.end(); ++c0) {
                std::vector<size_t> ihash0 = {c0->key(), i1->key()};
                const std::unique_ptr<double[]> idata0 = f1_->get_block(ihash0);

                std::vector<size_t> ihash1 = {i0->key(), c0->key(), i2->key(), i3->key()};
                std::vector<size_t> ihash2 = {i0->key(), i3->key(), i2->key(), c0->key()};
                const std::unique_ptr<double[]> idata1 = t2_->get_block(ihash1);
                const std::unique_ptr<double[]> idata2 = t2_->get_block(ihash2);

                std::unique_ptr<double[]> idata3(new double[t2_->get_size(ihash1)]);
                sort_indices<1,0,2,3,0,1,2,1>(idata1, idata3, i0->size(), c0->size(), i2->size(), i3->size());
                sort_indices<3,0,2,1,1,1,-1,1>(idata2, idata3, i0->size(), i3->size(), i2->size(), c0->size());

                const int common = c0->size();
                const int isize0 = i1->size();
                const int isize1 = i0->size() * i2->size() * i3->size();
                dgemm_("T", "N", isize0, isize1, common, 1.0, idata0, common, idata3, common, 0.0, idata4, isize0);

                sort_indices<1,0,2,3,1,1,1,1>(idata4, odata, i1->size(), i0->size(), i2->size(), i3->size());
              }
              r2_->put_block(ohash, odata);
            }
          }
        }
      }
    };

  public:
    Task1(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      if (t.size() != 3 || i.size() != 2) throw std::logic_error("Task1 error");
      r2_ =  t[0];
      f1_ =  t[1];
      t2_ =  t[2];
      closed_ = i[0];
      virt_   = i[1];
    };
    ~Task1() {};

};


template <typename T>
class Task2 : public Task<T> {
  protected:
    std::shared_ptr<Tensor<T> > r2_;
    std::shared_ptr<Tensor<T> > v2_;
    IndexRange closed_;
    IndexRange virt_;

    void compute_() {
      for (auto i3 = virt_.begin(); i3 != virt_.end(); ++i3) {
        for (auto i2 = closed_.begin(); i2 != closed_.end(); ++i2) {
          for (auto i1 = virt_.begin(); i1 != virt_.end(); ++i1) {
            for (auto i0 = closed_.begin(); i0 != closed_.end(); ++i0) {
              std::vector<size_t> h = {i0->key(), i1->key(), i2->key(), i3->key()};
              std::vector<size_t> g = {i0->key(), i3->key(), i2->key(), i1->key()};
              std::unique_ptr<double[]> odata = r2_->get_block(h);
              std::unique_ptr<double[]> data0 = v2_->get_block(h);
              const std::unique_ptr<double[]> data1 = v2_->get_block(g);
              sort_indices<0,3,2,1,2,1,-1,1>(data1, data0, i0->size(), i3->size(), i2->size(), i1->size());
              daxpy_(v2_->get_size(h), 1.0, data0, 1, odata, 1);
              r2_->put_block(h,odata);
            }
          }
        }
      }
    };

  public:
    Task2(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      if (t.size() != 2 || i.size() != 2) throw std::logic_error("Task2 error");
      r2_ =  t[0];
      v2_ =  t[1];
      closed_ = i[0];
      virt_   = i[1];
    };
    ~Task2() {};

};


template <typename T>
class Task3 : public Task<T> {
  protected:
    std::shared_ptr<Tensor<T> > r2_;
    IndexRange closed_;
    IndexRange virt_;

    void compute_() {
      r2_->zero();
    };

  public:
    Task3(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      if (t.size() != 1 || i.size() != 2) throw std::logic_error("Task3 error");
      r2_ =  t[0];
      closed_ = i[0];
      virt_   = i[1];
    };
    ~Task3() {};

};


template <typename T>
class Task4 : public Task<T> {
  protected:
    std::shared_ptr<Tensor<T> > r2_;
    IndexRange closed_;
    IndexRange virt_;

    void compute_() {
      *r2_ = *(r2_->add_dagger());
    };

  public:
    Task4(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : Task<T>() {
      if (t.size() != 1 || i.size() != 2) throw std::logic_error("Task4 error");
      r2_ =  t[0];
      closed_ = i[0];
      virt_   = i[1];
    };
    ~Task4() {};

};


template <typename T>
class Task5 : public EnergyTask<T> {
  protected:
    std::shared_ptr<Tensor<T> > t2_;
    std::shared_ptr<Tensor<T> > v2_;
    IndexRange closed_;
    IndexRange virt_;

    void compute_() {
      std::vector<IndexRange> o = t2_->indexrange();
      assert(o.size() == 4);
      this->energy_ = 0.0;
      for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
        for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
          for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
            for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
              std::vector<size_t> h(4); h[0] = i0->key(); h[1] = i1->key(); h[2] = i2->key(); h[3] = i3->key();
              std::vector<size_t> g(4); g[0] = i0->key(); g[1] = i3->key(); g[2] = i2->key(); g[3] = i1->key();
              const std::unique_ptr<double[]> data = v2_->get_block(h);
              std::unique_ptr<double[]> data0 = t2_->get_block(h);
              const std::unique_ptr<double[]> data1 = t2_->get_block(g);
              sort_indices<0,3,2,1,2,1,-1,1>(data1, data0, i0->size(), i3->size(), i2->size(), i1->size());
              this->energy_ += 0.5*ddot_(t2_->get_size(h), data, 1, data0, 1);
            }
          }
        }
      }
    };

  public:
    Task5(std::vector<std::shared_ptr<Tensor<T> > > t, std::vector<IndexRange> i) : EnergyTask<T>() {
      if (t.size() != 2 || i.size() != 2) throw std::logic_error("Task5 error");
      t2_ =  t[0];
      v2_ =  t[1];
      closed_ = i[0];
      virt_   = i[1];
    };
    ~Task5() {};

};

}
}

#endif
