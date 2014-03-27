//
// BAGEL - Parallel electron correlation program.
// Filename: MP2.h
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


#ifndef __SRC_SMITH_MP2_H
#define __SRC_SMITH_MP2_H

#include <src/smith/spinfreebase.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <iostream>
#include <iomanip>
#include <src/smith/queue.h>
#include <src/smith/MP2_tasks.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace MP2{

template <typename T>
class MP2 : public SpinFreeMethod<T> {
  protected:
    using SpinFreeMethod<T>::ref_;
    using SpinFreeMethod<T>::closed_;
    using SpinFreeMethod<T>::active_;
    using SpinFreeMethod<T>::virt_;

    std::shared_ptr<Tensor<T>> t2;
    std::shared_ptr<Tensor<T>> r;
    double e0_;

    std::pair<std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>> make_queue_() {
      std::shared_ptr<Queue<T>> queue_(new Queue<T>());
      std::vector<IndexRange> index = {closed_, active_, virt_};

      std::vector<std::shared_ptr<Tensor<T>>> tensor0 = {r};
      std::shared_ptr<Task0<T>> task0(new Task0<T>(tensor0, index));
      queue_->add_task(task0);

      std::vector<IndexRange> Gamma0_index;
      std::shared_ptr<Tensor<T>> Gamma0(new Tensor<T>(Gamma0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor1 = {Gamma0, this->rdm1_, this->f1_};
      std::shared_ptr<Task1<T>> task1(new Task1<T>(tensor1, index));
      task1->add_dep(task0);
      queue_->add_task(task1);

      std::vector<IndexRange> I0_index = {closed_, virt_, closed_, virt_};
      std::shared_ptr<Tensor<T>> I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor2 = {r, I0};
      std::shared_ptr<Task2<T>> task2(new Task2<T>(tensor2, index));
      task2->add_dep(task0);
      queue_->add_task(task2);


      std::vector<std::shared_ptr<Tensor<T>>> tensor3 = {I0, t2, this->v2_};
      std::shared_ptr<Task3<T>> task3(new Task3<T>(tensor3, index, this->e0_));
      task2->add_dep(task3);
      task3->add_dep(task0);
      queue_->add_task(task3);


      std::vector<IndexRange> I1_index;
      std::shared_ptr<Tensor<T>> I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor4 = {I0, t2, I1};
      std::shared_ptr<Task4<T>> task4(new Task4<T>(tensor4, index));
      task2->add_dep(task4);
      task4->add_dep(task0);
      queue_->add_task(task4);


      std::vector<std::shared_ptr<Tensor<T>>> tensor5 = {I1, Gamma0};
      std::shared_ptr<Task5<T>> task5(new Task5<T>(tensor5, index));
      task4->add_dep(task5);
      task5->add_dep(task0);
      queue_->add_task(task5);

      task5->add_dep(task1);

      std::vector<IndexRange> I3_index;
      std::shared_ptr<Tensor<T>> I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor6 = {I0, t2, I3};
      std::shared_ptr<Task6<T>> task6(new Task6<T>(tensor6, index));
      task2->add_dep(task6);
      task6->add_dep(task0);
      queue_->add_task(task6);


      std::vector<std::shared_ptr<Tensor<T>>> tensor7 = {I3, Gamma0};
      std::shared_ptr<Task7<T>> task7(new Task7<T>(tensor7, index));
      task6->add_dep(task7);
      task7->add_dep(task0);
      queue_->add_task(task7);

      task7->add_dep(task1);

      std::vector<IndexRange> I4_index = {closed_, virt_, virt_, closed_};
      std::shared_ptr<Tensor<T>> I4(new Tensor<T>(I4_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor8 = {r, I4};
      std::shared_ptr<Task8<T>> task8(new Task8<T>(tensor8, index));
      task8->add_dep(task0);
      queue_->add_task(task8);


      std::vector<IndexRange> I5_index = {closed_, virt_, closed_, virt_};
      std::shared_ptr<Tensor<T>> I5(new Tensor<T>(I5_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor9 = {I4, this->f1_, I5};
      std::shared_ptr<Task9<T>> task9(new Task9<T>(tensor9, index));
      task8->add_dep(task9);
      task9->add_dep(task0);
      queue_->add_task(task9);


      std::vector<std::shared_ptr<Tensor<T>>> tensor10 = {I5, t2};
      std::shared_ptr<Task10<T>> task10(new Task10<T>(tensor10, index));
      task9->add_dep(task10);
      task10->add_dep(task0);
      queue_->add_task(task10);


      std::vector<IndexRange> I9_index = {closed_, virt_, closed_, virt_};
      std::shared_ptr<Tensor<T>> I9(new Tensor<T>(I9_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor11 = {I4, this->f1_, I9};
      std::shared_ptr<Task11<T>> task11(new Task11<T>(tensor11, index));
      task8->add_dep(task11);
      task11->add_dep(task0);
      queue_->add_task(task11);


      std::vector<std::shared_ptr<Tensor<T>>> tensor12 = {I9, t2};
      std::shared_ptr<Task12<T>> task12(new Task12<T>(tensor12, index));
      task11->add_dep(task12);
      task12->add_dep(task0);
      queue_->add_task(task12);


      std::shared_ptr<Queue<T>> energy_(new Queue<T>());
      std::vector<IndexRange> I16_index;
      std::shared_ptr<Tensor<T>> I16(new Tensor<T>(I16_index, false));
      std::vector<IndexRange> I17_index = {closed_, virt_, closed_, virt_};
      std::shared_ptr<Tensor<T>> I17(new Tensor<T>(I17_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor13 = {I16, t2, I17};
      std::shared_ptr<Task13<T>> task13(new Task13<T>(tensor13, index));
      energy_->add_task(task13);


      std::vector<std::shared_ptr<Tensor<T>>> tensor14 = {I17, this->v2_};
      std::shared_ptr<Task14<T>> task14(new Task14<T>(tensor14, index));
      task13->add_dep(task14);
      energy_->add_task(task14);


      return make_pair(queue_, energy_);
    };

  public:
    MP2(std::shared_ptr<const SMITH_Info> ref) : SpinFreeMethod<T>(ref) {
      this->eig_ = this->f1_->diag();
      t2 = this->v2_->clone();
      e0_ = this->e0();
#if 1
      this->update_amplitude(t2, this->v2_, true);
      t2->scale(2.0);
#endif
      r = t2->clone();
    };
    ~MP2() {};

    void solve() {
      this->print_iteration();
      int iter = 0;
      for ( ; iter != ref_->maxiter(); ++iter) {
        std::pair<std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>> q = make_queue_();
        std::shared_ptr<Queue<T>> queue = q.first;
        std::shared_ptr<Queue<T>> energ = q.second;
        while (!queue->done())
          queue->next_compute();
//      *r = *(r->add_dagger());
        this->update_amplitude(t2, r);
        const double err = r->rms();
        r->zero();
        energy_ = energy(energ);
        this->print_iteration(iter, energy_, err);
        if (err < ref_->thresh()) break;
      }
      this->print_iteration(iter == ref_->maxiter());
    };

    double energy(std::shared_ptr<Queue<T>> energ) {
      double en = 0.0;
      while (!energ->done()) {
        std::shared_ptr<Task<T>> c = energ->next_compute();
        en += c->energy() * 0.25; // FIXME
      }
      return en;
    };
};

}
}
}
#endif

