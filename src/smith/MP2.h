//
// Newint - Parallel electron correlation program.
// Filename: MP2.h
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


#ifndef __SRC_SMITH_MP2_H 
#define __SRC_SMITH_MP2_H 

#include <src/smith/spinfreebase.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <iostream>
#include <iomanip>
#include <src/smith/queue.h>
#include <src/smith/MP2_tasks.h>
#include <src/smith/smith.h>

namespace SMITH {
namespace MP2{

template <typename T>
class MP2 : public SpinFreeMethod<T>, SMITH_info {
  protected:
    std::shared_ptr<Queue<T> > queue_;
    std::shared_ptr<Queue<T> > energy_;
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > r;

  public:
    MP2(std::shared_ptr<Reference> ref) : SpinFreeMethod<T>(ref), SMITH_info(), queue_(new Queue<T>()), energy_(new Queue<T>()) {
      this->eig_ = this->f1_->diag();
      t2 = this->v2_->clone();
      r = t2->clone();
      std::vector<IndexRange> index = vec(this->closed_, this->act_, this->virt_);

      std::vector<std::shared_ptr<Tensor<T> > > tensor0 = vec(r);
      std::shared_ptr<Task0<T> > task0(new Task0<T>(tensor0, index));
      queue_->add_task(task0);

#if 0
      std::vector<IndexRange> I0_index = vec(this->closed_, this->virt_, this->closed_, this->virt_);
      std::shared_ptr<Tensor<T> > I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor1 = vec(r, I0);
      std::shared_ptr<Task1<T> > task1(new Task1<T>(tensor1, index));
      task0->add_dep(task1);
      queue_->add_task(task1);

      std::vector<IndexRange> I1_index;
      std::shared_ptr<Tensor<T> > I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor2 = vec(I0, t2, I1);
      std::shared_ptr<Task2<T> > task2(new Task2<T>(tensor2, index));
      task1->add_dep(task2);
      queue_->add_task(task2);

      std::vector<IndexRange> I3_index;
      std::shared_ptr<Tensor<T> > I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor3 = vec(I0, t2, I3);
      std::shared_ptr<Task3<T> > task3(new Task3<T>(tensor3, index));
      task1->add_dep(task3);
      queue_->add_task(task3);
#endif

      std::vector<IndexRange> I4_index = vec(this->closed_, this->virt_, this->virt_, this->closed_);
      std::shared_ptr<Tensor<T> > I4(new Tensor<T>(I4_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor4 = vec(r, I4);
      std::shared_ptr<Task4<T> > task4(new Task4<T>(tensor4, index));
      task0->add_dep(task4);
      queue_->add_task(task4);

      std::vector<IndexRange> I5_index = vec(this->closed_, this->virt_, this->closed_, this->virt_);
      std::shared_ptr<Tensor<T> > I5(new Tensor<T>(I5_index, true));
      std::vector<std::shared_ptr<Tensor<T> > > tensor5 = vec(I4, this->f1_, I5);
      std::shared_ptr<Task5<T> > task5(new Task5<T>(tensor5, index));
      task4->add_dep(task5);
      queue_->add_task(task5);

      std::vector<IndexRange> I9_index = vec(this->closed_, this->virt_, this->closed_, this->virt_);
      std::shared_ptr<Tensor<T> > I9(new Tensor<T>(I9_index, true));
      std::vector<std::shared_ptr<Tensor<T> > > tensor6 = vec(I4, this->f1_, I9);
      std::shared_ptr<Task6<T> > task6(new Task6<T>(tensor6, index));
      task4->add_dep(task6);
      queue_->add_task(task6);

    };
    ~MP2() {}; 

    void solve() {
      t2->zero();
      this->print_iteration();
      int iter = 0;
      for ( ; iter != maxiter_; ++iter) {
        queue_->initialize();
        while (!queue_->done())
          queue_->next()->compute();
        update_amplitude(t2, r);
        const double en = energy();
        const double err = r->rms();
        this->print_iteration(iter, en, err);
        if (err < thresh_residual()) break;
      }
      this->print_iteration(iter == maxiter_);
    };

    double energy() {
      double en = 0.0;
      energy_->initialize();
      while (!energy_->done()) {
        std::shared_ptr<Task<T> > c = energy_->next();
        c->compute();
        en += c->energy();
      }   
      return en; 
    };  
};

}
}
#endif

