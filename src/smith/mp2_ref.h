//
// BAGEL - Parallel electron correlation program.
// Filename: mp2_ref.h
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


#ifndef __SRC_SMITH_MP2_REF_H
#define __SRC_SMITH_MP2_REF_H

#include <src/smith/spinfreebase.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <iostream>
#include <iomanip>
#include <src/smith/queue.h>
#include <src/smith/mp2_ref_task.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {

template <typename T>
class MP2_Ref : public SpinFreeMethod<T>, SMITH_info {
  protected:
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > r2;

    std::pair<std::shared_ptr<Queue<T> >, std::shared_ptr<Queue<T> > > make_queue_() {
      std::shared_ptr<Queue<T> > queue_(new Queue<T>());
      std::shared_ptr<Queue<T> > energy_(new Queue<T>());

      std::vector<std::shared_ptr<Tensor<T> > > tensor0 = vec(r2, this->f1_, t2);
      std::vector<std::shared_ptr<Tensor<T> > > tensor1 = vec(r2, this->v2_);
      std::vector<std::shared_ptr<Tensor<T> > > tensor2 = vec(r2);
      std::vector<IndexRange> index0 =  vec(this->closed_, this->virt_);

      std::shared_ptr<Task0<T> > t0(new Task0<T>(tensor0, index0));
      std::shared_ptr<Task1<T> > t1(new Task1<T>(tensor0, index0));
      std::shared_ptr<Task2<T> > tt2(new Task2<T>(tensor1, index0));
      std::shared_ptr<Task3<T> > t3(new Task3<T>(tensor2, index0));
      t0->add_dep(t3);
      t1->add_dep(t3);
      tt2->add_dep(t3);
      std::shared_ptr<Task4<T> > t4(new Task4<T>(tensor2, index0));
      t4->add_dep(t0);
      t4->add_dep(t1);
      t4->add_dep(tt2);
      t4->add_dep(t3);

      queue_->add_task(t3);
      queue_->add_task(t4);
      queue_->add_task(tt2);
      queue_->add_task(t1);
      queue_->add_task(t0);

      std::vector<std::shared_ptr<Tensor<T> > > tensor5 = vec(t2, this->v2_);
      std::shared_ptr<Task5<T> > t5(new Task5<T>(tensor5, index0));
      energy_->add_task(t5);
      return make_pair(queue_, energy_);
    };

  public:
    MP2_Ref(std::shared_ptr<Reference> r) : SpinFreeMethod<T>(r), SMITH_info() {
      this->eig_ = this->f1_->diag();
      t2 = this->v2_->clone();
      r2 = t2->clone();
    };

    ~MP2_Ref() {};

    void solve() {
      t2->zero();
      this->print_iteration();
      int iter;
      for (iter = 0; iter != maxiter_; ++iter) {

        std::pair<std::shared_ptr<Queue<T> >, std::shared_ptr<Queue<T> > >  q = make_queue_();
        std::shared_ptr<Queue<T> > queue = q.first;
        std::shared_ptr<Queue<T> > eng = q.second;
        while (!queue->done()) queue->next_compute();

        this->update_amplitude(t2, r2);
        const double err = r2->rms();
        const double en = energy(eng);

        this->print_iteration(iter, en, err);
        if (err < thresh_residual()) break;
      }
      this->print_iteration(iter == maxiter_);
    };

    double energy(std::shared_ptr<Queue<T> > eng) {
      double en = 0.0;
      eng->initialize();
      while (!eng->done()) {
        std::shared_ptr<Task<T> > c = eng->next_compute();
        en += c->energy();
      }
      return en;
    };

};

}
}

#endif
