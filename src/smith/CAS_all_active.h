//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_all_active.h
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


#ifndef __SRC_SMITH_CAS_all_active_H 
#define __SRC_SMITH_CAS_all_active_H 

#include <src/smith/spinfreebase.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <iostream>
#include <iomanip>
#include <src/smith/queue.h>
#include <src/smith/CAS_all_active_tasks.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace CAS_all_active{

template <typename T>
class CAS_all_active : public SpinFreeMethod<T>, SMITH_info {
  protected:
    std::shared_ptr<Tensor<T> > t2;
    std::shared_ptr<Tensor<T> > r;
    double e0_;

    std::pair<std::shared_ptr<Queue<T> >, std::shared_ptr<Queue<T> > > make_queue_() {
      std::shared_ptr<Queue<T> > queue_(new Queue<T>());
      std::vector<IndexRange> index = {this->closed_, this->active_, this->virt_};

      std::vector<std::shared_ptr<Tensor<T> > > tensor0 = {r};
      std::shared_ptr<Task0<T> > task0(new Task0<T>(tensor0, index));
      queue_->add_task(task0);

      std::vector<IndexRange> Gamma0_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > Gamma0(new Tensor<T>(Gamma0_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor1 = {Gamma0, this->rdm3_, this->f1_};
      std::shared_ptr<Task1<T> > task1(new Task1<T>(tensor1, index));
      task1->add_dep(task0);
      queue_->add_task(task1);

      std::vector<IndexRange> Gamma6_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > Gamma6(new Tensor<T>(Gamma6_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor2 = {Gamma6, this->rdm2_};
      std::shared_ptr<Task2<T> > task2(new Task2<T>(tensor2, index));
      task2->add_dep(task0);
      queue_->add_task(task2);

      std::vector<IndexRange> Gamma5_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > Gamma5(new Tensor<T>(Gamma5_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor3 = {Gamma5, this->rdm2_, this->rdm3_};
      std::shared_ptr<Task3<T> > task3(new Task3<T>(tensor3, index));
      task3->add_dep(task0);
      queue_->add_task(task3);

      std::vector<IndexRange> Gamma2_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > Gamma2(new Tensor<T>(Gamma2_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor4 = {Gamma2, this->rdm2_, this->rdm3_, this->rdm4_, this->f1_};
      std::shared_ptr<Task4<T> > task4(new Task4<T>(tensor4, index));
      task4->add_dep(task0);
      queue_->add_task(task4);

      std::vector<IndexRange> Gamma3_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > Gamma3(new Tensor<T>(Gamma3_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor5 = {Gamma3, this->rdm1_, this->rdm2_, this->rdm3_};
      std::shared_ptr<Task5<T> > task5(new Task5<T>(tensor5, index));
      task5->add_dep(task0);
      queue_->add_task(task5);

      std::vector<IndexRange> Gamma4_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > Gamma4(new Tensor<T>(Gamma4_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor6 = {Gamma4, this->rdm2_, this->rdm3_};
      std::shared_ptr<Task6<T> > task6(new Task6<T>(tensor6, index));
      task6->add_dep(task0);
      queue_->add_task(task6);

      std::vector<IndexRange> I0_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor7 = {r, I0};
      std::shared_ptr<Task7<T> > task7(new Task7<T>(tensor7, index));
      task7->add_dep(task0);
      queue_->add_task(task7);


      std::vector<IndexRange> I1_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor8 = {I0, t2, I1};
      std::shared_ptr<Task8<T> > task8(new Task8<T>(tensor8, index));
      task7->add_dep(task8);
      task8->add_dep(task0);
      queue_->add_task(task8);


      std::vector<std::shared_ptr<Tensor<T> > > tensor9 = {I1, Gamma0};
      std::shared_ptr<Task9<T> > task9(new Task9<T>(tensor9, index));
      task8->add_dep(task9);
      task9->add_dep(task0);
      queue_->add_task(task9);

      task9->add_dep(task1);

      std::vector<IndexRange> I17_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I17(new Tensor<T>(I17_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor10 = {I0, t2, I17};
      std::shared_ptr<Task10<T> > task10(new Task10<T>(tensor10, index));
      task7->add_dep(task10);
      task10->add_dep(task0);
      queue_->add_task(task10);


      std::vector<std::shared_ptr<Tensor<T> > > tensor11 = {I17, Gamma6};
      std::shared_ptr<Task11<T> > task11(new Task11<T>(tensor11, index, this->e0_));
      task10->add_dep(task11);
      task11->add_dep(task0);
      queue_->add_task(task11);

      task11->add_dep(task2);

      std::vector<IndexRange> I19_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I19(new Tensor<T>(I19_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor12 = {I0, this->v2_, I19};
      std::shared_ptr<Task12<T> > task12(new Task12<T>(tensor12, index));
      task7->add_dep(task12);
      task12->add_dep(task0);
      queue_->add_task(task12);


      std::vector<std::shared_ptr<Tensor<T> > > tensor13 = {I19, Gamma6};
      std::shared_ptr<Task13<T> > task13(new Task13<T>(tensor13, index));
      task12->add_dep(task13);
      task13->add_dep(task0);
      queue_->add_task(task13);

      task13->add_dep(task2);

      std::vector<IndexRange> I2_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I2(new Tensor<T>(I2_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor14 = {r, I2};
      std::shared_ptr<Task14<T> > task14(new Task14<T>(tensor14, index));
      task14->add_dep(task0);
      queue_->add_task(task14);


      std::vector<IndexRange> I3_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor15 = {I2, this->f1_, I3};
      std::shared_ptr<Task15<T> > task15(new Task15<T>(tensor15, index));
      task14->add_dep(task15);
      task15->add_dep(task0);
      queue_->add_task(task15);


      std::vector<IndexRange> I4_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I4(new Tensor<T>(I4_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor16 = {I3, t2, I4};
      std::shared_ptr<Task16<T> > task16(new Task16<T>(tensor16, index));
      task15->add_dep(task16);
      task16->add_dep(task0);
      queue_->add_task(task16);


      std::vector<std::shared_ptr<Tensor<T> > > tensor17 = {I4, Gamma6};
      std::shared_ptr<Task17<T> > task17(new Task17<T>(tensor17, index));
      task16->add_dep(task17);
      task17->add_dep(task0);
      queue_->add_task(task17);

      task17->add_dep(task2);

      std::vector<IndexRange> I14_index = {this->active_, this->active_, this->active_, this->virt_};
      std::shared_ptr<Tensor<T> > I14(new Tensor<T>(I14_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor18 = {I2, this->f1_, I14};
      std::shared_ptr<Task18<T> > task18(new Task18<T>(tensor18, index));
      task14->add_dep(task18);
      task18->add_dep(task0);
      queue_->add_task(task18);


      std::vector<IndexRange> I15_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I15(new Tensor<T>(I15_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor19 = {I14, t2, I15};
      std::shared_ptr<Task19<T> > task19(new Task19<T>(tensor19, index));
      task18->add_dep(task19);
      task19->add_dep(task0);
      queue_->add_task(task19);


      std::vector<std::shared_ptr<Tensor<T> > > tensor20 = {I15, Gamma5};
      std::shared_ptr<Task20<T> > task20(new Task20<T>(tensor20, index));
      task19->add_dep(task20);
      task20->add_dep(task0);
      queue_->add_task(task20);

      task20->add_dep(task3);

      std::vector<IndexRange> I5_index = {this->active_, this->active_, this->active_, this->virt_};
      std::shared_ptr<Tensor<T> > I5(new Tensor<T>(I5_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor21 = {r, I5};
      std::shared_ptr<Task21<T> > task21(new Task21<T>(tensor21, index));
      task21->add_dep(task0);
      queue_->add_task(task21);


      std::vector<IndexRange> I6_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I6(new Tensor<T>(I6_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor22 = {I5, t2, I6};
      std::shared_ptr<Task22<T> > task22(new Task22<T>(tensor22, index));
      task21->add_dep(task22);
      task22->add_dep(task0);
      queue_->add_task(task22);


      std::vector<std::shared_ptr<Tensor<T> > > tensor23 = {I6, Gamma2};
      std::shared_ptr<Task23<T> > task23(new Task23<T>(tensor23, index));
      task22->add_dep(task23);
      task23->add_dep(task0);
      queue_->add_task(task23);

      task23->add_dep(task4);

      std::vector<IndexRange> I8_index = {this->active_, this->active_, this->active_, this->virt_};
      std::shared_ptr<Tensor<T> > I8(new Tensor<T>(I8_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor24 = {I5, this->f1_, I8};
      std::shared_ptr<Task24<T> > task24(new Task24<T>(tensor24, index));
      task21->add_dep(task24);
      task24->add_dep(task0);
      queue_->add_task(task24);


      std::vector<IndexRange> I9_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I9(new Tensor<T>(I9_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor25 = {I8, t2, I9};
      std::shared_ptr<Task25<T> > task25(new Task25<T>(tensor25, index));
      task24->add_dep(task25);
      task25->add_dep(task0);
      queue_->add_task(task25);


      std::vector<std::shared_ptr<Tensor<T> > > tensor26 = {I9, Gamma3};
      std::shared_ptr<Task26<T> > task26(new Task26<T>(tensor26, index));
      task25->add_dep(task26);
      task26->add_dep(task0);
      queue_->add_task(task26);

      task26->add_dep(task5);

      std::vector<IndexRange> I11_index = {this->active_, this->active_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I11(new Tensor<T>(I11_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor27 = {I5, this->f1_, I11};
      std::shared_ptr<Task27<T> > task27(new Task27<T>(tensor27, index));
      task21->add_dep(task27);
      task27->add_dep(task0);
      queue_->add_task(task27);


      std::vector<IndexRange> I12_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I12(new Tensor<T>(I12_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor28 = {I11, t2, I12};
      std::shared_ptr<Task28<T> > task28(new Task28<T>(tensor28, index));
      task27->add_dep(task28);
      task28->add_dep(task0);
      queue_->add_task(task28);


      std::vector<std::shared_ptr<Tensor<T> > > tensor29 = {I12, Gamma4};
      std::shared_ptr<Task29<T> > task29(new Task29<T>(tensor29, index));
      task28->add_dep(task29);
      task29->add_dep(task0);
      queue_->add_task(task29);

      task29->add_dep(task6);

      std::shared_ptr<Queue<T> > energy_(new Queue<T>());
      std::vector<IndexRange> I20_index;
      std::shared_ptr<Tensor<T> > I20(new Tensor<T>(I20_index, false));
      std::vector<IndexRange> I21_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I21(new Tensor<T>(I21_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor30 = {I20, t2, I21};
      std::shared_ptr<Task30<T> > task30(new Task30<T>(tensor30, index));
      energy_->add_task(task30);


      std::vector<IndexRange> I22_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I22(new Tensor<T>(I22_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor31 = {I21, this->v2_, I22};
      std::shared_ptr<Task31<T> > task31(new Task31<T>(tensor31, index));
      task30->add_dep(task31);
      energy_->add_task(task31);


      std::vector<std::shared_ptr<Tensor<T> > > tensor32 = {I22, Gamma6};
      std::shared_ptr<Task32<T> > task32(new Task32<T>(tensor32, index));
      task31->add_dep(task32);
      energy_->add_task(task32);

      task32->add_dep(task2);

      std::vector<IndexRange> I25_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I25(new Tensor<T>(I25_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor33 = {I21, r, I25};
      std::shared_ptr<Task33<T> > task33(new Task33<T>(tensor33, index));
      task30->add_dep(task33);
      energy_->add_task(task33);


      std::vector<std::shared_ptr<Tensor<T> > > tensor34 = {I25, Gamma6};
      std::shared_ptr<Task34<T> > task34(new Task34<T>(tensor34, index));
      task33->add_dep(task34);
      energy_->add_task(task34);

      task34->add_dep(task2);

      return make_pair(queue_, energy_);
    };

  public:
    CAS_all_active(std::shared_ptr<const Reference> ref) : SpinFreeMethod<T>(ref), SMITH_info() {
      this->eig_ = this->f1_->diag();
      t2 = this->v2_->clone();
      e0_ = this->compute_e0();
#if 1
      this->update_amplitude_start(t2, this->v2_);
      t2->scale(2.0);
#endif
      r = t2->clone();
    };
    ~CAS_all_active() {}; 

    void solve() {
      this->print_iteration();
      int iter = 0;
      for ( ; iter != maxiter_; ++iter) {
        std::pair<std::shared_ptr<Queue<T> >, std::shared_ptr<Queue<T> > > q = make_queue_();
        std::shared_ptr<Queue<T> > queue = q.first;
        std::shared_ptr<Queue<T> > energ = q.second;
        while (!queue->done())
          queue->next_compute();
        r->scale(0.25); // FIXME
//      *r = *(r->add_dagger());
        this->update_amplitude(t2, r);
        const double err = r->rms();
        r->zero();
        const double en = energy(energ);
        this->print_iteration(iter, en, err);
        if (err < thresh_residual()) break;
      }
      this->print_iteration(iter == maxiter_);
    };

    double energy(std::shared_ptr<Queue<T> > energ) {
      double en = 0.0;
      while (!energ->done()) {
        std::shared_ptr<Task<T> > c = energ->next_compute();
        en += c->energy() * 0.25; // FIXME
      }   
      return en; 
    };  
};

}
}
}
#endif

