//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test.h
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


#ifndef __SRC_SMITH_CAS_test_H 
#define __SRC_SMITH_CAS_test_H 

#include <src/smith/spinfreebase.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <iostream>
#include <iomanip>
#include <src/smith/queue.h>
#include <src/smith/CAS_test_tasks.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace CAS_test{

template <typename T>
class CAS_test : public SpinFreeMethod<T>, SMITH_info {
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

      std::vector<IndexRange> Gamma0_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > Gamma0(new Tensor<T>(Gamma0_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor1 = {Gamma0, this->rdm2_, this->f1_};
      std::shared_ptr<Task1<T> > task1(new Task1<T>(tensor1, index));
      task1->add_dep(task0);
      queue_->add_task(task1);

      std::vector<IndexRange> Gamma2_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > Gamma2(new Tensor<T>(Gamma2_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor2 = {Gamma2, this->rdm1_};
      std::shared_ptr<Task2<T> > task2(new Task2<T>(tensor2, index));
      task2->add_dep(task0);
      queue_->add_task(task2);

      std::vector<IndexRange> I0_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor3 = {r, I0};
      std::shared_ptr<Task3<T> > task3(new Task3<T>(tensor3, index));
      task3->add_dep(task0);
      queue_->add_task(task3);


      std::vector<IndexRange> I1_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor4 = {I0, t2, I1};
      std::shared_ptr<Task4<T> > task4(new Task4<T>(tensor4, index));
      task3->add_dep(task4);
      task4->add_dep(task0);
      queue_->add_task(task4);


      std::vector<std::shared_ptr<Tensor<T> > > tensor5 = {I1, Gamma0};
      std::shared_ptr<Task5<T> > task5(new Task5<T>(tensor5, index));
      task4->add_dep(task5);
      task5->add_dep(task0);
      queue_->add_task(task5);

      task5->add_dep(task1);

      std::vector<IndexRange> I3_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor6 = {I0, t2, I3};
      std::shared_ptr<Task6<T> > task6(new Task6<T>(tensor6, index));
      task3->add_dep(task6);
      task6->add_dep(task0);
      queue_->add_task(task6);


      std::vector<std::shared_ptr<Tensor<T> > > tensor7 = {I3, Gamma0};
      std::shared_ptr<Task7<T> > task7(new Task7<T>(tensor7, index));
      task6->add_dep(task7);
      task7->add_dep(task0);
      queue_->add_task(task7);

      task7->add_dep(task1);

      std::vector<IndexRange> I5_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I5(new Tensor<T>(I5_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor8 = {I0, this->f1_, I5};
      std::shared_ptr<Task8<T> > task8(new Task8<T>(tensor8, index));
      task3->add_dep(task8);
      task8->add_dep(task0);
      queue_->add_task(task8);


      std::vector<IndexRange> I6_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I6(new Tensor<T>(I6_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor9 = {I5, t2, I6};
      std::shared_ptr<Task9<T> > task9(new Task9<T>(tensor9, index));
      task8->add_dep(task9);
      task9->add_dep(task0);
      queue_->add_task(task9);


      std::vector<std::shared_ptr<Tensor<T> > > tensor10 = {I6, Gamma2};
      std::shared_ptr<Task10<T> > task10(new Task10<T>(tensor10, index));
      task9->add_dep(task10);
      task10->add_dep(task0);
      queue_->add_task(task10);

      task10->add_dep(task2);

      std::vector<IndexRange> I9_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I9(new Tensor<T>(I9_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor11 = {I5, t2, I9};
      std::shared_ptr<Task11<T> > task11(new Task11<T>(tensor11, index));
      task8->add_dep(task11);
      task11->add_dep(task0);
      queue_->add_task(task11);


      std::vector<std::shared_ptr<Tensor<T> > > tensor12 = {I9, Gamma2};
      std::shared_ptr<Task12<T> > task12(new Task12<T>(tensor12, index));
      task11->add_dep(task12);
      task12->add_dep(task0);
      queue_->add_task(task12);

      task12->add_dep(task2);

      std::vector<IndexRange> I11_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I11(new Tensor<T>(I11_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor13 = {I0, this->f1_, I11};
      std::shared_ptr<Task13<T> > task13(new Task13<T>(tensor13, index));
      task3->add_dep(task13);
      task13->add_dep(task0);
      queue_->add_task(task13);


      std::vector<IndexRange> I12_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I12(new Tensor<T>(I12_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor14 = {I11, t2, I12};
      std::shared_ptr<Task14<T> > task14(new Task14<T>(tensor14, index));
      task13->add_dep(task14);
      task14->add_dep(task0);
      queue_->add_task(task14);


      std::vector<std::shared_ptr<Tensor<T> > > tensor15 = {I12, Gamma2};
      std::shared_ptr<Task15<T> > task15(new Task15<T>(tensor15, index));
      task14->add_dep(task15);
      task15->add_dep(task0);
      queue_->add_task(task15);

      task15->add_dep(task2);

      std::vector<IndexRange> I15_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I15(new Tensor<T>(I15_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor16 = {I11, t2, I15};
      std::shared_ptr<Task16<T> > task16(new Task16<T>(tensor16, index));
      task13->add_dep(task16);
      task16->add_dep(task0);
      queue_->add_task(task16);


      std::vector<std::shared_ptr<Tensor<T> > > tensor17 = {I15, Gamma2};
      std::shared_ptr<Task17<T> > task17(new Task17<T>(tensor17, index));
      task16->add_dep(task17);
      task17->add_dep(task0);
      queue_->add_task(task17);

      task17->add_dep(task2);

      std::vector<IndexRange> I17_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T> > I17(new Tensor<T>(I17_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor18 = {I0, this->f1_, I17};
      std::shared_ptr<Task18<T> > task18(new Task18<T>(tensor18, index));
      task3->add_dep(task18);
      task18->add_dep(task0);
      queue_->add_task(task18);


      std::vector<IndexRange> I18_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I18(new Tensor<T>(I18_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor19 = {I17, t2, I18};
      std::shared_ptr<Task19<T> > task19(new Task19<T>(tensor19, index));
      task18->add_dep(task19);
      task19->add_dep(task0);
      queue_->add_task(task19);


      std::vector<std::shared_ptr<Tensor<T> > > tensor20 = {I18, Gamma2};
      std::shared_ptr<Task20<T> > task20(new Task20<T>(tensor20, index));
      task19->add_dep(task20);
      task20->add_dep(task0);
      queue_->add_task(task20);

      task20->add_dep(task2);

      std::vector<IndexRange> I21_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I21(new Tensor<T>(I21_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor21 = {I17, t2, I21};
      std::shared_ptr<Task21<T> > task21(new Task21<T>(tensor21, index));
      task18->add_dep(task21);
      task21->add_dep(task0);
      queue_->add_task(task21);


      std::vector<std::shared_ptr<Tensor<T> > > tensor22 = {I21, Gamma2};
      std::shared_ptr<Task22<T> > task22(new Task22<T>(tensor22, index));
      task21->add_dep(task22);
      task22->add_dep(task0);
      queue_->add_task(task22);

      task22->add_dep(task2);

      std::vector<IndexRange> I23_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I23(new Tensor<T>(I23_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor23 = {I0, t2, I23};
      std::shared_ptr<Task23<T> > task23(new Task23<T>(tensor23, index));
      task3->add_dep(task23);
      task23->add_dep(task0);
      queue_->add_task(task23);


      std::vector<std::shared_ptr<Tensor<T> > > tensor24 = {I23, Gamma2};
      std::shared_ptr<Task24<T> > task24(new Task24<T>(tensor24, index, this->e0_));
      task23->add_dep(task24);
      task24->add_dep(task0);
      queue_->add_task(task24);

      task24->add_dep(task2);

      std::vector<IndexRange> I25_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I25(new Tensor<T>(I25_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor25 = {I0, t2, I25};
      std::shared_ptr<Task25<T> > task25(new Task25<T>(tensor25, index));
      task3->add_dep(task25);
      task25->add_dep(task0);
      queue_->add_task(task25);


      std::vector<std::shared_ptr<Tensor<T> > > tensor26 = {I25, Gamma2};
      std::shared_ptr<Task26<T> > task26(new Task26<T>(tensor26, index, this->e0_));
      task25->add_dep(task26);
      task26->add_dep(task0);
      queue_->add_task(task26);

      task26->add_dep(task2);

      std::vector<IndexRange> I27_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I27(new Tensor<T>(I27_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor27 = {I0, this->v2_, I27};
      std::shared_ptr<Task27<T> > task27(new Task27<T>(tensor27, index));
      task3->add_dep(task27);
      task27->add_dep(task0);
      queue_->add_task(task27);


      std::vector<std::shared_ptr<Tensor<T> > > tensor28 = {I27, Gamma2};
      std::shared_ptr<Task28<T> > task28(new Task28<T>(tensor28, index));
      task27->add_dep(task28);
      task28->add_dep(task0);
      queue_->add_task(task28);

      task28->add_dep(task2);

      std::vector<IndexRange> I29_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I29(new Tensor<T>(I29_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor29 = {I0, this->v2_, I29};
      std::shared_ptr<Task29<T> > task29(new Task29<T>(tensor29, index));
      task3->add_dep(task29);
      task29->add_dep(task0);
      queue_->add_task(task29);


      std::vector<std::shared_ptr<Tensor<T> > > tensor30 = {I29, Gamma2};
      std::shared_ptr<Task30<T> > task30(new Task30<T>(tensor30, index));
      task29->add_dep(task30);
      task30->add_dep(task0);
      queue_->add_task(task30);

      task30->add_dep(task2);

      std::shared_ptr<Queue<T> > energy_(new Queue<T>());
      std::vector<IndexRange> I30_index;
      std::shared_ptr<Tensor<T> > I30(new Tensor<T>(I30_index, false));
      std::vector<IndexRange> I31_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T> > I31(new Tensor<T>(I31_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor31 = {I30, t2, I31};
      std::shared_ptr<Task31<T> > task31(new Task31<T>(tensor31, index));
      energy_->add_task(task31);


      std::vector<IndexRange> I32_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I32(new Tensor<T>(I32_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor32 = {I31, this->v2_, I32};
      std::shared_ptr<Task32<T> > task32(new Task32<T>(tensor32, index));
      task31->add_dep(task32);
      energy_->add_task(task32);


      std::vector<std::shared_ptr<Tensor<T> > > tensor33 = {I32, Gamma2};
      std::shared_ptr<Task33<T> > task33(new Task33<T>(tensor33, index));
      task32->add_dep(task33);
      energy_->add_task(task33);

      task33->add_dep(task2);

      std::vector<IndexRange> I35_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I35(new Tensor<T>(I35_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor34 = {I31, this->v2_, I35};
      std::shared_ptr<Task34<T> > task34(new Task34<T>(tensor34, index));
      task31->add_dep(task34);
      energy_->add_task(task34);


      std::vector<std::shared_ptr<Tensor<T> > > tensor35 = {I35, Gamma2};
      std::shared_ptr<Task35<T> > task35(new Task35<T>(tensor35, index));
      task34->add_dep(task35);
      energy_->add_task(task35);

      task35->add_dep(task2);

      std::vector<IndexRange> I38_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I38(new Tensor<T>(I38_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor36 = {I31, r, I38};
      std::shared_ptr<Task36<T> > task36(new Task36<T>(tensor36, index));
      task31->add_dep(task36);
      energy_->add_task(task36);


      std::vector<std::shared_ptr<Tensor<T> > > tensor37 = {I38, Gamma2};
      std::shared_ptr<Task37<T> > task37(new Task37<T>(tensor37, index));
      task36->add_dep(task37);
      energy_->add_task(task37);

      task37->add_dep(task2);

      std::vector<IndexRange> I41_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T> > I41(new Tensor<T>(I41_index, false));
      std::vector<std::shared_ptr<Tensor<T> > > tensor38 = {I31, r, I41};
      std::shared_ptr<Task38<T> > task38(new Task38<T>(tensor38, index));
      task31->add_dep(task38);
      energy_->add_task(task38);


      std::vector<std::shared_ptr<Tensor<T> > > tensor39 = {I41, Gamma2};
      std::shared_ptr<Task39<T> > task39(new Task39<T>(tensor39, index));
      task38->add_dep(task39);
      energy_->add_task(task39);

      task39->add_dep(task2);

      return make_pair(queue_, energy_);
    };

  public:
    CAS_test(std::shared_ptr<const Reference> ref) : SpinFreeMethod<T>(ref), SMITH_info() {
      this->eig_ = this->f1_->diag();
      t2 = this->v2_->clone();
      e0_ = this->e0();
#if 1
      this->update_amplitude(t2, this->v2_, true);
      t2->scale(2.0);
#endif
      r = t2->clone();
    };
    ~CAS_test() {}; 

    void solve() {
      this->print_iteration();
      int iter = 0;
      for ( ; iter != maxiter_; ++iter) {
        std::pair<std::shared_ptr<Queue<T> >, std::shared_ptr<Queue<T> > > q = make_queue_();
        std::shared_ptr<Queue<T> > queue = q.first;
        std::shared_ptr<Queue<T> > energ = q.second;
        while (!queue->done())
          queue->next_compute();
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

