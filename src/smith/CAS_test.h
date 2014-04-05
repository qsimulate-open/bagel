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
#include <tuple>
#include <iomanip>
#include <src/smith/queue.h>
#include <src/smith/CAS_test_tasks.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace CAS_test{

template <typename T>
class CAS_test : public SpinFreeMethod<T>{
  protected:
    using SpinFreeMethod<T>::ref_;

    std::shared_ptr<Tensor<T>> t2;
    std::shared_ptr<Tensor<T>> r;
    double e0_;
    std::shared_ptr<Tensor<T>> sigma_;
    std::shared_ptr<Tensor<T>> den1;
    std::shared_ptr<Tensor<T>> den2;
    double correlated_norm;
    std::shared_ptr<Tensor<T>> deci;

    std::tuple<std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>> make_queue_() {
      std::shared_ptr<Queue<T>> queue_(new Queue<T>());
      std::array<std::shared_ptr<const IndexRange>,3> pindex = {{this->rclosed_, this->ractive_, this->rvirt_}};
      std::array<std::shared_ptr<const IndexRange>,4> cindex = {{this->rclosed_, this->ractive_, this->rvirt_, this->rci_}};

      std::vector<std::shared_ptr<Tensor<T>>> tensor0 = {r};
      std::shared_ptr<Task0<T>> task0(new Task0<T>(tensor0));
      queue_->add_task(task0);

      std::vector<IndexRange> Gamma0_index;
      std::shared_ptr<Tensor<T>> Gamma0(new Tensor<T>(Gamma0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor1 = {Gamma0, this->rdm1_, this->f1_};
      std::shared_ptr<Task1<T>> task1(new Task1<T>(tensor1, pindex));
      task1->add_dep(task0);
      queue_->add_task(task1);

      std::vector<IndexRange> Gamma2_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma2(new Tensor<T>(Gamma2_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor2 = {Gamma2, this->rdm1_};
      std::shared_ptr<Task2<T>> task2(new Task2<T>(tensor2, pindex));
      task2->add_dep(task0);
      queue_->add_task(task2);

      std::vector<IndexRange> Gamma4_index = {this->ci_};
      std::shared_ptr<Tensor<T>> Gamma4(new Tensor<T>(Gamma4_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor3 = {Gamma4, this->rdm1deriv_, this->f1_};
      std::shared_ptr<Task3<T>> task3(new Task3<T>(tensor3, cindex));
      task3->add_dep(task0);
      queue_->add_task(task3);

      std::vector<IndexRange> I0_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor4 = {r, I0};
      std::shared_ptr<Task4<T>> task4(new Task4<T>(tensor4, pindex));
      task4->add_dep(task0);
      queue_->add_task(task4);


      std::vector<std::shared_ptr<Tensor<T>>> tensor5 = {I0, t2, this->v2_};
      std::shared_ptr<Task5<T>> task5(new Task5<T>(tensor5, pindex, this->e0_));
      task4->add_dep(task5);
      task5->add_dep(task0);
      queue_->add_task(task5);


      std::vector<IndexRange> I1_index;
      std::shared_ptr<Tensor<T>> I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor6 = {I0, t2, I1};
      std::shared_ptr<Task6<T>> task6(new Task6<T>(tensor6, pindex));
      task4->add_dep(task6);
      task6->add_dep(task0);
      queue_->add_task(task6);


      std::vector<std::shared_ptr<Tensor<T>>> tensor7 = {I1, Gamma0};
      std::shared_ptr<Task7<T>> task7(new Task7<T>(tensor7, pindex));
      task6->add_dep(task7);
      task7->add_dep(task0);
      queue_->add_task(task7);

      task7->add_dep(task1);

      std::vector<IndexRange> I3_index;
      std::shared_ptr<Tensor<T>> I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor8 = {I0, t2, I3};
      std::shared_ptr<Task8<T>> task8(new Task8<T>(tensor8, pindex));
      task4->add_dep(task8);
      task8->add_dep(task0);
      queue_->add_task(task8);


      std::vector<std::shared_ptr<Tensor<T>>> tensor9 = {I3, Gamma0};
      std::shared_ptr<Task9<T>> task9(new Task9<T>(tensor9, pindex));
      task8->add_dep(task9);
      task9->add_dep(task0);
      queue_->add_task(task9);

      task9->add_dep(task1);

      std::vector<IndexRange> I4_index = {this->closed_, this->virt_, this->virt_, this->closed_};
      std::shared_ptr<Tensor<T>> I4(new Tensor<T>(I4_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor10 = {r, I4};
      std::shared_ptr<Task10<T>> task10(new Task10<T>(tensor10, pindex));
      task10->add_dep(task0);
      queue_->add_task(task10);


      std::vector<IndexRange> I5_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I5(new Tensor<T>(I5_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor11 = {I4, this->f1_, I5};
      std::shared_ptr<Task11<T>> task11(new Task11<T>(tensor11, pindex));
      task10->add_dep(task11);
      task11->add_dep(task0);
      queue_->add_task(task11);


      std::vector<std::shared_ptr<Tensor<T>>> tensor12 = {I5, t2};
      std::shared_ptr<Task12<T>> task12(new Task12<T>(tensor12, pindex));
      task11->add_dep(task12);
      task12->add_dep(task0);
      queue_->add_task(task12);


      std::vector<IndexRange> I9_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I9(new Tensor<T>(I9_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor13 = {I4, this->f1_, I9};
      std::shared_ptr<Task13<T>> task13(new Task13<T>(tensor13, pindex));
      task10->add_dep(task13);
      task13->add_dep(task0);
      queue_->add_task(task13);


      std::vector<std::shared_ptr<Tensor<T>>> tensor14 = {I9, t2};
      std::shared_ptr<Task14<T>> task14(new Task14<T>(tensor14, pindex));
      task13->add_dep(task14);
      task14->add_dep(task0);
      queue_->add_task(task14);


      std::shared_ptr<Queue<T>> energy_(new Queue<T>());
      std::vector<IndexRange> I16_index;
      std::shared_ptr<Tensor<T>> I16(new Tensor<T>(I16_index, false));
      std::vector<IndexRange> I17_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I17(new Tensor<T>(I17_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor15 = {I16, t2, I17};
      std::shared_ptr<Task15<T>> task15(new Task15<T>(tensor15, pindex));
      energy_->add_task(task15);


      std::vector<std::shared_ptr<Tensor<T>>> tensor16 = {I17, this->v2_};
      std::shared_ptr<Task16<T>> task16(new Task16<T>(tensor16, pindex));
      task15->add_dep(task16);
      energy_->add_task(task16);


      std::shared_ptr<Queue<T>> correction_(new Queue<T>());
      std::vector<IndexRange> I20_index;
      std::shared_ptr<Tensor<T>> I20(new Tensor<T>(I20_index, false));
      std::vector<IndexRange> I21_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I21(new Tensor<T>(I21_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor17 = {I20, t2, I21};
      std::shared_ptr<Task17<T>> task17(new Task17<T>(tensor17, pindex));
      correction_->add_task(task17);


      std::vector<std::shared_ptr<Tensor<T>>> tensor18 = {I21, t2};
      std::shared_ptr<Task18<T>> task18(new Task18<T>(tensor18, pindex));
      task17->add_dep(task18);
      correction_->add_task(task18);


      std::shared_ptr<Queue<T>> density_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor19 = {den1};
      std::shared_ptr<Task19<T>> task19(new Task19<T>(tensor19));
      density_->add_task(task19);

      std::vector<IndexRange> I24_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I24(new Tensor<T>(I24_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor20 = {den1, I24};
      std::shared_ptr<Task20<T>> task20(new Task20<T>(tensor20, pindex));
      task20->add_dep(task19);
      density_->add_task(task20);


      std::vector<IndexRange> I25_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I25(new Tensor<T>(I25_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor21 = {I24, t2, I25};
      std::shared_ptr<Task21<T>> task21(new Task21<T>(tensor21, pindex));
      task20->add_dep(task21);
      task21->add_dep(task19);
      density_->add_task(task21);


      std::vector<IndexRange> I26_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I26(new Tensor<T>(I26_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor22 = {I25, t2, I26};
      std::shared_ptr<Task22<T>> task22(new Task22<T>(tensor22, pindex));
      task21->add_dep(task22);
      task22->add_dep(task19);
      density_->add_task(task22);


      std::vector<std::shared_ptr<Tensor<T>>> tensor23 = {I26, Gamma2};
      std::shared_ptr<Task23<T>> task23(new Task23<T>(tensor23, pindex));
      task22->add_dep(task23);
      task23->add_dep(task19);
      density_->add_task(task23);

      task23->add_dep(task2);

      std::vector<IndexRange> I29_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I29(new Tensor<T>(I29_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor24 = {I25, t2, I29};
      std::shared_ptr<Task24<T>> task24(new Task24<T>(tensor24, pindex));
      task21->add_dep(task24);
      task24->add_dep(task19);
      density_->add_task(task24);


      std::vector<std::shared_ptr<Tensor<T>>> tensor25 = {I29, Gamma2};
      std::shared_ptr<Task25<T>> task25(new Task25<T>(tensor25, pindex));
      task24->add_dep(task25);
      task25->add_dep(task19);
      density_->add_task(task25);

      task25->add_dep(task2);

      std::vector<IndexRange> I30_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor<T>> I30(new Tensor<T>(I30_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor26 = {den1, I30};
      std::shared_ptr<Task26<T>> task26(new Task26<T>(tensor26, pindex));
      task26->add_dep(task19);
      density_->add_task(task26);


      std::vector<IndexRange> I31_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I31(new Tensor<T>(I31_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor27 = {I30, t2, I31};
      std::shared_ptr<Task27<T>> task27(new Task27<T>(tensor27, pindex));
      task26->add_dep(task27);
      task27->add_dep(task19);
      density_->add_task(task27);


      std::vector<std::shared_ptr<Tensor<T>>> tensor28 = {I31, t2};
      std::shared_ptr<Task28<T>> task28(new Task28<T>(tensor28, pindex));
      task27->add_dep(task28);
      task28->add_dep(task19);
      density_->add_task(task28);


      std::vector<IndexRange> I34_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I34(new Tensor<T>(I34_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor29 = {den1, I34};
      std::shared_ptr<Task29<T>> task29(new Task29<T>(tensor29, pindex));
      task29->add_dep(task19);
      density_->add_task(task29);


      std::vector<IndexRange> I35_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I35(new Tensor<T>(I35_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor30 = {I34, t2, I35};
      std::shared_ptr<Task30<T>> task30(new Task30<T>(tensor30, pindex));
      task29->add_dep(task30);
      task30->add_dep(task19);
      density_->add_task(task30);


      std::vector<std::shared_ptr<Tensor<T>>> tensor31 = {I35, t2};
      std::shared_ptr<Task31<T>> task31(new Task31<T>(tensor31, pindex));
      task30->add_dep(task31);
      task31->add_dep(task19);
      density_->add_task(task31);


      std::shared_ptr<Queue<T>> density2_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor32 = {den2};
      std::shared_ptr<Task32<T>> task32(new Task32<T>(tensor32));
      density2_->add_task(task32);

      std::vector<IndexRange> I38_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I38(new Tensor<T>(I38_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor33 = {den2, I38};
      std::shared_ptr<Task33<T>> task33(new Task33<T>(tensor33, pindex));
      task33->add_dep(task32);
      density2_->add_task(task33);


      std::vector<std::shared_ptr<Tensor<T>>> tensor34 = {I38, t2};
      std::shared_ptr<Task34<T>> task34(new Task34<T>(tensor34, pindex));
      task33->add_dep(task34);
      task34->add_dep(task32);
      density2_->add_task(task34);


      std::shared_ptr<Queue<T>> dedci_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor35 = {deci};
      std::shared_ptr<Task35<T>> task35(new Task35<T>(tensor35));
      dedci_->add_task(task35);

      std::vector<IndexRange> I40_index = {this->ci_};
      std::shared_ptr<Tensor<T>> I40(new Tensor<T>(I40_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor36 = {deci, I40};
      std::shared_ptr<Task36<T>> task36(new Task36<T>(tensor36, cindex));
      task36->add_dep(task35);
      dedci_->add_task(task36);


      std::vector<IndexRange> I41_index = {this->ci_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I41(new Tensor<T>(I41_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor37 = {I40, t2, I41};
      std::shared_ptr<Task37<T>> task37(new Task37<T>(tensor37, cindex));
      task36->add_dep(task37);
      task37->add_dep(task35);
      dedci_->add_task(task37);


      std::vector<IndexRange> I42_index = {this->ci_};
      std::shared_ptr<Tensor<T>> I42(new Tensor<T>(I42_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor38 = {I41, t2, I42};
      std::shared_ptr<Task38<T>> task38(new Task38<T>(tensor38, cindex));
      task37->add_dep(task38);
      task38->add_dep(task35);
      dedci_->add_task(task38);


      std::vector<std::shared_ptr<Tensor<T>>> tensor39 = {I42, Gamma4};
      std::shared_ptr<Task39<T>> task39(new Task39<T>(tensor39, cindex));
      task38->add_dep(task39);
      task39->add_dep(task35);
      dedci_->add_task(task39);

      task39->add_dep(task3);
      task39->add_dep(task3);

      std::vector<IndexRange> I45_index = {this->ci_};
      std::shared_ptr<Tensor<T>> I45(new Tensor<T>(I45_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor40 = {I41, t2, I45};
      std::shared_ptr<Task40<T>> task40(new Task40<T>(tensor40, cindex));
      task37->add_dep(task40);
      task40->add_dep(task35);
      dedci_->add_task(task40);


      std::vector<std::shared_ptr<Tensor<T>>> tensor41 = {I45, Gamma4};
      std::shared_ptr<Task41<T>> task41(new Task41<T>(tensor41, cindex));
      task40->add_dep(task41);
      task41->add_dep(task35);
      dedci_->add_task(task41);

      task41->add_dep(task3);
      task41->add_dep(task3);

      return make_tuple(queue_, energy_, correction_, density_, density2_, dedci_);
    };

  public:
    CAS_test(std::shared_ptr<const SMITH_Info> ref) : SpinFreeMethod<T>(ref) {
      this->eig_ = this->f1_->diag();
      t2 = this->v2_->clone();
      e0_ = this->e0();
      sigma_ = this->sigma();
      this->update_amplitude(t2, this->v2_, true);
      t2->scale(2.0);
      r = t2->clone();
      den1 = this->h1_->clone();
      den2 = this->v2_->clone();
      deci = this->rdm0deriv_->clone();
    };
    ~CAS_test() {};

    void solve() {
      this->print_iteration();
      int iter = 0;
      std::shared_ptr<Queue<T>> queue, energ, correct, dens, dens2, dec;
      for ( ; iter != ref_->maxiter(); ++iter) {
        std::tie(queue, energ, correct, dens, dens2, dec) = make_queue_();
        while (!queue->done())
          queue->next_compute();
        this->update_amplitude(t2, r);
        const double err = r->rms();
        r->zero();
        this->energy_ = energy(energ);
        this->print_iteration(iter, this->energy_, err);
        if (err < ref_->thresh()) break;
      }
      this->print_iteration(iter == ref_->maxiter());

      correlated_norm = correction(correct);
      std::cout << "Norm, correlated overlap: <1|1> = " << std::setprecision(10) << correlated_norm << std::endl;

      std::cout << " === Computing unrelaxed density matrix, dm1, <1|E_pq|1> + 2<0|E_pq|1> ===" << std::endl;
      while (!dens->done())
        dens->next_compute();
      std::cout << " === Computing unrelaxed density matrix, dm2, <0|E_pqrs|1>  ===" << std::endl;
      while (!dens2->done())
        dens2->next_compute();
      std::cout << " === Calculating cI derivative dE/dcI ===" << std::endl;
      while (!dec->done())
        dec->next_compute();
      deci->ax_plus_y(-correlated_norm, sigma_);
      std::cout << std::endl;

    };

    double energy(std::shared_ptr<Queue<T>> energ) {
      double en = 0.0;
      while (!energ->done()) {
        std::shared_ptr<Task<T>> c = energ->next_compute();
        en += c->energy();
      }
      return en;
    }

    double correction(std::shared_ptr<Queue<T>> correct) {
      double n = 0.0;
      while (!correct->done()) {
        std::shared_ptr<Task<T>> c = correct->next_compute();
        n += c->correction();
      }
      return n;
    }

    std::shared_ptr<const Matrix> rdm1() const { return den1->matrix(); }
    std::shared_ptr<const Matrix> rdm2() const { return den2->matrix2(); }

    double rdm1_correction() const { return correlated_norm; }

    std::shared_ptr<const Civec> ci_deriv() const { return deci->civec(this->det_); }

};

}
}
}
#endif

