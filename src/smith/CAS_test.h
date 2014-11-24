//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test.h
// Copyright (C) 2014 Shiozaki group
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
    std::shared_ptr<Tensor<T>> Den1;
    double correlated_norm;
    std::shared_ptr<Tensor<T>> deci;

    std::tuple<std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>,  std::shared_ptr<Queue<T>>,  std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>> make_queue_() {
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

      std::vector<IndexRange> Gamma6_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma6(new Tensor<T>(Gamma6_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor3 = {Gamma6, this->rdm2_, this->f1_};
      std::shared_ptr<Task3<T>> task3(new Task3<T>(tensor3, pindex));
      task3->add_dep(task0);
      queue_->add_task(task3);

      std::vector<IndexRange> Gamma14_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma14(new Tensor<T>(Gamma14_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor4 = {Gamma14, this->rdm2_};
      std::shared_ptr<Task4<T>> task4(new Task4<T>(tensor4, pindex));
      task4->add_dep(task0);
      queue_->add_task(task4);

      std::vector<IndexRange> Gamma16_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma16(new Tensor<T>(Gamma16_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor5 = {Gamma16, this->rdm3_, this->f1_};
      std::shared_ptr<Task5<T>> task5(new Task5<T>(tensor5, pindex));
      task5->add_dep(task0);
      queue_->add_task(task5);

      std::vector<IndexRange> Gamma67_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma67(new Tensor<T>(Gamma67_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor6 = {Gamma67, this->rdm3_};
      std::shared_ptr<Task6<T>> task6(new Task6<T>(tensor6, pindex));
      task6->add_dep(task0);
      queue_->add_task(task6);

      std::vector<IndexRange> Gamma72_index = {this->ci_};
      std::shared_ptr<Tensor<T>> Gamma72(new Tensor<T>(Gamma72_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor7 = {Gamma72, this->rdm1deriv_, this->f1_};
      std::shared_ptr<Task7<T>> task7(new Task7<T>(tensor7, cindex));
      task7->add_dep(task0);
      queue_->add_task(task7);

      std::vector<IndexRange> Gamma74_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma74(new Tensor<T>(Gamma74_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor8 = {Gamma74, this->rdm1deriv_};
      std::shared_ptr<Task8<T>> task8(new Task8<T>(tensor8, cindex));
      task8->add_dep(task0);
      queue_->add_task(task8);

      std::vector<IndexRange> Gamma78_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma78(new Tensor<T>(Gamma78_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor9 = {Gamma78, this->rdm2deriv_, this->f1_};
      std::shared_ptr<Task9<T>> task9(new Task9<T>(tensor9, cindex));
      task9->add_dep(task0);
      queue_->add_task(task9);

      std::vector<IndexRange> Gamma86_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma86(new Tensor<T>(Gamma86_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor10 = {Gamma86, this->rdm2deriv_};
      std::shared_ptr<Task10<T>> task10(new Task10<T>(tensor10, cindex));
      task10->add_dep(task0);
      queue_->add_task(task10);

      std::vector<IndexRange> Gamma88_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma88(new Tensor<T>(Gamma88_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor11 = {Gamma88, this->rdm3deriv_, this->f1_};
      std::shared_ptr<Task11<T>> task11(new Task11<T>(tensor11, cindex));
      task11->add_dep(task0);
      queue_->add_task(task11);

      std::vector<IndexRange> I0_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor12 = {r, I0};
      std::shared_ptr<Task12<T>> task12(new Task12<T>(tensor12, pindex));
      task12->add_dep(task0);
      queue_->add_task(task12);


      std::vector<std::shared_ptr<Tensor<T>>> tensor13 = {I0, t2, this->v2_};
      std::shared_ptr<Task13<T>> task13(new Task13<T>(tensor13, pindex, this->e0_));
      task12->add_dep(task13);
      task13->add_dep(task0);
      queue_->add_task(task13);


      std::vector<IndexRange> I1_index;
      std::shared_ptr<Tensor<T>> I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor14 = {I0, t2, I1};
      std::shared_ptr<Task14<T>> task14(new Task14<T>(tensor14, pindex));
      task12->add_dep(task14);
      task14->add_dep(task0);
      queue_->add_task(task14);


      std::vector<std::shared_ptr<Tensor<T>>> tensor15 = {I1, Gamma0};
      std::shared_ptr<Task15<T>> task15(new Task15<T>(tensor15, pindex));
      task14->add_dep(task15);
      task15->add_dep(task0);
      queue_->add_task(task15);

      task15->add_dep(task1);

      std::vector<IndexRange> I3_index;
      std::shared_ptr<Tensor<T>> I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor16 = {I0, t2, I3};
      std::shared_ptr<Task16<T>> task16(new Task16<T>(tensor16, pindex));
      task12->add_dep(task16);
      task16->add_dep(task0);
      queue_->add_task(task16);


      std::vector<std::shared_ptr<Tensor<T>>> tensor17 = {I3, Gamma0};
      std::shared_ptr<Task17<T>> task17(new Task17<T>(tensor17, pindex));
      task16->add_dep(task17);
      task17->add_dep(task0);
      queue_->add_task(task17);

      task17->add_dep(task1);

      std::vector<IndexRange> I4_index = {this->closed_, this->virt_, this->virt_, this->closed_};
      std::shared_ptr<Tensor<T>> I4(new Tensor<T>(I4_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor18 = {r, I4};
      std::shared_ptr<Task18<T>> task18(new Task18<T>(tensor18, pindex));
      task18->add_dep(task0);
      queue_->add_task(task18);


      std::vector<IndexRange> I5_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I5(new Tensor<T>(I5_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor19 = {I4, this->f1_, I5};
      std::shared_ptr<Task19<T>> task19(new Task19<T>(tensor19, pindex));
      task18->add_dep(task19);
      task19->add_dep(task0);
      queue_->add_task(task19);


      std::vector<std::shared_ptr<Tensor<T>>> tensor20 = {I5, t2};
      std::shared_ptr<Task20<T>> task20(new Task20<T>(tensor20, pindex));
      task19->add_dep(task20);
      task20->add_dep(task0);
      queue_->add_task(task20);


      std::vector<IndexRange> I9_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I9(new Tensor<T>(I9_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor21 = {I4, this->f1_, I9};
      std::shared_ptr<Task21<T>> task21(new Task21<T>(tensor21, pindex));
      task18->add_dep(task21);
      task21->add_dep(task0);
      queue_->add_task(task21);


      std::vector<std::shared_ptr<Tensor<T>>> tensor22 = {I9, t2};
      std::shared_ptr<Task22<T>> task22(new Task22<T>(tensor22, pindex));
      task21->add_dep(task22);
      task22->add_dep(task0);
      queue_->add_task(task22);


      std::vector<IndexRange> I13_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I13(new Tensor<T>(I13_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor23 = {I4, this->f1_, I13};
      std::shared_ptr<Task23<T>> task23(new Task23<T>(tensor23, pindex));
      task18->add_dep(task23);
      task23->add_dep(task0);
      queue_->add_task(task23);


      std::vector<IndexRange> I14_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I14(new Tensor<T>(I14_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor24 = {I13, t2, I14};
      std::shared_ptr<Task24<T>> task24(new Task24<T>(tensor24, pindex));
      task23->add_dep(task24);
      task24->add_dep(task0);
      queue_->add_task(task24);


      std::vector<std::shared_ptr<Tensor<T>>> tensor25 = {I14, Gamma2};
      std::shared_ptr<Task25<T>> task25(new Task25<T>(tensor25, pindex));
      task24->add_dep(task25);
      task25->add_dep(task0);
      queue_->add_task(task25);

      task25->add_dep(task2);

      std::vector<IndexRange> I17_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I17(new Tensor<T>(I17_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor26 = {I13, t2, I17};
      std::shared_ptr<Task26<T>> task26(new Task26<T>(tensor26, pindex));
      task23->add_dep(task26);
      task26->add_dep(task0);
      queue_->add_task(task26);


      std::vector<std::shared_ptr<Tensor<T>>> tensor27 = {I17, Gamma2};
      std::shared_ptr<Task27<T>> task27(new Task27<T>(tensor27, pindex));
      task26->add_dep(task27);
      task27->add_dep(task0);
      queue_->add_task(task27);

      task27->add_dep(task2);

      std::vector<IndexRange> I18_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I18(new Tensor<T>(I18_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor28 = {r, I18};
      std::shared_ptr<Task28<T>> task28(new Task28<T>(tensor28, pindex));
      task28->add_dep(task0);
      queue_->add_task(task28);


      std::vector<IndexRange> I19_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I19(new Tensor<T>(I19_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor29 = {I18, this->f1_, I19};
      std::shared_ptr<Task29<T>> task29(new Task29<T>(tensor29, pindex));
      task28->add_dep(task29);
      task29->add_dep(task0);
      queue_->add_task(task29);


      std::vector<IndexRange> I20_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I20(new Tensor<T>(I20_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor30 = {I19, t2, I20};
      std::shared_ptr<Task30<T>> task30(new Task30<T>(tensor30, pindex));
      task29->add_dep(task30);
      task30->add_dep(task0);
      queue_->add_task(task30);


      std::vector<std::shared_ptr<Tensor<T>>> tensor31 = {I20, Gamma2};
      std::shared_ptr<Task31<T>> task31(new Task31<T>(tensor31, pindex));
      task30->add_dep(task31);
      task31->add_dep(task0);
      queue_->add_task(task31);

      task31->add_dep(task2);

      std::vector<IndexRange> I23_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I23(new Tensor<T>(I23_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor32 = {I19, t2, I23};
      std::shared_ptr<Task32<T>> task32(new Task32<T>(tensor32, pindex));
      task29->add_dep(task32);
      task32->add_dep(task0);
      queue_->add_task(task32);


      std::vector<std::shared_ptr<Tensor<T>>> tensor33 = {I23, Gamma2};
      std::shared_ptr<Task33<T>> task33(new Task33<T>(tensor33, pindex));
      task32->add_dep(task33);
      task33->add_dep(task0);
      queue_->add_task(task33);

      task33->add_dep(task2);

      std::vector<IndexRange> I25_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I25(new Tensor<T>(I25_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor34 = {I18, t2, I25};
      std::shared_ptr<Task34<T>> task34(new Task34<T>(tensor34, pindex));
      task28->add_dep(task34);
      task34->add_dep(task0);
      queue_->add_task(task34);


      std::vector<std::shared_ptr<Tensor<T>>> tensor35 = {I25, Gamma6};
      std::shared_ptr<Task35<T>> task35(new Task35<T>(tensor35, pindex));
      task34->add_dep(task35);
      task35->add_dep(task0);
      queue_->add_task(task35);

      task35->add_dep(task3);

      std::vector<IndexRange> I27_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I27(new Tensor<T>(I27_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor36 = {I18, t2, I27};
      std::shared_ptr<Task36<T>> task36(new Task36<T>(tensor36, pindex));
      task28->add_dep(task36);
      task36->add_dep(task0);
      queue_->add_task(task36);


      std::vector<std::shared_ptr<Tensor<T>>> tensor37 = {I27, Gamma6};
      std::shared_ptr<Task37<T>> task37(new Task37<T>(tensor37, pindex));
      task36->add_dep(task37);
      task37->add_dep(task0);
      queue_->add_task(task37);

      task37->add_dep(task3);

      std::vector<IndexRange> I29_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I29(new Tensor<T>(I29_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor38 = {I18, this->f1_, I29};
      std::shared_ptr<Task38<T>> task38(new Task38<T>(tensor38, pindex));
      task28->add_dep(task38);
      task38->add_dep(task0);
      queue_->add_task(task38);


      std::vector<IndexRange> I30_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I30(new Tensor<T>(I30_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor39 = {I29, t2, I30};
      std::shared_ptr<Task39<T>> task39(new Task39<T>(tensor39, pindex));
      task38->add_dep(task39);
      task39->add_dep(task0);
      queue_->add_task(task39);


      std::vector<std::shared_ptr<Tensor<T>>> tensor40 = {I30, Gamma2};
      std::shared_ptr<Task40<T>> task40(new Task40<T>(tensor40, pindex));
      task39->add_dep(task40);
      task40->add_dep(task0);
      queue_->add_task(task40);

      task40->add_dep(task2);

      std::vector<IndexRange> I33_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I33(new Tensor<T>(I33_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor41 = {I29, t2, I33};
      std::shared_ptr<Task41<T>> task41(new Task41<T>(tensor41, pindex));
      task38->add_dep(task41);
      task41->add_dep(task0);
      queue_->add_task(task41);


      std::vector<std::shared_ptr<Tensor<T>>> tensor42 = {I33, Gamma2};
      std::shared_ptr<Task42<T>> task42(new Task42<T>(tensor42, pindex));
      task41->add_dep(task42);
      task42->add_dep(task0);
      queue_->add_task(task42);

      task42->add_dep(task2);

      std::vector<IndexRange> I35_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I35(new Tensor<T>(I35_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor43 = {I18, this->f1_, I35};
      std::shared_ptr<Task43<T>> task43(new Task43<T>(tensor43, pindex));
      task28->add_dep(task43);
      task43->add_dep(task0);
      queue_->add_task(task43);


      std::vector<IndexRange> I36_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I36(new Tensor<T>(I36_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor44 = {I35, t2, I36};
      std::shared_ptr<Task44<T>> task44(new Task44<T>(tensor44, pindex));
      task43->add_dep(task44);
      task44->add_dep(task0);
      queue_->add_task(task44);


      std::vector<std::shared_ptr<Tensor<T>>> tensor45 = {I36, Gamma2};
      std::shared_ptr<Task45<T>> task45(new Task45<T>(tensor45, pindex));
      task44->add_dep(task45);
      task45->add_dep(task0);
      queue_->add_task(task45);

      task45->add_dep(task2);

      std::vector<IndexRange> I39_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I39(new Tensor<T>(I39_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor46 = {I35, t2, I39};
      std::shared_ptr<Task46<T>> task46(new Task46<T>(tensor46, pindex));
      task43->add_dep(task46);
      task46->add_dep(task0);
      queue_->add_task(task46);


      std::vector<std::shared_ptr<Tensor<T>>> tensor47 = {I39, Gamma2};
      std::shared_ptr<Task47<T>> task47(new Task47<T>(tensor47, pindex));
      task46->add_dep(task47);
      task47->add_dep(task0);
      queue_->add_task(task47);

      task47->add_dep(task2);

      std::vector<IndexRange> I41_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I41(new Tensor<T>(I41_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor48 = {I18, this->f1_, I41};
      std::shared_ptr<Task48<T>> task48(new Task48<T>(tensor48, pindex));
      task28->add_dep(task48);
      task48->add_dep(task0);
      queue_->add_task(task48);


      std::vector<IndexRange> I42_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I42(new Tensor<T>(I42_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor49 = {I41, t2, I42};
      std::shared_ptr<Task49<T>> task49(new Task49<T>(tensor49, pindex));
      task48->add_dep(task49);
      task49->add_dep(task0);
      queue_->add_task(task49);


      std::vector<std::shared_ptr<Tensor<T>>> tensor50 = {I42, Gamma2};
      std::shared_ptr<Task50<T>> task50(new Task50<T>(tensor50, pindex));
      task49->add_dep(task50);
      task50->add_dep(task0);
      queue_->add_task(task50);

      task50->add_dep(task2);

      std::vector<IndexRange> I45_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I45(new Tensor<T>(I45_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor51 = {I41, t2, I45};
      std::shared_ptr<Task51<T>> task51(new Task51<T>(tensor51, pindex));
      task48->add_dep(task51);
      task51->add_dep(task0);
      queue_->add_task(task51);


      std::vector<std::shared_ptr<Tensor<T>>> tensor52 = {I45, Gamma2};
      std::shared_ptr<Task52<T>> task52(new Task52<T>(tensor52, pindex));
      task51->add_dep(task52);
      task52->add_dep(task0);
      queue_->add_task(task52);

      task52->add_dep(task2);

      std::vector<IndexRange> I47_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I47(new Tensor<T>(I47_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor53 = {I18, this->f1_, I47};
      std::shared_ptr<Task53<T>> task53(new Task53<T>(tensor53, pindex));
      task28->add_dep(task53);
      task53->add_dep(task0);
      queue_->add_task(task53);


      std::vector<IndexRange> I48_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I48(new Tensor<T>(I48_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor54 = {I47, t2, I48};
      std::shared_ptr<Task54<T>> task54(new Task54<T>(tensor54, pindex));
      task53->add_dep(task54);
      task54->add_dep(task0);
      queue_->add_task(task54);


      std::vector<std::shared_ptr<Tensor<T>>> tensor55 = {I48, Gamma14};
      std::shared_ptr<Task55<T>> task55(new Task55<T>(tensor55, pindex));
      task54->add_dep(task55);
      task55->add_dep(task0);
      queue_->add_task(task55);

      task55->add_dep(task4);

      std::vector<IndexRange> I60_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I60(new Tensor<T>(I60_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor56 = {I18, t2, I60};
      std::shared_ptr<Task56<T>> task56(new Task56<T>(tensor56, pindex));
      task28->add_dep(task56);
      task56->add_dep(task0);
      queue_->add_task(task56);


      std::vector<std::shared_ptr<Tensor<T>>> tensor57 = {I60, Gamma2};
      std::shared_ptr<Task57<T>> task57(new Task57<T>(tensor57, pindex, this->e0_));
      task56->add_dep(task57);
      task57->add_dep(task0);
      queue_->add_task(task57);

      task57->add_dep(task2);

      std::vector<IndexRange> I62_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I62(new Tensor<T>(I62_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor58 = {I18, t2, I62};
      std::shared_ptr<Task58<T>> task58(new Task58<T>(tensor58, pindex));
      task28->add_dep(task58);
      task58->add_dep(task0);
      queue_->add_task(task58);


      std::vector<std::shared_ptr<Tensor<T>>> tensor59 = {I62, Gamma2};
      std::shared_ptr<Task59<T>> task59(new Task59<T>(tensor59, pindex, this->e0_));
      task58->add_dep(task59);
      task59->add_dep(task0);
      queue_->add_task(task59);

      task59->add_dep(task2);

      std::vector<IndexRange> I68_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I68(new Tensor<T>(I68_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor60 = {I18, this->v2_, I68};
      std::shared_ptr<Task60<T>> task60(new Task60<T>(tensor60, pindex));
      task28->add_dep(task60);
      task60->add_dep(task0);
      queue_->add_task(task60);


      std::vector<std::shared_ptr<Tensor<T>>> tensor61 = {I68, Gamma2};
      std::shared_ptr<Task61<T>> task61(new Task61<T>(tensor61, pindex));
      task60->add_dep(task61);
      task61->add_dep(task0);
      queue_->add_task(task61);

      task61->add_dep(task2);

      std::vector<IndexRange> I70_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I70(new Tensor<T>(I70_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor62 = {I18, this->v2_, I70};
      std::shared_ptr<Task62<T>> task62(new Task62<T>(tensor62, pindex));
      task28->add_dep(task62);
      task62->add_dep(task0);
      queue_->add_task(task62);


      std::vector<std::shared_ptr<Tensor<T>>> tensor63 = {I70, Gamma2};
      std::shared_ptr<Task63<T>> task63(new Task63<T>(tensor63, pindex));
      task62->add_dep(task63);
      task63->add_dep(task0);
      queue_->add_task(task63);

      task63->add_dep(task2);

      std::vector<IndexRange> I49_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I49(new Tensor<T>(I49_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor64 = {r, I49};
      std::shared_ptr<Task64<T>> task64(new Task64<T>(tensor64, pindex));
      task64->add_dep(task0);
      queue_->add_task(task64);


      std::vector<IndexRange> I50_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I50(new Tensor<T>(I50_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor65 = {I49, this->f1_, I50};
      std::shared_ptr<Task65<T>> task65(new Task65<T>(tensor65, pindex));
      task64->add_dep(task65);
      task65->add_dep(task0);
      queue_->add_task(task65);


      std::vector<IndexRange> I51_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I51(new Tensor<T>(I51_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor66 = {I50, t2, I51};
      std::shared_ptr<Task66<T>> task66(new Task66<T>(tensor66, pindex));
      task65->add_dep(task66);
      task66->add_dep(task0);
      queue_->add_task(task66);


      std::vector<std::shared_ptr<Tensor<T>>> tensor67 = {I51, Gamma14};
      std::shared_ptr<Task67<T>> task67(new Task67<T>(tensor67, pindex));
      task66->add_dep(task67);
      task67->add_dep(task0);
      queue_->add_task(task67);

      task67->add_dep(task4);

      std::vector<IndexRange> I55_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I55(new Tensor<T>(I55_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor68 = {I49, this->f1_, I55};
      std::shared_ptr<Task68<T>> task68(new Task68<T>(tensor68, pindex));
      task64->add_dep(task68);
      task68->add_dep(task0);
      queue_->add_task(task68);


      std::vector<IndexRange> I56_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I56(new Tensor<T>(I56_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor69 = {I55, t2, I56};
      std::shared_ptr<Task69<T>> task69(new Task69<T>(tensor69, pindex));
      task68->add_dep(task69);
      task69->add_dep(task0);
      queue_->add_task(task69);


      std::vector<std::shared_ptr<Tensor<T>>> tensor70 = {I56, Gamma14};
      std::shared_ptr<Task70<T>> task70(new Task70<T>(tensor70, pindex));
      task69->add_dep(task70);
      task70->add_dep(task0);
      queue_->add_task(task70);

      task70->add_dep(task4);

      std::vector<IndexRange> I52_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I52(new Tensor<T>(I52_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor71 = {r, I52};
      std::shared_ptr<Task71<T>> task71(new Task71<T>(tensor71, pindex));
      task71->add_dep(task0);
      queue_->add_task(task71);


      std::vector<IndexRange> I53_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I53(new Tensor<T>(I53_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor72 = {I52, t2, I53};
      std::shared_ptr<Task72<T>> task72(new Task72<T>(tensor72, pindex));
      task71->add_dep(task72);
      task72->add_dep(task0);
      queue_->add_task(task72);


      std::vector<std::shared_ptr<Tensor<T>>> tensor73 = {I53, Gamma16};
      std::shared_ptr<Task73<T>> task73(new Task73<T>(tensor73, pindex));
      task72->add_dep(task73);
      task73->add_dep(task0);
      queue_->add_task(task73);

      task73->add_dep(task5);

      std::vector<IndexRange> I64_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I64(new Tensor<T>(I64_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor74 = {I52, t2, I64};
      std::shared_ptr<Task74<T>> task74(new Task74<T>(tensor74, pindex));
      task71->add_dep(task74);
      task74->add_dep(task0);
      queue_->add_task(task74);


      std::vector<std::shared_ptr<Tensor<T>>> tensor75 = {I64, Gamma14};
      std::shared_ptr<Task75<T>> task75(new Task75<T>(tensor75, pindex, this->e0_));
      task74->add_dep(task75);
      task75->add_dep(task0);
      queue_->add_task(task75);

      task75->add_dep(task4);

      std::vector<IndexRange> I72_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I72(new Tensor<T>(I72_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor76 = {I52, this->v2_, I72};
      std::shared_ptr<Task76<T>> task76(new Task76<T>(tensor76, pindex));
      task71->add_dep(task76);
      task76->add_dep(task0);
      queue_->add_task(task76);


      std::vector<std::shared_ptr<Tensor<T>>> tensor77 = {I72, Gamma14};
      std::shared_ptr<Task77<T>> task77(new Task77<T>(tensor77, pindex));
      task76->add_dep(task77);
      task77->add_dep(task0);
      queue_->add_task(task77);

      task77->add_dep(task4);

      std::shared_ptr<Queue<T>> energy_(new Queue<T>());
      std::vector<IndexRange> I73_index;
      std::shared_ptr<Tensor<T>> I73(new Tensor<T>(I73_index, false));
      std::vector<IndexRange> I74_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I74(new Tensor<T>(I74_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor78 = {I73, t2, I74};
      std::shared_ptr<Task78<T>> task78(new Task78<T>(tensor78, pindex));
      energy_->add_task(task78);


      std::vector<std::shared_ptr<Tensor<T>>> tensor79 = {I74, t2, this->v2_};
      std::shared_ptr<Task79<T>> task79(new Task79<T>(tensor79, pindex, this->e0_));
      task78->add_dep(task79);
      energy_->add_task(task79);


      std::vector<IndexRange> I75_index;
      std::shared_ptr<Tensor<T>> I75(new Tensor<T>(I75_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor80 = {I74, t2, I75};
      std::shared_ptr<Task80<T>> task80(new Task80<T>(tensor80, pindex));
      task78->add_dep(task80);
      energy_->add_task(task80);


      std::vector<std::shared_ptr<Tensor<T>>> tensor81 = {I75, Gamma0};
      std::shared_ptr<Task81<T>> task81(new Task81<T>(tensor81, pindex));
      task80->add_dep(task81);
      energy_->add_task(task81);

      task81->add_dep(task1);

      std::vector<IndexRange> I78_index;
      std::shared_ptr<Tensor<T>> I78(new Tensor<T>(I78_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor82 = {I74, t2, I78};
      std::shared_ptr<Task82<T>> task82(new Task82<T>(tensor82, pindex));
      task78->add_dep(task82);
      energy_->add_task(task82);


      std::vector<std::shared_ptr<Tensor<T>>> tensor83 = {I78, Gamma0};
      std::shared_ptr<Task83<T>> task83(new Task83<T>(tensor83, pindex));
      task82->add_dep(task83);
      energy_->add_task(task83);

      task83->add_dep(task1);

      std::vector<IndexRange> I81_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I81(new Tensor<T>(I81_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor84 = {I74, this->f1_, I81};
      std::shared_ptr<Task84<T>> task84(new Task84<T>(tensor84, pindex));
      task78->add_dep(task84);
      energy_->add_task(task84);


      std::vector<std::shared_ptr<Tensor<T>>> tensor85 = {I81, t2};
      std::shared_ptr<Task85<T>> task85(new Task85<T>(tensor85, pindex));
      task84->add_dep(task85);
      energy_->add_task(task85);


      std::vector<IndexRange> I87_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I87(new Tensor<T>(I87_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor86 = {I74, this->f1_, I87};
      std::shared_ptr<Task86<T>> task86(new Task86<T>(tensor86, pindex));
      task78->add_dep(task86);
      energy_->add_task(task86);


      std::vector<std::shared_ptr<Tensor<T>>> tensor87 = {I87, t2};
      std::shared_ptr<Task87<T>> task87(new Task87<T>(tensor87, pindex));
      task86->add_dep(task87);
      energy_->add_task(task87);


      std::vector<IndexRange> I93_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I93(new Tensor<T>(I93_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor88 = {I74, this->f1_, I93};
      std::shared_ptr<Task88<T>> task88(new Task88<T>(tensor88, pindex));
      task78->add_dep(task88);
      energy_->add_task(task88);


      std::vector<IndexRange> I94_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I94(new Tensor<T>(I94_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor89 = {I93, t2, I94};
      std::shared_ptr<Task89<T>> task89(new Task89<T>(tensor89, pindex));
      task88->add_dep(task89);
      energy_->add_task(task89);


      std::vector<std::shared_ptr<Tensor<T>>> tensor90 = {I94, Gamma2};
      std::shared_ptr<Task90<T>> task90(new Task90<T>(tensor90, pindex));
      task89->add_dep(task90);
      energy_->add_task(task90);

      task90->add_dep(task2);

      std::vector<IndexRange> I98_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I98(new Tensor<T>(I98_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor91 = {I93, t2, I98};
      std::shared_ptr<Task91<T>> task91(new Task91<T>(tensor91, pindex));
      task88->add_dep(task91);
      energy_->add_task(task91);


      std::vector<std::shared_ptr<Tensor<T>>> tensor92 = {I98, Gamma2};
      std::shared_ptr<Task92<T>> task92(new Task92<T>(tensor92, pindex));
      task91->add_dep(task92);
      energy_->add_task(task92);

      task92->add_dep(task2);

      std::vector<IndexRange> I100_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I100(new Tensor<T>(I100_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor93 = {I73, t2, I100};
      std::shared_ptr<Task93<T>> task93(new Task93<T>(tensor93, pindex));
      task78->add_dep(task93);
      energy_->add_task(task93);


      std::vector<IndexRange> I101_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I101(new Tensor<T>(I101_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor94 = {I100, this->f1_, I101};
      std::shared_ptr<Task94<T>> task94(new Task94<T>(tensor94, pindex));
      task93->add_dep(task94);
      energy_->add_task(task94);


      std::vector<IndexRange> I102_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I102(new Tensor<T>(I102_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor95 = {I101, t2, I102};
      std::shared_ptr<Task95<T>> task95(new Task95<T>(tensor95, pindex));
      task94->add_dep(task95);
      energy_->add_task(task95);


      std::vector<std::shared_ptr<Tensor<T>>> tensor96 = {I102, Gamma2};
      std::shared_ptr<Task96<T>> task96(new Task96<T>(tensor96, pindex));
      task95->add_dep(task96);
      energy_->add_task(task96);

      task96->add_dep(task2);

      std::vector<IndexRange> I106_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I106(new Tensor<T>(I106_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor97 = {I101, t2, I106};
      std::shared_ptr<Task97<T>> task97(new Task97<T>(tensor97, pindex));
      task94->add_dep(task97);
      energy_->add_task(task97);


      std::vector<std::shared_ptr<Tensor<T>>> tensor98 = {I106, Gamma2};
      std::shared_ptr<Task98<T>> task98(new Task98<T>(tensor98, pindex));
      task97->add_dep(task98);
      energy_->add_task(task98);

      task98->add_dep(task2);

      std::vector<IndexRange> I109_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I109(new Tensor<T>(I109_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor99 = {I100, t2, I109};
      std::shared_ptr<Task99<T>> task99(new Task99<T>(tensor99, pindex));
      task93->add_dep(task99);
      energy_->add_task(task99);


      std::vector<std::shared_ptr<Tensor<T>>> tensor100 = {I109, Gamma6};
      std::shared_ptr<Task100<T>> task100(new Task100<T>(tensor100, pindex));
      task99->add_dep(task100);
      energy_->add_task(task100);

      task100->add_dep(task3);

      std::vector<IndexRange> I112_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I112(new Tensor<T>(I112_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor101 = {I100, t2, I112};
      std::shared_ptr<Task101<T>> task101(new Task101<T>(tensor101, pindex));
      task93->add_dep(task101);
      energy_->add_task(task101);


      std::vector<std::shared_ptr<Tensor<T>>> tensor102 = {I112, Gamma6};
      std::shared_ptr<Task102<T>> task102(new Task102<T>(tensor102, pindex));
      task101->add_dep(task102);
      energy_->add_task(task102);

      task102->add_dep(task3);

      std::vector<IndexRange> I115_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I115(new Tensor<T>(I115_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor103 = {I100, this->f1_, I115};
      std::shared_ptr<Task103<T>> task103(new Task103<T>(tensor103, pindex));
      task93->add_dep(task103);
      energy_->add_task(task103);


      std::vector<IndexRange> I116_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I116(new Tensor<T>(I116_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor104 = {I115, t2, I116};
      std::shared_ptr<Task104<T>> task104(new Task104<T>(tensor104, pindex));
      task103->add_dep(task104);
      energy_->add_task(task104);


      std::vector<std::shared_ptr<Tensor<T>>> tensor105 = {I116, Gamma2};
      std::shared_ptr<Task105<T>> task105(new Task105<T>(tensor105, pindex));
      task104->add_dep(task105);
      energy_->add_task(task105);

      task105->add_dep(task2);

      std::vector<IndexRange> I120_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I120(new Tensor<T>(I120_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor106 = {I115, t2, I120};
      std::shared_ptr<Task106<T>> task106(new Task106<T>(tensor106, pindex));
      task103->add_dep(task106);
      energy_->add_task(task106);


      std::vector<std::shared_ptr<Tensor<T>>> tensor107 = {I120, Gamma2};
      std::shared_ptr<Task107<T>> task107(new Task107<T>(tensor107, pindex));
      task106->add_dep(task107);
      energy_->add_task(task107);

      task107->add_dep(task2);

      std::vector<IndexRange> I123_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I123(new Tensor<T>(I123_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor108 = {I100, this->f1_, I123};
      std::shared_ptr<Task108<T>> task108(new Task108<T>(tensor108, pindex));
      task93->add_dep(task108);
      energy_->add_task(task108);


      std::vector<IndexRange> I124_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I124(new Tensor<T>(I124_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor109 = {I123, t2, I124};
      std::shared_ptr<Task109<T>> task109(new Task109<T>(tensor109, pindex));
      task108->add_dep(task109);
      energy_->add_task(task109);


      std::vector<std::shared_ptr<Tensor<T>>> tensor110 = {I124, Gamma2};
      std::shared_ptr<Task110<T>> task110(new Task110<T>(tensor110, pindex));
      task109->add_dep(task110);
      energy_->add_task(task110);

      task110->add_dep(task2);

      std::vector<IndexRange> I128_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I128(new Tensor<T>(I128_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor111 = {I123, t2, I128};
      std::shared_ptr<Task111<T>> task111(new Task111<T>(tensor111, pindex));
      task108->add_dep(task111);
      energy_->add_task(task111);


      std::vector<std::shared_ptr<Tensor<T>>> tensor112 = {I128, Gamma2};
      std::shared_ptr<Task112<T>> task112(new Task112<T>(tensor112, pindex));
      task111->add_dep(task112);
      energy_->add_task(task112);

      task112->add_dep(task2);

      std::vector<IndexRange> I131_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I131(new Tensor<T>(I131_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor113 = {I100, this->f1_, I131};
      std::shared_ptr<Task113<T>> task113(new Task113<T>(tensor113, pindex));
      task93->add_dep(task113);
      energy_->add_task(task113);


      std::vector<IndexRange> I132_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I132(new Tensor<T>(I132_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor114 = {I131, t2, I132};
      std::shared_ptr<Task114<T>> task114(new Task114<T>(tensor114, pindex));
      task113->add_dep(task114);
      energy_->add_task(task114);


      std::vector<std::shared_ptr<Tensor<T>>> tensor115 = {I132, Gamma2};
      std::shared_ptr<Task115<T>> task115(new Task115<T>(tensor115, pindex));
      task114->add_dep(task115);
      energy_->add_task(task115);

      task115->add_dep(task2);

      std::vector<IndexRange> I136_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I136(new Tensor<T>(I136_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor116 = {I131, t2, I136};
      std::shared_ptr<Task116<T>> task116(new Task116<T>(tensor116, pindex));
      task113->add_dep(task116);
      energy_->add_task(task116);


      std::vector<std::shared_ptr<Tensor<T>>> tensor117 = {I136, Gamma2};
      std::shared_ptr<Task117<T>> task117(new Task117<T>(tensor117, pindex));
      task116->add_dep(task117);
      energy_->add_task(task117);

      task117->add_dep(task2);

      std::vector<IndexRange> I139_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I139(new Tensor<T>(I139_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor118 = {I100, this->f1_, I139};
      std::shared_ptr<Task118<T>> task118(new Task118<T>(tensor118, pindex));
      task93->add_dep(task118);
      energy_->add_task(task118);


      std::vector<IndexRange> I140_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I140(new Tensor<T>(I140_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor119 = {I139, t2, I140};
      std::shared_ptr<Task119<T>> task119(new Task119<T>(tensor119, pindex));
      task118->add_dep(task119);
      energy_->add_task(task119);


      std::vector<std::shared_ptr<Tensor<T>>> tensor120 = {I140, Gamma14};
      std::shared_ptr<Task120<T>> task120(new Task120<T>(tensor120, pindex));
      task119->add_dep(task120);
      energy_->add_task(task120);

      task120->add_dep(task4);

      std::vector<IndexRange> I158_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I158(new Tensor<T>(I158_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor121 = {I100, t2, I158};
      std::shared_ptr<Task121<T>> task121(new Task121<T>(tensor121, pindex));
      task93->add_dep(task121);
      energy_->add_task(task121);


      std::vector<std::shared_ptr<Tensor<T>>> tensor122 = {I158, Gamma2};
      std::shared_ptr<Task122<T>> task122(new Task122<T>(tensor122, pindex, this->e0_));
      task121->add_dep(task122);
      energy_->add_task(task122);

      task122->add_dep(task2);

      std::vector<IndexRange> I161_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I161(new Tensor<T>(I161_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor123 = {I100, t2, I161};
      std::shared_ptr<Task123<T>> task123(new Task123<T>(tensor123, pindex));
      task93->add_dep(task123);
      energy_->add_task(task123);


      std::vector<std::shared_ptr<Tensor<T>>> tensor124 = {I161, Gamma2};
      std::shared_ptr<Task124<T>> task124(new Task124<T>(tensor124, pindex, this->e0_));
      task123->add_dep(task124);
      energy_->add_task(task124);

      task124->add_dep(task2);

      std::vector<IndexRange> I171_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I171(new Tensor<T>(I171_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor125 = {I100, this->v2_, I171};
      std::shared_ptr<Task125<T>> task125(new Task125<T>(tensor125, pindex));
      task93->add_dep(task125);
      energy_->add_task(task125);


      std::vector<std::shared_ptr<Tensor<T>>> tensor126 = {I171, Gamma2};
      std::shared_ptr<Task126<T>> task126(new Task126<T>(tensor126, pindex));
      task125->add_dep(task126);
      energy_->add_task(task126);

      task126->add_dep(task2);

      std::vector<IndexRange> I174_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I174(new Tensor<T>(I174_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor127 = {I100, this->v2_, I174};
      std::shared_ptr<Task127<T>> task127(new Task127<T>(tensor127, pindex));
      task93->add_dep(task127);
      energy_->add_task(task127);


      std::vector<std::shared_ptr<Tensor<T>>> tensor128 = {I174, Gamma2};
      std::shared_ptr<Task128<T>> task128(new Task128<T>(tensor128, pindex));
      task127->add_dep(task128);
      energy_->add_task(task128);

      task128->add_dep(task2);

      std::vector<IndexRange> I142_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I142(new Tensor<T>(I142_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor129 = {I73, t2, I142};
      std::shared_ptr<Task129<T>> task129(new Task129<T>(tensor129, pindex));
      task78->add_dep(task129);
      energy_->add_task(task129);


      std::vector<IndexRange> I143_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I143(new Tensor<T>(I143_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor130 = {I142, this->f1_, I143};
      std::shared_ptr<Task130<T>> task130(new Task130<T>(tensor130, pindex));
      task129->add_dep(task130);
      energy_->add_task(task130);


      std::vector<IndexRange> I144_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I144(new Tensor<T>(I144_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor131 = {I143, t2, I144};
      std::shared_ptr<Task131<T>> task131(new Task131<T>(tensor131, pindex));
      task130->add_dep(task131);
      energy_->add_task(task131);


      std::vector<std::shared_ptr<Tensor<T>>> tensor132 = {I144, Gamma14};
      std::shared_ptr<Task132<T>> task132(new Task132<T>(tensor132, pindex));
      task131->add_dep(task132);
      energy_->add_task(task132);

      task132->add_dep(task4);

      std::vector<IndexRange> I147_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I147(new Tensor<T>(I147_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor133 = {I142, t2, I147};
      std::shared_ptr<Task133<T>> task133(new Task133<T>(tensor133, pindex));
      task129->add_dep(task133);
      energy_->add_task(task133);


      std::vector<std::shared_ptr<Tensor<T>>> tensor134 = {I147, Gamma16};
      std::shared_ptr<Task134<T>> task134(new Task134<T>(tensor134, pindex));
      task133->add_dep(task134);
      energy_->add_task(task134);

      task134->add_dep(task5);

      std::vector<IndexRange> I150_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I150(new Tensor<T>(I150_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor135 = {I142, this->f1_, I150};
      std::shared_ptr<Task135<T>> task135(new Task135<T>(tensor135, pindex));
      task129->add_dep(task135);
      energy_->add_task(task135);


      std::vector<IndexRange> I151_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I151(new Tensor<T>(I151_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor136 = {I150, t2, I151};
      std::shared_ptr<Task136<T>> task136(new Task136<T>(tensor136, pindex));
      task135->add_dep(task136);
      energy_->add_task(task136);


      std::vector<std::shared_ptr<Tensor<T>>> tensor137 = {I151, Gamma14};
      std::shared_ptr<Task137<T>> task137(new Task137<T>(tensor137, pindex));
      task136->add_dep(task137);
      energy_->add_task(task137);

      task137->add_dep(task4);

      std::vector<IndexRange> I164_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I164(new Tensor<T>(I164_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor138 = {I142, t2, I164};
      std::shared_ptr<Task138<T>> task138(new Task138<T>(tensor138, pindex));
      task129->add_dep(task138);
      energy_->add_task(task138);


      std::vector<std::shared_ptr<Tensor<T>>> tensor139 = {I164, Gamma14};
      std::shared_ptr<Task139<T>> task139(new Task139<T>(tensor139, pindex, this->e0_));
      task138->add_dep(task139);
      energy_->add_task(task139);

      task139->add_dep(task4);

      std::vector<IndexRange> I177_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I177(new Tensor<T>(I177_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor140 = {I142, this->v2_, I177};
      std::shared_ptr<Task140<T>> task140(new Task140<T>(tensor140, pindex));
      task129->add_dep(task140);
      energy_->add_task(task140);


      std::vector<std::shared_ptr<Tensor<T>>> tensor141 = {I177, Gamma14};
      std::shared_ptr<Task141<T>> task141(new Task141<T>(tensor141, pindex));
      task140->add_dep(task141);
      energy_->add_task(task141);

      task141->add_dep(task4);

      std::shared_ptr<Queue<T>> correction_(new Queue<T>());
      std::vector<IndexRange> I178_index;
      std::shared_ptr<Tensor<T>> I178(new Tensor<T>(I178_index, false));
      std::vector<IndexRange> I179_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I179(new Tensor<T>(I179_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor142 = {I178, t2, I179};
      std::shared_ptr<Task142<T>> task142(new Task142<T>(tensor142, pindex));
      correction_->add_task(task142);


      std::vector<std::shared_ptr<Tensor<T>>> tensor143 = {I179, t2};
      std::shared_ptr<Task143<T>> task143(new Task143<T>(tensor143, pindex));
      task142->add_dep(task143);
      correction_->add_task(task143);


      std::vector<IndexRange> I183_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I183(new Tensor<T>(I183_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor144 = {I178, t2, I183};
      std::shared_ptr<Task144<T>> task144(new Task144<T>(tensor144, pindex));
      task142->add_dep(task144);
      correction_->add_task(task144);


      std::vector<IndexRange> I184_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I184(new Tensor<T>(I184_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor145 = {I183, t2, I184};
      std::shared_ptr<Task145<T>> task145(new Task145<T>(tensor145, pindex));
      task144->add_dep(task145);
      correction_->add_task(task145);


      std::vector<std::shared_ptr<Tensor<T>>> tensor146 = {I184, Gamma2};
      std::shared_ptr<Task146<T>> task146(new Task146<T>(tensor146, pindex));
      task145->add_dep(task146);
      correction_->add_task(task146);

      task146->add_dep(task2);

      std::vector<IndexRange> I187_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I187(new Tensor<T>(I187_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor147 = {I183, t2, I187};
      std::shared_ptr<Task147<T>> task147(new Task147<T>(tensor147, pindex));
      task144->add_dep(task147);
      correction_->add_task(task147);


      std::vector<std::shared_ptr<Tensor<T>>> tensor148 = {I187, Gamma2};
      std::shared_ptr<Task148<T>> task148(new Task148<T>(tensor148, pindex));
      task147->add_dep(task148);
      correction_->add_task(task148);

      task148->add_dep(task2);

      std::vector<IndexRange> I189_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I189(new Tensor<T>(I189_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor149 = {I178, t2, I189};
      std::shared_ptr<Task149<T>> task149(new Task149<T>(tensor149, pindex));
      task142->add_dep(task149);
      correction_->add_task(task149);


      std::vector<IndexRange> I190_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I190(new Tensor<T>(I190_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor150 = {I189, t2, I190};
      std::shared_ptr<Task150<T>> task150(new Task150<T>(tensor150, pindex));
      task149->add_dep(task150);
      correction_->add_task(task150);


      std::vector<std::shared_ptr<Tensor<T>>> tensor151 = {I190, Gamma14};
      std::shared_ptr<Task151<T>> task151(new Task151<T>(tensor151, pindex));
      task150->add_dep(task151);
      correction_->add_task(task151);

      task151->add_dep(task4);

      std::shared_ptr<Queue<T>> density_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor152 = {den2};
      std::shared_ptr<Task152<T>> task152(new Task152<T>(tensor152));
      density_->add_task(task152);

      std::vector<IndexRange> I191_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I191(new Tensor<T>(I191_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor153 = {den2, I191};
      std::shared_ptr<Task153<T>> task153(new Task153<T>(tensor153, pindex));
      task153->add_dep(task152);
      density_->add_task(task153);


      std::vector<IndexRange> I192_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I192(new Tensor<T>(I192_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor154 = {I191, t2, I192};
      std::shared_ptr<Task154<T>> task154(new Task154<T>(tensor154, pindex));
      task153->add_dep(task154);
      task154->add_dep(task152);
      density_->add_task(task154);


      std::vector<IndexRange> I193_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I193(new Tensor<T>(I193_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor155 = {I192, t2, I193};
      std::shared_ptr<Task155<T>> task155(new Task155<T>(tensor155, pindex));
      task154->add_dep(task155);
      task155->add_dep(task152);
      density_->add_task(task155);


      std::vector<std::shared_ptr<Tensor<T>>> tensor156 = {I193, Gamma2};
      std::shared_ptr<Task156<T>> task156(new Task156<T>(tensor156, pindex));
      task155->add_dep(task156);
      task156->add_dep(task152);
      density_->add_task(task156);

      task156->add_dep(task2);

      std::vector<IndexRange> I196_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I196(new Tensor<T>(I196_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor157 = {I192, t2, I196};
      std::shared_ptr<Task157<T>> task157(new Task157<T>(tensor157, pindex));
      task154->add_dep(task157);
      task157->add_dep(task152);
      density_->add_task(task157);


      std::vector<std::shared_ptr<Tensor<T>>> tensor158 = {I196, Gamma2};
      std::shared_ptr<Task158<T>> task158(new Task158<T>(tensor158, pindex));
      task157->add_dep(task158);
      task158->add_dep(task152);
      density_->add_task(task158);

      task158->add_dep(task2);

      std::vector<IndexRange> I197_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor<T>> I197(new Tensor<T>(I197_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor159 = {den2, I197};
      std::shared_ptr<Task159<T>> task159(new Task159<T>(tensor159, pindex));
      task159->add_dep(task152);
      density_->add_task(task159);


      std::vector<IndexRange> I198_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I198(new Tensor<T>(I198_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor160 = {I197, t2, I198};
      std::shared_ptr<Task160<T>> task160(new Task160<T>(tensor160, pindex));
      task159->add_dep(task160);
      task160->add_dep(task152);
      density_->add_task(task160);


      std::vector<std::shared_ptr<Tensor<T>>> tensor161 = {I198, t2};
      std::shared_ptr<Task161<T>> task161(new Task161<T>(tensor161, pindex));
      task160->add_dep(task161);
      task161->add_dep(task152);
      density_->add_task(task161);


      std::vector<IndexRange> I201_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I201(new Tensor<T>(I201_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor162 = {den2, I201};
      std::shared_ptr<Task162<T>> task162(new Task162<T>(tensor162, pindex));
      task162->add_dep(task152);
      density_->add_task(task162);


      std::vector<IndexRange> I202_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I202(new Tensor<T>(I202_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor163 = {I201, t2, I202};
      std::shared_ptr<Task163<T>> task163(new Task163<T>(tensor163, pindex));
      task162->add_dep(task163);
      task163->add_dep(task152);
      density_->add_task(task163);


      std::vector<std::shared_ptr<Tensor<T>>> tensor164 = {I202, t2};
      std::shared_ptr<Task164<T>> task164(new Task164<T>(tensor164, pindex));
      task163->add_dep(task164);
      task164->add_dep(task152);
      density_->add_task(task164);


      std::vector<IndexRange> I205_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor<T>> I205(new Tensor<T>(I205_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor165 = {den2, I205};
      std::shared_ptr<Task165<T>> task165(new Task165<T>(tensor165, pindex));
      task165->add_dep(task152);
      density_->add_task(task165);


      std::vector<IndexRange> I206_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I206(new Tensor<T>(I206_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor166 = {I205, t2, I206};
      std::shared_ptr<Task166<T>> task166(new Task166<T>(tensor166, pindex));
      task165->add_dep(task166);
      task166->add_dep(task152);
      density_->add_task(task166);


      std::vector<IndexRange> I207_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I207(new Tensor<T>(I207_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor167 = {I206, t2, I207};
      std::shared_ptr<Task167<T>> task167(new Task167<T>(tensor167, pindex));
      task166->add_dep(task167);
      task167->add_dep(task152);
      density_->add_task(task167);


      std::vector<std::shared_ptr<Tensor<T>>> tensor168 = {I207, Gamma2};
      std::shared_ptr<Task168<T>> task168(new Task168<T>(tensor168, pindex));
      task167->add_dep(task168);
      task168->add_dep(task152);
      density_->add_task(task168);

      task168->add_dep(task2);

      std::vector<IndexRange> I210_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I210(new Tensor<T>(I210_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor169 = {I206, t2, I210};
      std::shared_ptr<Task169<T>> task169(new Task169<T>(tensor169, pindex));
      task166->add_dep(task169);
      task169->add_dep(task152);
      density_->add_task(task169);


      std::vector<std::shared_ptr<Tensor<T>>> tensor170 = {I210, Gamma2};
      std::shared_ptr<Task170<T>> task170(new Task170<T>(tensor170, pindex));
      task169->add_dep(task170);
      task170->add_dep(task152);
      density_->add_task(task170);

      task170->add_dep(task2);

      std::vector<IndexRange> I211_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor<T>> I211(new Tensor<T>(I211_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor171 = {den2, I211};
      std::shared_ptr<Task171<T>> task171(new Task171<T>(tensor171, pindex));
      task171->add_dep(task152);
      density_->add_task(task171);


      std::vector<IndexRange> I212_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I212(new Tensor<T>(I212_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor172 = {I211, t2, I212};
      std::shared_ptr<Task172<T>> task172(new Task172<T>(tensor172, pindex));
      task171->add_dep(task172);
      task172->add_dep(task152);
      density_->add_task(task172);


      std::vector<IndexRange> I213_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I213(new Tensor<T>(I213_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor173 = {I212, t2, I213};
      std::shared_ptr<Task173<T>> task173(new Task173<T>(tensor173, pindex));
      task172->add_dep(task173);
      task173->add_dep(task152);
      density_->add_task(task173);


      std::vector<std::shared_ptr<Tensor<T>>> tensor174 = {I213, Gamma2};
      std::shared_ptr<Task174<T>> task174(new Task174<T>(tensor174, pindex));
      task173->add_dep(task174);
      task174->add_dep(task152);
      density_->add_task(task174);

      task174->add_dep(task2);

      std::vector<IndexRange> I216_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I216(new Tensor<T>(I216_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor175 = {I212, t2, I216};
      std::shared_ptr<Task175<T>> task175(new Task175<T>(tensor175, pindex));
      task172->add_dep(task175);
      task175->add_dep(task152);
      density_->add_task(task175);


      std::vector<std::shared_ptr<Tensor<T>>> tensor176 = {I216, Gamma2};
      std::shared_ptr<Task176<T>> task176(new Task176<T>(tensor176, pindex));
      task175->add_dep(task176);
      task176->add_dep(task152);
      density_->add_task(task176);

      task176->add_dep(task2);

      std::vector<IndexRange> I217_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I217(new Tensor<T>(I217_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor177 = {den2, I217};
      std::shared_ptr<Task177<T>> task177(new Task177<T>(tensor177, pindex));
      task177->add_dep(task152);
      density_->add_task(task177);


      std::vector<IndexRange> I218_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I218(new Tensor<T>(I218_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor178 = {I217, t2, I218};
      std::shared_ptr<Task178<T>> task178(new Task178<T>(tensor178, pindex));
      task177->add_dep(task178);
      task178->add_dep(task152);
      density_->add_task(task178);


      std::vector<IndexRange> I219_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I219(new Tensor<T>(I219_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor179 = {I218, t2, I219};
      std::shared_ptr<Task179<T>> task179(new Task179<T>(tensor179, pindex));
      task178->add_dep(task179);
      task179->add_dep(task152);
      density_->add_task(task179);


      std::vector<std::shared_ptr<Tensor<T>>> tensor180 = {I219, Gamma14};
      std::shared_ptr<Task180<T>> task180(new Task180<T>(tensor180, pindex));
      task179->add_dep(task180);
      task180->add_dep(task152);
      density_->add_task(task180);

      task180->add_dep(task4);

      std::vector<IndexRange> I222_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I222(new Tensor<T>(I222_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor181 = {I218, t2, I222};
      std::shared_ptr<Task181<T>> task181(new Task181<T>(tensor181, pindex));
      task178->add_dep(task181);
      task181->add_dep(task152);
      density_->add_task(task181);


      std::vector<std::shared_ptr<Tensor<T>>> tensor182 = {I222, Gamma14};
      std::shared_ptr<Task182<T>> task182(new Task182<T>(tensor182, pindex));
      task181->add_dep(task182);
      task182->add_dep(task152);
      density_->add_task(task182);

      task182->add_dep(task4);

      std::vector<IndexRange> I223_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor<T>> I223(new Tensor<T>(I223_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor183 = {den2, I223};
      std::shared_ptr<Task183<T>> task183(new Task183<T>(tensor183, pindex));
      task183->add_dep(task152);
      density_->add_task(task183);


      std::vector<IndexRange> I224_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I224(new Tensor<T>(I224_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor184 = {I223, t2, I224};
      std::shared_ptr<Task184<T>> task184(new Task184<T>(tensor184, pindex));
      task183->add_dep(task184);
      task184->add_dep(task152);
      density_->add_task(task184);


      std::vector<IndexRange> I225_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I225(new Tensor<T>(I225_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor185 = {I224, t2, I225};
      std::shared_ptr<Task185<T>> task185(new Task185<T>(tensor185, pindex));
      task184->add_dep(task185);
      task185->add_dep(task152);
      density_->add_task(task185);


      std::vector<std::shared_ptr<Tensor<T>>> tensor186 = {I225, Gamma2};
      std::shared_ptr<Task186<T>> task186(new Task186<T>(tensor186, pindex));
      task185->add_dep(task186);
      task186->add_dep(task152);
      density_->add_task(task186);

      task186->add_dep(task2);

      std::vector<IndexRange> I228_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I228(new Tensor<T>(I228_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor187 = {I224, t2, I228};
      std::shared_ptr<Task187<T>> task187(new Task187<T>(tensor187, pindex));
      task184->add_dep(task187);
      task187->add_dep(task152);
      density_->add_task(task187);


      std::vector<std::shared_ptr<Tensor<T>>> tensor188 = {I228, Gamma2};
      std::shared_ptr<Task188<T>> task188(new Task188<T>(tensor188, pindex));
      task187->add_dep(task188);
      task188->add_dep(task152);
      density_->add_task(task188);

      task188->add_dep(task2);

      std::vector<IndexRange> I229_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I229(new Tensor<T>(I229_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor189 = {den2, I229};
      std::shared_ptr<Task189<T>> task189(new Task189<T>(tensor189, pindex));
      task189->add_dep(task152);
      density_->add_task(task189);


      std::vector<IndexRange> I230_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I230(new Tensor<T>(I230_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor190 = {I229, t2, I230};
      std::shared_ptr<Task190<T>> task190(new Task190<T>(tensor190, pindex));
      task189->add_dep(task190);
      task190->add_dep(task152);
      density_->add_task(task190);


      std::vector<IndexRange> I231_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I231(new Tensor<T>(I231_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor191 = {I230, t2, I231};
      std::shared_ptr<Task191<T>> task191(new Task191<T>(tensor191, pindex));
      task190->add_dep(task191);
      task191->add_dep(task152);
      density_->add_task(task191);


      std::vector<std::shared_ptr<Tensor<T>>> tensor192 = {I231, Gamma2};
      std::shared_ptr<Task192<T>> task192(new Task192<T>(tensor192, pindex));
      task191->add_dep(task192);
      task192->add_dep(task152);
      density_->add_task(task192);

      task192->add_dep(task2);

      std::vector<IndexRange> I234_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I234(new Tensor<T>(I234_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor193 = {I230, t2, I234};
      std::shared_ptr<Task193<T>> task193(new Task193<T>(tensor193, pindex));
      task190->add_dep(task193);
      task193->add_dep(task152);
      density_->add_task(task193);


      std::vector<std::shared_ptr<Tensor<T>>> tensor194 = {I234, Gamma2};
      std::shared_ptr<Task194<T>> task194(new Task194<T>(tensor194, pindex));
      task193->add_dep(task194);
      task194->add_dep(task152);
      density_->add_task(task194);

      task194->add_dep(task2);

      std::vector<IndexRange> I235_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I235(new Tensor<T>(I235_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor195 = {den2, I235};
      std::shared_ptr<Task195<T>> task195(new Task195<T>(tensor195, pindex));
      task195->add_dep(task152);
      density_->add_task(task195);


      std::vector<IndexRange> I236_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I236(new Tensor<T>(I236_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor196 = {I235, t2, I236};
      std::shared_ptr<Task196<T>> task196(new Task196<T>(tensor196, pindex));
      task195->add_dep(task196);
      task196->add_dep(task152);
      density_->add_task(task196);


      std::vector<IndexRange> I237_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I237(new Tensor<T>(I237_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor197 = {I236, t2, I237};
      std::shared_ptr<Task197<T>> task197(new Task197<T>(tensor197, pindex));
      task196->add_dep(task197);
      task197->add_dep(task152);
      density_->add_task(task197);


      std::vector<std::shared_ptr<Tensor<T>>> tensor198 = {I237, Gamma2};
      std::shared_ptr<Task198<T>> task198(new Task198<T>(tensor198, pindex));
      task197->add_dep(task198);
      task198->add_dep(task152);
      density_->add_task(task198);

      task198->add_dep(task2);

      std::vector<IndexRange> I240_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I240(new Tensor<T>(I240_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor199 = {I236, t2, I240};
      std::shared_ptr<Task199<T>> task199(new Task199<T>(tensor199, pindex));
      task196->add_dep(task199);
      task199->add_dep(task152);
      density_->add_task(task199);


      std::vector<std::shared_ptr<Tensor<T>>> tensor200 = {I240, Gamma2};
      std::shared_ptr<Task200<T>> task200(new Task200<T>(tensor200, pindex));
      task199->add_dep(task200);
      task200->add_dep(task152);
      density_->add_task(task200);

      task200->add_dep(task2);

      std::vector<IndexRange> I241_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor<T>> I241(new Tensor<T>(I241_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor201 = {den2, I241};
      std::shared_ptr<Task201<T>> task201(new Task201<T>(tensor201, pindex));
      task201->add_dep(task152);
      density_->add_task(task201);


      std::vector<IndexRange> I242_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I242(new Tensor<T>(I242_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor202 = {I241, t2, I242};
      std::shared_ptr<Task202<T>> task202(new Task202<T>(tensor202, pindex));
      task201->add_dep(task202);
      task202->add_dep(task152);
      density_->add_task(task202);


      std::vector<IndexRange> I243_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I243(new Tensor<T>(I243_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor203 = {I242, t2, I243};
      std::shared_ptr<Task203<T>> task203(new Task203<T>(tensor203, pindex));
      task202->add_dep(task203);
      task203->add_dep(task152);
      density_->add_task(task203);


      std::vector<std::shared_ptr<Tensor<T>>> tensor204 = {I243, Gamma14};
      std::shared_ptr<Task204<T>> task204(new Task204<T>(tensor204, pindex));
      task203->add_dep(task204);
      task204->add_dep(task152);
      density_->add_task(task204);

      task204->add_dep(task4);

      std::vector<IndexRange> I244_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor<T>> I244(new Tensor<T>(I244_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor205 = {den2, I244};
      std::shared_ptr<Task205<T>> task205(new Task205<T>(tensor205, pindex));
      task205->add_dep(task152);
      density_->add_task(task205);


      std::vector<IndexRange> I245_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I245(new Tensor<T>(I245_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor206 = {I244, t2, I245};
      std::shared_ptr<Task206<T>> task206(new Task206<T>(tensor206, pindex));
      task205->add_dep(task206);
      task206->add_dep(task152);
      density_->add_task(task206);


      std::vector<IndexRange> I246_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I246(new Tensor<T>(I246_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor207 = {I245, t2, I246};
      std::shared_ptr<Task207<T>> task207(new Task207<T>(tensor207, pindex));
      task206->add_dep(task207);
      task207->add_dep(task152);
      density_->add_task(task207);


      std::vector<std::shared_ptr<Tensor<T>>> tensor208 = {I246, Gamma14};
      std::shared_ptr<Task208<T>> task208(new Task208<T>(tensor208, pindex));
      task207->add_dep(task208);
      task208->add_dep(task152);
      density_->add_task(task208);

      task208->add_dep(task4);

      std::vector<IndexRange> I247_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I247(new Tensor<T>(I247_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor209 = {den2, I247};
      std::shared_ptr<Task209<T>> task209(new Task209<T>(tensor209, pindex));
      task209->add_dep(task152);
      density_->add_task(task209);


      std::vector<IndexRange> I248_index = {this->active_, this->active_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I248(new Tensor<T>(I248_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor210 = {I247, t2, I248};
      std::shared_ptr<Task210<T>> task210(new Task210<T>(tensor210, pindex));
      task209->add_dep(task210);
      task210->add_dep(task152);
      density_->add_task(task210);


      std::vector<IndexRange> I249_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I249(new Tensor<T>(I249_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor211 = {I248, t2, I249};
      std::shared_ptr<Task211<T>> task211(new Task211<T>(tensor211, pindex));
      task210->add_dep(task211);
      task211->add_dep(task152);
      density_->add_task(task211);


      std::vector<std::shared_ptr<Tensor<T>>> tensor212 = {I249, Gamma67};
      std::shared_ptr<Task212<T>> task212(new Task212<T>(tensor212, pindex));
      task211->add_dep(task212);
      task212->add_dep(task152);
      density_->add_task(task212);

      task212->add_dep(task6);

      std::vector<IndexRange> I250_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I250(new Tensor<T>(I250_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor213 = {den2, I250};
      std::shared_ptr<Task213<T>> task213(new Task213<T>(tensor213, pindex));
      task213->add_dep(task152);
      density_->add_task(task213);


      std::vector<IndexRange> I251_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I251(new Tensor<T>(I251_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor214 = {I250, t2, I251};
      std::shared_ptr<Task214<T>> task214(new Task214<T>(tensor214, pindex));
      task213->add_dep(task214);
      task214->add_dep(task152);
      density_->add_task(task214);


      std::vector<IndexRange> I252_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I252(new Tensor<T>(I252_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor215 = {I251, t2, I252};
      std::shared_ptr<Task215<T>> task215(new Task215<T>(tensor215, pindex));
      task214->add_dep(task215);
      task215->add_dep(task152);
      density_->add_task(task215);


      std::vector<std::shared_ptr<Tensor<T>>> tensor216 = {I252, Gamma14};
      std::shared_ptr<Task216<T>> task216(new Task216<T>(tensor216, pindex));
      task215->add_dep(task216);
      task216->add_dep(task152);
      density_->add_task(task216);

      task216->add_dep(task4);

      std::shared_ptr<Queue<T>> density1_(new Queue<T>());
      std::shared_ptr<Queue<T>> density2_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor217 = {Den1};
      std::shared_ptr<Task217<T>> task217(new Task217<T>(tensor217));
      density2_->add_task(task217);

      std::vector<IndexRange> I253_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I253(new Tensor<T>(I253_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor218 = {Den1, I253};
      std::shared_ptr<Task218<T>> task218(new Task218<T>(tensor218, pindex));
      task218->add_dep(task217);
      density2_->add_task(task218);


      std::vector<std::shared_ptr<Tensor<T>>> tensor219 = {I253, t2};
      std::shared_ptr<Task219<T>> task219(new Task219<T>(tensor219, pindex));
      task218->add_dep(task219);
      task219->add_dep(task217);
      density2_->add_task(task219);


      std::vector<IndexRange> I255_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I255(new Tensor<T>(I255_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor220 = {Den1, I255};
      std::shared_ptr<Task220<T>> task220(new Task220<T>(tensor220, pindex));
      task220->add_dep(task217);
      density2_->add_task(task220);


      std::vector<IndexRange> I256_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I256(new Tensor<T>(I256_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor221 = {I255, t2, I256};
      std::shared_ptr<Task221<T>> task221(new Task221<T>(tensor221, pindex));
      task220->add_dep(task221);
      task221->add_dep(task217);
      density2_->add_task(task221);


      std::vector<std::shared_ptr<Tensor<T>>> tensor222 = {I256, Gamma2};
      std::shared_ptr<Task222<T>> task222(new Task222<T>(tensor222, pindex));
      task221->add_dep(task222);
      task222->add_dep(task217);
      density2_->add_task(task222);

      task222->add_dep(task2);

      std::vector<IndexRange> I258_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I258(new Tensor<T>(I258_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor223 = {I255, t2, I258};
      std::shared_ptr<Task223<T>> task223(new Task223<T>(tensor223, pindex));
      task220->add_dep(task223);
      task223->add_dep(task217);
      density2_->add_task(task223);


      std::vector<std::shared_ptr<Tensor<T>>> tensor224 = {I258, Gamma2};
      std::shared_ptr<Task224<T>> task224(new Task224<T>(tensor224, pindex));
      task223->add_dep(task224);
      task224->add_dep(task217);
      density2_->add_task(task224);

      task224->add_dep(task2);

      std::vector<IndexRange> I259_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I259(new Tensor<T>(I259_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor225 = {Den1, I259};
      std::shared_ptr<Task225<T>> task225(new Task225<T>(tensor225, pindex));
      task225->add_dep(task217);
      density2_->add_task(task225);


      std::vector<IndexRange> I260_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I260(new Tensor<T>(I260_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor226 = {I259, t2, I260};
      std::shared_ptr<Task226<T>> task226(new Task226<T>(tensor226, pindex));
      task225->add_dep(task226);
      task226->add_dep(task217);
      density2_->add_task(task226);


      std::vector<std::shared_ptr<Tensor<T>>> tensor227 = {I260, Gamma14};
      std::shared_ptr<Task227<T>> task227(new Task227<T>(tensor227, pindex));
      task226->add_dep(task227);
      task227->add_dep(task217);
      density2_->add_task(task227);

      task227->add_dep(task4);

      std::shared_ptr<Queue<T>> dedci_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor228 = {deci};
      std::shared_ptr<Task228<T>> task228(new Task228<T>(tensor228));
      dedci_->add_task(task228);

      std::vector<IndexRange> I261_index = {this->ci_};
      std::shared_ptr<Tensor<T>> I261(new Tensor<T>(I261_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor229 = {deci, I261};
      std::shared_ptr<Task229<T>> task229(new Task229<T>(tensor229, cindex));
      task229->add_dep(task228);
      dedci_->add_task(task229);


      std::vector<IndexRange> I262_index = {this->ci_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I262(new Tensor<T>(I262_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor230 = {I261, t2, I262};
      std::shared_ptr<Task230<T>> task230(new Task230<T>(tensor230, cindex));
      task229->add_dep(task230);
      task230->add_dep(task228);
      dedci_->add_task(task230);


      std::vector<IndexRange> I263_index = {this->ci_};
      std::shared_ptr<Tensor<T>> I263(new Tensor<T>(I263_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor231 = {I262, t2, I263};
      std::shared_ptr<Task231<T>> task231(new Task231<T>(tensor231, cindex));
      task230->add_dep(task231);
      task231->add_dep(task228);
      dedci_->add_task(task231);


      std::vector<std::shared_ptr<Tensor<T>>> tensor232 = {I263, Gamma72};
      std::shared_ptr<Task232<T>> task232(new Task232<T>(tensor232, cindex));
      task231->add_dep(task232);
      task232->add_dep(task228);
      dedci_->add_task(task232);

      task232->add_dep(task7);
      task232->add_dep(task7);

      std::vector<IndexRange> I266_index = {this->ci_};
      std::shared_ptr<Tensor<T>> I266(new Tensor<T>(I266_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor233 = {I262, t2, I266};
      std::shared_ptr<Task233<T>> task233(new Task233<T>(tensor233, cindex));
      task230->add_dep(task233);
      task233->add_dep(task228);
      dedci_->add_task(task233);


      std::vector<std::shared_ptr<Tensor<T>>> tensor234 = {I266, Gamma72};
      std::shared_ptr<Task234<T>> task234(new Task234<T>(tensor234, cindex));
      task233->add_dep(task234);
      task234->add_dep(task228);
      dedci_->add_task(task234);

      task234->add_dep(task7);
      task234->add_dep(task7);

      std::vector<IndexRange> I269_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I269(new Tensor<T>(I269_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor235 = {I262, this->f1_, I269};
      std::shared_ptr<Task235<T>> task235(new Task235<T>(tensor235, cindex));
      task230->add_dep(task235);
      task235->add_dep(task228);
      dedci_->add_task(task235);


      std::vector<IndexRange> I270_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I270(new Tensor<T>(I270_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor236 = {I269, t2, I270};
      std::shared_ptr<Task236<T>> task236(new Task236<T>(tensor236, cindex));
      task235->add_dep(task236);
      task236->add_dep(task228);
      dedci_->add_task(task236);


      std::vector<std::shared_ptr<Tensor<T>>> tensor237 = {I270, Gamma74};
      std::shared_ptr<Task237<T>> task237(new Task237<T>(tensor237, cindex));
      task236->add_dep(task237);
      task237->add_dep(task228);
      dedci_->add_task(task237);

      task237->add_dep(task8);

      std::vector<IndexRange> I274_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I274(new Tensor<T>(I274_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor238 = {I269, t2, I274};
      std::shared_ptr<Task238<T>> task238(new Task238<T>(tensor238, cindex));
      task235->add_dep(task238);
      task238->add_dep(task228);
      dedci_->add_task(task238);


      std::vector<std::shared_ptr<Tensor<T>>> tensor239 = {I274, Gamma74};
      std::shared_ptr<Task239<T>> task239(new Task239<T>(tensor239, cindex));
      task238->add_dep(task239);
      task239->add_dep(task228);
      dedci_->add_task(task239);

      task239->add_dep(task8);

      std::vector<IndexRange> I336_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I336(new Tensor<T>(I336_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor240 = {I262, this->f1_, I336};
      std::shared_ptr<Task240<T>> task240(new Task240<T>(tensor240, cindex));
      task230->add_dep(task240);
      task240->add_dep(task228);
      dedci_->add_task(task240);


      std::vector<IndexRange> I337_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I337(new Tensor<T>(I337_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor241 = {I336, t2, I337};
      std::shared_ptr<Task241<T>> task241(new Task241<T>(tensor241, cindex));
      task240->add_dep(task241);
      task241->add_dep(task228);
      dedci_->add_task(task241);


      std::vector<std::shared_ptr<Tensor<T>>> tensor242 = {I337, Gamma74};
      std::shared_ptr<Task242<T>> task242(new Task242<T>(tensor242, cindex));
      task241->add_dep(task242);
      task242->add_dep(task228);
      dedci_->add_task(task242);

      task242->add_dep(task8);

      std::vector<IndexRange> I341_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I341(new Tensor<T>(I341_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor243 = {I336, t2, I341};
      std::shared_ptr<Task243<T>> task243(new Task243<T>(tensor243, cindex));
      task240->add_dep(task243);
      task243->add_dep(task228);
      dedci_->add_task(task243);


      std::vector<std::shared_ptr<Tensor<T>>> tensor244 = {I341, Gamma74};
      std::shared_ptr<Task244<T>> task244(new Task244<T>(tensor244, cindex));
      task243->add_dep(task244);
      task244->add_dep(task228);
      dedci_->add_task(task244);

      task244->add_dep(task8);

      std::vector<IndexRange> I276_index = {this->ci_, this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I276(new Tensor<T>(I276_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor245 = {I261, t2, I276};
      std::shared_ptr<Task245<T>> task245(new Task245<T>(tensor245, cindex));
      task229->add_dep(task245);
      task245->add_dep(task228);
      dedci_->add_task(task245);


      std::vector<IndexRange> I277_index = {this->ci_, this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I277(new Tensor<T>(I277_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor246 = {I276, this->f1_, I277};
      std::shared_ptr<Task246<T>> task246(new Task246<T>(tensor246, cindex));
      task245->add_dep(task246);
      task246->add_dep(task228);
      dedci_->add_task(task246);


      std::vector<IndexRange> I278_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I278(new Tensor<T>(I278_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor247 = {I277, t2, I278};
      std::shared_ptr<Task247<T>> task247(new Task247<T>(tensor247, cindex));
      task246->add_dep(task247);
      task247->add_dep(task228);
      dedci_->add_task(task247);


      std::vector<std::shared_ptr<Tensor<T>>> tensor248 = {I278, Gamma74};
      std::shared_ptr<Task248<T>> task248(new Task248<T>(tensor248, cindex));
      task247->add_dep(task248);
      task248->add_dep(task228);
      dedci_->add_task(task248);

      task248->add_dep(task8);

      std::vector<IndexRange> I282_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I282(new Tensor<T>(I282_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor249 = {I277, t2, I282};
      std::shared_ptr<Task249<T>> task249(new Task249<T>(tensor249, cindex));
      task246->add_dep(task249);
      task249->add_dep(task228);
      dedci_->add_task(task249);


      std::vector<std::shared_ptr<Tensor<T>>> tensor250 = {I282, Gamma74};
      std::shared_ptr<Task250<T>> task250(new Task250<T>(tensor250, cindex));
      task249->add_dep(task250);
      task250->add_dep(task228);
      dedci_->add_task(task250);

      task250->add_dep(task8);

      std::vector<IndexRange> I285_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I285(new Tensor<T>(I285_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor251 = {I276, t2, I285};
      std::shared_ptr<Task251<T>> task251(new Task251<T>(tensor251, cindex));
      task245->add_dep(task251);
      task251->add_dep(task228);
      dedci_->add_task(task251);


      std::vector<std::shared_ptr<Tensor<T>>> tensor252 = {I285, Gamma78};
      std::shared_ptr<Task252<T>> task252(new Task252<T>(tensor252, cindex));
      task251->add_dep(task252);
      task252->add_dep(task228);
      dedci_->add_task(task252);

      task252->add_dep(task9);

      std::vector<IndexRange> I288_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I288(new Tensor<T>(I288_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor253 = {I276, t2, I288};
      std::shared_ptr<Task253<T>> task253(new Task253<T>(tensor253, cindex));
      task245->add_dep(task253);
      task253->add_dep(task228);
      dedci_->add_task(task253);


      std::vector<std::shared_ptr<Tensor<T>>> tensor254 = {I288, Gamma78};
      std::shared_ptr<Task254<T>> task254(new Task254<T>(tensor254, cindex));
      task253->add_dep(task254);
      task254->add_dep(task228);
      dedci_->add_task(task254);

      task254->add_dep(task9);

      std::vector<IndexRange> I291_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I291(new Tensor<T>(I291_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor255 = {I276, this->f1_, I291};
      std::shared_ptr<Task255<T>> task255(new Task255<T>(tensor255, cindex));
      task245->add_dep(task255);
      task255->add_dep(task228);
      dedci_->add_task(task255);


      std::vector<IndexRange> I292_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I292(new Tensor<T>(I292_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor256 = {I291, t2, I292};
      std::shared_ptr<Task256<T>> task256(new Task256<T>(tensor256, cindex));
      task255->add_dep(task256);
      task256->add_dep(task228);
      dedci_->add_task(task256);


      std::vector<std::shared_ptr<Tensor<T>>> tensor257 = {I292, Gamma74};
      std::shared_ptr<Task257<T>> task257(new Task257<T>(tensor257, cindex));
      task256->add_dep(task257);
      task257->add_dep(task228);
      dedci_->add_task(task257);

      task257->add_dep(task8);

      std::vector<IndexRange> I296_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I296(new Tensor<T>(I296_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor258 = {I291, t2, I296};
      std::shared_ptr<Task258<T>> task258(new Task258<T>(tensor258, cindex));
      task255->add_dep(task258);
      task258->add_dep(task228);
      dedci_->add_task(task258);


      std::vector<std::shared_ptr<Tensor<T>>> tensor259 = {I296, Gamma74};
      std::shared_ptr<Task259<T>> task259(new Task259<T>(tensor259, cindex));
      task258->add_dep(task259);
      task259->add_dep(task228);
      dedci_->add_task(task259);

      task259->add_dep(task8);

      std::vector<IndexRange> I299_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I299(new Tensor<T>(I299_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor260 = {I276, this->f1_, I299};
      std::shared_ptr<Task260<T>> task260(new Task260<T>(tensor260, cindex));
      task245->add_dep(task260);
      task260->add_dep(task228);
      dedci_->add_task(task260);


      std::vector<IndexRange> I300_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I300(new Tensor<T>(I300_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor261 = {I299, t2, I300};
      std::shared_ptr<Task261<T>> task261(new Task261<T>(tensor261, cindex));
      task260->add_dep(task261);
      task261->add_dep(task228);
      dedci_->add_task(task261);


      std::vector<std::shared_ptr<Tensor<T>>> tensor262 = {I300, Gamma74};
      std::shared_ptr<Task262<T>> task262(new Task262<T>(tensor262, cindex));
      task261->add_dep(task262);
      task262->add_dep(task228);
      dedci_->add_task(task262);

      task262->add_dep(task8);

      std::vector<IndexRange> I304_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I304(new Tensor<T>(I304_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor263 = {I299, t2, I304};
      std::shared_ptr<Task263<T>> task263(new Task263<T>(tensor263, cindex));
      task260->add_dep(task263);
      task263->add_dep(task228);
      dedci_->add_task(task263);


      std::vector<std::shared_ptr<Tensor<T>>> tensor264 = {I304, Gamma74};
      std::shared_ptr<Task264<T>> task264(new Task264<T>(tensor264, cindex));
      task263->add_dep(task264);
      task264->add_dep(task228);
      dedci_->add_task(task264);

      task264->add_dep(task8);

      std::vector<IndexRange> I307_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I307(new Tensor<T>(I307_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor265 = {I276, this->f1_, I307};
      std::shared_ptr<Task265<T>> task265(new Task265<T>(tensor265, cindex));
      task245->add_dep(task265);
      task265->add_dep(task228);
      dedci_->add_task(task265);


      std::vector<IndexRange> I308_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I308(new Tensor<T>(I308_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor266 = {I307, t2, I308};
      std::shared_ptr<Task266<T>> task266(new Task266<T>(tensor266, cindex));
      task265->add_dep(task266);
      task266->add_dep(task228);
      dedci_->add_task(task266);


      std::vector<std::shared_ptr<Tensor<T>>> tensor267 = {I308, Gamma74};
      std::shared_ptr<Task267<T>> task267(new Task267<T>(tensor267, cindex));
      task266->add_dep(task267);
      task267->add_dep(task228);
      dedci_->add_task(task267);

      task267->add_dep(task8);

      std::vector<IndexRange> I312_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I312(new Tensor<T>(I312_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor268 = {I307, t2, I312};
      std::shared_ptr<Task268<T>> task268(new Task268<T>(tensor268, cindex));
      task265->add_dep(task268);
      task268->add_dep(task228);
      dedci_->add_task(task268);


      std::vector<std::shared_ptr<Tensor<T>>> tensor269 = {I312, Gamma74};
      std::shared_ptr<Task269<T>> task269(new Task269<T>(tensor269, cindex));
      task268->add_dep(task269);
      task269->add_dep(task228);
      dedci_->add_task(task269);

      task269->add_dep(task8);

      std::vector<IndexRange> I315_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I315(new Tensor<T>(I315_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor270 = {I276, this->f1_, I315};
      std::shared_ptr<Task270<T>> task270(new Task270<T>(tensor270, cindex));
      task245->add_dep(task270);
      task270->add_dep(task228);
      dedci_->add_task(task270);


      std::vector<IndexRange> I316_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I316(new Tensor<T>(I316_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor271 = {I315, t2, I316};
      std::shared_ptr<Task271<T>> task271(new Task271<T>(tensor271, cindex));
      task270->add_dep(task271);
      task271->add_dep(task228);
      dedci_->add_task(task271);


      std::vector<std::shared_ptr<Tensor<T>>> tensor272 = {I316, Gamma86};
      std::shared_ptr<Task272<T>> task272(new Task272<T>(tensor272, cindex));
      task271->add_dep(task272);
      task272->add_dep(task228);
      dedci_->add_task(task272);

      task272->add_dep(task10);

      std::vector<IndexRange> I397_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I397(new Tensor<T>(I397_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor273 = {I276, t2, I397};
      std::shared_ptr<Task273<T>> task273(new Task273<T>(tensor273, cindex));
      task245->add_dep(task273);
      task273->add_dep(task228);
      dedci_->add_task(task273);


      std::vector<std::shared_ptr<Tensor<T>>> tensor274 = {I397, Gamma74};
      std::shared_ptr<Task274<T>> task274(new Task274<T>(tensor274, cindex, this->e0_));
      task273->add_dep(task274);
      task274->add_dep(task228);
      dedci_->add_task(task274);

      task274->add_dep(task8);

      std::vector<IndexRange> I400_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I400(new Tensor<T>(I400_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor275 = {I276, t2, I400};
      std::shared_ptr<Task275<T>> task275(new Task275<T>(tensor275, cindex));
      task245->add_dep(task275);
      task275->add_dep(task228);
      dedci_->add_task(task275);


      std::vector<std::shared_ptr<Tensor<T>>> tensor276 = {I400, Gamma74};
      std::shared_ptr<Task276<T>> task276(new Task276<T>(tensor276, cindex, this->e0_));
      task275->add_dep(task276);
      task276->add_dep(task228);
      dedci_->add_task(task276);

      task276->add_dep(task8);

      std::vector<IndexRange> I415_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I415(new Tensor<T>(I415_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor277 = {I276, this->v2_, I415};
      std::shared_ptr<Task277<T>> task277(new Task277<T>(tensor277, cindex));
      task245->add_dep(task277);
      task277->add_dep(task228);
      dedci_->add_task(task277);


      std::vector<std::shared_ptr<Tensor<T>>> tensor278 = {I415, Gamma74};
      std::shared_ptr<Task278<T>> task278(new Task278<T>(tensor278, cindex));
      task277->add_dep(task278);
      task278->add_dep(task228);
      dedci_->add_task(task278);

      task278->add_dep(task8);

      std::vector<IndexRange> I418_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I418(new Tensor<T>(I418_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor279 = {I276, this->v2_, I418};
      std::shared_ptr<Task279<T>> task279(new Task279<T>(tensor279, cindex));
      task245->add_dep(task279);
      task279->add_dep(task228);
      dedci_->add_task(task279);


      std::vector<std::shared_ptr<Tensor<T>>> tensor280 = {I418, Gamma74};
      std::shared_ptr<Task280<T>> task280(new Task280<T>(tensor280, cindex));
      task279->add_dep(task280);
      task280->add_dep(task228);
      dedci_->add_task(task280);

      task280->add_dep(task8);

      std::vector<IndexRange> I318_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I318(new Tensor<T>(I318_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor281 = {I261, t2, I318};
      std::shared_ptr<Task281<T>> task281(new Task281<T>(tensor281, cindex));
      task229->add_dep(task281);
      task281->add_dep(task228);
      dedci_->add_task(task281);


      std::vector<IndexRange> I319_index = {this->ci_, this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I319(new Tensor<T>(I319_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor282 = {I318, this->f1_, I319};
      std::shared_ptr<Task282<T>> task282(new Task282<T>(tensor282, cindex));
      task281->add_dep(task282);
      task282->add_dep(task228);
      dedci_->add_task(task282);


      std::vector<IndexRange> I320_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I320(new Tensor<T>(I320_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor283 = {I319, t2, I320};
      std::shared_ptr<Task283<T>> task283(new Task283<T>(tensor283, cindex));
      task282->add_dep(task283);
      task283->add_dep(task228);
      dedci_->add_task(task283);


      std::vector<std::shared_ptr<Tensor<T>>> tensor284 = {I320, Gamma86};
      std::shared_ptr<Task284<T>> task284(new Task284<T>(tensor284, cindex));
      task283->add_dep(task284);
      task284->add_dep(task228);
      dedci_->add_task(task284);

      task284->add_dep(task10);

      std::vector<IndexRange> I323_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I323(new Tensor<T>(I323_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor285 = {I318, t2, I323};
      std::shared_ptr<Task285<T>> task285(new Task285<T>(tensor285, cindex));
      task281->add_dep(task285);
      task285->add_dep(task228);
      dedci_->add_task(task285);


      std::vector<std::shared_ptr<Tensor<T>>> tensor286 = {I323, Gamma88};
      std::shared_ptr<Task286<T>> task286(new Task286<T>(tensor286, cindex));
      task285->add_dep(task286);
      task286->add_dep(task228);
      dedci_->add_task(task286);

      task286->add_dep(task11);

      std::vector<IndexRange> I326_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I326(new Tensor<T>(I326_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor287 = {I318, this->f1_, I326};
      std::shared_ptr<Task287<T>> task287(new Task287<T>(tensor287, cindex));
      task281->add_dep(task287);
      task287->add_dep(task228);
      dedci_->add_task(task287);


      std::vector<IndexRange> I327_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I327(new Tensor<T>(I327_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor288 = {I326, t2, I327};
      std::shared_ptr<Task288<T>> task288(new Task288<T>(tensor288, cindex));
      task287->add_dep(task288);
      task288->add_dep(task228);
      dedci_->add_task(task288);


      std::vector<std::shared_ptr<Tensor<T>>> tensor289 = {I327, Gamma86};
      std::shared_ptr<Task289<T>> task289(new Task289<T>(tensor289, cindex));
      task288->add_dep(task289);
      task289->add_dep(task228);
      dedci_->add_task(task289);

      task289->add_dep(task10);

      std::vector<IndexRange> I403_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I403(new Tensor<T>(I403_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor290 = {I318, t2, I403};
      std::shared_ptr<Task290<T>> task290(new Task290<T>(tensor290, cindex));
      task281->add_dep(task290);
      task290->add_dep(task228);
      dedci_->add_task(task290);


      std::vector<std::shared_ptr<Tensor<T>>> tensor291 = {I403, Gamma86};
      std::shared_ptr<Task291<T>> task291(new Task291<T>(tensor291, cindex, this->e0_));
      task290->add_dep(task291);
      task291->add_dep(task228);
      dedci_->add_task(task291);

      task291->add_dep(task10);

      std::vector<IndexRange> I421_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I421(new Tensor<T>(I421_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor292 = {I318, this->v2_, I421};
      std::shared_ptr<Task292<T>> task292(new Task292<T>(tensor292, cindex));
      task281->add_dep(task292);
      task292->add_dep(task228);
      dedci_->add_task(task292);


      std::vector<std::shared_ptr<Tensor<T>>> tensor293 = {I421, Gamma86};
      std::shared_ptr<Task293<T>> task293(new Task293<T>(tensor293, cindex));
      task292->add_dep(task293);
      task293->add_dep(task228);
      dedci_->add_task(task293);

      task293->add_dep(task10);

      std::vector<IndexRange> I343_index = {this->ci_, this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I343(new Tensor<T>(I343_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor294 = {I261, t2, I343};
      std::shared_ptr<Task294<T>> task294(new Task294<T>(tensor294, cindex));
      task229->add_dep(task294);
      task294->add_dep(task228);
      dedci_->add_task(task294);


      std::vector<IndexRange> I344_index = {this->ci_, this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I344(new Tensor<T>(I344_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor295 = {I343, this->f1_, I344};
      std::shared_ptr<Task295<T>> task295(new Task295<T>(tensor295, cindex));
      task294->add_dep(task295);
      task295->add_dep(task228);
      dedci_->add_task(task295);


      std::vector<IndexRange> I345_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I345(new Tensor<T>(I345_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor296 = {I344, t2, I345};
      std::shared_ptr<Task296<T>> task296(new Task296<T>(tensor296, cindex));
      task295->add_dep(task296);
      task296->add_dep(task228);
      dedci_->add_task(task296);


      std::vector<std::shared_ptr<Tensor<T>>> tensor297 = {I345, Gamma74};
      std::shared_ptr<Task297<T>> task297(new Task297<T>(tensor297, cindex));
      task296->add_dep(task297);
      task297->add_dep(task228);
      dedci_->add_task(task297);

      task297->add_dep(task8);

      std::vector<IndexRange> I349_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I349(new Tensor<T>(I349_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor298 = {I344, t2, I349};
      std::shared_ptr<Task298<T>> task298(new Task298<T>(tensor298, cindex));
      task295->add_dep(task298);
      task298->add_dep(task228);
      dedci_->add_task(task298);


      std::vector<std::shared_ptr<Tensor<T>>> tensor299 = {I349, Gamma74};
      std::shared_ptr<Task299<T>> task299(new Task299<T>(tensor299, cindex));
      task298->add_dep(task299);
      task299->add_dep(task228);
      dedci_->add_task(task299);

      task299->add_dep(task8);

      std::vector<IndexRange> I358_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I358(new Tensor<T>(I358_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor300 = {I343, this->f1_, I358};
      std::shared_ptr<Task300<T>> task300(new Task300<T>(tensor300, cindex));
      task294->add_dep(task300);
      task300->add_dep(task228);
      dedci_->add_task(task300);


      std::vector<IndexRange> I359_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I359(new Tensor<T>(I359_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor301 = {I358, t2, I359};
      std::shared_ptr<Task301<T>> task301(new Task301<T>(tensor301, cindex));
      task300->add_dep(task301);
      task301->add_dep(task228);
      dedci_->add_task(task301);


      std::vector<std::shared_ptr<Tensor<T>>> tensor302 = {I359, Gamma74};
      std::shared_ptr<Task302<T>> task302(new Task302<T>(tensor302, cindex));
      task301->add_dep(task302);
      task302->add_dep(task228);
      dedci_->add_task(task302);

      task302->add_dep(task8);

      std::vector<IndexRange> I363_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I363(new Tensor<T>(I363_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor303 = {I358, t2, I363};
      std::shared_ptr<Task303<T>> task303(new Task303<T>(tensor303, cindex));
      task300->add_dep(task303);
      task303->add_dep(task228);
      dedci_->add_task(task303);


      std::vector<std::shared_ptr<Tensor<T>>> tensor304 = {I363, Gamma74};
      std::shared_ptr<Task304<T>> task304(new Task304<T>(tensor304, cindex));
      task303->add_dep(task304);
      task304->add_dep(task228);
      dedci_->add_task(task304);

      task304->add_dep(task8);

      std::vector<IndexRange> I366_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I366(new Tensor<T>(I366_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor305 = {I343, this->f1_, I366};
      std::shared_ptr<Task305<T>> task305(new Task305<T>(tensor305, cindex));
      task294->add_dep(task305);
      task305->add_dep(task228);
      dedci_->add_task(task305);


      std::vector<IndexRange> I367_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I367(new Tensor<T>(I367_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor306 = {I366, t2, I367};
      std::shared_ptr<Task306<T>> task306(new Task306<T>(tensor306, cindex));
      task305->add_dep(task306);
      task306->add_dep(task228);
      dedci_->add_task(task306);


      std::vector<std::shared_ptr<Tensor<T>>> tensor307 = {I367, Gamma74};
      std::shared_ptr<Task307<T>> task307(new Task307<T>(tensor307, cindex));
      task306->add_dep(task307);
      task307->add_dep(task228);
      dedci_->add_task(task307);

      task307->add_dep(task8);

      std::vector<IndexRange> I371_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I371(new Tensor<T>(I371_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor308 = {I366, t2, I371};
      std::shared_ptr<Task308<T>> task308(new Task308<T>(tensor308, cindex));
      task305->add_dep(task308);
      task308->add_dep(task228);
      dedci_->add_task(task308);


      std::vector<std::shared_ptr<Tensor<T>>> tensor309 = {I371, Gamma74};
      std::shared_ptr<Task309<T>> task309(new Task309<T>(tensor309, cindex));
      task308->add_dep(task309);
      task309->add_dep(task228);
      dedci_->add_task(task309);

      task309->add_dep(task8);

      std::vector<IndexRange> I374_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I374(new Tensor<T>(I374_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor310 = {I343, this->f1_, I374};
      std::shared_ptr<Task310<T>> task310(new Task310<T>(tensor310, cindex));
      task294->add_dep(task310);
      task310->add_dep(task228);
      dedci_->add_task(task310);


      std::vector<IndexRange> I375_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I375(new Tensor<T>(I375_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor311 = {I374, t2, I375};
      std::shared_ptr<Task311<T>> task311(new Task311<T>(tensor311, cindex));
      task310->add_dep(task311);
      task311->add_dep(task228);
      dedci_->add_task(task311);


      std::vector<std::shared_ptr<Tensor<T>>> tensor312 = {I375, Gamma74};
      std::shared_ptr<Task312<T>> task312(new Task312<T>(tensor312, cindex));
      task311->add_dep(task312);
      task312->add_dep(task228);
      dedci_->add_task(task312);

      task312->add_dep(task8);

      std::vector<IndexRange> I379_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I379(new Tensor<T>(I379_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor313 = {I374, t2, I379};
      std::shared_ptr<Task313<T>> task313(new Task313<T>(tensor313, cindex));
      task310->add_dep(task313);
      task313->add_dep(task228);
      dedci_->add_task(task313);


      std::vector<std::shared_ptr<Tensor<T>>> tensor314 = {I379, Gamma74};
      std::shared_ptr<Task314<T>> task314(new Task314<T>(tensor314, cindex));
      task313->add_dep(task314);
      task314->add_dep(task228);
      dedci_->add_task(task314);

      task314->add_dep(task8);

      std::vector<IndexRange> I406_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I406(new Tensor<T>(I406_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor315 = {I343, t2, I406};
      std::shared_ptr<Task315<T>> task315(new Task315<T>(tensor315, cindex));
      task294->add_dep(task315);
      task315->add_dep(task228);
      dedci_->add_task(task315);


      std::vector<std::shared_ptr<Tensor<T>>> tensor316 = {I406, Gamma74};
      std::shared_ptr<Task316<T>> task316(new Task316<T>(tensor316, cindex, this->e0_));
      task315->add_dep(task316);
      task316->add_dep(task228);
      dedci_->add_task(task316);

      task316->add_dep(task8);

      std::vector<IndexRange> I409_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I409(new Tensor<T>(I409_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor317 = {I343, t2, I409};
      std::shared_ptr<Task317<T>> task317(new Task317<T>(tensor317, cindex));
      task294->add_dep(task317);
      task317->add_dep(task228);
      dedci_->add_task(task317);


      std::vector<std::shared_ptr<Tensor<T>>> tensor318 = {I409, Gamma74};
      std::shared_ptr<Task318<T>> task318(new Task318<T>(tensor318, cindex, this->e0_));
      task317->add_dep(task318);
      task318->add_dep(task228);
      dedci_->add_task(task318);

      task318->add_dep(task8);

      std::vector<IndexRange> I424_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I424(new Tensor<T>(I424_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor319 = {I343, this->v2_, I424};
      std::shared_ptr<Task319<T>> task319(new Task319<T>(tensor319, cindex));
      task294->add_dep(task319);
      task319->add_dep(task228);
      dedci_->add_task(task319);


      std::vector<std::shared_ptr<Tensor<T>>> tensor320 = {I424, Gamma74};
      std::shared_ptr<Task320<T>> task320(new Task320<T>(tensor320, cindex));
      task319->add_dep(task320);
      task320->add_dep(task228);
      dedci_->add_task(task320);

      task320->add_dep(task8);

      std::vector<IndexRange> I427_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I427(new Tensor<T>(I427_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor321 = {I343, this->v2_, I427};
      std::shared_ptr<Task321<T>> task321(new Task321<T>(tensor321, cindex));
      task294->add_dep(task321);
      task321->add_dep(task228);
      dedci_->add_task(task321);


      std::vector<std::shared_ptr<Tensor<T>>> tensor322 = {I427, Gamma74};
      std::shared_ptr<Task322<T>> task322(new Task322<T>(tensor322, cindex));
      task321->add_dep(task322);
      task322->add_dep(task228);
      dedci_->add_task(task322);

      task322->add_dep(task8);

      std::vector<IndexRange> I351_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I351(new Tensor<T>(I351_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor323 = {I261, t2, I351};
      std::shared_ptr<Task323<T>> task323(new Task323<T>(tensor323, cindex));
      task229->add_dep(task323);
      task323->add_dep(task228);
      dedci_->add_task(task323);


      std::vector<IndexRange> I352_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I352(new Tensor<T>(I352_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor324 = {I351, t2, I352};
      std::shared_ptr<Task324<T>> task324(new Task324<T>(tensor324, cindex));
      task323->add_dep(task324);
      task324->add_dep(task228);
      dedci_->add_task(task324);


      std::vector<std::shared_ptr<Tensor<T>>> tensor325 = {I352, Gamma78};
      std::shared_ptr<Task325<T>> task325(new Task325<T>(tensor325, cindex));
      task324->add_dep(task325);
      task325->add_dep(task228);
      dedci_->add_task(task325);

      task325->add_dep(task9);

      std::vector<IndexRange> I355_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I355(new Tensor<T>(I355_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor326 = {I351, t2, I355};
      std::shared_ptr<Task326<T>> task326(new Task326<T>(tensor326, cindex));
      task323->add_dep(task326);
      task326->add_dep(task228);
      dedci_->add_task(task326);


      std::vector<std::shared_ptr<Tensor<T>>> tensor327 = {I355, Gamma78};
      std::shared_ptr<Task327<T>> task327(new Task327<T>(tensor327, cindex));
      task326->add_dep(task327);
      task327->add_dep(task228);
      dedci_->add_task(task327);

      task327->add_dep(task9);

      std::vector<IndexRange> I382_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I382(new Tensor<T>(I382_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor328 = {I351, this->f1_, I382};
      std::shared_ptr<Task328<T>> task328(new Task328<T>(tensor328, cindex));
      task323->add_dep(task328);
      task328->add_dep(task228);
      dedci_->add_task(task328);


      std::vector<IndexRange> I383_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I383(new Tensor<T>(I383_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor329 = {I382, t2, I383};
      std::shared_ptr<Task329<T>> task329(new Task329<T>(tensor329, cindex));
      task328->add_dep(task329);
      task329->add_dep(task228);
      dedci_->add_task(task329);


      std::vector<std::shared_ptr<Tensor<T>>> tensor330 = {I383, Gamma86};
      std::shared_ptr<Task330<T>> task330(new Task330<T>(tensor330, cindex));
      task329->add_dep(task330);
      task330->add_dep(task228);
      dedci_->add_task(task330);

      task330->add_dep(task10);

      std::vector<IndexRange> I385_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I385(new Tensor<T>(I385_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor331 = {I261, t2, I385};
      std::shared_ptr<Task331<T>> task331(new Task331<T>(tensor331, cindex));
      task229->add_dep(task331);
      task331->add_dep(task228);
      dedci_->add_task(task331);


      std::vector<IndexRange> I386_index = {this->ci_, this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I386(new Tensor<T>(I386_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor332 = {I385, this->f1_, I386};
      std::shared_ptr<Task332<T>> task332(new Task332<T>(tensor332, cindex));
      task331->add_dep(task332);
      task332->add_dep(task228);
      dedci_->add_task(task332);


      std::vector<IndexRange> I387_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I387(new Tensor<T>(I387_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor333 = {I386, t2, I387};
      std::shared_ptr<Task333<T>> task333(new Task333<T>(tensor333, cindex));
      task332->add_dep(task333);
      task333->add_dep(task228);
      dedci_->add_task(task333);


      std::vector<std::shared_ptr<Tensor<T>>> tensor334 = {I387, Gamma86};
      std::shared_ptr<Task334<T>> task334(new Task334<T>(tensor334, cindex));
      task333->add_dep(task334);
      task334->add_dep(task228);
      dedci_->add_task(task334);

      task334->add_dep(task10);

      std::vector<IndexRange> I393_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I393(new Tensor<T>(I393_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor335 = {I385, this->f1_, I393};
      std::shared_ptr<Task335<T>> task335(new Task335<T>(tensor335, cindex));
      task331->add_dep(task335);
      task335->add_dep(task228);
      dedci_->add_task(task335);


      std::vector<IndexRange> I394_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I394(new Tensor<T>(I394_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor336 = {I393, t2, I394};
      std::shared_ptr<Task336<T>> task336(new Task336<T>(tensor336, cindex));
      task335->add_dep(task336);
      task336->add_dep(task228);
      dedci_->add_task(task336);


      std::vector<std::shared_ptr<Tensor<T>>> tensor337 = {I394, Gamma86};
      std::shared_ptr<Task337<T>> task337(new Task337<T>(tensor337, cindex));
      task336->add_dep(task337);
      task337->add_dep(task228);
      dedci_->add_task(task337);

      task337->add_dep(task10);

      std::vector<IndexRange> I412_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I412(new Tensor<T>(I412_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor338 = {I385, t2, I412};
      std::shared_ptr<Task338<T>> task338(new Task338<T>(tensor338, cindex));
      task331->add_dep(task338);
      task338->add_dep(task228);
      dedci_->add_task(task338);


      std::vector<std::shared_ptr<Tensor<T>>> tensor339 = {I412, Gamma86};
      std::shared_ptr<Task339<T>> task339(new Task339<T>(tensor339, cindex, this->e0_));
      task338->add_dep(task339);
      task339->add_dep(task228);
      dedci_->add_task(task339);

      task339->add_dep(task10);

      std::vector<IndexRange> I430_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I430(new Tensor<T>(I430_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor340 = {I385, this->v2_, I430};
      std::shared_ptr<Task340<T>> task340(new Task340<T>(tensor340, cindex));
      task331->add_dep(task340);
      task340->add_dep(task228);
      dedci_->add_task(task340);


      std::vector<std::shared_ptr<Tensor<T>>> tensor341 = {I430, Gamma86};
      std::shared_ptr<Task341<T>> task341(new Task341<T>(tensor341, cindex));
      task340->add_dep(task341);
      task341->add_dep(task228);
      dedci_->add_task(task341);

      task341->add_dep(task10);

      std::vector<IndexRange> I389_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I389(new Tensor<T>(I389_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor342 = {I261, t2, I389};
      std::shared_ptr<Task342<T>> task342(new Task342<T>(tensor342, cindex));
      task229->add_dep(task342);
      task342->add_dep(task228);
      dedci_->add_task(task342);


      std::vector<IndexRange> I390_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I390(new Tensor<T>(I390_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor343 = {I389, t2, I390};
      std::shared_ptr<Task343<T>> task343(new Task343<T>(tensor343, cindex));
      task342->add_dep(task343);
      task343->add_dep(task228);
      dedci_->add_task(task343);


      std::vector<std::shared_ptr<Tensor<T>>> tensor344 = {I390, Gamma88};
      std::shared_ptr<Task344<T>> task344(new Task344<T>(tensor344, cindex));
      task343->add_dep(task344);
      task344->add_dep(task228);
      dedci_->add_task(task344);

      task344->add_dep(task11);

      return make_tuple(queue_, energy_, correction_, density_, density1_, density2_, dedci_);
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
      den2 = this->h1_->clone();
      Den1 = this->v2_->clone();
      deci = this->rdm0deriv_->clone();
    };
    ~CAS_test() {};

    void solve() {
      this->print_iteration();
      int iter = 0;
      std::shared_ptr<Queue<T>> queue, energ, correct, dens2, dens1, Dens1, dec;
      for ( ; iter != ref_->maxiter(); ++iter) {
        std::tie(queue, energ, correct, dens2, dens1, Dens1, dec) = make_queue_();
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

      std::cout << " === Computing correlated overlap, <1|1> ===" << std::endl;
      correlated_norm = correction(correct);
      std::cout << std::endl;
      std::cout << "      Norm  = " << std::setprecision(10) << correlated_norm << std::endl;
      std::cout << std::endl;

      std::cout << " === Computing unrelaxed one-body density matrix, dm2, <1|E_pq|1>  ===" << std::endl;
      while (!dens2->done())
        dens2->next_compute();
#if 0
      den2->print2("smith d1 correlated one-body density matrix dm2 ", 1.0e-5);
#endif
      std::cout << std::endl;
      std::cout << " === Computing unrelaxed one-body density matrix, dm1, 2<0|E_pq|1> ===" << std::endl;
      while (!dens1->done())
        dens1->next_compute();
#if 0
      den1->print2("smith d1 correlated one-body density matrix dm1", 1.0e-5);
#endif
      std::cout << std::endl;
      std::cout << " === Computing unrelaxed two-body density matrix, D1, <0|E_pqrs|1>  ===" << std::endl;
      while (!Dens1->done())
        Dens1->next_compute();
#if 0
      Den1->print4("smith d2 correlated two-body density matrix D1", 1.0e-5);
#endif
      std::cout << std::endl;

      std::cout << " === Computing cI derivative dE/dcI ===" << std::endl;
      while (!dec->done())
        dec->next_compute();
      deci->print1("cI derivative tensor: ", 1.0e-15);
      std::cout << std::endl;
      std::cout << "      cI derivative * cI    = " << std::setprecision(10) <<  deci->dot_product(this->rdm0deriv_) << std::endl;
      std::cout << "      Expecting 2E          = " << std::setprecision(10) <<  2.0*this->energy_ << std::endl;
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

    std::shared_ptr<const Matrix> rdm11() const { return den1->matrix(); }
    std::shared_ptr<const Matrix> rdm12() const { return den2->matrix(); }
    std::shared_ptr<const Matrix> rdm21() const { return Den1->matrix2(); }

    double rdm1_correction() const { return correlated_norm; }

    std::shared_ptr<const Civec> ci_deriv() const { return deci->civec(this->det_); }

};

}
}
}
#endif

