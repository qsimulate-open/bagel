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
class CAS_test : public SpinFreeMethod<T>, SMITH_info {
  protected:
    std::shared_ptr<Tensor<T>> t2;
    std::shared_ptr<Tensor<T>> r;
    double e0_;

    std::tuple<std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>> make_queue_() {
      std::shared_ptr<Queue<T>> queue_(new Queue<T>());
      std::array<std::shared_ptr<const IndexRange>,3> pindex = {{this->rclosed_, this->ractive_, this->rvirt_}};

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

      std::vector<IndexRange> Gamma27_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma27(new Tensor<T>(Gamma27_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor6 = {Gamma27, this->rdm1_};
      std::shared_ptr<Task6<T>> task6(new Task6<T>(tensor6, pindex));
      task6->add_dep(task0);
      queue_->add_task(task6);

      std::vector<IndexRange> Gamma29_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma29(new Tensor<T>(Gamma29_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor7 = {Gamma29, this->rdm2_};
      std::shared_ptr<Task7<T>> task7(new Task7<T>(tensor7, pindex));
      task7->add_dep(task0);
      queue_->add_task(task7);

      std::vector<IndexRange> Gamma52_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma52(new Tensor<T>(Gamma52_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor8 = {Gamma52, this->rdm3_};
      std::shared_ptr<Task8<T>> task8(new Task8<T>(tensor8, pindex));
      task8->add_dep(task0);
      queue_->add_task(task8);

      std::vector<IndexRange> I0_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor9 = {r, I0};
      std::shared_ptr<Task9<T>> task9(new Task9<T>(tensor9, pindex));
      task9->add_dep(task0);
      queue_->add_task(task9);


      std::vector<std::shared_ptr<Tensor<T>>> tensor10 = {I0, t2, this->v2_};
      std::shared_ptr<Task10<T>> task10(new Task10<T>(tensor10, pindex, this->e0_));
      task9->add_dep(task10);
      task10->add_dep(task0);
      queue_->add_task(task10);


      std::vector<IndexRange> I1_index;
      std::shared_ptr<Tensor<T>> I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor11 = {I0, t2, I1};
      std::shared_ptr<Task11<T>> task11(new Task11<T>(tensor11, pindex));
      task9->add_dep(task11);
      task11->add_dep(task0);
      queue_->add_task(task11);


      std::vector<std::shared_ptr<Tensor<T>>> tensor12 = {I1, Gamma0};
      std::shared_ptr<Task12<T>> task12(new Task12<T>(tensor12, pindex));
      task11->add_dep(task12);
      task12->add_dep(task0);
      queue_->add_task(task12);

      task12->add_dep(task1);

      std::vector<IndexRange> I3_index;
      std::shared_ptr<Tensor<T>> I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor13 = {I0, t2, I3};
      std::shared_ptr<Task13<T>> task13(new Task13<T>(tensor13, pindex));
      task9->add_dep(task13);
      task13->add_dep(task0);
      queue_->add_task(task13);


      std::vector<std::shared_ptr<Tensor<T>>> tensor14 = {I3, Gamma0};
      std::shared_ptr<Task14<T>> task14(new Task14<T>(tensor14, pindex));
      task13->add_dep(task14);
      task14->add_dep(task0);
      queue_->add_task(task14);

      task14->add_dep(task1);

      std::vector<IndexRange> I4_index = {this->closed_, this->virt_, this->virt_, this->closed_};
      std::shared_ptr<Tensor<T>> I4(new Tensor<T>(I4_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor15 = {r, I4};
      std::shared_ptr<Task15<T>> task15(new Task15<T>(tensor15, pindex));
      task15->add_dep(task0);
      queue_->add_task(task15);


      std::vector<IndexRange> I5_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I5(new Tensor<T>(I5_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor16 = {I4, this->f1_, I5};
      std::shared_ptr<Task16<T>> task16(new Task16<T>(tensor16, pindex));
      task15->add_dep(task16);
      task16->add_dep(task0);
      queue_->add_task(task16);


      std::vector<std::shared_ptr<Tensor<T>>> tensor17 = {I5, t2};
      std::shared_ptr<Task17<T>> task17(new Task17<T>(tensor17, pindex));
      task16->add_dep(task17);
      task17->add_dep(task0);
      queue_->add_task(task17);


      std::vector<IndexRange> I9_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I9(new Tensor<T>(I9_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor18 = {I4, this->f1_, I9};
      std::shared_ptr<Task18<T>> task18(new Task18<T>(tensor18, pindex));
      task15->add_dep(task18);
      task18->add_dep(task0);
      queue_->add_task(task18);


      std::vector<std::shared_ptr<Tensor<T>>> tensor19 = {I9, t2};
      std::shared_ptr<Task19<T>> task19(new Task19<T>(tensor19, pindex));
      task18->add_dep(task19);
      task19->add_dep(task0);
      queue_->add_task(task19);


      std::vector<IndexRange> I13_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I13(new Tensor<T>(I13_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor20 = {I4, this->f1_, I13};
      std::shared_ptr<Task20<T>> task20(new Task20<T>(tensor20, pindex));
      task15->add_dep(task20);
      task20->add_dep(task0);
      queue_->add_task(task20);


      std::vector<IndexRange> I14_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I14(new Tensor<T>(I14_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor21 = {I13, t2, I14};
      std::shared_ptr<Task21<T>> task21(new Task21<T>(tensor21, pindex));
      task20->add_dep(task21);
      task21->add_dep(task0);
      queue_->add_task(task21);


      std::vector<std::shared_ptr<Tensor<T>>> tensor22 = {I14, Gamma2};
      std::shared_ptr<Task22<T>> task22(new Task22<T>(tensor22, pindex));
      task21->add_dep(task22);
      task22->add_dep(task0);
      queue_->add_task(task22);

      task22->add_dep(task2);

      std::vector<IndexRange> I17_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I17(new Tensor<T>(I17_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor23 = {I13, t2, I17};
      std::shared_ptr<Task23<T>> task23(new Task23<T>(tensor23, pindex));
      task20->add_dep(task23);
      task23->add_dep(task0);
      queue_->add_task(task23);


      std::vector<std::shared_ptr<Tensor<T>>> tensor24 = {I17, Gamma2};
      std::shared_ptr<Task24<T>> task24(new Task24<T>(tensor24, pindex));
      task23->add_dep(task24);
      task24->add_dep(task0);
      queue_->add_task(task24);

      task24->add_dep(task2);

      std::vector<IndexRange> I18_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I18(new Tensor<T>(I18_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor25 = {r, I18};
      std::shared_ptr<Task25<T>> task25(new Task25<T>(tensor25, pindex));
      task25->add_dep(task0);
      queue_->add_task(task25);


      std::vector<IndexRange> I19_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I19(new Tensor<T>(I19_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor26 = {I18, this->f1_, I19};
      std::shared_ptr<Task26<T>> task26(new Task26<T>(tensor26, pindex));
      task25->add_dep(task26);
      task26->add_dep(task0);
      queue_->add_task(task26);


      std::vector<IndexRange> I20_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I20(new Tensor<T>(I20_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor27 = {I19, t2, I20};
      std::shared_ptr<Task27<T>> task27(new Task27<T>(tensor27, pindex));
      task26->add_dep(task27);
      task27->add_dep(task0);
      queue_->add_task(task27);


      std::vector<std::shared_ptr<Tensor<T>>> tensor28 = {I20, Gamma2};
      std::shared_ptr<Task28<T>> task28(new Task28<T>(tensor28, pindex));
      task27->add_dep(task28);
      task28->add_dep(task0);
      queue_->add_task(task28);

      task28->add_dep(task2);

      std::vector<IndexRange> I23_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I23(new Tensor<T>(I23_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor29 = {I19, t2, I23};
      std::shared_ptr<Task29<T>> task29(new Task29<T>(tensor29, pindex));
      task26->add_dep(task29);
      task29->add_dep(task0);
      queue_->add_task(task29);


      std::vector<std::shared_ptr<Tensor<T>>> tensor30 = {I23, Gamma2};
      std::shared_ptr<Task30<T>> task30(new Task30<T>(tensor30, pindex));
      task29->add_dep(task30);
      task30->add_dep(task0);
      queue_->add_task(task30);

      task30->add_dep(task2);

      std::vector<IndexRange> I25_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I25(new Tensor<T>(I25_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor31 = {I18, t2, I25};
      std::shared_ptr<Task31<T>> task31(new Task31<T>(tensor31, pindex));
      task25->add_dep(task31);
      task31->add_dep(task0);
      queue_->add_task(task31);


      std::vector<std::shared_ptr<Tensor<T>>> tensor32 = {I25, Gamma6};
      std::shared_ptr<Task32<T>> task32(new Task32<T>(tensor32, pindex));
      task31->add_dep(task32);
      task32->add_dep(task0);
      queue_->add_task(task32);

      task32->add_dep(task3);

      std::vector<IndexRange> I27_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I27(new Tensor<T>(I27_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor33 = {I18, t2, I27};
      std::shared_ptr<Task33<T>> task33(new Task33<T>(tensor33, pindex));
      task25->add_dep(task33);
      task33->add_dep(task0);
      queue_->add_task(task33);


      std::vector<std::shared_ptr<Tensor<T>>> tensor34 = {I27, Gamma6};
      std::shared_ptr<Task34<T>> task34(new Task34<T>(tensor34, pindex));
      task33->add_dep(task34);
      task34->add_dep(task0);
      queue_->add_task(task34);

      task34->add_dep(task3);

      std::vector<IndexRange> I29_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I29(new Tensor<T>(I29_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor35 = {I18, this->f1_, I29};
      std::shared_ptr<Task35<T>> task35(new Task35<T>(tensor35, pindex));
      task25->add_dep(task35);
      task35->add_dep(task0);
      queue_->add_task(task35);


      std::vector<IndexRange> I30_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I30(new Tensor<T>(I30_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor36 = {I29, t2, I30};
      std::shared_ptr<Task36<T>> task36(new Task36<T>(tensor36, pindex));
      task35->add_dep(task36);
      task36->add_dep(task0);
      queue_->add_task(task36);


      std::vector<std::shared_ptr<Tensor<T>>> tensor37 = {I30, Gamma2};
      std::shared_ptr<Task37<T>> task37(new Task37<T>(tensor37, pindex));
      task36->add_dep(task37);
      task37->add_dep(task0);
      queue_->add_task(task37);

      task37->add_dep(task2);

      std::vector<IndexRange> I33_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I33(new Tensor<T>(I33_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor38 = {I29, t2, I33};
      std::shared_ptr<Task38<T>> task38(new Task38<T>(tensor38, pindex));
      task35->add_dep(task38);
      task38->add_dep(task0);
      queue_->add_task(task38);


      std::vector<std::shared_ptr<Tensor<T>>> tensor39 = {I33, Gamma2};
      std::shared_ptr<Task39<T>> task39(new Task39<T>(tensor39, pindex));
      task38->add_dep(task39);
      task39->add_dep(task0);
      queue_->add_task(task39);

      task39->add_dep(task2);

      std::vector<IndexRange> I35_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I35(new Tensor<T>(I35_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor40 = {I18, this->f1_, I35};
      std::shared_ptr<Task40<T>> task40(new Task40<T>(tensor40, pindex));
      task25->add_dep(task40);
      task40->add_dep(task0);
      queue_->add_task(task40);


      std::vector<IndexRange> I36_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I36(new Tensor<T>(I36_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor41 = {I35, t2, I36};
      std::shared_ptr<Task41<T>> task41(new Task41<T>(tensor41, pindex));
      task40->add_dep(task41);
      task41->add_dep(task0);
      queue_->add_task(task41);


      std::vector<std::shared_ptr<Tensor<T>>> tensor42 = {I36, Gamma2};
      std::shared_ptr<Task42<T>> task42(new Task42<T>(tensor42, pindex));
      task41->add_dep(task42);
      task42->add_dep(task0);
      queue_->add_task(task42);

      task42->add_dep(task2);

      std::vector<IndexRange> I39_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I39(new Tensor<T>(I39_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor43 = {I35, t2, I39};
      std::shared_ptr<Task43<T>> task43(new Task43<T>(tensor43, pindex));
      task40->add_dep(task43);
      task43->add_dep(task0);
      queue_->add_task(task43);


      std::vector<std::shared_ptr<Tensor<T>>> tensor44 = {I39, Gamma2};
      std::shared_ptr<Task44<T>> task44(new Task44<T>(tensor44, pindex));
      task43->add_dep(task44);
      task44->add_dep(task0);
      queue_->add_task(task44);

      task44->add_dep(task2);

      std::vector<IndexRange> I41_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I41(new Tensor<T>(I41_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor45 = {I18, this->f1_, I41};
      std::shared_ptr<Task45<T>> task45(new Task45<T>(tensor45, pindex));
      task25->add_dep(task45);
      task45->add_dep(task0);
      queue_->add_task(task45);


      std::vector<IndexRange> I42_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I42(new Tensor<T>(I42_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor46 = {I41, t2, I42};
      std::shared_ptr<Task46<T>> task46(new Task46<T>(tensor46, pindex));
      task45->add_dep(task46);
      task46->add_dep(task0);
      queue_->add_task(task46);


      std::vector<std::shared_ptr<Tensor<T>>> tensor47 = {I42, Gamma2};
      std::shared_ptr<Task47<T>> task47(new Task47<T>(tensor47, pindex));
      task46->add_dep(task47);
      task47->add_dep(task0);
      queue_->add_task(task47);

      task47->add_dep(task2);

      std::vector<IndexRange> I45_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I45(new Tensor<T>(I45_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor48 = {I41, t2, I45};
      std::shared_ptr<Task48<T>> task48(new Task48<T>(tensor48, pindex));
      task45->add_dep(task48);
      task48->add_dep(task0);
      queue_->add_task(task48);


      std::vector<std::shared_ptr<Tensor<T>>> tensor49 = {I45, Gamma2};
      std::shared_ptr<Task49<T>> task49(new Task49<T>(tensor49, pindex));
      task48->add_dep(task49);
      task49->add_dep(task0);
      queue_->add_task(task49);

      task49->add_dep(task2);

      std::vector<IndexRange> I47_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I47(new Tensor<T>(I47_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor50 = {I18, this->f1_, I47};
      std::shared_ptr<Task50<T>> task50(new Task50<T>(tensor50, pindex));
      task25->add_dep(task50);
      task50->add_dep(task0);
      queue_->add_task(task50);


      std::vector<IndexRange> I48_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I48(new Tensor<T>(I48_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor51 = {I47, t2, I48};
      std::shared_ptr<Task51<T>> task51(new Task51<T>(tensor51, pindex));
      task50->add_dep(task51);
      task51->add_dep(task0);
      queue_->add_task(task51);


      std::vector<std::shared_ptr<Tensor<T>>> tensor52 = {I48, Gamma14};
      std::shared_ptr<Task52<T>> task52(new Task52<T>(tensor52, pindex));
      task51->add_dep(task52);
      task52->add_dep(task0);
      queue_->add_task(task52);

      task52->add_dep(task4);

      std::vector<IndexRange> I60_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I60(new Tensor<T>(I60_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor53 = {I18, t2, I60};
      std::shared_ptr<Task53<T>> task53(new Task53<T>(tensor53, pindex));
      task25->add_dep(task53);
      task53->add_dep(task0);
      queue_->add_task(task53);


      std::vector<std::shared_ptr<Tensor<T>>> tensor54 = {I60, Gamma2};
      std::shared_ptr<Task54<T>> task54(new Task54<T>(tensor54, pindex, this->e0_));
      task53->add_dep(task54);
      task54->add_dep(task0);
      queue_->add_task(task54);

      task54->add_dep(task2);

      std::vector<IndexRange> I62_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I62(new Tensor<T>(I62_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor55 = {I18, t2, I62};
      std::shared_ptr<Task55<T>> task55(new Task55<T>(tensor55, pindex));
      task25->add_dep(task55);
      task55->add_dep(task0);
      queue_->add_task(task55);


      std::vector<std::shared_ptr<Tensor<T>>> tensor56 = {I62, Gamma2};
      std::shared_ptr<Task56<T>> task56(new Task56<T>(tensor56, pindex, this->e0_));
      task55->add_dep(task56);
      task56->add_dep(task0);
      queue_->add_task(task56);

      task56->add_dep(task2);

      std::vector<IndexRange> I68_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I68(new Tensor<T>(I68_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor57 = {I18, this->v2_, I68};
      std::shared_ptr<Task57<T>> task57(new Task57<T>(tensor57, pindex));
      task25->add_dep(task57);
      task57->add_dep(task0);
      queue_->add_task(task57);


      std::vector<std::shared_ptr<Tensor<T>>> tensor58 = {I68, Gamma2};
      std::shared_ptr<Task58<T>> task58(new Task58<T>(tensor58, pindex));
      task57->add_dep(task58);
      task58->add_dep(task0);
      queue_->add_task(task58);

      task58->add_dep(task2);

      std::vector<IndexRange> I70_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I70(new Tensor<T>(I70_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor59 = {I18, this->v2_, I70};
      std::shared_ptr<Task59<T>> task59(new Task59<T>(tensor59, pindex));
      task25->add_dep(task59);
      task59->add_dep(task0);
      queue_->add_task(task59);


      std::vector<std::shared_ptr<Tensor<T>>> tensor60 = {I70, Gamma2};
      std::shared_ptr<Task60<T>> task60(new Task60<T>(tensor60, pindex));
      task59->add_dep(task60);
      task60->add_dep(task0);
      queue_->add_task(task60);

      task60->add_dep(task2);

      std::vector<IndexRange> I49_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I49(new Tensor<T>(I49_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor61 = {r, I49};
      std::shared_ptr<Task61<T>> task61(new Task61<T>(tensor61, pindex));
      task61->add_dep(task0);
      queue_->add_task(task61);


      std::vector<IndexRange> I50_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I50(new Tensor<T>(I50_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor62 = {I49, this->f1_, I50};
      std::shared_ptr<Task62<T>> task62(new Task62<T>(tensor62, pindex));
      task61->add_dep(task62);
      task62->add_dep(task0);
      queue_->add_task(task62);


      std::vector<IndexRange> I51_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I51(new Tensor<T>(I51_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor63 = {I50, t2, I51};
      std::shared_ptr<Task63<T>> task63(new Task63<T>(tensor63, pindex));
      task62->add_dep(task63);
      task63->add_dep(task0);
      queue_->add_task(task63);


      std::vector<std::shared_ptr<Tensor<T>>> tensor64 = {I51, Gamma14};
      std::shared_ptr<Task64<T>> task64(new Task64<T>(tensor64, pindex));
      task63->add_dep(task64);
      task64->add_dep(task0);
      queue_->add_task(task64);

      task64->add_dep(task4);

      std::vector<IndexRange> I55_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I55(new Tensor<T>(I55_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor65 = {I49, this->f1_, I55};
      std::shared_ptr<Task65<T>> task65(new Task65<T>(tensor65, pindex));
      task61->add_dep(task65);
      task65->add_dep(task0);
      queue_->add_task(task65);


      std::vector<IndexRange> I56_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I56(new Tensor<T>(I56_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor66 = {I55, t2, I56};
      std::shared_ptr<Task66<T>> task66(new Task66<T>(tensor66, pindex));
      task65->add_dep(task66);
      task66->add_dep(task0);
      queue_->add_task(task66);


      std::vector<std::shared_ptr<Tensor<T>>> tensor67 = {I56, Gamma14};
      std::shared_ptr<Task67<T>> task67(new Task67<T>(tensor67, pindex));
      task66->add_dep(task67);
      task67->add_dep(task0);
      queue_->add_task(task67);

      task67->add_dep(task4);

      std::vector<IndexRange> I52_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I52(new Tensor<T>(I52_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor68 = {r, I52};
      std::shared_ptr<Task68<T>> task68(new Task68<T>(tensor68, pindex));
      task68->add_dep(task0);
      queue_->add_task(task68);


      std::vector<IndexRange> I53_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I53(new Tensor<T>(I53_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor69 = {I52, t2, I53};
      std::shared_ptr<Task69<T>> task69(new Task69<T>(tensor69, pindex));
      task68->add_dep(task69);
      task69->add_dep(task0);
      queue_->add_task(task69);


      std::vector<std::shared_ptr<Tensor<T>>> tensor70 = {I53, Gamma16};
      std::shared_ptr<Task70<T>> task70(new Task70<T>(tensor70, pindex));
      task69->add_dep(task70);
      task70->add_dep(task0);
      queue_->add_task(task70);

      task70->add_dep(task5);

      std::vector<IndexRange> I64_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I64(new Tensor<T>(I64_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor71 = {I52, t2, I64};
      std::shared_ptr<Task71<T>> task71(new Task71<T>(tensor71, pindex));
      task68->add_dep(task71);
      task71->add_dep(task0);
      queue_->add_task(task71);


      std::vector<std::shared_ptr<Tensor<T>>> tensor72 = {I64, Gamma14};
      std::shared_ptr<Task72<T>> task72(new Task72<T>(tensor72, pindex, this->e0_));
      task71->add_dep(task72);
      task72->add_dep(task0);
      queue_->add_task(task72);

      task72->add_dep(task4);

      std::vector<IndexRange> I72_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I72(new Tensor<T>(I72_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor73 = {I52, this->v2_, I72};
      std::shared_ptr<Task73<T>> task73(new Task73<T>(tensor73, pindex));
      task68->add_dep(task73);
      task73->add_dep(task0);
      queue_->add_task(task73);


      std::vector<std::shared_ptr<Tensor<T>>> tensor74 = {I72, Gamma14};
      std::shared_ptr<Task74<T>> task74(new Task74<T>(tensor74, pindex));
      task73->add_dep(task74);
      task74->add_dep(task0);
      queue_->add_task(task74);

      task74->add_dep(task4);

      std::shared_ptr<Queue<T>> energy_(new Queue<T>());
      std::vector<IndexRange> I73_index;
      std::shared_ptr<Tensor<T>> I73(new Tensor<T>(I73_index, false));
      std::vector<IndexRange> I74_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I74(new Tensor<T>(I74_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor75 = {I73, t2, I74};
      std::shared_ptr<Task75<T>> task75(new Task75<T>(tensor75, pindex));
      energy_->add_task(task75);


      std::vector<std::shared_ptr<Tensor<T>>> tensor76 = {I74, this->v2_};
      std::shared_ptr<Task76<T>> task76(new Task76<T>(tensor76, pindex));
      task75->add_dep(task76);
      energy_->add_task(task76);


      std::vector<IndexRange> I78_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I78(new Tensor<T>(I78_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor77 = {I73, t2, I78};
      std::shared_ptr<Task77<T>> task77(new Task77<T>(tensor77, pindex));
      task75->add_dep(task77);
      energy_->add_task(task77);


      std::vector<IndexRange> I79_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I79(new Tensor<T>(I79_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor78 = {I78, this->v2_, I79};
      std::shared_ptr<Task78<T>> task78(new Task78<T>(tensor78, pindex));
      task77->add_dep(task78);
      energy_->add_task(task78);


      std::vector<std::shared_ptr<Tensor<T>>> tensor79 = {I79, Gamma2};
      std::shared_ptr<Task79<T>> task79(new Task79<T>(tensor79, pindex));
      task78->add_dep(task79);
      energy_->add_task(task79);

      task79->add_dep(task2);

      std::vector<IndexRange> I82_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I82(new Tensor<T>(I82_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor80 = {I78, this->v2_, I82};
      std::shared_ptr<Task80<T>> task80(new Task80<T>(tensor80, pindex));
      task77->add_dep(task80);
      energy_->add_task(task80);


      std::vector<std::shared_ptr<Tensor<T>>> tensor81 = {I82, Gamma2};
      std::shared_ptr<Task81<T>> task81(new Task81<T>(tensor81, pindex));
      task80->add_dep(task81);
      energy_->add_task(task81);

      task81->add_dep(task2);

      std::vector<IndexRange> I84_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I84(new Tensor<T>(I84_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor82 = {I73, t2, I84};
      std::shared_ptr<Task82<T>> task82(new Task82<T>(tensor82, pindex));
      task75->add_dep(task82);
      energy_->add_task(task82);


      std::vector<IndexRange> I85_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I85(new Tensor<T>(I85_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor83 = {I84, this->v2_, I85};
      std::shared_ptr<Task83<T>> task83(new Task83<T>(tensor83, pindex));
      task82->add_dep(task83);
      energy_->add_task(task83);


      std::vector<std::shared_ptr<Tensor<T>>> tensor84 = {I85, Gamma14};
      std::shared_ptr<Task84<T>> task84(new Task84<T>(tensor84, pindex));
      task83->add_dep(task84);
      energy_->add_task(task84);

      task84->add_dep(task4);

      std::shared_ptr<Queue<T>> dedci_(new Queue<T>());
      std::vector<IndexRange> I86_index;
      std::shared_ptr<Tensor<T>> I86(new Tensor<T>(I86_index, false));
      std::vector<IndexRange> I87_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I87(new Tensor<T>(I87_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor85 = {I86, t2, I87};
      std::shared_ptr<Task85<T>> task85(new Task85<T>(tensor85, pindex));
      dedci_->add_task(task85);


      std::vector<std::shared_ptr<Tensor<T>>> tensor86 = {I87, this->v2_};
      std::shared_ptr<Task86<T>> task86(new Task86<T>(tensor86, pindex));
      task85->add_dep(task86);
      dedci_->add_task(task86);


      std::vector<IndexRange> I91_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I91(new Tensor<T>(I91_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor87 = {I86, t2, I91};
      std::shared_ptr<Task87<T>> task87(new Task87<T>(tensor87, pindex));
      task85->add_dep(task87);
      dedci_->add_task(task87);


      std::vector<IndexRange> I92_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I92(new Tensor<T>(I92_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor88 = {I91, this->v2_, I92};
      std::shared_ptr<Task88<T>> task88(new Task88<T>(tensor88, pindex));
      task87->add_dep(task88);
      dedci_->add_task(task88);


      std::vector<std::shared_ptr<Tensor<T>>> tensor89 = {I92, Gamma27};
      std::shared_ptr<Task89<T>> task89(new Task89<T>(tensor89, pindex));
      task88->add_dep(task89);
      dedci_->add_task(task89);

      task89->add_dep(task6);

      std::vector<IndexRange> I95_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I95(new Tensor<T>(I95_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor90 = {I91, this->v2_, I95};
      std::shared_ptr<Task90<T>> task90(new Task90<T>(tensor90, pindex));
      task87->add_dep(task90);
      dedci_->add_task(task90);


      std::vector<std::shared_ptr<Tensor<T>>> tensor91 = {I95, Gamma27};
      std::shared_ptr<Task91<T>> task91(new Task91<T>(tensor91, pindex));
      task90->add_dep(task91);
      dedci_->add_task(task91);

      task91->add_dep(task6);

      std::vector<IndexRange> I97_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I97(new Tensor<T>(I97_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor92 = {I86, t2, I97};
      std::shared_ptr<Task92<T>> task92(new Task92<T>(tensor92, pindex));
      task85->add_dep(task92);
      dedci_->add_task(task92);


      std::vector<IndexRange> I98_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I98(new Tensor<T>(I98_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor93 = {I97, this->v2_, I98};
      std::shared_ptr<Task93<T>> task93(new Task93<T>(tensor93, pindex));
      task92->add_dep(task93);
      dedci_->add_task(task93);


      std::vector<std::shared_ptr<Tensor<T>>> tensor94 = {I98, Gamma29};
      std::shared_ptr<Task94<T>> task94(new Task94<T>(tensor94, pindex));
      task93->add_dep(task94);
      dedci_->add_task(task94);

      task94->add_dep(task7);

      std::vector<IndexRange> I104_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I104(new Tensor<T>(I104_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor95 = {I86, t2, I104};
      std::shared_ptr<Task95<T>> task95(new Task95<T>(tensor95, pindex));
      task85->add_dep(task95);
      dedci_->add_task(task95);


      std::vector<IndexRange> I105_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I105(new Tensor<T>(I105_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor96 = {I104, this->v2_, I105};
      std::shared_ptr<Task96<T>> task96(new Task96<T>(tensor96, pindex));
      task95->add_dep(task96);
      dedci_->add_task(task96);


      std::vector<std::shared_ptr<Tensor<T>>> tensor97 = {I105, Gamma27};
      std::shared_ptr<Task97<T>> task97(new Task97<T>(tensor97, pindex));
      task96->add_dep(task97);
      dedci_->add_task(task97);

      task97->add_dep(task6);

      std::vector<IndexRange> I108_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I108(new Tensor<T>(I108_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor98 = {I104, this->v2_, I108};
      std::shared_ptr<Task98<T>> task98(new Task98<T>(tensor98, pindex));
      task95->add_dep(task98);
      dedci_->add_task(task98);


      std::vector<std::shared_ptr<Tensor<T>>> tensor99 = {I108, Gamma27};
      std::shared_ptr<Task99<T>> task99(new Task99<T>(tensor99, pindex));
      task98->add_dep(task99);
      dedci_->add_task(task99);

      task99->add_dep(task6);

      std::vector<IndexRange> I110_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I110(new Tensor<T>(I110_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor100 = {I86, t2, I110};
      std::shared_ptr<Task100<T>> task100(new Task100<T>(tensor100, pindex));
      task85->add_dep(task100);
      dedci_->add_task(task100);


      std::vector<IndexRange> I111_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I111(new Tensor<T>(I111_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor101 = {I110, this->v2_, I111};
      std::shared_ptr<Task101<T>> task101(new Task101<T>(tensor101, pindex));
      task100->add_dep(task101);
      dedci_->add_task(task101);


      std::vector<std::shared_ptr<Tensor<T>>> tensor102 = {I111, Gamma29};
      std::shared_ptr<Task102<T>> task102(new Task102<T>(tensor102, pindex));
      task101->add_dep(task102);
      dedci_->add_task(task102);

      task102->add_dep(task7);

      std::shared_ptr<Queue<T>> correction_(new Queue<T>());
      std::vector<IndexRange> I112_index;
      std::shared_ptr<Tensor<T>> I112(new Tensor<T>(I112_index, false));
      std::vector<IndexRange> I113_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I113(new Tensor<T>(I113_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor103 = {I112, t2, I113};
      std::shared_ptr<Task103<T>> task103(new Task103<T>(tensor103, pindex));
      correction_->add_task(task103);


      std::vector<std::shared_ptr<Tensor<T>>> tensor104 = {I113, t2};
      std::shared_ptr<Task104<T>> task104(new Task104<T>(tensor104, pindex));
      task103->add_dep(task104);
      correction_->add_task(task104);


      std::vector<IndexRange> I117_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I117(new Tensor<T>(I117_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor105 = {I112, t2, I117};
      std::shared_ptr<Task105<T>> task105(new Task105<T>(tensor105, pindex));
      task103->add_dep(task105);
      correction_->add_task(task105);


      std::vector<IndexRange> I118_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I118(new Tensor<T>(I118_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor106 = {I117, t2, I118};
      std::shared_ptr<Task106<T>> task106(new Task106<T>(tensor106, pindex));
      task105->add_dep(task106);
      correction_->add_task(task106);


      std::vector<std::shared_ptr<Tensor<T>>> tensor107 = {I118, Gamma2};
      std::shared_ptr<Task107<T>> task107(new Task107<T>(tensor107, pindex));
      task106->add_dep(task107);
      correction_->add_task(task107);

      task107->add_dep(task2);

      std::vector<IndexRange> I121_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I121(new Tensor<T>(I121_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor108 = {I117, t2, I121};
      std::shared_ptr<Task108<T>> task108(new Task108<T>(tensor108, pindex));
      task105->add_dep(task108);
      correction_->add_task(task108);


      std::vector<std::shared_ptr<Tensor<T>>> tensor109 = {I121, Gamma2};
      std::shared_ptr<Task109<T>> task109(new Task109<T>(tensor109, pindex));
      task108->add_dep(task109);
      correction_->add_task(task109);

      task109->add_dep(task2);

      std::vector<IndexRange> I123_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I123(new Tensor<T>(I123_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor110 = {I112, t2, I123};
      std::shared_ptr<Task110<T>> task110(new Task110<T>(tensor110, pindex));
      task103->add_dep(task110);
      correction_->add_task(task110);


      std::vector<IndexRange> I124_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I124(new Tensor<T>(I124_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor111 = {I123, t2, I124};
      std::shared_ptr<Task111<T>> task111(new Task111<T>(tensor111, pindex));
      task110->add_dep(task111);
      correction_->add_task(task111);


      std::vector<std::shared_ptr<Tensor<T>>> tensor112 = {I124, Gamma14};
      std::shared_ptr<Task112<T>> task112(new Task112<T>(tensor112, pindex));
      task111->add_dep(task112);
      correction_->add_task(task112);

      task112->add_dep(task4);

      std::shared_ptr<Queue<T>> density_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor113 = {this->den1_};
      std::shared_ptr<Task113<T>> task113(new Task113<T>(tensor113));
      density_->add_task(task113);

      std::vector<IndexRange> I125_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I125(new Tensor<T>(I125_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor114 = {this->den1_, I125};
      std::shared_ptr<Task114<T>> task114(new Task114<T>(tensor114, pindex));
      task114->add_dep(task113);
      density_->add_task(task114);


      std::vector<IndexRange> I126_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I126(new Tensor<T>(I126_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor115 = {I125, t2, I126};
      std::shared_ptr<Task115<T>> task115(new Task115<T>(tensor115, pindex));
      task114->add_dep(task115);
      task115->add_dep(task113);
      density_->add_task(task115);


      std::vector<IndexRange> I127_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I127(new Tensor<T>(I127_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor116 = {I126, t2, I127};
      std::shared_ptr<Task116<T>> task116(new Task116<T>(tensor116, pindex));
      task115->add_dep(task116);
      task116->add_dep(task113);
      density_->add_task(task116);


      std::vector<std::shared_ptr<Tensor<T>>> tensor117 = {I127, Gamma2};
      std::shared_ptr<Task117<T>> task117(new Task117<T>(tensor117, pindex));
      task116->add_dep(task117);
      task117->add_dep(task113);
      density_->add_task(task117);

      task117->add_dep(task2);

      std::vector<IndexRange> I130_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I130(new Tensor<T>(I130_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor118 = {I126, t2, I130};
      std::shared_ptr<Task118<T>> task118(new Task118<T>(tensor118, pindex));
      task115->add_dep(task118);
      task118->add_dep(task113);
      density_->add_task(task118);


      std::vector<std::shared_ptr<Tensor<T>>> tensor119 = {I130, Gamma2};
      std::shared_ptr<Task119<T>> task119(new Task119<T>(tensor119, pindex));
      task118->add_dep(task119);
      task119->add_dep(task113);
      density_->add_task(task119);

      task119->add_dep(task2);

      std::vector<IndexRange> I131_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor<T>> I131(new Tensor<T>(I131_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor120 = {this->den1_, I131};
      std::shared_ptr<Task120<T>> task120(new Task120<T>(tensor120, pindex));
      task120->add_dep(task113);
      density_->add_task(task120);


      std::vector<IndexRange> I132_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I132(new Tensor<T>(I132_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor121 = {I131, t2, I132};
      std::shared_ptr<Task121<T>> task121(new Task121<T>(tensor121, pindex));
      task120->add_dep(task121);
      task121->add_dep(task113);
      density_->add_task(task121);


      std::vector<std::shared_ptr<Tensor<T>>> tensor122 = {I132, t2};
      std::shared_ptr<Task122<T>> task122(new Task122<T>(tensor122, pindex));
      task121->add_dep(task122);
      task122->add_dep(task113);
      density_->add_task(task122);


      std::vector<IndexRange> I135_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I135(new Tensor<T>(I135_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor123 = {this->den1_, I135};
      std::shared_ptr<Task123<T>> task123(new Task123<T>(tensor123, pindex));
      task123->add_dep(task113);
      density_->add_task(task123);


      std::vector<IndexRange> I136_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I136(new Tensor<T>(I136_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor124 = {I135, t2, I136};
      std::shared_ptr<Task124<T>> task124(new Task124<T>(tensor124, pindex));
      task123->add_dep(task124);
      task124->add_dep(task113);
      density_->add_task(task124);


      std::vector<std::shared_ptr<Tensor<T>>> tensor125 = {I136, t2};
      std::shared_ptr<Task125<T>> task125(new Task125<T>(tensor125, pindex));
      task124->add_dep(task125);
      task125->add_dep(task113);
      density_->add_task(task125);


      std::vector<IndexRange> I139_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor<T>> I139(new Tensor<T>(I139_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor126 = {this->den1_, I139};
      std::shared_ptr<Task126<T>> task126(new Task126<T>(tensor126, pindex));
      task126->add_dep(task113);
      density_->add_task(task126);


      std::vector<IndexRange> I140_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I140(new Tensor<T>(I140_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor127 = {I139, t2, I140};
      std::shared_ptr<Task127<T>> task127(new Task127<T>(tensor127, pindex));
      task126->add_dep(task127);
      task127->add_dep(task113);
      density_->add_task(task127);


      std::vector<IndexRange> I141_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I141(new Tensor<T>(I141_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor128 = {I140, t2, I141};
      std::shared_ptr<Task128<T>> task128(new Task128<T>(tensor128, pindex));
      task127->add_dep(task128);
      task128->add_dep(task113);
      density_->add_task(task128);


      std::vector<std::shared_ptr<Tensor<T>>> tensor129 = {I141, Gamma2};
      std::shared_ptr<Task129<T>> task129(new Task129<T>(tensor129, pindex));
      task128->add_dep(task129);
      task129->add_dep(task113);
      density_->add_task(task129);

      task129->add_dep(task2);

      std::vector<IndexRange> I144_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I144(new Tensor<T>(I144_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor130 = {I140, t2, I144};
      std::shared_ptr<Task130<T>> task130(new Task130<T>(tensor130, pindex));
      task127->add_dep(task130);
      task130->add_dep(task113);
      density_->add_task(task130);


      std::vector<std::shared_ptr<Tensor<T>>> tensor131 = {I144, Gamma2};
      std::shared_ptr<Task131<T>> task131(new Task131<T>(tensor131, pindex));
      task130->add_dep(task131);
      task131->add_dep(task113);
      density_->add_task(task131);

      task131->add_dep(task2);

      std::vector<IndexRange> I145_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor<T>> I145(new Tensor<T>(I145_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor132 = {this->den1_, I145};
      std::shared_ptr<Task132<T>> task132(new Task132<T>(tensor132, pindex));
      task132->add_dep(task113);
      density_->add_task(task132);


      std::vector<IndexRange> I146_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I146(new Tensor<T>(I146_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor133 = {I145, t2, I146};
      std::shared_ptr<Task133<T>> task133(new Task133<T>(tensor133, pindex));
      task132->add_dep(task133);
      task133->add_dep(task113);
      density_->add_task(task133);


      std::vector<IndexRange> I147_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I147(new Tensor<T>(I147_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor134 = {I146, t2, I147};
      std::shared_ptr<Task134<T>> task134(new Task134<T>(tensor134, pindex));
      task133->add_dep(task134);
      task134->add_dep(task113);
      density_->add_task(task134);


      std::vector<std::shared_ptr<Tensor<T>>> tensor135 = {I147, Gamma2};
      std::shared_ptr<Task135<T>> task135(new Task135<T>(tensor135, pindex));
      task134->add_dep(task135);
      task135->add_dep(task113);
      density_->add_task(task135);

      task135->add_dep(task2);

      std::vector<IndexRange> I150_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I150(new Tensor<T>(I150_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor136 = {I146, t2, I150};
      std::shared_ptr<Task136<T>> task136(new Task136<T>(tensor136, pindex));
      task133->add_dep(task136);
      task136->add_dep(task113);
      density_->add_task(task136);


      std::vector<std::shared_ptr<Tensor<T>>> tensor137 = {I150, Gamma2};
      std::shared_ptr<Task137<T>> task137(new Task137<T>(tensor137, pindex));
      task136->add_dep(task137);
      task137->add_dep(task113);
      density_->add_task(task137);

      task137->add_dep(task2);

      std::vector<IndexRange> I151_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I151(new Tensor<T>(I151_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor138 = {this->den1_, I151};
      std::shared_ptr<Task138<T>> task138(new Task138<T>(tensor138, pindex));
      task138->add_dep(task113);
      density_->add_task(task138);


      std::vector<IndexRange> I152_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I152(new Tensor<T>(I152_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor139 = {I151, t2, I152};
      std::shared_ptr<Task139<T>> task139(new Task139<T>(tensor139, pindex));
      task138->add_dep(task139);
      task139->add_dep(task113);
      density_->add_task(task139);


      std::vector<IndexRange> I153_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I153(new Tensor<T>(I153_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor140 = {I152, t2, I153};
      std::shared_ptr<Task140<T>> task140(new Task140<T>(tensor140, pindex));
      task139->add_dep(task140);
      task140->add_dep(task113);
      density_->add_task(task140);


      std::vector<std::shared_ptr<Tensor<T>>> tensor141 = {I153, Gamma14};
      std::shared_ptr<Task141<T>> task141(new Task141<T>(tensor141, pindex));
      task140->add_dep(task141);
      task141->add_dep(task113);
      density_->add_task(task141);

      task141->add_dep(task4);

      std::vector<IndexRange> I156_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I156(new Tensor<T>(I156_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor142 = {I152, t2, I156};
      std::shared_ptr<Task142<T>> task142(new Task142<T>(tensor142, pindex));
      task139->add_dep(task142);
      task142->add_dep(task113);
      density_->add_task(task142);


      std::vector<std::shared_ptr<Tensor<T>>> tensor143 = {I156, Gamma14};
      std::shared_ptr<Task143<T>> task143(new Task143<T>(tensor143, pindex));
      task142->add_dep(task143);
      task143->add_dep(task113);
      density_->add_task(task143);

      task143->add_dep(task4);

      std::vector<IndexRange> I157_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor<T>> I157(new Tensor<T>(I157_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor144 = {this->den1_, I157};
      std::shared_ptr<Task144<T>> task144(new Task144<T>(tensor144, pindex));
      task144->add_dep(task113);
      density_->add_task(task144);


      std::vector<IndexRange> I158_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I158(new Tensor<T>(I158_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor145 = {I157, t2, I158};
      std::shared_ptr<Task145<T>> task145(new Task145<T>(tensor145, pindex));
      task144->add_dep(task145);
      task145->add_dep(task113);
      density_->add_task(task145);


      std::vector<IndexRange> I159_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I159(new Tensor<T>(I159_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor146 = {I158, t2, I159};
      std::shared_ptr<Task146<T>> task146(new Task146<T>(tensor146, pindex));
      task145->add_dep(task146);
      task146->add_dep(task113);
      density_->add_task(task146);


      std::vector<std::shared_ptr<Tensor<T>>> tensor147 = {I159, Gamma2};
      std::shared_ptr<Task147<T>> task147(new Task147<T>(tensor147, pindex));
      task146->add_dep(task147);
      task147->add_dep(task113);
      density_->add_task(task147);

      task147->add_dep(task2);

      std::vector<IndexRange> I162_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I162(new Tensor<T>(I162_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor148 = {I158, t2, I162};
      std::shared_ptr<Task148<T>> task148(new Task148<T>(tensor148, pindex));
      task145->add_dep(task148);
      task148->add_dep(task113);
      density_->add_task(task148);


      std::vector<std::shared_ptr<Tensor<T>>> tensor149 = {I162, Gamma2};
      std::shared_ptr<Task149<T>> task149(new Task149<T>(tensor149, pindex));
      task148->add_dep(task149);
      task149->add_dep(task113);
      density_->add_task(task149);

      task149->add_dep(task2);

      std::vector<IndexRange> I163_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I163(new Tensor<T>(I163_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor150 = {this->den1_, I163};
      std::shared_ptr<Task150<T>> task150(new Task150<T>(tensor150, pindex));
      task150->add_dep(task113);
      density_->add_task(task150);


      std::vector<IndexRange> I164_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I164(new Tensor<T>(I164_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor151 = {I163, t2, I164};
      std::shared_ptr<Task151<T>> task151(new Task151<T>(tensor151, pindex));
      task150->add_dep(task151);
      task151->add_dep(task113);
      density_->add_task(task151);


      std::vector<IndexRange> I165_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I165(new Tensor<T>(I165_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor152 = {I164, t2, I165};
      std::shared_ptr<Task152<T>> task152(new Task152<T>(tensor152, pindex));
      task151->add_dep(task152);
      task152->add_dep(task113);
      density_->add_task(task152);


      std::vector<std::shared_ptr<Tensor<T>>> tensor153 = {I165, Gamma2};
      std::shared_ptr<Task153<T>> task153(new Task153<T>(tensor153, pindex));
      task152->add_dep(task153);
      task153->add_dep(task113);
      density_->add_task(task153);

      task153->add_dep(task2);

      std::vector<IndexRange> I168_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I168(new Tensor<T>(I168_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor154 = {I164, t2, I168};
      std::shared_ptr<Task154<T>> task154(new Task154<T>(tensor154, pindex));
      task151->add_dep(task154);
      task154->add_dep(task113);
      density_->add_task(task154);


      std::vector<std::shared_ptr<Tensor<T>>> tensor155 = {I168, Gamma2};
      std::shared_ptr<Task155<T>> task155(new Task155<T>(tensor155, pindex));
      task154->add_dep(task155);
      task155->add_dep(task113);
      density_->add_task(task155);

      task155->add_dep(task2);

      std::vector<IndexRange> I169_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I169(new Tensor<T>(I169_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor156 = {this->den1_, I169};
      std::shared_ptr<Task156<T>> task156(new Task156<T>(tensor156, pindex));
      task156->add_dep(task113);
      density_->add_task(task156);


      std::vector<IndexRange> I170_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I170(new Tensor<T>(I170_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor157 = {I169, t2, I170};
      std::shared_ptr<Task157<T>> task157(new Task157<T>(tensor157, pindex));
      task156->add_dep(task157);
      task157->add_dep(task113);
      density_->add_task(task157);


      std::vector<IndexRange> I171_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I171(new Tensor<T>(I171_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor158 = {I170, t2, I171};
      std::shared_ptr<Task158<T>> task158(new Task158<T>(tensor158, pindex));
      task157->add_dep(task158);
      task158->add_dep(task113);
      density_->add_task(task158);


      std::vector<std::shared_ptr<Tensor<T>>> tensor159 = {I171, Gamma2};
      std::shared_ptr<Task159<T>> task159(new Task159<T>(tensor159, pindex));
      task158->add_dep(task159);
      task159->add_dep(task113);
      density_->add_task(task159);

      task159->add_dep(task2);

      std::vector<IndexRange> I174_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I174(new Tensor<T>(I174_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor160 = {I170, t2, I174};
      std::shared_ptr<Task160<T>> task160(new Task160<T>(tensor160, pindex));
      task157->add_dep(task160);
      task160->add_dep(task113);
      density_->add_task(task160);


      std::vector<std::shared_ptr<Tensor<T>>> tensor161 = {I174, Gamma2};
      std::shared_ptr<Task161<T>> task161(new Task161<T>(tensor161, pindex));
      task160->add_dep(task161);
      task161->add_dep(task113);
      density_->add_task(task161);

      task161->add_dep(task2);

      std::vector<IndexRange> I175_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor<T>> I175(new Tensor<T>(I175_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor162 = {this->den1_, I175};
      std::shared_ptr<Task162<T>> task162(new Task162<T>(tensor162, pindex));
      task162->add_dep(task113);
      density_->add_task(task162);


      std::vector<IndexRange> I176_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I176(new Tensor<T>(I176_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor163 = {I175, t2, I176};
      std::shared_ptr<Task163<T>> task163(new Task163<T>(tensor163, pindex));
      task162->add_dep(task163);
      task163->add_dep(task113);
      density_->add_task(task163);


      std::vector<IndexRange> I177_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I177(new Tensor<T>(I177_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor164 = {I176, t2, I177};
      std::shared_ptr<Task164<T>> task164(new Task164<T>(tensor164, pindex));
      task163->add_dep(task164);
      task164->add_dep(task113);
      density_->add_task(task164);


      std::vector<std::shared_ptr<Tensor<T>>> tensor165 = {I177, Gamma14};
      std::shared_ptr<Task165<T>> task165(new Task165<T>(tensor165, pindex));
      task164->add_dep(task165);
      task165->add_dep(task113);
      density_->add_task(task165);

      task165->add_dep(task4);

      std::vector<IndexRange> I178_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor<T>> I178(new Tensor<T>(I178_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor166 = {this->den1_, I178};
      std::shared_ptr<Task166<T>> task166(new Task166<T>(tensor166, pindex));
      task166->add_dep(task113);
      density_->add_task(task166);


      std::vector<IndexRange> I179_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I179(new Tensor<T>(I179_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor167 = {I178, t2, I179};
      std::shared_ptr<Task167<T>> task167(new Task167<T>(tensor167, pindex));
      task166->add_dep(task167);
      task167->add_dep(task113);
      density_->add_task(task167);


      std::vector<IndexRange> I180_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I180(new Tensor<T>(I180_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor168 = {I179, t2, I180};
      std::shared_ptr<Task168<T>> task168(new Task168<T>(tensor168, pindex));
      task167->add_dep(task168);
      task168->add_dep(task113);
      density_->add_task(task168);


      std::vector<std::shared_ptr<Tensor<T>>> tensor169 = {I180, Gamma14};
      std::shared_ptr<Task169<T>> task169(new Task169<T>(tensor169, pindex));
      task168->add_dep(task169);
      task169->add_dep(task113);
      density_->add_task(task169);

      task169->add_dep(task4);

      std::vector<IndexRange> I181_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I181(new Tensor<T>(I181_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor170 = {this->den1_, I181};
      std::shared_ptr<Task170<T>> task170(new Task170<T>(tensor170, pindex));
      task170->add_dep(task113);
      density_->add_task(task170);


      std::vector<IndexRange> I182_index = {this->active_, this->active_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I182(new Tensor<T>(I182_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor171 = {I181, t2, I182};
      std::shared_ptr<Task171<T>> task171(new Task171<T>(tensor171, pindex));
      task170->add_dep(task171);
      task171->add_dep(task113);
      density_->add_task(task171);


      std::vector<IndexRange> I183_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I183(new Tensor<T>(I183_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor172 = {I182, t2, I183};
      std::shared_ptr<Task172<T>> task172(new Task172<T>(tensor172, pindex));
      task171->add_dep(task172);
      task172->add_dep(task113);
      density_->add_task(task172);


      std::vector<std::shared_ptr<Tensor<T>>> tensor173 = {I183, Gamma52};
      std::shared_ptr<Task173<T>> task173(new Task173<T>(tensor173, pindex));
      task172->add_dep(task173);
      task173->add_dep(task113);
      density_->add_task(task173);

      task173->add_dep(task8);

      std::vector<IndexRange> I184_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I184(new Tensor<T>(I184_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor174 = {this->den1_, I184};
      std::shared_ptr<Task174<T>> task174(new Task174<T>(tensor174, pindex));
      task174->add_dep(task113);
      density_->add_task(task174);


      std::vector<IndexRange> I185_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I185(new Tensor<T>(I185_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor175 = {I184, t2, I185};
      std::shared_ptr<Task175<T>> task175(new Task175<T>(tensor175, pindex));
      task174->add_dep(task175);
      task175->add_dep(task113);
      density_->add_task(task175);


      std::vector<IndexRange> I186_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I186(new Tensor<T>(I186_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor176 = {I185, t2, I186};
      std::shared_ptr<Task176<T>> task176(new Task176<T>(tensor176, pindex));
      task175->add_dep(task176);
      task176->add_dep(task113);
      density_->add_task(task176);


      std::vector<std::shared_ptr<Tensor<T>>> tensor177 = {I186, Gamma14};
      std::shared_ptr<Task177<T>> task177(new Task177<T>(tensor177, pindex));
      task176->add_dep(task177);
      task177->add_dep(task113);
      density_->add_task(task177);

      task177->add_dep(task4);

      std::shared_ptr<Queue<T>> density2_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor178 = {this->den2_};
      std::shared_ptr<Task178<T>> task178(new Task178<T>(tensor178));
      density2_->add_task(task178);

      std::vector<IndexRange> I187_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I187(new Tensor<T>(I187_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor179 = {this->den2_, I187};
      std::shared_ptr<Task179<T>> task179(new Task179<T>(tensor179, pindex));
      task179->add_dep(task178);
      density2_->add_task(task179);


      std::vector<std::shared_ptr<Tensor<T>>> tensor180 = {I187, t2};
      std::shared_ptr<Task180<T>> task180(new Task180<T>(tensor180, pindex));
      task179->add_dep(task180);
      task180->add_dep(task178);
      density2_->add_task(task180);


      std::vector<IndexRange> I189_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I189(new Tensor<T>(I189_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor181 = {this->den2_, I189};
      std::shared_ptr<Task181<T>> task181(new Task181<T>(tensor181, pindex));
      task181->add_dep(task178);
      density2_->add_task(task181);


      std::vector<IndexRange> I190_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I190(new Tensor<T>(I190_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor182 = {I189, t2, I190};
      std::shared_ptr<Task182<T>> task182(new Task182<T>(tensor182, pindex));
      task181->add_dep(task182);
      task182->add_dep(task178);
      density2_->add_task(task182);


      std::vector<std::shared_ptr<Tensor<T>>> tensor183 = {I190, Gamma2};
      std::shared_ptr<Task183<T>> task183(new Task183<T>(tensor183, pindex));
      task182->add_dep(task183);
      task183->add_dep(task178);
      density2_->add_task(task183);

      task183->add_dep(task2);

      std::vector<IndexRange> I192_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I192(new Tensor<T>(I192_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor184 = {I189, t2, I192};
      std::shared_ptr<Task184<T>> task184(new Task184<T>(tensor184, pindex));
      task181->add_dep(task184);
      task184->add_dep(task178);
      density2_->add_task(task184);


      std::vector<std::shared_ptr<Tensor<T>>> tensor185 = {I192, Gamma2};
      std::shared_ptr<Task185<T>> task185(new Task185<T>(tensor185, pindex));
      task184->add_dep(task185);
      task185->add_dep(task178);
      density2_->add_task(task185);

      task185->add_dep(task2);

      std::vector<IndexRange> I193_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I193(new Tensor<T>(I193_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor186 = {this->den2_, I193};
      std::shared_ptr<Task186<T>> task186(new Task186<T>(tensor186, pindex));
      task186->add_dep(task178);
      density2_->add_task(task186);


      std::vector<IndexRange> I194_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I194(new Tensor<T>(I194_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor187 = {I193, t2, I194};
      std::shared_ptr<Task187<T>> task187(new Task187<T>(tensor187, pindex));
      task186->add_dep(task187);
      task187->add_dep(task178);
      density2_->add_task(task187);


      std::vector<std::shared_ptr<Tensor<T>>> tensor188 = {I194, Gamma14};
      std::shared_ptr<Task188<T>> task188(new Task188<T>(tensor188, pindex));
      task187->add_dep(task188);
      task188->add_dep(task178);
      density2_->add_task(task188);

      task188->add_dep(task4);

      return make_tuple(queue_, energy_, dedci_, density_, correction_, density2_);
    };

  public:
    CAS_test(std::shared_ptr<const Reference> ref) : SpinFreeMethod<T>(ref), SMITH_info() {
      this->eig_ = this->f1_->diag();
      t2 = this->v2_->clone();
      e0_ = this->e0();
      this->update_amplitude(t2, this->v2_, true);
      t2->scale(2.0);
      r = t2->clone();
      this->den1_ = this->h1_->clone();
      this->den2_ = this->v2_->clone();
    };
    ~CAS_test() {}; 

    void solve() {
      this->print_iteration();
      int iter = 0;
      std::shared_ptr<Queue<T>> queue, energ, dec, dens, correct, dens2;
      for ( ; iter != maxiter_; ++iter) {
        std::tie(queue, energ, dec, dens, correct, dens2) = make_queue_();
        while (!queue->done())
          queue->next_compute();
        this->update_amplitude(t2, r);
        const double err = r->rms();
        r->zero();
        const double en = energy(energ);
        this->print_iteration(iter, en, err);
        if (err < thresh_residual()) break;
      }
      this->print_iteration(iter == maxiter_);
      const double de = dedci(dec);
      std::cout << "CI derivative energy (testing) = " << std::setprecision(10) << de  << std::endl;
      std::cout << std::endl;

      std::cout << " === Unrelaxed density matrix, dm1, <1|E_pq|1> + 2<0|E_pq|1> ===" << std::endl; 
      while (!dens->done())
        dens->next_compute();
#if 0
      this->den1_->print2("density matrix", 1.0e-5);
#endif
      this->correct_den1_ = correction(correct);
      std::cout << "Unlinked correction term, <1|1>*rdm1 = " << std::setprecision(10) << this->correct_den1_ << "*rdm1" << std::endl;
#if 0
      this->rdm1_->print2("rdm1", 1.0e-5);
#endif
      std::cout << " === Unrelaxed density matrix, dm2, <0|E_pqrs|1>  ===" << std::endl; 
      while (!dens2->done())
        dens2->next_compute();
#if 0
      this->den2_->print4("density matrix", 1.0e-5);
#endif
    };

    double energy(std::shared_ptr<Queue<T>> energ) {
      double en = 0.0;
      while (!energ->done()) {
        std::shared_ptr<Task<T>> c = energ->next_compute();
        en += c->energy();
      }   
      return en; 
    };  

    double dedci(std::shared_ptr<Queue<T>> dec) {
      double de = 0.0;
      while (!dec->done()) {
        std::shared_ptr<Task<T>> d = dec->next_compute();
        de += d->dedci();
      }   
      return de; 
    };  

    double correction(std::shared_ptr<Queue<T>> correct) {
      double n = 0.0;
      while (!correct->done()) {
        std::shared_ptr<Task<T>> c = correct->next_compute();
        n += c->correction();
      }   
      return n; 
    };  

};

}
}
}
#endif

