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
    std::shared_ptr<Tensor<T>> den1;
    std::shared_ptr<Tensor<T>> den2;
    double correct_den1;
    std::shared_ptr<Tensor<T>> deci;

    std::tuple<std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>, std::shared_ptr<Queue<T>>> make_queue_() {
      std::shared_ptr<Queue<T>> queue_(new Queue<T>());
      std::array<std::shared_ptr<const IndexRange>,3> pindex = {{this->rclosed_, this->ractive_, this->rvirt_}};
      std::array<std::shared_ptr<const IndexRange>,4> cindex = {{this->rclosed_, this->ractive_, this->rvirt_, this->rci_}};

      std::vector<std::shared_ptr<Tensor<T>>> tensor0 = {r};
      std::shared_ptr<Task0<T>> task0(new Task0<T>(tensor0));
      queue_->add_task(task0);

      std::vector<IndexRange> Gamma0_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma0(new Tensor<T>(Gamma0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor1 = {Gamma0, this->rdm2_, this->f1_};
      std::shared_ptr<Task1<T>> task1(new Task1<T>(tensor1, pindex));
      task1->add_dep(task0);
      queue_->add_task(task1);

      std::vector<IndexRange> Gamma2_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma2(new Tensor<T>(Gamma2_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor2 = {Gamma2, this->rdm1_};
      std::shared_ptr<Task2<T>> task2(new Task2<T>(tensor2, pindex));
      task2->add_dep(task0);
      queue_->add_task(task2);

      std::vector<IndexRange> Gamma14_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma14(new Tensor<T>(Gamma14_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor3 = {Gamma14, this->rdm1deriv_};
      std::shared_ptr<Task3<T>> task3(new Task3<T>(tensor3, cindex));
      task3->add_dep(task0);
      queue_->add_task(task3);

      std::vector<IndexRange> Gamma20_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> Gamma20(new Tensor<T>(Gamma20_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor4 = {Gamma20, this->rdm2_};
      std::shared_ptr<Task4<T>> task4(new Task4<T>(tensor4, pindex));
      task4->add_dep(task0);
      queue_->add_task(task4);

      std::vector<IndexRange> I0_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I0(new Tensor<T>(I0_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor5 = {r, I0};
      std::shared_ptr<Task5<T>> task5(new Task5<T>(tensor5, pindex));
      task5->add_dep(task0);
      queue_->add_task(task5);


      std::vector<IndexRange> I1_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I1(new Tensor<T>(I1_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor6 = {I0, t2, I1};
      std::shared_ptr<Task6<T>> task6(new Task6<T>(tensor6, pindex));
      task5->add_dep(task6);
      task6->add_dep(task0);
      queue_->add_task(task6);


      std::vector<std::shared_ptr<Tensor<T>>> tensor7 = {I1, Gamma0};
      std::shared_ptr<Task7<T>> task7(new Task7<T>(tensor7, pindex));
      task6->add_dep(task7);
      task7->add_dep(task0);
      queue_->add_task(task7);

      task7->add_dep(task1);

      std::vector<IndexRange> I3_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I3(new Tensor<T>(I3_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor8 = {I0, t2, I3};
      std::shared_ptr<Task8<T>> task8(new Task8<T>(tensor8, pindex));
      task5->add_dep(task8);
      task8->add_dep(task0);
      queue_->add_task(task8);


      std::vector<std::shared_ptr<Tensor<T>>> tensor9 = {I3, Gamma0};
      std::shared_ptr<Task9<T>> task9(new Task9<T>(tensor9, pindex));
      task8->add_dep(task9);
      task9->add_dep(task0);
      queue_->add_task(task9);

      task9->add_dep(task1);

      std::vector<IndexRange> I5_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I5(new Tensor<T>(I5_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor10 = {I0, this->f1_, I5};
      std::shared_ptr<Task10<T>> task10(new Task10<T>(tensor10, pindex));
      task5->add_dep(task10);
      task10->add_dep(task0);
      queue_->add_task(task10);


      std::vector<IndexRange> I6_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I6(new Tensor<T>(I6_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor11 = {I5, t2, I6};
      std::shared_ptr<Task11<T>> task11(new Task11<T>(tensor11, pindex));
      task10->add_dep(task11);
      task11->add_dep(task0);
      queue_->add_task(task11);


      std::vector<std::shared_ptr<Tensor<T>>> tensor12 = {I6, Gamma2};
      std::shared_ptr<Task12<T>> task12(new Task12<T>(tensor12, pindex));
      task11->add_dep(task12);
      task12->add_dep(task0);
      queue_->add_task(task12);

      task12->add_dep(task2);

      std::vector<IndexRange> I9_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I9(new Tensor<T>(I9_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor13 = {I5, t2, I9};
      std::shared_ptr<Task13<T>> task13(new Task13<T>(tensor13, pindex));
      task10->add_dep(task13);
      task13->add_dep(task0);
      queue_->add_task(task13);


      std::vector<std::shared_ptr<Tensor<T>>> tensor14 = {I9, Gamma2};
      std::shared_ptr<Task14<T>> task14(new Task14<T>(tensor14, pindex));
      task13->add_dep(task14);
      task14->add_dep(task0);
      queue_->add_task(task14);

      task14->add_dep(task2);

      std::vector<IndexRange> I11_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I11(new Tensor<T>(I11_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor15 = {I0, this->f1_, I11};
      std::shared_ptr<Task15<T>> task15(new Task15<T>(tensor15, pindex));
      task5->add_dep(task15);
      task15->add_dep(task0);
      queue_->add_task(task15);


      std::vector<IndexRange> I12_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I12(new Tensor<T>(I12_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor16 = {I11, t2, I12};
      std::shared_ptr<Task16<T>> task16(new Task16<T>(tensor16, pindex));
      task15->add_dep(task16);
      task16->add_dep(task0);
      queue_->add_task(task16);


      std::vector<std::shared_ptr<Tensor<T>>> tensor17 = {I12, Gamma2};
      std::shared_ptr<Task17<T>> task17(new Task17<T>(tensor17, pindex));
      task16->add_dep(task17);
      task17->add_dep(task0);
      queue_->add_task(task17);

      task17->add_dep(task2);

      std::vector<IndexRange> I15_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I15(new Tensor<T>(I15_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor18 = {I11, t2, I15};
      std::shared_ptr<Task18<T>> task18(new Task18<T>(tensor18, pindex));
      task15->add_dep(task18);
      task18->add_dep(task0);
      queue_->add_task(task18);


      std::vector<std::shared_ptr<Tensor<T>>> tensor19 = {I15, Gamma2};
      std::shared_ptr<Task19<T>> task19(new Task19<T>(tensor19, pindex));
      task18->add_dep(task19);
      task19->add_dep(task0);
      queue_->add_task(task19);

      task19->add_dep(task2);

      std::vector<IndexRange> I17_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I17(new Tensor<T>(I17_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor20 = {I0, this->f1_, I17};
      std::shared_ptr<Task20<T>> task20(new Task20<T>(tensor20, pindex));
      task5->add_dep(task20);
      task20->add_dep(task0);
      queue_->add_task(task20);


      std::vector<IndexRange> I18_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I18(new Tensor<T>(I18_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor21 = {I17, t2, I18};
      std::shared_ptr<Task21<T>> task21(new Task21<T>(tensor21, pindex));
      task20->add_dep(task21);
      task21->add_dep(task0);
      queue_->add_task(task21);


      std::vector<std::shared_ptr<Tensor<T>>> tensor22 = {I18, Gamma2};
      std::shared_ptr<Task22<T>> task22(new Task22<T>(tensor22, pindex));
      task21->add_dep(task22);
      task22->add_dep(task0);
      queue_->add_task(task22);

      task22->add_dep(task2);

      std::vector<IndexRange> I21_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I21(new Tensor<T>(I21_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor23 = {I17, t2, I21};
      std::shared_ptr<Task23<T>> task23(new Task23<T>(tensor23, pindex));
      task20->add_dep(task23);
      task23->add_dep(task0);
      queue_->add_task(task23);


      std::vector<std::shared_ptr<Tensor<T>>> tensor24 = {I21, Gamma2};
      std::shared_ptr<Task24<T>> task24(new Task24<T>(tensor24, pindex));
      task23->add_dep(task24);
      task24->add_dep(task0);
      queue_->add_task(task24);

      task24->add_dep(task2);

      std::vector<IndexRange> I23_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I23(new Tensor<T>(I23_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor25 = {I0, t2, I23};
      std::shared_ptr<Task25<T>> task25(new Task25<T>(tensor25, pindex));
      task5->add_dep(task25);
      task25->add_dep(task0);
      queue_->add_task(task25);


      std::vector<std::shared_ptr<Tensor<T>>> tensor26 = {I23, Gamma2};
      std::shared_ptr<Task26<T>> task26(new Task26<T>(tensor26, pindex, this->e0_));
      task25->add_dep(task26);
      task26->add_dep(task0);
      queue_->add_task(task26);

      task26->add_dep(task2);

      std::vector<IndexRange> I25_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I25(new Tensor<T>(I25_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor27 = {I0, t2, I25};
      std::shared_ptr<Task27<T>> task27(new Task27<T>(tensor27, pindex));
      task5->add_dep(task27);
      task27->add_dep(task0);
      queue_->add_task(task27);


      std::vector<std::shared_ptr<Tensor<T>>> tensor28 = {I25, Gamma2};
      std::shared_ptr<Task28<T>> task28(new Task28<T>(tensor28, pindex, this->e0_));
      task27->add_dep(task28);
      task28->add_dep(task0);
      queue_->add_task(task28);

      task28->add_dep(task2);

      std::vector<IndexRange> I27_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I27(new Tensor<T>(I27_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor29 = {I0, this->v2_, I27};
      std::shared_ptr<Task29<T>> task29(new Task29<T>(tensor29, pindex));
      task5->add_dep(task29);
      task29->add_dep(task0);
      queue_->add_task(task29);


      std::vector<std::shared_ptr<Tensor<T>>> tensor30 = {I27, Gamma2};
      std::shared_ptr<Task30<T>> task30(new Task30<T>(tensor30, pindex));
      task29->add_dep(task30);
      task30->add_dep(task0);
      queue_->add_task(task30);

      task30->add_dep(task2);

      std::vector<IndexRange> I29_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I29(new Tensor<T>(I29_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor31 = {I0, this->v2_, I29};
      std::shared_ptr<Task31<T>> task31(new Task31<T>(tensor31, pindex));
      task5->add_dep(task31);
      task31->add_dep(task0);
      queue_->add_task(task31);


      std::vector<std::shared_ptr<Tensor<T>>> tensor32 = {I29, Gamma2};
      std::shared_ptr<Task32<T>> task32(new Task32<T>(tensor32, pindex));
      task31->add_dep(task32);
      task32->add_dep(task0);
      queue_->add_task(task32);

      task32->add_dep(task2);

      std::shared_ptr<Queue<T>> energy_(new Queue<T>());
      std::vector<IndexRange> I30_index;
      std::shared_ptr<Tensor<T>> I30(new Tensor<T>(I30_index, false));
      std::vector<IndexRange> I31_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I31(new Tensor<T>(I31_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor33 = {I30, t2, I31};
      std::shared_ptr<Task33<T>> task33(new Task33<T>(tensor33, pindex));
      energy_->add_task(task33);


      std::vector<IndexRange> I32_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I32(new Tensor<T>(I32_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor34 = {I31, this->v2_, I32};
      std::shared_ptr<Task34<T>> task34(new Task34<T>(tensor34, pindex));
      task33->add_dep(task34);
      energy_->add_task(task34);


      std::vector<std::shared_ptr<Tensor<T>>> tensor35 = {I32, Gamma2};
      std::shared_ptr<Task35<T>> task35(new Task35<T>(tensor35, pindex));
      task34->add_dep(task35);
      energy_->add_task(task35);

      task35->add_dep(task2);

      std::vector<IndexRange> I35_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I35(new Tensor<T>(I35_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor36 = {I31, this->v2_, I35};
      std::shared_ptr<Task36<T>> task36(new Task36<T>(tensor36, pindex));
      task33->add_dep(task36);
      energy_->add_task(task36);


      std::vector<std::shared_ptr<Tensor<T>>> tensor37 = {I35, Gamma2};
      std::shared_ptr<Task37<T>> task37(new Task37<T>(tensor37, pindex));
      task36->add_dep(task37);
      energy_->add_task(task37);

      task37->add_dep(task2);

      std::shared_ptr<Queue<T>> dedci_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor38 = {deci};
      std::shared_ptr<Task38<T>> task38(new Task38<T>(tensor38));
      dedci_->add_task(task38);

      std::vector<IndexRange> I36_index = {this->ci_};
      std::shared_ptr<Tensor<T>> I36(new Tensor<T>(I36_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor39 = {deci, I36};
      std::shared_ptr<Task39<T>> task39(new Task39<T>(tensor39, cindex));
      task39->add_dep(task38);
      dedci_->add_task(task39);


      std::vector<IndexRange> I37_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I37(new Tensor<T>(I37_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor40 = {I36, t2, I37};
      std::shared_ptr<Task40<T>> task40(new Task40<T>(tensor40, cindex));
      task39->add_dep(task40);
      task40->add_dep(task38);
      dedci_->add_task(task40);


      std::vector<IndexRange> I38_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I38(new Tensor<T>(I38_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor41 = {I37, this->v2_, I38};
      std::shared_ptr<Task41<T>> task41(new Task41<T>(tensor41, cindex));
      task40->add_dep(task41);
      task41->add_dep(task38);
      dedci_->add_task(task41);


      std::vector<std::shared_ptr<Tensor<T>>> tensor42 = {I38, Gamma14};
      std::shared_ptr<Task42<T>> task42(new Task42<T>(tensor42, cindex));
      task41->add_dep(task42);
      task42->add_dep(task38);
      dedci_->add_task(task42);

      task42->add_dep(task3);

      std::vector<IndexRange> I41_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I41(new Tensor<T>(I41_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor43 = {I37, this->v2_, I41};
      std::shared_ptr<Task43<T>> task43(new Task43<T>(tensor43, cindex));
      task40->add_dep(task43);
      task43->add_dep(task38);
      dedci_->add_task(task43);


      std::vector<std::shared_ptr<Tensor<T>>> tensor44 = {I41, Gamma14};
      std::shared_ptr<Task44<T>> task44(new Task44<T>(tensor44, cindex));
      task43->add_dep(task44);
      task44->add_dep(task38);
      dedci_->add_task(task44);

      task44->add_dep(task3);

      std::vector<IndexRange> I43_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I43(new Tensor<T>(I43_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor45 = {I36, t2, I43};
      std::shared_ptr<Task45<T>> task45(new Task45<T>(tensor45, cindex));
      task39->add_dep(task45);
      task45->add_dep(task38);
      dedci_->add_task(task45);


      std::vector<IndexRange> I44_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I44(new Tensor<T>(I44_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor46 = {I43, this->v2_, I44};
      std::shared_ptr<Task46<T>> task46(new Task46<T>(tensor46, cindex));
      task45->add_dep(task46);
      task46->add_dep(task38);
      dedci_->add_task(task46);


      std::vector<std::shared_ptr<Tensor<T>>> tensor47 = {I44, Gamma14};
      std::shared_ptr<Task47<T>> task47(new Task47<T>(tensor47, cindex));
      task46->add_dep(task47);
      task47->add_dep(task38);
      dedci_->add_task(task47);

      task47->add_dep(task3);

      std::vector<IndexRange> I47_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I47(new Tensor<T>(I47_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor48 = {I43, this->v2_, I47};
      std::shared_ptr<Task48<T>> task48(new Task48<T>(tensor48, cindex));
      task45->add_dep(task48);
      task48->add_dep(task38);
      dedci_->add_task(task48);


      std::vector<std::shared_ptr<Tensor<T>>> tensor49 = {I47, Gamma14};
      std::shared_ptr<Task49<T>> task49(new Task49<T>(tensor49, cindex));
      task48->add_dep(task49);
      task49->add_dep(task38);
      dedci_->add_task(task49);

      task49->add_dep(task3);

      std::shared_ptr<Queue<T>> correction_(new Queue<T>());
      std::vector<IndexRange> I48_index;
      std::shared_ptr<Tensor<T>> I48(new Tensor<T>(I48_index, false));
      std::vector<IndexRange> I49_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I49(new Tensor<T>(I49_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor50 = {I48, t2, I49};
      std::shared_ptr<Task50<T>> task50(new Task50<T>(tensor50, pindex));
      correction_->add_task(task50);


      std::vector<IndexRange> I50_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I50(new Tensor<T>(I50_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor51 = {I49, t2, I50};
      std::shared_ptr<Task51<T>> task51(new Task51<T>(tensor51, pindex));
      task50->add_dep(task51);
      correction_->add_task(task51);


      std::vector<std::shared_ptr<Tensor<T>>> tensor52 = {I50, Gamma2};
      std::shared_ptr<Task52<T>> task52(new Task52<T>(tensor52, pindex));
      task51->add_dep(task52);
      correction_->add_task(task52);

      task52->add_dep(task2);

      std::vector<IndexRange> I53_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I53(new Tensor<T>(I53_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor53 = {I49, t2, I53};
      std::shared_ptr<Task53<T>> task53(new Task53<T>(tensor53, pindex));
      task50->add_dep(task53);
      correction_->add_task(task53);


      std::vector<std::shared_ptr<Tensor<T>>> tensor54 = {I53, Gamma2};
      std::shared_ptr<Task54<T>> task54(new Task54<T>(tensor54, pindex));
      task53->add_dep(task54);
      correction_->add_task(task54);

      task54->add_dep(task2);

      std::shared_ptr<Queue<T>> density_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor55 = {den1};
      std::shared_ptr<Task55<T>> task55(new Task55<T>(tensor55));
      density_->add_task(task55);

      std::vector<IndexRange> I54_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I54(new Tensor<T>(I54_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor56 = {den1, I54};
      std::shared_ptr<Task56<T>> task56(new Task56<T>(tensor56, pindex));
      task56->add_dep(task55);
      density_->add_task(task56);


      std::vector<IndexRange> I55_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I55(new Tensor<T>(I55_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor57 = {I54, t2, I55};
      std::shared_ptr<Task57<T>> task57(new Task57<T>(tensor57, pindex));
      task56->add_dep(task57);
      task57->add_dep(task55);
      density_->add_task(task57);


      std::vector<IndexRange> I56_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I56(new Tensor<T>(I56_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor58 = {I55, t2, I56};
      std::shared_ptr<Task58<T>> task58(new Task58<T>(tensor58, pindex));
      task57->add_dep(task58);
      task58->add_dep(task55);
      density_->add_task(task58);


      std::vector<std::shared_ptr<Tensor<T>>> tensor59 = {I56, Gamma20};
      std::shared_ptr<Task59<T>> task59(new Task59<T>(tensor59, pindex));
      task58->add_dep(task59);
      task59->add_dep(task55);
      density_->add_task(task59);

      task59->add_dep(task4);

      std::vector<IndexRange> I59_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I59(new Tensor<T>(I59_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor60 = {I55, t2, I59};
      std::shared_ptr<Task60<T>> task60(new Task60<T>(tensor60, pindex));
      task57->add_dep(task60);
      task60->add_dep(task55);
      density_->add_task(task60);


      std::vector<std::shared_ptr<Tensor<T>>> tensor61 = {I59, Gamma20};
      std::shared_ptr<Task61<T>> task61(new Task61<T>(tensor61, pindex));
      task60->add_dep(task61);
      task61->add_dep(task55);
      density_->add_task(task61);

      task61->add_dep(task4);

      std::vector<IndexRange> I60_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor<T>> I60(new Tensor<T>(I60_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor62 = {den1, I60};
      std::shared_ptr<Task62<T>> task62(new Task62<T>(tensor62, pindex));
      task62->add_dep(task55);
      density_->add_task(task62);


      std::vector<IndexRange> I61_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I61(new Tensor<T>(I61_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor63 = {I60, t2, I61};
      std::shared_ptr<Task63<T>> task63(new Task63<T>(tensor63, pindex));
      task62->add_dep(task63);
      task63->add_dep(task55);
      density_->add_task(task63);


      std::vector<IndexRange> I62_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I62(new Tensor<T>(I62_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor64 = {I61, t2, I62};
      std::shared_ptr<Task64<T>> task64(new Task64<T>(tensor64, pindex));
      task63->add_dep(task64);
      task64->add_dep(task55);
      density_->add_task(task64);


      std::vector<std::shared_ptr<Tensor<T>>> tensor65 = {I62, Gamma2};
      std::shared_ptr<Task65<T>> task65(new Task65<T>(tensor65, pindex));
      task64->add_dep(task65);
      task65->add_dep(task55);
      density_->add_task(task65);

      task65->add_dep(task2);

      std::vector<IndexRange> I65_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I65(new Tensor<T>(I65_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor66 = {I61, t2, I65};
      std::shared_ptr<Task66<T>> task66(new Task66<T>(tensor66, pindex));
      task63->add_dep(task66);
      task66->add_dep(task55);
      density_->add_task(task66);


      std::vector<std::shared_ptr<Tensor<T>>> tensor67 = {I65, Gamma2};
      std::shared_ptr<Task67<T>> task67(new Task67<T>(tensor67, pindex));
      task66->add_dep(task67);
      task67->add_dep(task55);
      density_->add_task(task67);

      task67->add_dep(task2);

      std::vector<IndexRange> I66_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I66(new Tensor<T>(I66_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor68 = {den1, I66};
      std::shared_ptr<Task68<T>> task68(new Task68<T>(tensor68, pindex));
      task68->add_dep(task55);
      density_->add_task(task68);


      std::vector<IndexRange> I67_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I67(new Tensor<T>(I67_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor69 = {I66, t2, I67};
      std::shared_ptr<Task69<T>> task69(new Task69<T>(tensor69, pindex));
      task68->add_dep(task69);
      task69->add_dep(task55);
      density_->add_task(task69);


      std::vector<IndexRange> I68_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I68(new Tensor<T>(I68_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor70 = {I67, t2, I68};
      std::shared_ptr<Task70<T>> task70(new Task70<T>(tensor70, pindex));
      task69->add_dep(task70);
      task70->add_dep(task55);
      density_->add_task(task70);


      std::vector<std::shared_ptr<Tensor<T>>> tensor71 = {I68, Gamma2};
      std::shared_ptr<Task71<T>> task71(new Task71<T>(tensor71, pindex));
      task70->add_dep(task71);
      task71->add_dep(task55);
      density_->add_task(task71);

      task71->add_dep(task2);

      std::vector<IndexRange> I71_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I71(new Tensor<T>(I71_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor72 = {I67, t2, I71};
      std::shared_ptr<Task72<T>> task72(new Task72<T>(tensor72, pindex));
      task69->add_dep(task72);
      task72->add_dep(task55);
      density_->add_task(task72);


      std::vector<std::shared_ptr<Tensor<T>>> tensor73 = {I71, Gamma2};
      std::shared_ptr<Task73<T>> task73(new Task73<T>(tensor73, pindex));
      task72->add_dep(task73);
      task73->add_dep(task55);
      density_->add_task(task73);

      task73->add_dep(task2);

      std::vector<IndexRange> I72_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor<T>> I72(new Tensor<T>(I72_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor74 = {den1, I72};
      std::shared_ptr<Task74<T>> task74(new Task74<T>(tensor74, pindex));
      task74->add_dep(task55);
      density_->add_task(task74);


      std::vector<IndexRange> I73_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I73(new Tensor<T>(I73_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor75 = {I72, t2, I73};
      std::shared_ptr<Task75<T>> task75(new Task75<T>(tensor75, pindex));
      task74->add_dep(task75);
      task75->add_dep(task55);
      density_->add_task(task75);


      std::vector<IndexRange> I74_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I74(new Tensor<T>(I74_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor76 = {I73, t2, I74};
      std::shared_ptr<Task76<T>> task76(new Task76<T>(tensor76, pindex));
      task75->add_dep(task76);
      task76->add_dep(task55);
      density_->add_task(task76);


      std::vector<std::shared_ptr<Tensor<T>>> tensor77 = {I74, Gamma2};
      std::shared_ptr<Task77<T>> task77(new Task77<T>(tensor77, pindex));
      task76->add_dep(task77);
      task77->add_dep(task55);
      density_->add_task(task77);

      task77->add_dep(task2);

      std::vector<IndexRange> I77_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I77(new Tensor<T>(I77_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor78 = {I73, t2, I77};
      std::shared_ptr<Task78<T>> task78(new Task78<T>(tensor78, pindex));
      task75->add_dep(task78);
      task78->add_dep(task55);
      density_->add_task(task78);


      std::vector<std::shared_ptr<Tensor<T>>> tensor79 = {I77, Gamma2};
      std::shared_ptr<Task79<T>> task79(new Task79<T>(tensor79, pindex));
      task78->add_dep(task79);
      task79->add_dep(task55);
      density_->add_task(task79);

      task79->add_dep(task2);

      std::shared_ptr<Queue<T>> density2_(new Queue<T>());
      std::vector<std::shared_ptr<Tensor<T>>> tensor80 = {den2};
      std::shared_ptr<Task80<T>> task80(new Task80<T>(tensor80));
      density2_->add_task(task80);

      std::vector<IndexRange> I78_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor<T>> I78(new Tensor<T>(I78_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor81 = {den2, I78};
      std::shared_ptr<Task81<T>> task81(new Task81<T>(tensor81, pindex));
      task81->add_dep(task80);
      density2_->add_task(task81);


      std::vector<IndexRange> I79_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I79(new Tensor<T>(I79_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor82 = {I78, t2, I79};
      std::shared_ptr<Task82<T>> task82(new Task82<T>(tensor82, pindex));
      task81->add_dep(task82);
      task82->add_dep(task80);
      density2_->add_task(task82);


      std::vector<std::shared_ptr<Tensor<T>>> tensor83 = {I79, Gamma2};
      std::shared_ptr<Task83<T>> task83(new Task83<T>(tensor83, pindex));
      task82->add_dep(task83);
      task83->add_dep(task80);
      density2_->add_task(task83);

      task83->add_dep(task2);

      std::vector<IndexRange> I81_index = {this->active_, this->active_};
      std::shared_ptr<Tensor<T>> I81(new Tensor<T>(I81_index, false));
      std::vector<std::shared_ptr<Tensor<T>>> tensor84 = {I78, t2, I81};
      std::shared_ptr<Task84<T>> task84(new Task84<T>(tensor84, pindex));
      task81->add_dep(task84);
      task84->add_dep(task80);
      density2_->add_task(task84);


      std::vector<std::shared_ptr<Tensor<T>>> tensor85 = {I81, Gamma2};
      std::shared_ptr<Task85<T>> task85(new Task85<T>(tensor85, pindex));
      task84->add_dep(task85);
      task85->add_dep(task80);
      density2_->add_task(task85);

      task85->add_dep(task2);

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
      den1 = this->h1_->clone();
      den2 = this->v2_->clone();
      deci = this->rdm0deriv_->clone();
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
      std::cout << " === Calculating CI derivative dE/dcI ===" << std::endl;
      while (!dec->done())
        dec->next_compute();
      deci->print1("CI derivative tensor: ", 1.0e-15);
      std::cout << std::endl;
      std::cout << "CI derivative * cI  = " << std::setprecision(10) <<  deci->dot_product(this->rdm0deriv_) << std::endl;
      std::cout << std::endl;

      std::cout << " === Unrelaxed density matrix, dm1, <1|E_pq|1> + 2<0|E_pq|1> ===" << std::endl;
      while (!dens->done())
        dens->next_compute();
#if 0
      den1->print2("density matrix", 1.0e-5);
#endif
      correct_den1 = correction(correct);
      std::cout << "Unlinked correction term, <1|1>*rdm1 = " << std::setprecision(10) << correct_den1 << "*rdm1" << std::endl;
      std::cout << " === Unrelaxed density matrix, dm2, <0|E_pqrs|1>  ===" << std::endl;
      while (!dens2->done())
        dens2->next_compute();
#if 0
      den2->print4("density matrix", 1.0e-5);
#endif
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

    double rdm1_correction() const { return correct_den1; }

    std::shared_ptr<const Civec> ci_deriv() const { return deci->civec(this->det_); }

};

}
}
}
#endif

