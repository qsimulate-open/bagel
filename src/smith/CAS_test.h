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

#include <iostream>
#include <tuple>
#include <iomanip>
#include <src/smith/spinfreebase.h>
#include <src/scf/hf/fock.h>
#include <src/util/f77.h>
#include <src/smith/queue.h>
#include <src/smith/CAS_test_tasks.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace CAS_test{

class CAS_test : public SpinFreeMethod {
  protected:
    using SpinFreeMethod::ref_;

    std::shared_ptr<Tensor> t2;
    std::shared_ptr<Tensor> r;
    double e0_;
    std::shared_ptr<Tensor> sigma_;
    std::shared_ptr<Tensor> den1;
    std::shared_ptr<Tensor> den2;
    std::shared_ptr<Tensor> Den1;
    double correlated_norm;
    std::shared_ptr<Tensor> deci;

    std::tuple<std::shared_ptr<Queue>, std::shared_ptr<Queue>, std::shared_ptr<Queue>,  std::shared_ptr<Queue>,  std::shared_ptr<Queue>, std::shared_ptr<Queue>, std::shared_ptr<Queue>> make_queue_() {
      auto queue_ = std::make_shared<Queue>();
      std::array<std::shared_ptr<const IndexRange>,3> pindex = {{this->rclosed_, this->ractive_, this->rvirt_}};
      std::array<std::shared_ptr<const IndexRange>,4> cindex = {{this->rclosed_, this->ractive_, this->rvirt_, this->rci_}};

      std::vector<std::shared_ptr<Tensor>> tensor0 = {r};
      auto task0 = std::make_shared<Task0>(tensor0);
      queue_->add_task(task0);

      std::vector<IndexRange> Gamma0_index;
      std::shared_ptr<Tensor> Gamma0(new Tensor(Gamma0_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor1 = {Gamma0, this->rdm1_, this->f1_};
      auto task1 = std::make_shared<Task1>(tensor1, pindex);
      task1->add_dep(task0);
      queue_->add_task(task1);

      std::vector<IndexRange> Gamma2_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma2(new Tensor(Gamma2_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor2 = {Gamma2, this->rdm1_};
      auto task2 = std::make_shared<Task2>(tensor2, pindex);
      task2->add_dep(task0);
      queue_->add_task(task2);

      std::vector<IndexRange> Gamma6_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma6(new Tensor(Gamma6_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor3 = {Gamma6, this->rdm2_, this->f1_};
      auto task3 = std::make_shared<Task3>(tensor3, pindex);
      task3->add_dep(task0);
      queue_->add_task(task3);

      std::vector<IndexRange> Gamma14_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma14(new Tensor(Gamma14_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor4 = {Gamma14, this->rdm2_};
      auto task4 = std::make_shared<Task4>(tensor4, pindex);
      task4->add_dep(task0);
      queue_->add_task(task4);

      std::vector<IndexRange> Gamma16_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma16(new Tensor(Gamma16_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor5 = {Gamma16, this->rdm3_, this->f1_};
      auto task5 = std::make_shared<Task5>(tensor5, pindex);
      task5->add_dep(task0);
      queue_->add_task(task5);

      std::vector<IndexRange> Gamma67_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma67(new Tensor(Gamma67_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor6 = {Gamma67, this->rdm3_};
      auto task6 = std::make_shared<Task6>(tensor6, pindex);
      task6->add_dep(task0);
      queue_->add_task(task6);

      std::vector<IndexRange> Gamma72_index = {this->ci_};
      std::shared_ptr<Tensor> Gamma72(new Tensor(Gamma72_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor7 = {Gamma72, this->rdm1deriv_, this->f1_};
      auto task7 = std::make_shared<Task7>(tensor7, cindex);
      task7->add_dep(task0);
      queue_->add_task(task7);

      std::vector<IndexRange> Gamma74_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma74(new Tensor(Gamma74_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor8 = {Gamma74, this->rdm1deriv_};
      auto task8 = std::make_shared<Task8>(tensor8, cindex);
      task8->add_dep(task0);
      queue_->add_task(task8);

      std::vector<IndexRange> Gamma78_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma78(new Tensor(Gamma78_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor9 = {Gamma78, this->rdm2deriv_, this->f1_};
      auto task9 = std::make_shared<Task9>(tensor9, cindex);
      task9->add_dep(task0);
      queue_->add_task(task9);

      std::vector<IndexRange> Gamma86_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma86(new Tensor(Gamma86_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor10 = {Gamma86, this->rdm2deriv_};
      auto task10 = std::make_shared<Task10>(tensor10, cindex);
      task10->add_dep(task0);
      queue_->add_task(task10);

      std::vector<IndexRange> Gamma88_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> Gamma88(new Tensor(Gamma88_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor11 = {Gamma88, this->rdm3deriv_, this->f1_};
      auto task11 = std::make_shared<Task11>(tensor11, cindex);
      task11->add_dep(task0);
      queue_->add_task(task11);

      std::vector<IndexRange> I0_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I0(new Tensor(I0_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor12 = {r, I0};
      auto task12 = std::make_shared<Task12>(tensor12, pindex);
      task12->add_dep(task0);
      queue_->add_task(task12);


      std::vector<std::shared_ptr<Tensor>> tensor13 = {I0, t2, this->v2_};
      auto task13 = std::make_shared<Task13>(tensor13, pindex, this->e0_);
      task12->add_dep(task13);
      task13->add_dep(task0);
      queue_->add_task(task13);


      std::vector<IndexRange> I1_index;
      std::shared_ptr<Tensor> I1(new Tensor(I1_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor14 = {I0, t2, I1};
      auto task14 = std::make_shared<Task14>(tensor14, pindex);
      task12->add_dep(task14);
      task14->add_dep(task0);
      queue_->add_task(task14);


      std::vector<std::shared_ptr<Tensor>> tensor15 = {I1, Gamma0};
      auto task15 = std::make_shared<Task15>(tensor15, pindex);
      task14->add_dep(task15);
      task15->add_dep(task0);
      queue_->add_task(task15);

      task15->add_dep(task1);

      std::vector<IndexRange> I3_index;
      std::shared_ptr<Tensor> I3(new Tensor(I3_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor16 = {I0, t2, I3};
      auto task16 = std::make_shared<Task16>(tensor16, pindex);
      task12->add_dep(task16);
      task16->add_dep(task0);
      queue_->add_task(task16);


      std::vector<std::shared_ptr<Tensor>> tensor17 = {I3, Gamma0};
      auto task17 = std::make_shared<Task17>(tensor17, pindex);
      task16->add_dep(task17);
      task17->add_dep(task0);
      queue_->add_task(task17);

      task17->add_dep(task1);

      std::vector<IndexRange> I4_index = {this->closed_, this->virt_, this->virt_, this->closed_};
      std::shared_ptr<Tensor> I4(new Tensor(I4_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor18 = {r, I4};
      auto task18 = std::make_shared<Task18>(tensor18, pindex);
      task18->add_dep(task0);
      queue_->add_task(task18);


      std::vector<IndexRange> I5_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I5(new Tensor(I5_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor19 = {I4, this->f1_, I5};
      auto task19 = std::make_shared<Task19>(tensor19, pindex);
      task18->add_dep(task19);
      task19->add_dep(task0);
      queue_->add_task(task19);


      std::vector<std::shared_ptr<Tensor>> tensor20 = {I5, t2};
      auto task20 = std::make_shared<Task20>(tensor20, pindex);
      task19->add_dep(task20);
      task20->add_dep(task0);
      queue_->add_task(task20);


      std::vector<IndexRange> I9_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I9(new Tensor(I9_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor21 = {I4, this->f1_, I9};
      auto task21 = std::make_shared<Task21>(tensor21, pindex);
      task18->add_dep(task21);
      task21->add_dep(task0);
      queue_->add_task(task21);


      std::vector<std::shared_ptr<Tensor>> tensor22 = {I9, t2};
      auto task22 = std::make_shared<Task22>(tensor22, pindex);
      task21->add_dep(task22);
      task22->add_dep(task0);
      queue_->add_task(task22);


      std::vector<IndexRange> I13_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I13(new Tensor(I13_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor23 = {I4, this->f1_, I13};
      auto task23 = std::make_shared<Task23>(tensor23, pindex);
      task18->add_dep(task23);
      task23->add_dep(task0);
      queue_->add_task(task23);


      std::vector<IndexRange> I14_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I14(new Tensor(I14_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor24 = {I13, t2, I14};
      auto task24 = std::make_shared<Task24>(tensor24, pindex);
      task23->add_dep(task24);
      task24->add_dep(task0);
      queue_->add_task(task24);


      std::vector<std::shared_ptr<Tensor>> tensor25 = {I14, Gamma2};
      auto task25 = std::make_shared<Task25>(tensor25, pindex);
      task24->add_dep(task25);
      task25->add_dep(task0);
      queue_->add_task(task25);

      task25->add_dep(task2);

      std::vector<IndexRange> I17_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I17(new Tensor(I17_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor26 = {I13, t2, I17};
      auto task26 = std::make_shared<Task26>(tensor26, pindex);
      task23->add_dep(task26);
      task26->add_dep(task0);
      queue_->add_task(task26);


      std::vector<std::shared_ptr<Tensor>> tensor27 = {I17, Gamma2};
      auto task27 = std::make_shared<Task27>(tensor27, pindex);
      task26->add_dep(task27);
      task27->add_dep(task0);
      queue_->add_task(task27);

      task27->add_dep(task2);

      std::vector<IndexRange> I18_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I18(new Tensor(I18_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor28 = {r, I18};
      auto task28 = std::make_shared<Task28>(tensor28, pindex);
      task28->add_dep(task0);
      queue_->add_task(task28);


      std::vector<IndexRange> I19_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I19(new Tensor(I19_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor29 = {I18, this->f1_, I19};
      auto task29 = std::make_shared<Task29>(tensor29, pindex);
      task28->add_dep(task29);
      task29->add_dep(task0);
      queue_->add_task(task29);


      std::vector<IndexRange> I20_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I20(new Tensor(I20_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor30 = {I19, t2, I20};
      auto task30 = std::make_shared<Task30>(tensor30, pindex);
      task29->add_dep(task30);
      task30->add_dep(task0);
      queue_->add_task(task30);


      std::vector<std::shared_ptr<Tensor>> tensor31 = {I20, Gamma2};
      auto task31 = std::make_shared<Task31>(tensor31, pindex);
      task30->add_dep(task31);
      task31->add_dep(task0);
      queue_->add_task(task31);

      task31->add_dep(task2);

      std::vector<IndexRange> I23_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I23(new Tensor(I23_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor32 = {I19, t2, I23};
      auto task32 = std::make_shared<Task32>(tensor32, pindex);
      task29->add_dep(task32);
      task32->add_dep(task0);
      queue_->add_task(task32);


      std::vector<std::shared_ptr<Tensor>> tensor33 = {I23, Gamma2};
      auto task33 = std::make_shared<Task33>(tensor33, pindex);
      task32->add_dep(task33);
      task33->add_dep(task0);
      queue_->add_task(task33);

      task33->add_dep(task2);

      std::vector<IndexRange> I25_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I25(new Tensor(I25_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor34 = {I18, t2, I25};
      auto task34 = std::make_shared<Task34>(tensor34, pindex);
      task28->add_dep(task34);
      task34->add_dep(task0);
      queue_->add_task(task34);


      std::vector<std::shared_ptr<Tensor>> tensor35 = {I25, Gamma6};
      auto task35 = std::make_shared<Task35>(tensor35, pindex);
      task34->add_dep(task35);
      task35->add_dep(task0);
      queue_->add_task(task35);

      task35->add_dep(task3);

      std::vector<IndexRange> I27_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I27(new Tensor(I27_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor36 = {I18, t2, I27};
      auto task36 = std::make_shared<Task36>(tensor36, pindex);
      task28->add_dep(task36);
      task36->add_dep(task0);
      queue_->add_task(task36);


      std::vector<std::shared_ptr<Tensor>> tensor37 = {I27, Gamma6};
      auto task37 = std::make_shared<Task37>(tensor37, pindex);
      task36->add_dep(task37);
      task37->add_dep(task0);
      queue_->add_task(task37);

      task37->add_dep(task3);

      std::vector<IndexRange> I29_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I29(new Tensor(I29_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor38 = {I18, this->f1_, I29};
      auto task38 = std::make_shared<Task38>(tensor38, pindex);
      task28->add_dep(task38);
      task38->add_dep(task0);
      queue_->add_task(task38);


      std::vector<IndexRange> I30_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I30(new Tensor(I30_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor39 = {I29, t2, I30};
      auto task39 = std::make_shared<Task39>(tensor39, pindex);
      task38->add_dep(task39);
      task39->add_dep(task0);
      queue_->add_task(task39);


      std::vector<std::shared_ptr<Tensor>> tensor40 = {I30, Gamma2};
      auto task40 = std::make_shared<Task40>(tensor40, pindex);
      task39->add_dep(task40);
      task40->add_dep(task0);
      queue_->add_task(task40);

      task40->add_dep(task2);

      std::vector<IndexRange> I33_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I33(new Tensor(I33_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor41 = {I29, t2, I33};
      auto task41 = std::make_shared<Task41>(tensor41, pindex);
      task38->add_dep(task41);
      task41->add_dep(task0);
      queue_->add_task(task41);


      std::vector<std::shared_ptr<Tensor>> tensor42 = {I33, Gamma2};
      auto task42 = std::make_shared<Task42>(tensor42, pindex);
      task41->add_dep(task42);
      task42->add_dep(task0);
      queue_->add_task(task42);

      task42->add_dep(task2);

      std::vector<IndexRange> I35_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I35(new Tensor(I35_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor43 = {I18, this->f1_, I35};
      auto task43 = std::make_shared<Task43>(tensor43, pindex);
      task28->add_dep(task43);
      task43->add_dep(task0);
      queue_->add_task(task43);


      std::vector<IndexRange> I36_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I36(new Tensor(I36_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor44 = {I35, t2, I36};
      auto task44 = std::make_shared<Task44>(tensor44, pindex);
      task43->add_dep(task44);
      task44->add_dep(task0);
      queue_->add_task(task44);


      std::vector<std::shared_ptr<Tensor>> tensor45 = {I36, Gamma2};
      auto task45 = std::make_shared<Task45>(tensor45, pindex);
      task44->add_dep(task45);
      task45->add_dep(task0);
      queue_->add_task(task45);

      task45->add_dep(task2);

      std::vector<IndexRange> I39_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I39(new Tensor(I39_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor46 = {I35, t2, I39};
      auto task46 = std::make_shared<Task46>(tensor46, pindex);
      task43->add_dep(task46);
      task46->add_dep(task0);
      queue_->add_task(task46);


      std::vector<std::shared_ptr<Tensor>> tensor47 = {I39, Gamma2};
      auto task47 = std::make_shared<Task47>(tensor47, pindex);
      task46->add_dep(task47);
      task47->add_dep(task0);
      queue_->add_task(task47);

      task47->add_dep(task2);

      std::vector<IndexRange> I41_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I41(new Tensor(I41_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor48 = {I18, this->f1_, I41};
      auto task48 = std::make_shared<Task48>(tensor48, pindex);
      task28->add_dep(task48);
      task48->add_dep(task0);
      queue_->add_task(task48);


      std::vector<IndexRange> I42_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I42(new Tensor(I42_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor49 = {I41, t2, I42};
      auto task49 = std::make_shared<Task49>(tensor49, pindex);
      task48->add_dep(task49);
      task49->add_dep(task0);
      queue_->add_task(task49);


      std::vector<std::shared_ptr<Tensor>> tensor50 = {I42, Gamma2};
      auto task50 = std::make_shared<Task50>(tensor50, pindex);
      task49->add_dep(task50);
      task50->add_dep(task0);
      queue_->add_task(task50);

      task50->add_dep(task2);

      std::vector<IndexRange> I45_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I45(new Tensor(I45_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor51 = {I41, t2, I45};
      auto task51 = std::make_shared<Task51>(tensor51, pindex);
      task48->add_dep(task51);
      task51->add_dep(task0);
      queue_->add_task(task51);


      std::vector<std::shared_ptr<Tensor>> tensor52 = {I45, Gamma2};
      auto task52 = std::make_shared<Task52>(tensor52, pindex);
      task51->add_dep(task52);
      task52->add_dep(task0);
      queue_->add_task(task52);

      task52->add_dep(task2);

      std::vector<IndexRange> I47_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I47(new Tensor(I47_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor53 = {I18, this->f1_, I47};
      auto task53 = std::make_shared<Task53>(tensor53, pindex);
      task28->add_dep(task53);
      task53->add_dep(task0);
      queue_->add_task(task53);


      std::vector<IndexRange> I48_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I48(new Tensor(I48_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor54 = {I47, t2, I48};
      auto task54 = std::make_shared<Task54>(tensor54, pindex);
      task53->add_dep(task54);
      task54->add_dep(task0);
      queue_->add_task(task54);


      std::vector<std::shared_ptr<Tensor>> tensor55 = {I48, Gamma14};
      auto task55 = std::make_shared<Task55>(tensor55, pindex);
      task54->add_dep(task55);
      task55->add_dep(task0);
      queue_->add_task(task55);

      task55->add_dep(task4);

      std::vector<IndexRange> I60_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I60(new Tensor(I60_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor56 = {I18, t2, I60};
      auto task56 = std::make_shared<Task56>(tensor56, pindex);
      task28->add_dep(task56);
      task56->add_dep(task0);
      queue_->add_task(task56);


      std::vector<std::shared_ptr<Tensor>> tensor57 = {I60, Gamma2};
      auto task57 = std::make_shared<Task57>(tensor57, pindex, this->e0_);
      task56->add_dep(task57);
      task57->add_dep(task0);
      queue_->add_task(task57);

      task57->add_dep(task2);

      std::vector<IndexRange> I62_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I62(new Tensor(I62_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor58 = {I18, t2, I62};
      auto task58 = std::make_shared<Task58>(tensor58, pindex);
      task28->add_dep(task58);
      task58->add_dep(task0);
      queue_->add_task(task58);


      std::vector<std::shared_ptr<Tensor>> tensor59 = {I62, Gamma2};
      auto task59 = std::make_shared<Task59>(tensor59, pindex, this->e0_);
      task58->add_dep(task59);
      task59->add_dep(task0);
      queue_->add_task(task59);

      task59->add_dep(task2);

      std::vector<IndexRange> I68_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I68(new Tensor(I68_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor60 = {I18, this->v2_, I68};
      auto task60 = std::make_shared<Task60>(tensor60, pindex);
      task28->add_dep(task60);
      task60->add_dep(task0);
      queue_->add_task(task60);


      std::vector<std::shared_ptr<Tensor>> tensor61 = {I68, Gamma2};
      auto task61 = std::make_shared<Task61>(tensor61, pindex);
      task60->add_dep(task61);
      task61->add_dep(task0);
      queue_->add_task(task61);

      task61->add_dep(task2);

      std::vector<IndexRange> I70_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I70(new Tensor(I70_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor62 = {I18, this->v2_, I70};
      auto task62 = std::make_shared<Task62>(tensor62, pindex);
      task28->add_dep(task62);
      task62->add_dep(task0);
      queue_->add_task(task62);


      std::vector<std::shared_ptr<Tensor>> tensor63 = {I70, Gamma2};
      auto task63 = std::make_shared<Task63>(tensor63, pindex);
      task62->add_dep(task63);
      task63->add_dep(task0);
      queue_->add_task(task63);

      task63->add_dep(task2);

      std::vector<IndexRange> I49_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I49(new Tensor(I49_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor64 = {r, I49};
      auto task64 = std::make_shared<Task64>(tensor64, pindex);
      task64->add_dep(task0);
      queue_->add_task(task64);


      std::vector<IndexRange> I50_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I50(new Tensor(I50_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor65 = {I49, this->f1_, I50};
      auto task65 = std::make_shared<Task65>(tensor65, pindex);
      task64->add_dep(task65);
      task65->add_dep(task0);
      queue_->add_task(task65);


      std::vector<IndexRange> I51_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I51(new Tensor(I51_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor66 = {I50, t2, I51};
      auto task66 = std::make_shared<Task66>(tensor66, pindex);
      task65->add_dep(task66);
      task66->add_dep(task0);
      queue_->add_task(task66);


      std::vector<std::shared_ptr<Tensor>> tensor67 = {I51, Gamma14};
      auto task67 = std::make_shared<Task67>(tensor67, pindex);
      task66->add_dep(task67);
      task67->add_dep(task0);
      queue_->add_task(task67);

      task67->add_dep(task4);

      std::vector<IndexRange> I55_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I55(new Tensor(I55_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor68 = {I49, this->f1_, I55};
      auto task68 = std::make_shared<Task68>(tensor68, pindex);
      task64->add_dep(task68);
      task68->add_dep(task0);
      queue_->add_task(task68);


      std::vector<IndexRange> I56_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I56(new Tensor(I56_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor69 = {I55, t2, I56};
      auto task69 = std::make_shared<Task69>(tensor69, pindex);
      task68->add_dep(task69);
      task69->add_dep(task0);
      queue_->add_task(task69);


      std::vector<std::shared_ptr<Tensor>> tensor70 = {I56, Gamma14};
      auto task70 = std::make_shared<Task70>(tensor70, pindex);
      task69->add_dep(task70);
      task70->add_dep(task0);
      queue_->add_task(task70);

      task70->add_dep(task4);

      std::vector<IndexRange> I52_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I52(new Tensor(I52_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor71 = {r, I52};
      auto task71 = std::make_shared<Task71>(tensor71, pindex);
      task71->add_dep(task0);
      queue_->add_task(task71);


      std::vector<IndexRange> I53_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I53(new Tensor(I53_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor72 = {I52, t2, I53};
      auto task72 = std::make_shared<Task72>(tensor72, pindex);
      task71->add_dep(task72);
      task72->add_dep(task0);
      queue_->add_task(task72);


      std::vector<std::shared_ptr<Tensor>> tensor73 = {I53, Gamma16};
      auto task73 = std::make_shared<Task73>(tensor73, pindex);
      task72->add_dep(task73);
      task73->add_dep(task0);
      queue_->add_task(task73);

      task73->add_dep(task5);

      std::vector<IndexRange> I64_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I64(new Tensor(I64_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor74 = {I52, t2, I64};
      auto task74 = std::make_shared<Task74>(tensor74, pindex);
      task71->add_dep(task74);
      task74->add_dep(task0);
      queue_->add_task(task74);


      std::vector<std::shared_ptr<Tensor>> tensor75 = {I64, Gamma14};
      auto task75 = std::make_shared<Task75>(tensor75, pindex, this->e0_);
      task74->add_dep(task75);
      task75->add_dep(task0);
      queue_->add_task(task75);

      task75->add_dep(task4);

      std::vector<IndexRange> I72_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I72(new Tensor(I72_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor76 = {I52, this->v2_, I72};
      auto task76 = std::make_shared<Task76>(tensor76, pindex);
      task71->add_dep(task76);
      task76->add_dep(task0);
      queue_->add_task(task76);


      std::vector<std::shared_ptr<Tensor>> tensor77 = {I72, Gamma14};
      auto task77 = std::make_shared<Task77>(tensor77, pindex);
      task76->add_dep(task77);
      task77->add_dep(task0);
      queue_->add_task(task77);

      task77->add_dep(task4);

      auto energy_ = std::make_shared<Queue>();
      std::vector<IndexRange> I73_index;
      std::shared_ptr<Tensor> I73(new Tensor(I73_index, false));
      std::vector<IndexRange> I74_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I74(new Tensor(I74_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor78 = {I73, t2, I74};
      auto task78 = std::make_shared<Task78>(tensor78, pindex);
      energy_->add_task(task78);


      std::vector<std::shared_ptr<Tensor>> tensor79 = {I74, t2, this->v2_};
      auto task79 = std::make_shared<Task79>(tensor79, pindex, this->e0_);
      task78->add_dep(task79);
      energy_->add_task(task79);


      std::vector<IndexRange> I75_index;
      std::shared_ptr<Tensor> I75(new Tensor(I75_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor80 = {I74, t2, I75};
      auto task80 = std::make_shared<Task80>(tensor80, pindex);
      task78->add_dep(task80);
      energy_->add_task(task80);


      std::vector<std::shared_ptr<Tensor>> tensor81 = {I75, Gamma0};
      auto task81 = std::make_shared<Task81>(tensor81, pindex);
      task80->add_dep(task81);
      energy_->add_task(task81);

      task81->add_dep(task1);

      std::vector<IndexRange> I78_index;
      std::shared_ptr<Tensor> I78(new Tensor(I78_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor82 = {I74, t2, I78};
      auto task82 = std::make_shared<Task82>(tensor82, pindex);
      task78->add_dep(task82);
      energy_->add_task(task82);


      std::vector<std::shared_ptr<Tensor>> tensor83 = {I78, Gamma0};
      auto task83 = std::make_shared<Task83>(tensor83, pindex);
      task82->add_dep(task83);
      energy_->add_task(task83);

      task83->add_dep(task1);

      std::vector<IndexRange> I81_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I81(new Tensor(I81_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor84 = {I74, this->f1_, I81};
      auto task84 = std::make_shared<Task84>(tensor84, pindex);
      task78->add_dep(task84);
      energy_->add_task(task84);


      std::vector<std::shared_ptr<Tensor>> tensor85 = {I81, t2};
      auto task85 = std::make_shared<Task85>(tensor85, pindex);
      task84->add_dep(task85);
      energy_->add_task(task85);


      std::vector<IndexRange> I87_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I87(new Tensor(I87_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor86 = {I74, this->f1_, I87};
      auto task86 = std::make_shared<Task86>(tensor86, pindex);
      task78->add_dep(task86);
      energy_->add_task(task86);


      std::vector<std::shared_ptr<Tensor>> tensor87 = {I87, t2};
      auto task87 = std::make_shared<Task87>(tensor87, pindex);
      task86->add_dep(task87);
      energy_->add_task(task87);


      std::vector<IndexRange> I93_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I93(new Tensor(I93_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor88 = {I74, this->f1_, I93};
      auto task88 = std::make_shared<Task88>(tensor88, pindex);
      task78->add_dep(task88);
      energy_->add_task(task88);


      std::vector<IndexRange> I94_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I94(new Tensor(I94_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor89 = {I93, t2, I94};
      auto task89 = std::make_shared<Task89>(tensor89, pindex);
      task88->add_dep(task89);
      energy_->add_task(task89);


      std::vector<std::shared_ptr<Tensor>> tensor90 = {I94, Gamma2};
      auto task90 = std::make_shared<Task90>(tensor90, pindex);
      task89->add_dep(task90);
      energy_->add_task(task90);

      task90->add_dep(task2);

      std::vector<IndexRange> I98_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I98(new Tensor(I98_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor91 = {I93, t2, I98};
      auto task91 = std::make_shared<Task91>(tensor91, pindex);
      task88->add_dep(task91);
      energy_->add_task(task91);


      std::vector<std::shared_ptr<Tensor>> tensor92 = {I98, Gamma2};
      auto task92 = std::make_shared<Task92>(tensor92, pindex);
      task91->add_dep(task92);
      energy_->add_task(task92);

      task92->add_dep(task2);

      std::vector<IndexRange> I100_index = {this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I100(new Tensor(I100_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor93 = {I73, t2, I100};
      auto task93 = std::make_shared<Task93>(tensor93, pindex);
      task78->add_dep(task93);
      energy_->add_task(task93);


      std::vector<IndexRange> I101_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I101(new Tensor(I101_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor94 = {I100, this->f1_, I101};
      auto task94 = std::make_shared<Task94>(tensor94, pindex);
      task93->add_dep(task94);
      energy_->add_task(task94);


      std::vector<IndexRange> I102_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I102(new Tensor(I102_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor95 = {I101, t2, I102};
      auto task95 = std::make_shared<Task95>(tensor95, pindex);
      task94->add_dep(task95);
      energy_->add_task(task95);


      std::vector<std::shared_ptr<Tensor>> tensor96 = {I102, Gamma2};
      auto task96 = std::make_shared<Task96>(tensor96, pindex);
      task95->add_dep(task96);
      energy_->add_task(task96);

      task96->add_dep(task2);

      std::vector<IndexRange> I106_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I106(new Tensor(I106_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor97 = {I101, t2, I106};
      auto task97 = std::make_shared<Task97>(tensor97, pindex);
      task94->add_dep(task97);
      energy_->add_task(task97);


      std::vector<std::shared_ptr<Tensor>> tensor98 = {I106, Gamma2};
      auto task98 = std::make_shared<Task98>(tensor98, pindex);
      task97->add_dep(task98);
      energy_->add_task(task98);

      task98->add_dep(task2);

      std::vector<IndexRange> I109_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I109(new Tensor(I109_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor99 = {I100, t2, I109};
      auto task99 = std::make_shared<Task99>(tensor99, pindex);
      task93->add_dep(task99);
      energy_->add_task(task99);


      std::vector<std::shared_ptr<Tensor>> tensor100 = {I109, Gamma6};
      auto task100 = std::make_shared<Task100>(tensor100, pindex);
      task99->add_dep(task100);
      energy_->add_task(task100);

      task100->add_dep(task3);

      std::vector<IndexRange> I112_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I112(new Tensor(I112_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor101 = {I100, t2, I112};
      auto task101 = std::make_shared<Task101>(tensor101, pindex);
      task93->add_dep(task101);
      energy_->add_task(task101);


      std::vector<std::shared_ptr<Tensor>> tensor102 = {I112, Gamma6};
      auto task102 = std::make_shared<Task102>(tensor102, pindex);
      task101->add_dep(task102);
      energy_->add_task(task102);

      task102->add_dep(task3);

      std::vector<IndexRange> I115_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I115(new Tensor(I115_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor103 = {I100, this->f1_, I115};
      auto task103 = std::make_shared<Task103>(tensor103, pindex);
      task93->add_dep(task103);
      energy_->add_task(task103);


      std::vector<IndexRange> I116_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I116(new Tensor(I116_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor104 = {I115, t2, I116};
      auto task104 = std::make_shared<Task104>(tensor104, pindex);
      task103->add_dep(task104);
      energy_->add_task(task104);


      std::vector<std::shared_ptr<Tensor>> tensor105 = {I116, Gamma2};
      auto task105 = std::make_shared<Task105>(tensor105, pindex);
      task104->add_dep(task105);
      energy_->add_task(task105);

      task105->add_dep(task2);

      std::vector<IndexRange> I120_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I120(new Tensor(I120_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor106 = {I115, t2, I120};
      auto task106 = std::make_shared<Task106>(tensor106, pindex);
      task103->add_dep(task106);
      energy_->add_task(task106);


      std::vector<std::shared_ptr<Tensor>> tensor107 = {I120, Gamma2};
      auto task107 = std::make_shared<Task107>(tensor107, pindex);
      task106->add_dep(task107);
      energy_->add_task(task107);

      task107->add_dep(task2);

      std::vector<IndexRange> I123_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I123(new Tensor(I123_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor108 = {I100, this->f1_, I123};
      auto task108 = std::make_shared<Task108>(tensor108, pindex);
      task93->add_dep(task108);
      energy_->add_task(task108);


      std::vector<IndexRange> I124_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I124(new Tensor(I124_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor109 = {I123, t2, I124};
      auto task109 = std::make_shared<Task109>(tensor109, pindex);
      task108->add_dep(task109);
      energy_->add_task(task109);


      std::vector<std::shared_ptr<Tensor>> tensor110 = {I124, Gamma2};
      auto task110 = std::make_shared<Task110>(tensor110, pindex);
      task109->add_dep(task110);
      energy_->add_task(task110);

      task110->add_dep(task2);

      std::vector<IndexRange> I128_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I128(new Tensor(I128_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor111 = {I123, t2, I128};
      auto task111 = std::make_shared<Task111>(tensor111, pindex);
      task108->add_dep(task111);
      energy_->add_task(task111);


      std::vector<std::shared_ptr<Tensor>> tensor112 = {I128, Gamma2};
      auto task112 = std::make_shared<Task112>(tensor112, pindex);
      task111->add_dep(task112);
      energy_->add_task(task112);

      task112->add_dep(task2);

      std::vector<IndexRange> I131_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I131(new Tensor(I131_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor113 = {I100, this->f1_, I131};
      auto task113 = std::make_shared<Task113>(tensor113, pindex);
      task93->add_dep(task113);
      energy_->add_task(task113);


      std::vector<IndexRange> I132_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I132(new Tensor(I132_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor114 = {I131, t2, I132};
      auto task114 = std::make_shared<Task114>(tensor114, pindex);
      task113->add_dep(task114);
      energy_->add_task(task114);


      std::vector<std::shared_ptr<Tensor>> tensor115 = {I132, Gamma2};
      auto task115 = std::make_shared<Task115>(tensor115, pindex);
      task114->add_dep(task115);
      energy_->add_task(task115);

      task115->add_dep(task2);

      std::vector<IndexRange> I136_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I136(new Tensor(I136_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor116 = {I131, t2, I136};
      auto task116 = std::make_shared<Task116>(tensor116, pindex);
      task113->add_dep(task116);
      energy_->add_task(task116);


      std::vector<std::shared_ptr<Tensor>> tensor117 = {I136, Gamma2};
      auto task117 = std::make_shared<Task117>(tensor117, pindex);
      task116->add_dep(task117);
      energy_->add_task(task117);

      task117->add_dep(task2);

      std::vector<IndexRange> I139_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I139(new Tensor(I139_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor118 = {I100, this->f1_, I139};
      auto task118 = std::make_shared<Task118>(tensor118, pindex);
      task93->add_dep(task118);
      energy_->add_task(task118);


      std::vector<IndexRange> I140_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I140(new Tensor(I140_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor119 = {I139, t2, I140};
      auto task119 = std::make_shared<Task119>(tensor119, pindex);
      task118->add_dep(task119);
      energy_->add_task(task119);


      std::vector<std::shared_ptr<Tensor>> tensor120 = {I140, Gamma14};
      auto task120 = std::make_shared<Task120>(tensor120, pindex);
      task119->add_dep(task120);
      energy_->add_task(task120);

      task120->add_dep(task4);

      std::vector<IndexRange> I158_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I158(new Tensor(I158_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor121 = {I100, t2, I158};
      auto task121 = std::make_shared<Task121>(tensor121, pindex);
      task93->add_dep(task121);
      energy_->add_task(task121);


      std::vector<std::shared_ptr<Tensor>> tensor122 = {I158, Gamma2};
      auto task122 = std::make_shared<Task122>(tensor122, pindex, this->e0_);
      task121->add_dep(task122);
      energy_->add_task(task122);

      task122->add_dep(task2);

      std::vector<IndexRange> I161_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I161(new Tensor(I161_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor123 = {I100, t2, I161};
      auto task123 = std::make_shared<Task123>(tensor123, pindex);
      task93->add_dep(task123);
      energy_->add_task(task123);


      std::vector<std::shared_ptr<Tensor>> tensor124 = {I161, Gamma2};
      auto task124 = std::make_shared<Task124>(tensor124, pindex, this->e0_);
      task123->add_dep(task124);
      energy_->add_task(task124);

      task124->add_dep(task2);

      std::vector<IndexRange> I171_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I171(new Tensor(I171_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor125 = {I100, this->v2_, I171};
      auto task125 = std::make_shared<Task125>(tensor125, pindex);
      task93->add_dep(task125);
      energy_->add_task(task125);


      std::vector<std::shared_ptr<Tensor>> tensor126 = {I171, Gamma2};
      auto task126 = std::make_shared<Task126>(tensor126, pindex);
      task125->add_dep(task126);
      energy_->add_task(task126);

      task126->add_dep(task2);

      std::vector<IndexRange> I174_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I174(new Tensor(I174_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor127 = {I100, this->v2_, I174};
      auto task127 = std::make_shared<Task127>(tensor127, pindex);
      task93->add_dep(task127);
      energy_->add_task(task127);


      std::vector<std::shared_ptr<Tensor>> tensor128 = {I174, Gamma2};
      auto task128 = std::make_shared<Task128>(tensor128, pindex);
      task127->add_dep(task128);
      energy_->add_task(task128);

      task128->add_dep(task2);

      std::vector<IndexRange> I142_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I142(new Tensor(I142_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor129 = {I73, t2, I142};
      auto task129 = std::make_shared<Task129>(tensor129, pindex);
      task78->add_dep(task129);
      energy_->add_task(task129);


      std::vector<IndexRange> I143_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I143(new Tensor(I143_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor130 = {I142, this->f1_, I143};
      auto task130 = std::make_shared<Task130>(tensor130, pindex);
      task129->add_dep(task130);
      energy_->add_task(task130);


      std::vector<IndexRange> I144_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I144(new Tensor(I144_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor131 = {I143, t2, I144};
      auto task131 = std::make_shared<Task131>(tensor131, pindex);
      task130->add_dep(task131);
      energy_->add_task(task131);


      std::vector<std::shared_ptr<Tensor>> tensor132 = {I144, Gamma14};
      auto task132 = std::make_shared<Task132>(tensor132, pindex);
      task131->add_dep(task132);
      energy_->add_task(task132);

      task132->add_dep(task4);

      std::vector<IndexRange> I147_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I147(new Tensor(I147_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor133 = {I142, t2, I147};
      auto task133 = std::make_shared<Task133>(tensor133, pindex);
      task129->add_dep(task133);
      energy_->add_task(task133);


      std::vector<std::shared_ptr<Tensor>> tensor134 = {I147, Gamma16};
      auto task134 = std::make_shared<Task134>(tensor134, pindex);
      task133->add_dep(task134);
      energy_->add_task(task134);

      task134->add_dep(task5);

      std::vector<IndexRange> I150_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I150(new Tensor(I150_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor135 = {I142, this->f1_, I150};
      auto task135 = std::make_shared<Task135>(tensor135, pindex);
      task129->add_dep(task135);
      energy_->add_task(task135);


      std::vector<IndexRange> I151_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I151(new Tensor(I151_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor136 = {I150, t2, I151};
      auto task136 = std::make_shared<Task136>(tensor136, pindex);
      task135->add_dep(task136);
      energy_->add_task(task136);


      std::vector<std::shared_ptr<Tensor>> tensor137 = {I151, Gamma14};
      auto task137 = std::make_shared<Task137>(tensor137, pindex);
      task136->add_dep(task137);
      energy_->add_task(task137);

      task137->add_dep(task4);

      std::vector<IndexRange> I164_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I164(new Tensor(I164_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor138 = {I142, t2, I164};
      auto task138 = std::make_shared<Task138>(tensor138, pindex);
      task129->add_dep(task138);
      energy_->add_task(task138);


      std::vector<std::shared_ptr<Tensor>> tensor139 = {I164, Gamma14};
      auto task139 = std::make_shared<Task139>(tensor139, pindex, this->e0_);
      task138->add_dep(task139);
      energy_->add_task(task139);

      task139->add_dep(task4);

      std::vector<IndexRange> I177_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I177(new Tensor(I177_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor140 = {I142, this->v2_, I177};
      auto task140 = std::make_shared<Task140>(tensor140, pindex);
      task129->add_dep(task140);
      energy_->add_task(task140);


      std::vector<std::shared_ptr<Tensor>> tensor141 = {I177, Gamma14};
      auto task141 = std::make_shared<Task141>(tensor141, pindex);
      task140->add_dep(task141);
      energy_->add_task(task141);

      task141->add_dep(task4);

      auto correction_ = std::make_shared<Queue>();
      std::vector<IndexRange> I178_index;
      std::shared_ptr<Tensor> I178(new Tensor(I178_index, false));
      std::vector<IndexRange> I179_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I179(new Tensor(I179_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor142 = {I178, t2, I179};
      auto task142 = std::make_shared<Task142>(tensor142, pindex);
      correction_->add_task(task142);


      std::vector<std::shared_ptr<Tensor>> tensor143 = {I179, t2};
      auto task143 = std::make_shared<Task143>(tensor143, pindex);
      task142->add_dep(task143);
      correction_->add_task(task143);


      std::vector<IndexRange> I183_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I183(new Tensor(I183_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor144 = {I178, t2, I183};
      auto task144 = std::make_shared<Task144>(tensor144, pindex);
      task142->add_dep(task144);
      correction_->add_task(task144);


      std::vector<IndexRange> I184_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I184(new Tensor(I184_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor145 = {I183, t2, I184};
      auto task145 = std::make_shared<Task145>(tensor145, pindex);
      task144->add_dep(task145);
      correction_->add_task(task145);


      std::vector<std::shared_ptr<Tensor>> tensor146 = {I184, Gamma2};
      auto task146 = std::make_shared<Task146>(tensor146, pindex);
      task145->add_dep(task146);
      correction_->add_task(task146);

      task146->add_dep(task2);

      std::vector<IndexRange> I187_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I187(new Tensor(I187_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor147 = {I183, t2, I187};
      auto task147 = std::make_shared<Task147>(tensor147, pindex);
      task144->add_dep(task147);
      correction_->add_task(task147);


      std::vector<std::shared_ptr<Tensor>> tensor148 = {I187, Gamma2};
      auto task148 = std::make_shared<Task148>(tensor148, pindex);
      task147->add_dep(task148);
      correction_->add_task(task148);

      task148->add_dep(task2);

      std::vector<IndexRange> I189_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I189(new Tensor(I189_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor149 = {I178, t2, I189};
      auto task149 = std::make_shared<Task149>(tensor149, pindex);
      task142->add_dep(task149);
      correction_->add_task(task149);


      std::vector<IndexRange> I190_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I190(new Tensor(I190_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor150 = {I189, t2, I190};
      auto task150 = std::make_shared<Task150>(tensor150, pindex);
      task149->add_dep(task150);
      correction_->add_task(task150);


      std::vector<std::shared_ptr<Tensor>> tensor151 = {I190, Gamma14};
      auto task151 = std::make_shared<Task151>(tensor151, pindex);
      task150->add_dep(task151);
      correction_->add_task(task151);

      task151->add_dep(task4);

      auto density_ = std::make_shared<Queue>();
      std::vector<std::shared_ptr<Tensor>> tensor152 = {den2};
      auto task152 = std::make_shared<Task152>(tensor152);
      density_->add_task(task152);

      std::vector<IndexRange> I191_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I191(new Tensor(I191_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor153 = {den2, I191};
      auto task153 = std::make_shared<Task153>(tensor153, pindex);
      task153->add_dep(task152);
      density_->add_task(task153);


      std::vector<IndexRange> I192_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I192(new Tensor(I192_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor154 = {I191, t2, I192};
      auto task154 = std::make_shared<Task154>(tensor154, pindex);
      task153->add_dep(task154);
      task154->add_dep(task152);
      density_->add_task(task154);


      std::vector<IndexRange> I193_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I193(new Tensor(I193_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor155 = {I192, t2, I193};
      auto task155 = std::make_shared<Task155>(tensor155, pindex);
      task154->add_dep(task155);
      task155->add_dep(task152);
      density_->add_task(task155);


      std::vector<std::shared_ptr<Tensor>> tensor156 = {I193, Gamma2};
      auto task156 = std::make_shared<Task156>(tensor156, pindex);
      task155->add_dep(task156);
      task156->add_dep(task152);
      density_->add_task(task156);

      task156->add_dep(task2);

      std::vector<IndexRange> I196_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I196(new Tensor(I196_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor157 = {I192, t2, I196};
      auto task157 = std::make_shared<Task157>(tensor157, pindex);
      task154->add_dep(task157);
      task157->add_dep(task152);
      density_->add_task(task157);


      std::vector<std::shared_ptr<Tensor>> tensor158 = {I196, Gamma2};
      auto task158 = std::make_shared<Task158>(tensor158, pindex);
      task157->add_dep(task158);
      task158->add_dep(task152);
      density_->add_task(task158);

      task158->add_dep(task2);

      std::vector<IndexRange> I197_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor> I197(new Tensor(I197_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor159 = {den2, I197};
      auto task159 = std::make_shared<Task159>(tensor159, pindex);
      task159->add_dep(task152);
      density_->add_task(task159);


      std::vector<IndexRange> I198_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I198(new Tensor(I198_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor160 = {I197, t2, I198};
      auto task160 = std::make_shared<Task160>(tensor160, pindex);
      task159->add_dep(task160);
      task160->add_dep(task152);
      density_->add_task(task160);


      std::vector<std::shared_ptr<Tensor>> tensor161 = {I198, t2};
      auto task161 = std::make_shared<Task161>(tensor161, pindex);
      task160->add_dep(task161);
      task161->add_dep(task152);
      density_->add_task(task161);


      std::vector<IndexRange> I201_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor> I201(new Tensor(I201_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor162 = {den2, I201};
      auto task162 = std::make_shared<Task162>(tensor162, pindex);
      task162->add_dep(task152);
      density_->add_task(task162);


      std::vector<IndexRange> I202_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I202(new Tensor(I202_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor163 = {I201, t2, I202};
      auto task163 = std::make_shared<Task163>(tensor163, pindex);
      task162->add_dep(task163);
      task163->add_dep(task152);
      density_->add_task(task163);


      std::vector<std::shared_ptr<Tensor>> tensor164 = {I202, t2};
      auto task164 = std::make_shared<Task164>(tensor164, pindex);
      task163->add_dep(task164);
      task164->add_dep(task152);
      density_->add_task(task164);


      std::vector<IndexRange> I205_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor> I205(new Tensor(I205_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor165 = {den2, I205};
      auto task165 = std::make_shared<Task165>(tensor165, pindex);
      task165->add_dep(task152);
      density_->add_task(task165);


      std::vector<IndexRange> I206_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I206(new Tensor(I206_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor166 = {I205, t2, I206};
      auto task166 = std::make_shared<Task166>(tensor166, pindex);
      task165->add_dep(task166);
      task166->add_dep(task152);
      density_->add_task(task166);


      std::vector<IndexRange> I207_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I207(new Tensor(I207_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor167 = {I206, t2, I207};
      auto task167 = std::make_shared<Task167>(tensor167, pindex);
      task166->add_dep(task167);
      task167->add_dep(task152);
      density_->add_task(task167);


      std::vector<std::shared_ptr<Tensor>> tensor168 = {I207, Gamma2};
      auto task168 = std::make_shared<Task168>(tensor168, pindex);
      task167->add_dep(task168);
      task168->add_dep(task152);
      density_->add_task(task168);

      task168->add_dep(task2);

      std::vector<IndexRange> I210_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I210(new Tensor(I210_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor169 = {I206, t2, I210};
      auto task169 = std::make_shared<Task169>(tensor169, pindex);
      task166->add_dep(task169);
      task169->add_dep(task152);
      density_->add_task(task169);


      std::vector<std::shared_ptr<Tensor>> tensor170 = {I210, Gamma2};
      auto task170 = std::make_shared<Task170>(tensor170, pindex);
      task169->add_dep(task170);
      task170->add_dep(task152);
      density_->add_task(task170);

      task170->add_dep(task2);

      std::vector<IndexRange> I211_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor> I211(new Tensor(I211_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor171 = {den2, I211};
      auto task171 = std::make_shared<Task171>(tensor171, pindex);
      task171->add_dep(task152);
      density_->add_task(task171);


      std::vector<IndexRange> I212_index = {this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I212(new Tensor(I212_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor172 = {I211, t2, I212};
      auto task172 = std::make_shared<Task172>(tensor172, pindex);
      task171->add_dep(task172);
      task172->add_dep(task152);
      density_->add_task(task172);


      std::vector<IndexRange> I213_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I213(new Tensor(I213_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor173 = {I212, t2, I213};
      auto task173 = std::make_shared<Task173>(tensor173, pindex);
      task172->add_dep(task173);
      task173->add_dep(task152);
      density_->add_task(task173);


      std::vector<std::shared_ptr<Tensor>> tensor174 = {I213, Gamma2};
      auto task174 = std::make_shared<Task174>(tensor174, pindex);
      task173->add_dep(task174);
      task174->add_dep(task152);
      density_->add_task(task174);

      task174->add_dep(task2);

      std::vector<IndexRange> I216_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I216(new Tensor(I216_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor175 = {I212, t2, I216};
      auto task175 = std::make_shared<Task175>(tensor175, pindex);
      task172->add_dep(task175);
      task175->add_dep(task152);
      density_->add_task(task175);


      std::vector<std::shared_ptr<Tensor>> tensor176 = {I216, Gamma2};
      auto task176 = std::make_shared<Task176>(tensor176, pindex);
      task175->add_dep(task176);
      task176->add_dep(task152);
      density_->add_task(task176);

      task176->add_dep(task2);

      std::vector<IndexRange> I217_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I217(new Tensor(I217_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor177 = {den2, I217};
      auto task177 = std::make_shared<Task177>(tensor177, pindex);
      task177->add_dep(task152);
      density_->add_task(task177);


      std::vector<IndexRange> I218_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I218(new Tensor(I218_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor178 = {I217, t2, I218};
      auto task178 = std::make_shared<Task178>(tensor178, pindex);
      task177->add_dep(task178);
      task178->add_dep(task152);
      density_->add_task(task178);


      std::vector<IndexRange> I219_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I219(new Tensor(I219_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor179 = {I218, t2, I219};
      auto task179 = std::make_shared<Task179>(tensor179, pindex);
      task178->add_dep(task179);
      task179->add_dep(task152);
      density_->add_task(task179);


      std::vector<std::shared_ptr<Tensor>> tensor180 = {I219, Gamma14};
      auto task180 = std::make_shared<Task180>(tensor180, pindex);
      task179->add_dep(task180);
      task180->add_dep(task152);
      density_->add_task(task180);

      task180->add_dep(task4);

      std::vector<IndexRange> I222_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I222(new Tensor(I222_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor181 = {I218, t2, I222};
      auto task181 = std::make_shared<Task181>(tensor181, pindex);
      task178->add_dep(task181);
      task181->add_dep(task152);
      density_->add_task(task181);


      std::vector<std::shared_ptr<Tensor>> tensor182 = {I222, Gamma14};
      auto task182 = std::make_shared<Task182>(tensor182, pindex);
      task181->add_dep(task182);
      task182->add_dep(task152);
      density_->add_task(task182);

      task182->add_dep(task4);

      std::vector<IndexRange> I223_index = {this->closed_, this->closed_};
      std::shared_ptr<Tensor> I223(new Tensor(I223_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor183 = {den2, I223};
      auto task183 = std::make_shared<Task183>(tensor183, pindex);
      task183->add_dep(task152);
      density_->add_task(task183);


      std::vector<IndexRange> I224_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I224(new Tensor(I224_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor184 = {I223, t2, I224};
      auto task184 = std::make_shared<Task184>(tensor184, pindex);
      task183->add_dep(task184);
      task184->add_dep(task152);
      density_->add_task(task184);


      std::vector<IndexRange> I225_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I225(new Tensor(I225_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor185 = {I224, t2, I225};
      auto task185 = std::make_shared<Task185>(tensor185, pindex);
      task184->add_dep(task185);
      task185->add_dep(task152);
      density_->add_task(task185);


      std::vector<std::shared_ptr<Tensor>> tensor186 = {I225, Gamma2};
      auto task186 = std::make_shared<Task186>(tensor186, pindex);
      task185->add_dep(task186);
      task186->add_dep(task152);
      density_->add_task(task186);

      task186->add_dep(task2);

      std::vector<IndexRange> I228_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I228(new Tensor(I228_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor187 = {I224, t2, I228};
      auto task187 = std::make_shared<Task187>(tensor187, pindex);
      task184->add_dep(task187);
      task187->add_dep(task152);
      density_->add_task(task187);


      std::vector<std::shared_ptr<Tensor>> tensor188 = {I228, Gamma2};
      auto task188 = std::make_shared<Task188>(tensor188, pindex);
      task187->add_dep(task188);
      task188->add_dep(task152);
      density_->add_task(task188);

      task188->add_dep(task2);

      std::vector<IndexRange> I229_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor> I229(new Tensor(I229_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor189 = {den2, I229};
      auto task189 = std::make_shared<Task189>(tensor189, pindex);
      task189->add_dep(task152);
      density_->add_task(task189);


      std::vector<IndexRange> I230_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I230(new Tensor(I230_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor190 = {I229, t2, I230};
      auto task190 = std::make_shared<Task190>(tensor190, pindex);
      task189->add_dep(task190);
      task190->add_dep(task152);
      density_->add_task(task190);


      std::vector<IndexRange> I231_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I231(new Tensor(I231_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor191 = {I230, t2, I231};
      auto task191 = std::make_shared<Task191>(tensor191, pindex);
      task190->add_dep(task191);
      task191->add_dep(task152);
      density_->add_task(task191);


      std::vector<std::shared_ptr<Tensor>> tensor192 = {I231, Gamma2};
      auto task192 = std::make_shared<Task192>(tensor192, pindex);
      task191->add_dep(task192);
      task192->add_dep(task152);
      density_->add_task(task192);

      task192->add_dep(task2);

      std::vector<IndexRange> I234_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I234(new Tensor(I234_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor193 = {I230, t2, I234};
      auto task193 = std::make_shared<Task193>(tensor193, pindex);
      task190->add_dep(task193);
      task193->add_dep(task152);
      density_->add_task(task193);


      std::vector<std::shared_ptr<Tensor>> tensor194 = {I234, Gamma2};
      auto task194 = std::make_shared<Task194>(tensor194, pindex);
      task193->add_dep(task194);
      task194->add_dep(task152);
      density_->add_task(task194);

      task194->add_dep(task2);

      std::vector<IndexRange> I235_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor> I235(new Tensor(I235_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor195 = {den2, I235};
      auto task195 = std::make_shared<Task195>(tensor195, pindex);
      task195->add_dep(task152);
      density_->add_task(task195);


      std::vector<IndexRange> I236_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I236(new Tensor(I236_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor196 = {I235, t2, I236};
      auto task196 = std::make_shared<Task196>(tensor196, pindex);
      task195->add_dep(task196);
      task196->add_dep(task152);
      density_->add_task(task196);


      std::vector<IndexRange> I237_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I237(new Tensor(I237_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor197 = {I236, t2, I237};
      auto task197 = std::make_shared<Task197>(tensor197, pindex);
      task196->add_dep(task197);
      task197->add_dep(task152);
      density_->add_task(task197);


      std::vector<std::shared_ptr<Tensor>> tensor198 = {I237, Gamma2};
      auto task198 = std::make_shared<Task198>(tensor198, pindex);
      task197->add_dep(task198);
      task198->add_dep(task152);
      density_->add_task(task198);

      task198->add_dep(task2);

      std::vector<IndexRange> I240_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I240(new Tensor(I240_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor199 = {I236, t2, I240};
      auto task199 = std::make_shared<Task199>(tensor199, pindex);
      task196->add_dep(task199);
      task199->add_dep(task152);
      density_->add_task(task199);


      std::vector<std::shared_ptr<Tensor>> tensor200 = {I240, Gamma2};
      auto task200 = std::make_shared<Task200>(tensor200, pindex);
      task199->add_dep(task200);
      task200->add_dep(task152);
      density_->add_task(task200);

      task200->add_dep(task2);

      std::vector<IndexRange> I241_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor> I241(new Tensor(I241_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor201 = {den2, I241};
      auto task201 = std::make_shared<Task201>(tensor201, pindex);
      task201->add_dep(task152);
      density_->add_task(task201);


      std::vector<IndexRange> I242_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I242(new Tensor(I242_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor202 = {I241, t2, I242};
      auto task202 = std::make_shared<Task202>(tensor202, pindex);
      task201->add_dep(task202);
      task202->add_dep(task152);
      density_->add_task(task202);


      std::vector<IndexRange> I243_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I243(new Tensor(I243_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor203 = {I242, t2, I243};
      auto task203 = std::make_shared<Task203>(tensor203, pindex);
      task202->add_dep(task203);
      task203->add_dep(task152);
      density_->add_task(task203);


      std::vector<std::shared_ptr<Tensor>> tensor204 = {I243, Gamma14};
      auto task204 = std::make_shared<Task204>(tensor204, pindex);
      task203->add_dep(task204);
      task204->add_dep(task152);
      density_->add_task(task204);

      task204->add_dep(task4);

      std::vector<IndexRange> I244_index = {this->active_, this->closed_};
      std::shared_ptr<Tensor> I244(new Tensor(I244_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor205 = {den2, I244};
      auto task205 = std::make_shared<Task205>(tensor205, pindex);
      task205->add_dep(task152);
      density_->add_task(task205);


      std::vector<IndexRange> I245_index = {this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I245(new Tensor(I245_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor206 = {I244, t2, I245};
      auto task206 = std::make_shared<Task206>(tensor206, pindex);
      task205->add_dep(task206);
      task206->add_dep(task152);
      density_->add_task(task206);


      std::vector<IndexRange> I246_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I246(new Tensor(I246_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor207 = {I245, t2, I246};
      auto task207 = std::make_shared<Task207>(tensor207, pindex);
      task206->add_dep(task207);
      task207->add_dep(task152);
      density_->add_task(task207);


      std::vector<std::shared_ptr<Tensor>> tensor208 = {I246, Gamma14};
      auto task208 = std::make_shared<Task208>(tensor208, pindex);
      task207->add_dep(task208);
      task208->add_dep(task152);
      density_->add_task(task208);

      task208->add_dep(task4);

      std::vector<IndexRange> I247_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I247(new Tensor(I247_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor209 = {den2, I247};
      auto task209 = std::make_shared<Task209>(tensor209, pindex);
      task209->add_dep(task152);
      density_->add_task(task209);


      std::vector<IndexRange> I248_index = {this->active_, this->active_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I248(new Tensor(I248_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor210 = {I247, t2, I248};
      auto task210 = std::make_shared<Task210>(tensor210, pindex);
      task209->add_dep(task210);
      task210->add_dep(task152);
      density_->add_task(task210);


      std::vector<IndexRange> I249_index = {this->active_, this->active_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I249(new Tensor(I249_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor211 = {I248, t2, I249};
      auto task211 = std::make_shared<Task211>(tensor211, pindex);
      task210->add_dep(task211);
      task211->add_dep(task152);
      density_->add_task(task211);


      std::vector<std::shared_ptr<Tensor>> tensor212 = {I249, Gamma67};
      auto task212 = std::make_shared<Task212>(tensor212, pindex);
      task211->add_dep(task212);
      task212->add_dep(task152);
      density_->add_task(task212);

      task212->add_dep(task6);

      std::vector<IndexRange> I250_index = {this->virt_, this->virt_};
      std::shared_ptr<Tensor> I250(new Tensor(I250_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor213 = {den2, I250};
      auto task213 = std::make_shared<Task213>(tensor213, pindex);
      task213->add_dep(task152);
      density_->add_task(task213);


      std::vector<IndexRange> I251_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I251(new Tensor(I251_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor214 = {I250, t2, I251};
      auto task214 = std::make_shared<Task214>(tensor214, pindex);
      task213->add_dep(task214);
      task214->add_dep(task152);
      density_->add_task(task214);


      std::vector<IndexRange> I252_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I252(new Tensor(I252_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor215 = {I251, t2, I252};
      auto task215 = std::make_shared<Task215>(tensor215, pindex);
      task214->add_dep(task215);
      task215->add_dep(task152);
      density_->add_task(task215);


      std::vector<std::shared_ptr<Tensor>> tensor216 = {I252, Gamma14};
      auto task216 = std::make_shared<Task216>(tensor216, pindex);
      task215->add_dep(task216);
      task216->add_dep(task152);
      density_->add_task(task216);

      task216->add_dep(task4);

      auto density1_ = std::make_shared<Queue>();
      auto density2_ = std::make_shared<Queue>();
      std::vector<std::shared_ptr<Tensor>> tensor217 = {Den1};
      auto task217 = std::make_shared<Task217>(tensor217);
      density2_->add_task(task217);

      std::vector<IndexRange> I253_index = {this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I253(new Tensor(I253_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor218 = {Den1, I253};
      auto task218 = std::make_shared<Task218>(tensor218, pindex);
      task218->add_dep(task217);
      density2_->add_task(task218);


      std::vector<std::shared_ptr<Tensor>> tensor219 = {I253, t2};
      auto task219 = std::make_shared<Task219>(tensor219, pindex);
      task218->add_dep(task219);
      task219->add_dep(task217);
      density2_->add_task(task219);


      std::vector<IndexRange> I255_index = {this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I255(new Tensor(I255_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor220 = {Den1, I255};
      auto task220 = std::make_shared<Task220>(tensor220, pindex);
      task220->add_dep(task217);
      density2_->add_task(task220);


      std::vector<IndexRange> I256_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I256(new Tensor(I256_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor221 = {I255, t2, I256};
      auto task221 = std::make_shared<Task221>(tensor221, pindex);
      task220->add_dep(task221);
      task221->add_dep(task217);
      density2_->add_task(task221);


      std::vector<std::shared_ptr<Tensor>> tensor222 = {I256, Gamma2};
      auto task222 = std::make_shared<Task222>(tensor222, pindex);
      task221->add_dep(task222);
      task222->add_dep(task217);
      density2_->add_task(task222);

      task222->add_dep(task2);

      std::vector<IndexRange> I258_index = {this->active_, this->active_};
      std::shared_ptr<Tensor> I258(new Tensor(I258_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor223 = {I255, t2, I258};
      auto task223 = std::make_shared<Task223>(tensor223, pindex);
      task220->add_dep(task223);
      task223->add_dep(task217);
      density2_->add_task(task223);


      std::vector<std::shared_ptr<Tensor>> tensor224 = {I258, Gamma2};
      auto task224 = std::make_shared<Task224>(tensor224, pindex);
      task223->add_dep(task224);
      task224->add_dep(task217);
      density2_->add_task(task224);

      task224->add_dep(task2);

      std::vector<IndexRange> I259_index = {this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I259(new Tensor(I259_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor225 = {Den1, I259};
      auto task225 = std::make_shared<Task225>(tensor225, pindex);
      task225->add_dep(task217);
      density2_->add_task(task225);


      std::vector<IndexRange> I260_index = {this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I260(new Tensor(I260_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor226 = {I259, t2, I260};
      auto task226 = std::make_shared<Task226>(tensor226, pindex);
      task225->add_dep(task226);
      task226->add_dep(task217);
      density2_->add_task(task226);


      std::vector<std::shared_ptr<Tensor>> tensor227 = {I260, Gamma14};
      auto task227 = std::make_shared<Task227>(tensor227, pindex);
      task226->add_dep(task227);
      task227->add_dep(task217);
      density2_->add_task(task227);

      task227->add_dep(task4);

      auto dedci_ = std::make_shared<Queue>();
      std::vector<std::shared_ptr<Tensor>> tensor228 = {deci};
      auto task228 = std::make_shared<Task228>(tensor228);
      dedci_->add_task(task228);

      std::vector<IndexRange> I261_index = {this->ci_};
      std::shared_ptr<Tensor> I261(new Tensor(I261_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor229 = {deci, I261};
      auto task229 = std::make_shared<Task229>(tensor229, cindex);
      task229->add_dep(task228);
      dedci_->add_task(task229);


      std::vector<IndexRange> I262_index = {this->ci_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I262(new Tensor(I262_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor230 = {I261, t2, I262};
      auto task230 = std::make_shared<Task230>(tensor230, cindex);
      task229->add_dep(task230);
      task230->add_dep(task228);
      dedci_->add_task(task230);


      std::vector<IndexRange> I263_index = {this->ci_};
      std::shared_ptr<Tensor> I263(new Tensor(I263_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor231 = {I262, t2, I263};
      auto task231 = std::make_shared<Task231>(tensor231, cindex);
      task230->add_dep(task231);
      task231->add_dep(task228);
      dedci_->add_task(task231);


      std::vector<std::shared_ptr<Tensor>> tensor232 = {I263, Gamma72};
      auto task232 = std::make_shared<Task232>(tensor232, cindex);
      task231->add_dep(task232);
      task232->add_dep(task228);
      dedci_->add_task(task232);

      task232->add_dep(task7);
      task232->add_dep(task7);

      std::vector<IndexRange> I266_index = {this->ci_};
      std::shared_ptr<Tensor> I266(new Tensor(I266_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor233 = {I262, t2, I266};
      auto task233 = std::make_shared<Task233>(tensor233, cindex);
      task230->add_dep(task233);
      task233->add_dep(task228);
      dedci_->add_task(task233);


      std::vector<std::shared_ptr<Tensor>> tensor234 = {I266, Gamma72};
      auto task234 = std::make_shared<Task234>(tensor234, cindex);
      task233->add_dep(task234);
      task234->add_dep(task228);
      dedci_->add_task(task234);

      task234->add_dep(task7);
      task234->add_dep(task7);

      std::vector<IndexRange> I269_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I269(new Tensor(I269_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor235 = {I262, this->f1_, I269};
      auto task235 = std::make_shared<Task235>(tensor235, cindex);
      task230->add_dep(task235);
      task235->add_dep(task228);
      dedci_->add_task(task235);


      std::vector<IndexRange> I270_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I270(new Tensor(I270_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor236 = {I269, t2, I270};
      auto task236 = std::make_shared<Task236>(tensor236, cindex);
      task235->add_dep(task236);
      task236->add_dep(task228);
      dedci_->add_task(task236);


      std::vector<std::shared_ptr<Tensor>> tensor237 = {I270, Gamma74};
      auto task237 = std::make_shared<Task237>(tensor237, cindex);
      task236->add_dep(task237);
      task237->add_dep(task228);
      dedci_->add_task(task237);

      task237->add_dep(task8);

      std::vector<IndexRange> I274_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I274(new Tensor(I274_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor238 = {I269, t2, I274};
      auto task238 = std::make_shared<Task238>(tensor238, cindex);
      task235->add_dep(task238);
      task238->add_dep(task228);
      dedci_->add_task(task238);


      std::vector<std::shared_ptr<Tensor>> tensor239 = {I274, Gamma74};
      auto task239 = std::make_shared<Task239>(tensor239, cindex);
      task238->add_dep(task239);
      task239->add_dep(task228);
      dedci_->add_task(task239);

      task239->add_dep(task8);

      std::vector<IndexRange> I336_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I336(new Tensor(I336_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor240 = {I262, this->f1_, I336};
      auto task240 = std::make_shared<Task240>(tensor240, cindex);
      task230->add_dep(task240);
      task240->add_dep(task228);
      dedci_->add_task(task240);


      std::vector<IndexRange> I337_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I337(new Tensor(I337_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor241 = {I336, t2, I337};
      auto task241 = std::make_shared<Task241>(tensor241, cindex);
      task240->add_dep(task241);
      task241->add_dep(task228);
      dedci_->add_task(task241);


      std::vector<std::shared_ptr<Tensor>> tensor242 = {I337, Gamma74};
      auto task242 = std::make_shared<Task242>(tensor242, cindex);
      task241->add_dep(task242);
      task242->add_dep(task228);
      dedci_->add_task(task242);

      task242->add_dep(task8);

      std::vector<IndexRange> I341_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I341(new Tensor(I341_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor243 = {I336, t2, I341};
      auto task243 = std::make_shared<Task243>(tensor243, cindex);
      task240->add_dep(task243);
      task243->add_dep(task228);
      dedci_->add_task(task243);


      std::vector<std::shared_ptr<Tensor>> tensor244 = {I341, Gamma74};
      auto task244 = std::make_shared<Task244>(tensor244, cindex);
      task243->add_dep(task244);
      task244->add_dep(task228);
      dedci_->add_task(task244);

      task244->add_dep(task8);

      std::vector<IndexRange> I276_index = {this->ci_, this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I276(new Tensor(I276_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor245 = {I261, t2, I276};
      auto task245 = std::make_shared<Task245>(tensor245, cindex);
      task229->add_dep(task245);
      task245->add_dep(task228);
      dedci_->add_task(task245);


      std::vector<IndexRange> I277_index = {this->ci_, this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I277(new Tensor(I277_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor246 = {I276, this->f1_, I277};
      auto task246 = std::make_shared<Task246>(tensor246, cindex);
      task245->add_dep(task246);
      task246->add_dep(task228);
      dedci_->add_task(task246);


      std::vector<IndexRange> I278_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I278(new Tensor(I278_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor247 = {I277, t2, I278};
      auto task247 = std::make_shared<Task247>(tensor247, cindex);
      task246->add_dep(task247);
      task247->add_dep(task228);
      dedci_->add_task(task247);


      std::vector<std::shared_ptr<Tensor>> tensor248 = {I278, Gamma74};
      auto task248 = std::make_shared<Task248>(tensor248, cindex);
      task247->add_dep(task248);
      task248->add_dep(task228);
      dedci_->add_task(task248);

      task248->add_dep(task8);

      std::vector<IndexRange> I282_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I282(new Tensor(I282_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor249 = {I277, t2, I282};
      auto task249 = std::make_shared<Task249>(tensor249, cindex);
      task246->add_dep(task249);
      task249->add_dep(task228);
      dedci_->add_task(task249);


      std::vector<std::shared_ptr<Tensor>> tensor250 = {I282, Gamma74};
      auto task250 = std::make_shared<Task250>(tensor250, cindex);
      task249->add_dep(task250);
      task250->add_dep(task228);
      dedci_->add_task(task250);

      task250->add_dep(task8);

      std::vector<IndexRange> I285_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I285(new Tensor(I285_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor251 = {I276, t2, I285};
      auto task251 = std::make_shared<Task251>(tensor251, cindex);
      task245->add_dep(task251);
      task251->add_dep(task228);
      dedci_->add_task(task251);


      std::vector<std::shared_ptr<Tensor>> tensor252 = {I285, Gamma78};
      auto task252 = std::make_shared<Task252>(tensor252, cindex);
      task251->add_dep(task252);
      task252->add_dep(task228);
      dedci_->add_task(task252);

      task252->add_dep(task9);

      std::vector<IndexRange> I288_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I288(new Tensor(I288_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor253 = {I276, t2, I288};
      auto task253 = std::make_shared<Task253>(tensor253, cindex);
      task245->add_dep(task253);
      task253->add_dep(task228);
      dedci_->add_task(task253);


      std::vector<std::shared_ptr<Tensor>> tensor254 = {I288, Gamma78};
      auto task254 = std::make_shared<Task254>(tensor254, cindex);
      task253->add_dep(task254);
      task254->add_dep(task228);
      dedci_->add_task(task254);

      task254->add_dep(task9);

      std::vector<IndexRange> I291_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I291(new Tensor(I291_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor255 = {I276, this->f1_, I291};
      auto task255 = std::make_shared<Task255>(tensor255, cindex);
      task245->add_dep(task255);
      task255->add_dep(task228);
      dedci_->add_task(task255);


      std::vector<IndexRange> I292_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I292(new Tensor(I292_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor256 = {I291, t2, I292};
      auto task256 = std::make_shared<Task256>(tensor256, cindex);
      task255->add_dep(task256);
      task256->add_dep(task228);
      dedci_->add_task(task256);


      std::vector<std::shared_ptr<Tensor>> tensor257 = {I292, Gamma74};
      auto task257 = std::make_shared<Task257>(tensor257, cindex);
      task256->add_dep(task257);
      task257->add_dep(task228);
      dedci_->add_task(task257);

      task257->add_dep(task8);

      std::vector<IndexRange> I296_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I296(new Tensor(I296_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor258 = {I291, t2, I296};
      auto task258 = std::make_shared<Task258>(tensor258, cindex);
      task255->add_dep(task258);
      task258->add_dep(task228);
      dedci_->add_task(task258);


      std::vector<std::shared_ptr<Tensor>> tensor259 = {I296, Gamma74};
      auto task259 = std::make_shared<Task259>(tensor259, cindex);
      task258->add_dep(task259);
      task259->add_dep(task228);
      dedci_->add_task(task259);

      task259->add_dep(task8);

      std::vector<IndexRange> I299_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I299(new Tensor(I299_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor260 = {I276, this->f1_, I299};
      auto task260 = std::make_shared<Task260>(tensor260, cindex);
      task245->add_dep(task260);
      task260->add_dep(task228);
      dedci_->add_task(task260);


      std::vector<IndexRange> I300_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I300(new Tensor(I300_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor261 = {I299, t2, I300};
      auto task261 = std::make_shared<Task261>(tensor261, cindex);
      task260->add_dep(task261);
      task261->add_dep(task228);
      dedci_->add_task(task261);


      std::vector<std::shared_ptr<Tensor>> tensor262 = {I300, Gamma74};
      auto task262 = std::make_shared<Task262>(tensor262, cindex);
      task261->add_dep(task262);
      task262->add_dep(task228);
      dedci_->add_task(task262);

      task262->add_dep(task8);

      std::vector<IndexRange> I304_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I304(new Tensor(I304_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor263 = {I299, t2, I304};
      auto task263 = std::make_shared<Task263>(tensor263, cindex);
      task260->add_dep(task263);
      task263->add_dep(task228);
      dedci_->add_task(task263);


      std::vector<std::shared_ptr<Tensor>> tensor264 = {I304, Gamma74};
      auto task264 = std::make_shared<Task264>(tensor264, cindex);
      task263->add_dep(task264);
      task264->add_dep(task228);
      dedci_->add_task(task264);

      task264->add_dep(task8);

      std::vector<IndexRange> I307_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I307(new Tensor(I307_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor265 = {I276, this->f1_, I307};
      auto task265 = std::make_shared<Task265>(tensor265, cindex);
      task245->add_dep(task265);
      task265->add_dep(task228);
      dedci_->add_task(task265);


      std::vector<IndexRange> I308_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I308(new Tensor(I308_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor266 = {I307, t2, I308};
      auto task266 = std::make_shared<Task266>(tensor266, cindex);
      task265->add_dep(task266);
      task266->add_dep(task228);
      dedci_->add_task(task266);


      std::vector<std::shared_ptr<Tensor>> tensor267 = {I308, Gamma74};
      auto task267 = std::make_shared<Task267>(tensor267, cindex);
      task266->add_dep(task267);
      task267->add_dep(task228);
      dedci_->add_task(task267);

      task267->add_dep(task8);

      std::vector<IndexRange> I312_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I312(new Tensor(I312_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor268 = {I307, t2, I312};
      auto task268 = std::make_shared<Task268>(tensor268, cindex);
      task265->add_dep(task268);
      task268->add_dep(task228);
      dedci_->add_task(task268);


      std::vector<std::shared_ptr<Tensor>> tensor269 = {I312, Gamma74};
      auto task269 = std::make_shared<Task269>(tensor269, cindex);
      task268->add_dep(task269);
      task269->add_dep(task228);
      dedci_->add_task(task269);

      task269->add_dep(task8);

      std::vector<IndexRange> I315_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I315(new Tensor(I315_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor270 = {I276, this->f1_, I315};
      auto task270 = std::make_shared<Task270>(tensor270, cindex);
      task245->add_dep(task270);
      task270->add_dep(task228);
      dedci_->add_task(task270);


      std::vector<IndexRange> I316_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I316(new Tensor(I316_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor271 = {I315, t2, I316};
      auto task271 = std::make_shared<Task271>(tensor271, cindex);
      task270->add_dep(task271);
      task271->add_dep(task228);
      dedci_->add_task(task271);


      std::vector<std::shared_ptr<Tensor>> tensor272 = {I316, Gamma86};
      auto task272 = std::make_shared<Task272>(tensor272, cindex);
      task271->add_dep(task272);
      task272->add_dep(task228);
      dedci_->add_task(task272);

      task272->add_dep(task10);

      std::vector<IndexRange> I397_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I397(new Tensor(I397_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor273 = {I276, t2, I397};
      auto task273 = std::make_shared<Task273>(tensor273, cindex);
      task245->add_dep(task273);
      task273->add_dep(task228);
      dedci_->add_task(task273);


      std::vector<std::shared_ptr<Tensor>> tensor274 = {I397, Gamma74};
      auto task274 = std::make_shared<Task274>(tensor274, cindex, this->e0_);
      task273->add_dep(task274);
      task274->add_dep(task228);
      dedci_->add_task(task274);

      task274->add_dep(task8);

      std::vector<IndexRange> I400_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I400(new Tensor(I400_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor275 = {I276, t2, I400};
      auto task275 = std::make_shared<Task275>(tensor275, cindex);
      task245->add_dep(task275);
      task275->add_dep(task228);
      dedci_->add_task(task275);


      std::vector<std::shared_ptr<Tensor>> tensor276 = {I400, Gamma74};
      auto task276 = std::make_shared<Task276>(tensor276, cindex, this->e0_);
      task275->add_dep(task276);
      task276->add_dep(task228);
      dedci_->add_task(task276);

      task276->add_dep(task8);

      std::vector<IndexRange> I415_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I415(new Tensor(I415_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor277 = {I276, this->v2_, I415};
      auto task277 = std::make_shared<Task277>(tensor277, cindex);
      task245->add_dep(task277);
      task277->add_dep(task228);
      dedci_->add_task(task277);


      std::vector<std::shared_ptr<Tensor>> tensor278 = {I415, Gamma74};
      auto task278 = std::make_shared<Task278>(tensor278, cindex);
      task277->add_dep(task278);
      task278->add_dep(task228);
      dedci_->add_task(task278);

      task278->add_dep(task8);

      std::vector<IndexRange> I418_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I418(new Tensor(I418_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor279 = {I276, this->v2_, I418};
      auto task279 = std::make_shared<Task279>(tensor279, cindex);
      task245->add_dep(task279);
      task279->add_dep(task228);
      dedci_->add_task(task279);


      std::vector<std::shared_ptr<Tensor>> tensor280 = {I418, Gamma74};
      auto task280 = std::make_shared<Task280>(tensor280, cindex);
      task279->add_dep(task280);
      task280->add_dep(task228);
      dedci_->add_task(task280);

      task280->add_dep(task8);

      std::vector<IndexRange> I318_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I318(new Tensor(I318_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor281 = {I261, t2, I318};
      auto task281 = std::make_shared<Task281>(tensor281, cindex);
      task229->add_dep(task281);
      task281->add_dep(task228);
      dedci_->add_task(task281);


      std::vector<IndexRange> I319_index = {this->ci_, this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I319(new Tensor(I319_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor282 = {I318, this->f1_, I319};
      auto task282 = std::make_shared<Task282>(tensor282, cindex);
      task281->add_dep(task282);
      task282->add_dep(task228);
      dedci_->add_task(task282);


      std::vector<IndexRange> I320_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I320(new Tensor(I320_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor283 = {I319, t2, I320};
      auto task283 = std::make_shared<Task283>(tensor283, cindex);
      task282->add_dep(task283);
      task283->add_dep(task228);
      dedci_->add_task(task283);


      std::vector<std::shared_ptr<Tensor>> tensor284 = {I320, Gamma86};
      auto task284 = std::make_shared<Task284>(tensor284, cindex);
      task283->add_dep(task284);
      task284->add_dep(task228);
      dedci_->add_task(task284);

      task284->add_dep(task10);

      std::vector<IndexRange> I323_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I323(new Tensor(I323_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor285 = {I318, t2, I323};
      auto task285 = std::make_shared<Task285>(tensor285, cindex);
      task281->add_dep(task285);
      task285->add_dep(task228);
      dedci_->add_task(task285);


      std::vector<std::shared_ptr<Tensor>> tensor286 = {I323, Gamma88};
      auto task286 = std::make_shared<Task286>(tensor286, cindex);
      task285->add_dep(task286);
      task286->add_dep(task228);
      dedci_->add_task(task286);

      task286->add_dep(task11);

      std::vector<IndexRange> I326_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I326(new Tensor(I326_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor287 = {I318, this->f1_, I326};
      auto task287 = std::make_shared<Task287>(tensor287, cindex);
      task281->add_dep(task287);
      task287->add_dep(task228);
      dedci_->add_task(task287);


      std::vector<IndexRange> I327_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I327(new Tensor(I327_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor288 = {I326, t2, I327};
      auto task288 = std::make_shared<Task288>(tensor288, cindex);
      task287->add_dep(task288);
      task288->add_dep(task228);
      dedci_->add_task(task288);


      std::vector<std::shared_ptr<Tensor>> tensor289 = {I327, Gamma86};
      auto task289 = std::make_shared<Task289>(tensor289, cindex);
      task288->add_dep(task289);
      task289->add_dep(task228);
      dedci_->add_task(task289);

      task289->add_dep(task10);

      std::vector<IndexRange> I403_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I403(new Tensor(I403_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor290 = {I318, t2, I403};
      auto task290 = std::make_shared<Task290>(tensor290, cindex);
      task281->add_dep(task290);
      task290->add_dep(task228);
      dedci_->add_task(task290);


      std::vector<std::shared_ptr<Tensor>> tensor291 = {I403, Gamma86};
      auto task291 = std::make_shared<Task291>(tensor291, cindex, this->e0_);
      task290->add_dep(task291);
      task291->add_dep(task228);
      dedci_->add_task(task291);

      task291->add_dep(task10);

      std::vector<IndexRange> I421_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I421(new Tensor(I421_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor292 = {I318, this->v2_, I421};
      auto task292 = std::make_shared<Task292>(tensor292, cindex);
      task281->add_dep(task292);
      task292->add_dep(task228);
      dedci_->add_task(task292);


      std::vector<std::shared_ptr<Tensor>> tensor293 = {I421, Gamma86};
      auto task293 = std::make_shared<Task293>(tensor293, cindex);
      task292->add_dep(task293);
      task293->add_dep(task228);
      dedci_->add_task(task293);

      task293->add_dep(task10);

      std::vector<IndexRange> I343_index = {this->ci_, this->active_, this->closed_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I343(new Tensor(I343_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor294 = {I261, t2, I343};
      auto task294 = std::make_shared<Task294>(tensor294, cindex);
      task229->add_dep(task294);
      task294->add_dep(task228);
      dedci_->add_task(task294);


      std::vector<IndexRange> I344_index = {this->ci_, this->active_, this->active_, this->closed_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I344(new Tensor(I344_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor295 = {I343, this->f1_, I344};
      auto task295 = std::make_shared<Task295>(tensor295, cindex);
      task294->add_dep(task295);
      task295->add_dep(task228);
      dedci_->add_task(task295);


      std::vector<IndexRange> I345_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I345(new Tensor(I345_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor296 = {I344, t2, I345};
      auto task296 = std::make_shared<Task296>(tensor296, cindex);
      task295->add_dep(task296);
      task296->add_dep(task228);
      dedci_->add_task(task296);


      std::vector<std::shared_ptr<Tensor>> tensor297 = {I345, Gamma74};
      auto task297 = std::make_shared<Task297>(tensor297, cindex);
      task296->add_dep(task297);
      task297->add_dep(task228);
      dedci_->add_task(task297);

      task297->add_dep(task8);

      std::vector<IndexRange> I349_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I349(new Tensor(I349_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor298 = {I344, t2, I349};
      auto task298 = std::make_shared<Task298>(tensor298, cindex);
      task295->add_dep(task298);
      task298->add_dep(task228);
      dedci_->add_task(task298);


      std::vector<std::shared_ptr<Tensor>> tensor299 = {I349, Gamma74};
      auto task299 = std::make_shared<Task299>(tensor299, cindex);
      task298->add_dep(task299);
      task299->add_dep(task228);
      dedci_->add_task(task299);

      task299->add_dep(task8);

      std::vector<IndexRange> I358_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I358(new Tensor(I358_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor300 = {I343, this->f1_, I358};
      auto task300 = std::make_shared<Task300>(tensor300, cindex);
      task294->add_dep(task300);
      task300->add_dep(task228);
      dedci_->add_task(task300);


      std::vector<IndexRange> I359_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I359(new Tensor(I359_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor301 = {I358, t2, I359};
      auto task301 = std::make_shared<Task301>(tensor301, cindex);
      task300->add_dep(task301);
      task301->add_dep(task228);
      dedci_->add_task(task301);


      std::vector<std::shared_ptr<Tensor>> tensor302 = {I359, Gamma74};
      auto task302 = std::make_shared<Task302>(tensor302, cindex);
      task301->add_dep(task302);
      task302->add_dep(task228);
      dedci_->add_task(task302);

      task302->add_dep(task8);

      std::vector<IndexRange> I363_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I363(new Tensor(I363_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor303 = {I358, t2, I363};
      auto task303 = std::make_shared<Task303>(tensor303, cindex);
      task300->add_dep(task303);
      task303->add_dep(task228);
      dedci_->add_task(task303);


      std::vector<std::shared_ptr<Tensor>> tensor304 = {I363, Gamma74};
      auto task304 = std::make_shared<Task304>(tensor304, cindex);
      task303->add_dep(task304);
      task304->add_dep(task228);
      dedci_->add_task(task304);

      task304->add_dep(task8);

      std::vector<IndexRange> I366_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I366(new Tensor(I366_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor305 = {I343, this->f1_, I366};
      auto task305 = std::make_shared<Task305>(tensor305, cindex);
      task294->add_dep(task305);
      task305->add_dep(task228);
      dedci_->add_task(task305);


      std::vector<IndexRange> I367_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I367(new Tensor(I367_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor306 = {I366, t2, I367};
      auto task306 = std::make_shared<Task306>(tensor306, cindex);
      task305->add_dep(task306);
      task306->add_dep(task228);
      dedci_->add_task(task306);


      std::vector<std::shared_ptr<Tensor>> tensor307 = {I367, Gamma74};
      auto task307 = std::make_shared<Task307>(tensor307, cindex);
      task306->add_dep(task307);
      task307->add_dep(task228);
      dedci_->add_task(task307);

      task307->add_dep(task8);

      std::vector<IndexRange> I371_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I371(new Tensor(I371_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor308 = {I366, t2, I371};
      auto task308 = std::make_shared<Task308>(tensor308, cindex);
      task305->add_dep(task308);
      task308->add_dep(task228);
      dedci_->add_task(task308);


      std::vector<std::shared_ptr<Tensor>> tensor309 = {I371, Gamma74};
      auto task309 = std::make_shared<Task309>(tensor309, cindex);
      task308->add_dep(task309);
      task309->add_dep(task228);
      dedci_->add_task(task309);

      task309->add_dep(task8);

      std::vector<IndexRange> I374_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I374(new Tensor(I374_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor310 = {I343, this->f1_, I374};
      auto task310 = std::make_shared<Task310>(tensor310, cindex);
      task294->add_dep(task310);
      task310->add_dep(task228);
      dedci_->add_task(task310);


      std::vector<IndexRange> I375_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I375(new Tensor(I375_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor311 = {I374, t2, I375};
      auto task311 = std::make_shared<Task311>(tensor311, cindex);
      task310->add_dep(task311);
      task311->add_dep(task228);
      dedci_->add_task(task311);


      std::vector<std::shared_ptr<Tensor>> tensor312 = {I375, Gamma74};
      auto task312 = std::make_shared<Task312>(tensor312, cindex);
      task311->add_dep(task312);
      task312->add_dep(task228);
      dedci_->add_task(task312);

      task312->add_dep(task8);

      std::vector<IndexRange> I379_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I379(new Tensor(I379_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor313 = {I374, t2, I379};
      auto task313 = std::make_shared<Task313>(tensor313, cindex);
      task310->add_dep(task313);
      task313->add_dep(task228);
      dedci_->add_task(task313);


      std::vector<std::shared_ptr<Tensor>> tensor314 = {I379, Gamma74};
      auto task314 = std::make_shared<Task314>(tensor314, cindex);
      task313->add_dep(task314);
      task314->add_dep(task228);
      dedci_->add_task(task314);

      task314->add_dep(task8);

      std::vector<IndexRange> I406_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I406(new Tensor(I406_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor315 = {I343, t2, I406};
      auto task315 = std::make_shared<Task315>(tensor315, cindex);
      task294->add_dep(task315);
      task315->add_dep(task228);
      dedci_->add_task(task315);


      std::vector<std::shared_ptr<Tensor>> tensor316 = {I406, Gamma74};
      auto task316 = std::make_shared<Task316>(tensor316, cindex, this->e0_);
      task315->add_dep(task316);
      task316->add_dep(task228);
      dedci_->add_task(task316);

      task316->add_dep(task8);

      std::vector<IndexRange> I409_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I409(new Tensor(I409_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor317 = {I343, t2, I409};
      auto task317 = std::make_shared<Task317>(tensor317, cindex);
      task294->add_dep(task317);
      task317->add_dep(task228);
      dedci_->add_task(task317);


      std::vector<std::shared_ptr<Tensor>> tensor318 = {I409, Gamma74};
      auto task318 = std::make_shared<Task318>(tensor318, cindex, this->e0_);
      task317->add_dep(task318);
      task318->add_dep(task228);
      dedci_->add_task(task318);

      task318->add_dep(task8);

      std::vector<IndexRange> I424_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I424(new Tensor(I424_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor319 = {I343, this->v2_, I424};
      auto task319 = std::make_shared<Task319>(tensor319, cindex);
      task294->add_dep(task319);
      task319->add_dep(task228);
      dedci_->add_task(task319);


      std::vector<std::shared_ptr<Tensor>> tensor320 = {I424, Gamma74};
      auto task320 = std::make_shared<Task320>(tensor320, cindex);
      task319->add_dep(task320);
      task320->add_dep(task228);
      dedci_->add_task(task320);

      task320->add_dep(task8);

      std::vector<IndexRange> I427_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I427(new Tensor(I427_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor321 = {I343, this->v2_, I427};
      auto task321 = std::make_shared<Task321>(tensor321, cindex);
      task294->add_dep(task321);
      task321->add_dep(task228);
      dedci_->add_task(task321);


      std::vector<std::shared_ptr<Tensor>> tensor322 = {I427, Gamma74};
      auto task322 = std::make_shared<Task322>(tensor322, cindex);
      task321->add_dep(task322);
      task322->add_dep(task228);
      dedci_->add_task(task322);

      task322->add_dep(task8);

      std::vector<IndexRange> I351_index = {this->ci_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I351(new Tensor(I351_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor323 = {I261, t2, I351};
      auto task323 = std::make_shared<Task323>(tensor323, cindex);
      task229->add_dep(task323);
      task323->add_dep(task228);
      dedci_->add_task(task323);


      std::vector<IndexRange> I352_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I352(new Tensor(I352_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor324 = {I351, t2, I352};
      auto task324 = std::make_shared<Task324>(tensor324, cindex);
      task323->add_dep(task324);
      task324->add_dep(task228);
      dedci_->add_task(task324);


      std::vector<std::shared_ptr<Tensor>> tensor325 = {I352, Gamma78};
      auto task325 = std::make_shared<Task325>(tensor325, cindex);
      task324->add_dep(task325);
      task325->add_dep(task228);
      dedci_->add_task(task325);

      task325->add_dep(task9);

      std::vector<IndexRange> I355_index = {this->ci_, this->active_, this->active_};
      std::shared_ptr<Tensor> I355(new Tensor(I355_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor326 = {I351, t2, I355};
      auto task326 = std::make_shared<Task326>(tensor326, cindex);
      task323->add_dep(task326);
      task326->add_dep(task228);
      dedci_->add_task(task326);


      std::vector<std::shared_ptr<Tensor>> tensor327 = {I355, Gamma78};
      auto task327 = std::make_shared<Task327>(tensor327, cindex);
      task326->add_dep(task327);
      task327->add_dep(task228);
      dedci_->add_task(task327);

      task327->add_dep(task9);

      std::vector<IndexRange> I382_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I382(new Tensor(I382_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor328 = {I351, this->f1_, I382};
      auto task328 = std::make_shared<Task328>(tensor328, cindex);
      task323->add_dep(task328);
      task328->add_dep(task228);
      dedci_->add_task(task328);


      std::vector<IndexRange> I383_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I383(new Tensor(I383_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor329 = {I382, t2, I383};
      auto task329 = std::make_shared<Task329>(tensor329, cindex);
      task328->add_dep(task329);
      task329->add_dep(task228);
      dedci_->add_task(task329);


      std::vector<std::shared_ptr<Tensor>> tensor330 = {I383, Gamma86};
      auto task330 = std::make_shared<Task330>(tensor330, cindex);
      task329->add_dep(task330);
      task330->add_dep(task228);
      dedci_->add_task(task330);

      task330->add_dep(task10);

      std::vector<IndexRange> I385_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I385(new Tensor(I385_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor331 = {I261, t2, I385};
      auto task331 = std::make_shared<Task331>(tensor331, cindex);
      task229->add_dep(task331);
      task331->add_dep(task228);
      dedci_->add_task(task331);


      std::vector<IndexRange> I386_index = {this->ci_, this->active_, this->active_, this->active_, this->virt_, this->closed_, this->virt_};
      std::shared_ptr<Tensor> I386(new Tensor(I386_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor332 = {I385, this->f1_, I386};
      auto task332 = std::make_shared<Task332>(tensor332, cindex);
      task331->add_dep(task332);
      task332->add_dep(task228);
      dedci_->add_task(task332);


      std::vector<IndexRange> I387_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I387(new Tensor(I387_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor333 = {I386, t2, I387};
      auto task333 = std::make_shared<Task333>(tensor333, cindex);
      task332->add_dep(task333);
      task333->add_dep(task228);
      dedci_->add_task(task333);


      std::vector<std::shared_ptr<Tensor>> tensor334 = {I387, Gamma86};
      auto task334 = std::make_shared<Task334>(tensor334, cindex);
      task333->add_dep(task334);
      task334->add_dep(task228);
      dedci_->add_task(task334);

      task334->add_dep(task10);

      std::vector<IndexRange> I393_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I393(new Tensor(I393_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor335 = {I385, this->f1_, I393};
      auto task335 = std::make_shared<Task335>(tensor335, cindex);
      task331->add_dep(task335);
      task335->add_dep(task228);
      dedci_->add_task(task335);


      std::vector<IndexRange> I394_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I394(new Tensor(I394_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor336 = {I393, t2, I394};
      auto task336 = std::make_shared<Task336>(tensor336, cindex);
      task335->add_dep(task336);
      task336->add_dep(task228);
      dedci_->add_task(task336);


      std::vector<std::shared_ptr<Tensor>> tensor337 = {I394, Gamma86};
      auto task337 = std::make_shared<Task337>(tensor337, cindex);
      task336->add_dep(task337);
      task337->add_dep(task228);
      dedci_->add_task(task337);

      task337->add_dep(task10);

      std::vector<IndexRange> I412_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I412(new Tensor(I412_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor338 = {I385, t2, I412};
      auto task338 = std::make_shared<Task338>(tensor338, cindex);
      task331->add_dep(task338);
      task338->add_dep(task228);
      dedci_->add_task(task338);


      std::vector<std::shared_ptr<Tensor>> tensor339 = {I412, Gamma86};
      auto task339 = std::make_shared<Task339>(tensor339, cindex, this->e0_);
      task338->add_dep(task339);
      task339->add_dep(task228);
      dedci_->add_task(task339);

      task339->add_dep(task10);

      std::vector<IndexRange> I430_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I430(new Tensor(I430_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor340 = {I385, this->v2_, I430};
      auto task340 = std::make_shared<Task340>(tensor340, cindex);
      task331->add_dep(task340);
      task340->add_dep(task228);
      dedci_->add_task(task340);


      std::vector<std::shared_ptr<Tensor>> tensor341 = {I430, Gamma86};
      auto task341 = std::make_shared<Task341>(tensor341, cindex);
      task340->add_dep(task341);
      task341->add_dep(task228);
      dedci_->add_task(task341);

      task341->add_dep(task10);

      std::vector<IndexRange> I389_index = {this->ci_, this->active_, this->active_, this->virt_, this->virt_};
      std::shared_ptr<Tensor> I389(new Tensor(I389_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor342 = {I261, t2, I389};
      auto task342 = std::make_shared<Task342>(tensor342, cindex);
      task229->add_dep(task342);
      task342->add_dep(task228);
      dedci_->add_task(task342);


      std::vector<IndexRange> I390_index = {this->ci_, this->active_, this->active_, this->active_, this->active_};
      std::shared_ptr<Tensor> I390(new Tensor(I390_index, false));
      std::vector<std::shared_ptr<Tensor>> tensor343 = {I389, t2, I390};
      auto task343 = std::make_shared<Task343>(tensor343, cindex);
      task342->add_dep(task343);
      task343->add_dep(task228);
      dedci_->add_task(task343);


      std::vector<std::shared_ptr<Tensor>> tensor344 = {I390, Gamma88};
      auto task344 = std::make_shared<Task344>(tensor344, cindex);
      task343->add_dep(task344);
      task344->add_dep(task228);
      dedci_->add_task(task344);

      task344->add_dep(task11);

      return make_tuple(queue_, energy_, correction_, density_, density1_, density2_, dedci_);
    };

  public:
    CAS_test(std::shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {
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
      std::shared_ptr<Queue> queue, energ, correct, dens2, dens1, Dens1, dec;
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

    double energy(std::shared_ptr<Queue> energ) {
      double en = 0.0;
      while (!energ->done()) {
        std::shared_ptr<Task> c = energ->next_compute();
        en += c->energy();
      }
      return en;
    }

    double correction(std::shared_ptr<Queue> correct) {
      double n = 0.0;
      while (!correct->done()) {
        std::shared_ptr<Task> c = correct->next_compute();
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

