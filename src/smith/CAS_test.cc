//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test.cc
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


#include <src/smith/CAS_test.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

tuple<shared_ptr<Queue>, shared_ptr<Queue>, shared_ptr<Queue>,  shared_ptr<Queue>,  shared_ptr<Queue>, shared_ptr<Queue>, shared_ptr<Queue>>
  CAS_test::CAS_test::make_queue_() {

  auto queue_ = make_shared<Queue>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  vector<shared_ptr<Tensor>> tensor0 = {r};
  auto task0 = make_shared<Task0>(tensor0);
  queue_->add_task(task0);

  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index, false);
  vector<shared_ptr<Tensor>> tensor1 = {Gamma0, rdm1_, rdm2_, rdm3_, rdm4_, f1_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  task1->add_dep(task0);
  queue_->add_task(task1);

  vector<IndexRange> Gamma1_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma1 = make_shared<Tensor>(Gamma1_index, false);
  vector<shared_ptr<Tensor>> tensor2 = {Gamma1, rdm1_, rdm2_, rdm3_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  task2->add_dep(task0);
  queue_->add_task(task2);

  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index, false);
  vector<shared_ptr<Tensor>> tensor3 = {Gamma3, rdm1_, rdm2_, rdm3_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  task3->add_dep(task0);
  queue_->add_task(task3);

  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index, false);
  vector<shared_ptr<Tensor>> tensor4 = {Gamma5, rdm1_, rdm2_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  task4->add_dep(task0);
  queue_->add_task(task4);

  vector<IndexRange> Gamma13_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma13 = make_shared<Tensor>(Gamma13_index, false);
  vector<shared_ptr<Tensor>> tensor5 = {Gamma13, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  task5->add_dep(task0);
  queue_->add_task(task5);

  vector<IndexRange> Gamma15_index = {active_, active_, active_, active_};
  auto Gamma15 = make_shared<Tensor>(Gamma15_index, false);
  vector<shared_ptr<Tensor>> tensor6 = {Gamma15, rdm1_, rdm2_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  task6->add_dep(task0);
  queue_->add_task(task6);

  vector<IndexRange> Gamma17_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma17 = make_shared<Tensor>(Gamma17_index, false);
  vector<shared_ptr<Tensor>> tensor7 = {Gamma17, rdm1deriv_, rdm2deriv_, rdm3deriv_, rdm4deriv_, f1_};
  auto task7 = make_shared<Task7>(tensor7, cindex);
  task7->add_dep(task0);
  queue_->add_task(task7);

  vector<IndexRange> Gamma18_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma18 = make_shared<Tensor>(Gamma18_index, false);
  vector<shared_ptr<Tensor>> tensor8 = {Gamma18, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task8 = make_shared<Task8>(tensor8, cindex);
  task8->add_dep(task0);
  queue_->add_task(task8);

  vector<IndexRange> Gamma23_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma23 = make_shared<Tensor>(Gamma23_index, false);
  vector<shared_ptr<Tensor>> tensor9 = {Gamma23, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task9 = make_shared<Task9>(tensor9, cindex);
  task9->add_dep(task0);
  queue_->add_task(task9);

  vector<IndexRange> Gamma27_index = {ci_, active_, active_, active_, active_};
  auto Gamma27 = make_shared<Tensor>(Gamma27_index, false);
  vector<shared_ptr<Tensor>> tensor10 = {Gamma27, rdm1deriv_, rdm2deriv_};
  auto task10 = make_shared<Task10>(tensor10, cindex);
  task10->add_dep(task0);
  queue_->add_task(task10);

  vector<IndexRange> Gamma25_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma25 = make_shared<Tensor>(Gamma25_index, false);
  vector<shared_ptr<Tensor>> tensor11 = {Gamma25, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task11 = make_shared<Task11>(tensor11, cindex);
  task11->add_dep(task0);
  queue_->add_task(task11);

  vector<IndexRange> Gamma28_index = {ci_, active_, active_, active_, active_};
  auto Gamma28 = make_shared<Tensor>(Gamma28_index, false);
  vector<shared_ptr<Tensor>> tensor12 = {Gamma28, rdm1deriv_, rdm2deriv_};
  auto task12 = make_shared<Task12>(tensor12, cindex);
  task12->add_dep(task0);
  queue_->add_task(task12);

  vector<IndexRange> I0_index = {active_, active_, active_, closed_};
  auto I0 = make_shared<Tensor>(I0_index, false);
  vector<shared_ptr<Tensor>> tensor13 = {r, I0};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  task13->add_dep(task0);
  queue_->add_task(task13);


  vector<IndexRange> I1_index = {active_, active_, active_, active_, active_, active_};
  auto I1 = make_shared<Tensor>(I1_index, false);
  vector<shared_ptr<Tensor>> tensor14 = {I0, t2, I1};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  task13->add_dep(task14);
  task14->add_dep(task0);
  queue_->add_task(task14);


  vector<shared_ptr<Tensor>> tensor15 = {I1, Gamma0};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  task14->add_dep(task15);
  task15->add_dep(task0);
  queue_->add_task(task15);

  task15->add_dep(task1);

  vector<IndexRange> I3_index = {active_, active_, active_, closed_};
  auto I3 = make_shared<Tensor>(I3_index, false);
  vector<shared_ptr<Tensor>> tensor16 = {I0, f1_, I3};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  task13->add_dep(task16);
  task16->add_dep(task0);
  queue_->add_task(task16);


  vector<IndexRange> I4_index = {active_, active_, active_, active_, active_, active_};
  auto I4 = make_shared<Tensor>(I4_index, false);
  vector<shared_ptr<Tensor>> tensor17 = {I3, t2, I4};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  task16->add_dep(task17);
  task17->add_dep(task0);
  queue_->add_task(task17);


  vector<shared_ptr<Tensor>> tensor18 = {I4, Gamma1};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  task17->add_dep(task18);
  task18->add_dep(task0);
  queue_->add_task(task18);

  task18->add_dep(task2);

  vector<IndexRange> I6_index = {active_, active_, active_, active_, active_, active_};
  auto I6 = make_shared<Tensor>(I6_index, false);
  vector<shared_ptr<Tensor>> tensor19 = {I0, t2, I6};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  task13->add_dep(task19);
  task19->add_dep(task0);
  queue_->add_task(task19);


  vector<shared_ptr<Tensor>> tensor20 = {I6, Gamma1};
  auto task20 = make_shared<Task20>(tensor20, pindex, this->e0_);
  task19->add_dep(task20);
  task20->add_dep(task0);
  queue_->add_task(task20);

  task20->add_dep(task2);

  vector<IndexRange> I8_index = {active_, active_, active_, active_, active_, active_};
  auto I8 = make_shared<Tensor>(I8_index, false);
  vector<shared_ptr<Tensor>> tensor21 = {I0, v2_, I8};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  task13->add_dep(task21);
  task21->add_dep(task0);
  queue_->add_task(task21);


  vector<shared_ptr<Tensor>> tensor22 = {I8, Gamma3};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  task21->add_dep(task22);
  task22->add_dep(task0);
  queue_->add_task(task22);

  task22->add_dep(task3);

  vector<IndexRange> I10_index = {active_, active_, active_, active_, active_, active_};
  auto I10 = make_shared<Tensor>(I10_index, false);
  vector<shared_ptr<Tensor>> tensor23 = {I0, v2_, I10};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  task13->add_dep(task23);
  task23->add_dep(task0);
  queue_->add_task(task23);


  vector<shared_ptr<Tensor>> tensor24 = {I10, Gamma1};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  task23->add_dep(task24);
  task24->add_dep(task0);
  queue_->add_task(task24);

  task24->add_dep(task2);

  vector<IndexRange> I12_index = {active_, active_, active_, active_};
  auto I12 = make_shared<Tensor>(I12_index, false);
  vector<shared_ptr<Tensor>> tensor25 = {I0, h1_, I12};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  task13->add_dep(task25);
  task25->add_dep(task0);
  queue_->add_task(task25);


  vector<shared_ptr<Tensor>> tensor26 = {I12, Gamma5};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  task25->add_dep(task26);
  task26->add_dep(task0);
  queue_->add_task(task26);

  task26->add_dep(task4);

  auto energy_ = make_shared<Queue>();
  vector<IndexRange> I13_index;
  auto I13 = make_shared<Tensor>(I13_index, false);
  vector<IndexRange> I14_index = {active_, active_, active_, closed_};
  auto I14 = make_shared<Tensor>(I14_index, false);
  vector<shared_ptr<Tensor>> tensor27 = {I13, t2, I14};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  energy_->add_task(task27);


  vector<IndexRange> I15_index = {active_, active_, active_, active_, active_, active_};
  auto I15 = make_shared<Tensor>(I15_index, false);
  vector<shared_ptr<Tensor>> tensor28 = {I14, t2, I15};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  task27->add_dep(task28);
  energy_->add_task(task28);


  vector<shared_ptr<Tensor>> tensor29 = {I15, Gamma0};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  task28->add_dep(task29);
  energy_->add_task(task29);

  task29->add_dep(task1);

  vector<IndexRange> I18_index = {active_, active_, active_, closed_};
  auto I18 = make_shared<Tensor>(I18_index, false);
  vector<shared_ptr<Tensor>> tensor30 = {I14, f1_, I18};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  task27->add_dep(task30);
  energy_->add_task(task30);


  vector<IndexRange> I19_index = {active_, active_, active_, active_, active_, active_};
  auto I19 = make_shared<Tensor>(I19_index, false);
  vector<shared_ptr<Tensor>> tensor31 = {I18, t2, I19};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  task30->add_dep(task31);
  energy_->add_task(task31);


  vector<shared_ptr<Tensor>> tensor32 = {I19, Gamma1};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  task31->add_dep(task32);
  energy_->add_task(task32);

  task32->add_dep(task2);

  vector<IndexRange> I22_index = {active_, active_, active_, active_, active_, active_};
  auto I22 = make_shared<Tensor>(I22_index, false);
  vector<shared_ptr<Tensor>> tensor33 = {I14, t2, I22};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  task27->add_dep(task33);
  energy_->add_task(task33);


  vector<shared_ptr<Tensor>> tensor34 = {I22, Gamma1};
  auto task34 = make_shared<Task34>(tensor34, pindex, this->e0_);
  task33->add_dep(task34);
  energy_->add_task(task34);

  task34->add_dep(task2);

  vector<IndexRange> I25_index = {active_, active_, active_, active_, active_, active_};
  auto I25 = make_shared<Tensor>(I25_index, false);
  vector<shared_ptr<Tensor>> tensor35 = {I14, v2_, I25};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  task27->add_dep(task35);
  energy_->add_task(task35);


  vector<shared_ptr<Tensor>> tensor36 = {I25, Gamma3};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  task35->add_dep(task36);
  energy_->add_task(task36);

  task36->add_dep(task3);

  vector<IndexRange> I28_index = {active_, active_, active_, active_, active_, active_};
  auto I28 = make_shared<Tensor>(I28_index, false);
  vector<shared_ptr<Tensor>> tensor37 = {I14, v2_, I28};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  task27->add_dep(task37);
  energy_->add_task(task37);


  vector<shared_ptr<Tensor>> tensor38 = {I28, Gamma1};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  task37->add_dep(task38);
  energy_->add_task(task38);

  task38->add_dep(task2);

  vector<IndexRange> I31_index = {active_, active_, active_, active_};
  auto I31 = make_shared<Tensor>(I31_index, false);
  vector<shared_ptr<Tensor>> tensor39 = {I14, h1_, I31};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  task27->add_dep(task39);
  energy_->add_task(task39);


  vector<shared_ptr<Tensor>> tensor40 = {I31, Gamma5};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  task39->add_dep(task40);
  energy_->add_task(task40);

  task40->add_dep(task4);

  auto correction_ = make_shared<Queue>();
  vector<IndexRange> I32_index;
  auto I32 = make_shared<Tensor>(I32_index, false);
  vector<IndexRange> I33_index = {active_, active_, active_, closed_};
  auto I33 = make_shared<Tensor>(I33_index, false);
  vector<shared_ptr<Tensor>> tensor41 = {I32, t2, I33};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  correction_->add_task(task41);


  vector<IndexRange> I34_index = {active_, active_, active_, active_, active_, active_};
  auto I34 = make_shared<Tensor>(I34_index, false);
  vector<shared_ptr<Tensor>> tensor42 = {I33, t2, I34};
  auto task42 = make_shared<Task42>(tensor42, pindex);
  task41->add_dep(task42);
  correction_->add_task(task42);


  vector<shared_ptr<Tensor>> tensor43 = {I34, Gamma1};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  task42->add_dep(task43);
  correction_->add_task(task43);

  task43->add_dep(task2);

  auto density_ = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor44 = {den2};
  auto task44 = make_shared<Task44>(tensor44);
  density_->add_task(task44);

  vector<IndexRange> I35_index = {active_, active_};
  auto I35 = make_shared<Tensor>(I35_index, false);
  vector<shared_ptr<Tensor>> tensor45 = {den2, I35};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  task45->add_dep(task44);
  density_->add_task(task45);


  vector<IndexRange> I36_index = {active_, active_, active_, active_, active_, closed_};
  auto I36 = make_shared<Tensor>(I36_index, false);
  vector<shared_ptr<Tensor>> tensor46 = {I35, t2, I36};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  task45->add_dep(task46);
  task46->add_dep(task44);
  density_->add_task(task46);


  vector<IndexRange> I37_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto I37 = make_shared<Tensor>(I37_index, false);
  vector<shared_ptr<Tensor>> tensor47 = {I36, t2, I37};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  task46->add_dep(task47);
  task47->add_dep(task44);
  density_->add_task(task47);


  vector<shared_ptr<Tensor>> tensor48 = {I37, Gamma13};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  task47->add_dep(task48);
  task48->add_dep(task44);
  density_->add_task(task48);

  task48->add_dep(task5);

  vector<IndexRange> I38_index = {closed_, closed_};
  auto I38 = make_shared<Tensor>(I38_index, false);
  vector<shared_ptr<Tensor>> tensor49 = {den2, I38};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  task49->add_dep(task44);
  density_->add_task(task49);


  vector<IndexRange> I39_index = {active_, active_, active_, closed_};
  auto I39 = make_shared<Tensor>(I39_index, false);
  vector<shared_ptr<Tensor>> tensor50 = {I38, t2, I39};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  task49->add_dep(task50);
  task50->add_dep(task44);
  density_->add_task(task50);


  vector<IndexRange> I40_index = {active_, active_, active_, active_, active_, active_};
  auto I40 = make_shared<Tensor>(I40_index, false);
  vector<shared_ptr<Tensor>> tensor51 = {I39, t2, I40};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  task50->add_dep(task51);
  task51->add_dep(task44);
  density_->add_task(task51);


  vector<shared_ptr<Tensor>> tensor52 = {I40, Gamma1};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  task51->add_dep(task52);
  task52->add_dep(task44);
  density_->add_task(task52);

  task52->add_dep(task2);

  auto density1_ = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor53 = {den1};
  auto task53 = make_shared<Task53>(tensor53);
  density1_->add_task(task53);

  vector<IndexRange> I41_index = {active_, closed_};
  auto I41 = make_shared<Tensor>(I41_index, false);
  vector<shared_ptr<Tensor>> tensor54 = {den1, I41};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  task54->add_dep(task53);
  density1_->add_task(task54);


  vector<IndexRange> I42_index = {active_, active_, active_, active_};
  auto I42 = make_shared<Tensor>(I42_index, false);
  vector<shared_ptr<Tensor>> tensor55 = {I41, t2, I42};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  task54->add_dep(task55);
  task55->add_dep(task53);
  density1_->add_task(task55);


  vector<shared_ptr<Tensor>> tensor56 = {I42, Gamma15};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  task55->add_dep(task56);
  task56->add_dep(task53);
  density1_->add_task(task56);

  task56->add_dep(task6);

  auto density2_ = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor57 = {Den1};
  auto task57 = make_shared<Task57>(tensor57);
  density2_->add_task(task57);

  vector<IndexRange> I43_index = {active_, active_, active_, closed_};
  auto I43 = make_shared<Tensor>(I43_index, false);
  vector<shared_ptr<Tensor>> tensor58 = {Den1, I43};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  task58->add_dep(task57);
  density2_->add_task(task58);


  vector<IndexRange> I44_index = {active_, active_, active_, active_, active_, active_};
  auto I44 = make_shared<Tensor>(I44_index, false);
  vector<shared_ptr<Tensor>> tensor59 = {I43, t2, I44};
  auto task59 = make_shared<Task59>(tensor59, pindex);
  task58->add_dep(task59);
  task59->add_dep(task57);
  density2_->add_task(task59);


  vector<shared_ptr<Tensor>> tensor60 = {I44, Gamma1};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  task59->add_dep(task60);
  task60->add_dep(task57);
  density2_->add_task(task60);

  task60->add_dep(task2);

  auto dedci_ = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor61 = {deci};
  auto task61 = make_shared<Task61>(tensor61);
  dedci_->add_task(task61);

  vector<IndexRange> I45_index = {ci_};
  auto I45 = make_shared<Tensor>(I45_index, false);
  vector<shared_ptr<Tensor>> tensor62 = {deci, I45};
  auto task62 = make_shared<Task62>(tensor62, cindex);
  task62->add_dep(task61);
  dedci_->add_task(task62);


  vector<IndexRange> I46_index = {ci_, active_, active_, active_, closed_};
  auto I46 = make_shared<Tensor>(I46_index, false);
  vector<shared_ptr<Tensor>> tensor63 = {I45, t2, I46};
  auto task63 = make_shared<Task63>(tensor63, cindex);
  task62->add_dep(task63);
  task63->add_dep(task61);
  dedci_->add_task(task63);


  vector<IndexRange> I47_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I47 = make_shared<Tensor>(I47_index, false);
  vector<shared_ptr<Tensor>> tensor64 = {I46, t2, I47};
  auto task64 = make_shared<Task64>(tensor64, cindex);
  task63->add_dep(task64);
  task64->add_dep(task61);
  dedci_->add_task(task64);


  vector<shared_ptr<Tensor>> tensor65 = {I47, Gamma17};
  auto task65 = make_shared<Task65>(tensor65, cindex);
  task64->add_dep(task65);
  task65->add_dep(task61);
  dedci_->add_task(task65);

  task65->add_dep(task7);

  vector<IndexRange> I50_index = {ci_, active_, active_, active_, closed_};
  auto I50 = make_shared<Tensor>(I50_index, false);
  vector<shared_ptr<Tensor>> tensor66 = {I46, f1_, I50};
  auto task66 = make_shared<Task66>(tensor66, cindex);
  task63->add_dep(task66);
  task66->add_dep(task61);
  dedci_->add_task(task66);


  vector<IndexRange> I51_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I51 = make_shared<Tensor>(I51_index, false);
  vector<shared_ptr<Tensor>> tensor67 = {I50, t2, I51};
  auto task67 = make_shared<Task67>(tensor67, cindex);
  task66->add_dep(task67);
  task67->add_dep(task61);
  dedci_->add_task(task67);


  vector<shared_ptr<Tensor>> tensor68 = {I51, Gamma18};
  auto task68 = make_shared<Task68>(tensor68, cindex);
  task67->add_dep(task68);
  task68->add_dep(task61);
  dedci_->add_task(task68);

  task68->add_dep(task8);

  vector<IndexRange> I61_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I61 = make_shared<Tensor>(I61_index, false);
  vector<shared_ptr<Tensor>> tensor69 = {I46, t2, I61};
  auto task69 = make_shared<Task69>(tensor69, cindex);
  task63->add_dep(task69);
  task69->add_dep(task61);
  dedci_->add_task(task69);


  vector<shared_ptr<Tensor>> tensor70 = {I61, Gamma18};
  auto task70 = make_shared<Task70>(tensor70, cindex, this->e0_);
  task69->add_dep(task70);
  task70->add_dep(task61);
  dedci_->add_task(task70);

  task70->add_dep(task8);

  vector<IndexRange> I67_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I67 = make_shared<Tensor>(I67_index, false);
  vector<shared_ptr<Tensor>> tensor71 = {I46, v2_, I67};
  auto task71 = make_shared<Task71>(tensor71, cindex);
  task63->add_dep(task71);
  task71->add_dep(task61);
  dedci_->add_task(task71);


  vector<shared_ptr<Tensor>> tensor72 = {I67, Gamma23};
  auto task72 = make_shared<Task72>(tensor72, cindex);
  task71->add_dep(task72);
  task72->add_dep(task61);
  dedci_->add_task(task72);

  task72->add_dep(task9);

  vector<IndexRange> I70_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I70 = make_shared<Tensor>(I70_index, false);
  vector<shared_ptr<Tensor>> tensor73 = {I46, v2_, I70};
  auto task73 = make_shared<Task73>(tensor73, cindex);
  task63->add_dep(task73);
  task73->add_dep(task61);
  dedci_->add_task(task73);


  vector<shared_ptr<Tensor>> tensor74 = {I70, Gamma18};
  auto task74 = make_shared<Task74>(tensor74, cindex);
  task73->add_dep(task74);
  task74->add_dep(task61);
  dedci_->add_task(task74);

  task74->add_dep(task8);

  vector<IndexRange> I79_index = {ci_, active_, active_, active_, active_};
  auto I79 = make_shared<Tensor>(I79_index, false);
  vector<shared_ptr<Tensor>> tensor75 = {I46, h1_, I79};
  auto task75 = make_shared<Task75>(tensor75, cindex);
  task63->add_dep(task75);
  task75->add_dep(task61);
  dedci_->add_task(task75);


  vector<shared_ptr<Tensor>> tensor76 = {I79, Gamma27};
  auto task76 = make_shared<Task76>(tensor76, cindex);
  task75->add_dep(task76);
  task76->add_dep(task61);
  dedci_->add_task(task76);

  task76->add_dep(task10);

  vector<IndexRange> I53_index = {ci_, active_, active_, active_, closed_};
  auto I53 = make_shared<Tensor>(I53_index, false);
  vector<shared_ptr<Tensor>> tensor77 = {I45, t2, I53};
  auto task77 = make_shared<Task77>(tensor77, cindex);
  task62->add_dep(task77);
  task77->add_dep(task61);
  dedci_->add_task(task77);


  vector<IndexRange> I54_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I54 = make_shared<Tensor>(I54_index, false);
  vector<shared_ptr<Tensor>> tensor78 = {I53, t2, I54};
  auto task78 = make_shared<Task78>(tensor78, cindex);
  task77->add_dep(task78);
  task78->add_dep(task61);
  dedci_->add_task(task78);


  vector<shared_ptr<Tensor>> tensor79 = {I54, Gamma17};
  auto task79 = make_shared<Task79>(tensor79, cindex);
  task78->add_dep(task79);
  task79->add_dep(task61);
  dedci_->add_task(task79);

  task79->add_dep(task7);

  vector<IndexRange> I56_index = {ci_, active_, active_, active_, closed_};
  auto I56 = make_shared<Tensor>(I56_index, false);
  vector<shared_ptr<Tensor>> tensor80 = {I45, t2, I56};
  auto task80 = make_shared<Task80>(tensor80, cindex);
  task62->add_dep(task80);
  task80->add_dep(task61);
  dedci_->add_task(task80);


  vector<IndexRange> I57_index = {ci_, active_, active_, active_, closed_};
  auto I57 = make_shared<Tensor>(I57_index, false);
  vector<shared_ptr<Tensor>> tensor81 = {I56, f1_, I57};
  auto task81 = make_shared<Task81>(tensor81, cindex);
  task80->add_dep(task81);
  task81->add_dep(task61);
  dedci_->add_task(task81);


  vector<IndexRange> I58_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I58 = make_shared<Tensor>(I58_index, false);
  vector<shared_ptr<Tensor>> tensor82 = {I57, t2, I58};
  auto task82 = make_shared<Task82>(tensor82, cindex);
  task81->add_dep(task82);
  task82->add_dep(task61);
  dedci_->add_task(task82);


  vector<shared_ptr<Tensor>> tensor83 = {I58, Gamma18};
  auto task83 = make_shared<Task83>(tensor83, cindex);
  task82->add_dep(task83);
  task83->add_dep(task61);
  dedci_->add_task(task83);

  task83->add_dep(task8);

  vector<IndexRange> I64_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I64 = make_shared<Tensor>(I64_index, false);
  vector<shared_ptr<Tensor>> tensor84 = {I56, t2, I64};
  auto task84 = make_shared<Task84>(tensor84, cindex);
  task80->add_dep(task84);
  task84->add_dep(task61);
  dedci_->add_task(task84);


  vector<shared_ptr<Tensor>> tensor85 = {I64, Gamma18};
  auto task85 = make_shared<Task85>(tensor85, cindex, this->e0_);
  task84->add_dep(task85);
  task85->add_dep(task61);
  dedci_->add_task(task85);

  task85->add_dep(task8);

  vector<IndexRange> I73_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I73 = make_shared<Tensor>(I73_index, false);
  vector<shared_ptr<Tensor>> tensor86 = {I56, v2_, I73};
  auto task86 = make_shared<Task86>(tensor86, cindex);
  task80->add_dep(task86);
  task86->add_dep(task61);
  dedci_->add_task(task86);


  vector<shared_ptr<Tensor>> tensor87 = {I73, Gamma25};
  auto task87 = make_shared<Task87>(tensor87, cindex);
  task86->add_dep(task87);
  task87->add_dep(task61);
  dedci_->add_task(task87);

  task87->add_dep(task11);

  vector<IndexRange> I76_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto I76 = make_shared<Tensor>(I76_index, false);
  vector<shared_ptr<Tensor>> tensor88 = {I56, v2_, I76};
  auto task88 = make_shared<Task88>(tensor88, cindex);
  task80->add_dep(task88);
  task88->add_dep(task61);
  dedci_->add_task(task88);


  vector<shared_ptr<Tensor>> tensor89 = {I76, Gamma18};
  auto task89 = make_shared<Task89>(tensor89, cindex);
  task88->add_dep(task89);
  task89->add_dep(task61);
  dedci_->add_task(task89);

  task89->add_dep(task8);

  vector<IndexRange> I81_index = {ci_, active_, active_, active_, closed_};
  auto I81 = make_shared<Tensor>(I81_index, false);
  vector<shared_ptr<Tensor>> tensor90 = {I45, t2, I81};
  auto task90 = make_shared<Task90>(tensor90, cindex);
  task62->add_dep(task90);
  task90->add_dep(task61);
  dedci_->add_task(task90);


  vector<IndexRange> I82_index = {ci_, active_, active_, active_, active_};
  auto I82 = make_shared<Tensor>(I82_index, false);
  vector<shared_ptr<Tensor>> tensor91 = {I81, h1_, I82};
  auto task91 = make_shared<Task91>(tensor91, cindex);
  task90->add_dep(task91);
  task91->add_dep(task61);
  dedci_->add_task(task91);


  vector<shared_ptr<Tensor>> tensor92 = {I82, Gamma28};
  auto task92 = make_shared<Task92>(tensor92, cindex);
  task91->add_dep(task92);
  task92->add_dep(task61);
  dedci_->add_task(task92);

  task92->add_dep(task12);

  return make_tuple(queue_, energy_, correction_, density_, density1_, density2_, dedci_);
}

CAS_test::CAS_test::CAS_test(shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {
  this->eig_ = f1_->diag();
  t2 = v2_->clone();
  e0_ = this->e0();
  sigma_ = this->sigma();
  this->update_amplitude(t2, v2_, true);
  t2->scale(2.0);
  r = t2->clone();
  den1 = h1_->clone();
  den2 = h1_->clone();
  Den1 = v2_->clone();
  deci = rdm0deriv_->clone();
}

void SMITH::CAS_test::CAS_test::solve() {
  Timer timer;
  this->print_iteration();
  int iter = 0;
  shared_ptr<Queue> queue, energ, correct, dens2, dens1, Dens1, dec;
  for ( ; iter != ref_->maxiter(); ++iter) {
    tie(queue, energ, correct, dens2, dens1, Dens1, dec) = make_queue_();
    while (!queue->done())
      queue->next_compute();
    this->update_amplitude(t2, r);
    const double err = r->rms();
    r->zero();
    this->energy_ = accumulate(energ);
    this->print_iteration(iter, this->energy_, err);
    if (err < ref_->thresh()) break;
  }
  this->print_iteration(iter == ref_->maxiter());
  timer.tick_print("CASPT2 energy evaluation");

  correlated_norm_ = accumulate(correct);
  timer.tick_print("T1 norm evaluation");

  while (!dens2->done())
    dens2->next_compute();
  while (!dens1->done())
    dens1->next_compute();
  while (!Dens1->done())
    Dens1->next_compute();
  timer.tick_print("Correlated density matrix evaluation");

  while (!dec->done())
    dec->next_compute();
  timer.tick_print("CI derivative evaluation");
  cout << endl;

}
