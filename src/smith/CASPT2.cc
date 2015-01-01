//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

tuple<shared_ptr<Queue>, shared_ptr<Queue>, shared_ptr<Queue>,  shared_ptr<Queue>,  shared_ptr<Queue>, shared_ptr<Queue>, shared_ptr<Queue>>
  CASPT2::CASPT2::make_queue_() {

  auto queue_ = make_shared<Queue>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  vector<shared_ptr<Tensor>> tensor0 = {r};
  auto task0 = make_shared<Task0>(tensor0);
  queue_->add_task(task0);

  vector<IndexRange> Gamma0_index;
  auto Gamma0 = make_shared<Tensor>(Gamma0_index, false);
  vector<shared_ptr<Tensor>> tensor1 = {Gamma0, rdm1_, f1_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  task1->add_dep(task0);
  queue_->add_task(task1);

  vector<IndexRange> Gamma4_index = {active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index, false);
  vector<shared_ptr<Tensor>> tensor2 = {Gamma4, rdm1_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  task2->add_dep(task0);
  queue_->add_task(task2);

  vector<IndexRange> Gamma6_index = {ci_};
  auto Gamma6 = make_shared<Tensor>(Gamma6_index, false);
  vector<shared_ptr<Tensor>> tensor3 = {Gamma6, rdm1deriv_, f1_};
  auto task3 = make_shared<Task3>(tensor3, cindex);
  task3->add_dep(task0);
  queue_->add_task(task3);

  vector<IndexRange> I0_index = {closed_, virt_, closed_, virt_};
  auto I0 = make_shared<Tensor>(I0_index, false);
  vector<shared_ptr<Tensor>> tensor4 = {r, I0};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  task4->add_dep(task0);
  queue_->add_task(task4);


  vector<shared_ptr<Tensor>> tensor5 = {I0, t2, v2_};
  auto task5 = make_shared<Task5>(tensor5, pindex, this->e0_);
  task4->add_dep(task5);
  task5->add_dep(task0);
  queue_->add_task(task5);


  vector<IndexRange> I1_index = {closed_, virt_, closed_, virt_};
  auto I1 = make_shared<Tensor>(I1_index, false);
  vector<shared_ptr<Tensor>> tensor6 = {I0, Gamma0, I1};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  task4->add_dep(task6);
  task6->add_dep(task0);
  queue_->add_task(task6);

  task6->add_dep(task1);

  vector<shared_ptr<Tensor>> tensor7 = {I1, t2};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  task6->add_dep(task7);
  task7->add_dep(task0);
  queue_->add_task(task7);


  vector<IndexRange> I3_index;
  auto I3 = make_shared<Tensor>(I3_index, false);
  vector<shared_ptr<Tensor>> tensor8 = {I0, t2, I3};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  task4->add_dep(task8);
  task8->add_dep(task0);
  queue_->add_task(task8);


  vector<shared_ptr<Tensor>> tensor9 = {I3, Gamma0};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  task8->add_dep(task9);
  task9->add_dep(task0);
  queue_->add_task(task9);

  task9->add_dep(task1);

  vector<IndexRange> I4_index = {closed_, virt_, virt_, closed_};
  auto I4 = make_shared<Tensor>(I4_index, false);
  vector<shared_ptr<Tensor>> tensor10 = {r, I4};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  task10->add_dep(task0);
  queue_->add_task(task10);


  vector<IndexRange> I5_index = {closed_, virt_, closed_, virt_};
  auto I5 = make_shared<Tensor>(I5_index, false);
  vector<shared_ptr<Tensor>> tensor11 = {I4, f1_, I5};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  task10->add_dep(task11);
  task11->add_dep(task0);
  queue_->add_task(task11);


  vector<shared_ptr<Tensor>> tensor12 = {I5, t2};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  task11->add_dep(task12);
  task12->add_dep(task0);
  queue_->add_task(task12);


  vector<IndexRange> I7_index = {closed_, closed_};
  auto I7 = make_shared<Tensor>(I7_index, false);
  vector<shared_ptr<Tensor>> tensor13 = {I4, t2, I7};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  task10->add_dep(task13);
  task13->add_dep(task0);
  queue_->add_task(task13);


  vector<shared_ptr<Tensor>> tensor14 = {I7, f1_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  task13->add_dep(task14);
  task14->add_dep(task0);
  queue_->add_task(task14);


  vector<IndexRange> I9_index = {closed_, virt_, closed_, virt_};
  auto I9 = make_shared<Tensor>(I9_index, false);
  vector<shared_ptr<Tensor>> tensor15 = {I4, f1_, I9};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  task10->add_dep(task15);
  task15->add_dep(task0);
  queue_->add_task(task15);


  vector<shared_ptr<Tensor>> tensor16 = {I9, t2};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  task15->add_dep(task16);
  task16->add_dep(task0);
  queue_->add_task(task16);


  vector<IndexRange> I11_index = {virt_, virt_};
  auto I11 = make_shared<Tensor>(I11_index, false);
  vector<shared_ptr<Tensor>> tensor17 = {I4, t2, I11};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  task10->add_dep(task17);
  task17->add_dep(task0);
  queue_->add_task(task17);


  vector<shared_ptr<Tensor>> tensor18 = {I11, f1_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  task17->add_dep(task18);
  task18->add_dep(task0);
  queue_->add_task(task18);


  auto energy_ = make_shared<Queue>();
  vector<IndexRange> I16_index;
  auto I16 = make_shared<Tensor>(I16_index, false);
  vector<IndexRange> I17_index;
  auto I17 = make_shared<Tensor>(I17_index, false);
  vector<shared_ptr<Tensor>> tensor19 = {I16, Gamma0, I17};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  energy_->add_task(task19);

  task19->add_dep(task1);

  vector<IndexRange> I18_index = {virt_, closed_, virt_, closed_};
  auto I18 = make_shared<Tensor>(I18_index, false);
  vector<shared_ptr<Tensor>> tensor20 = {I17, t2, I18};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  task19->add_dep(task20);
  energy_->add_task(task20);


  vector<shared_ptr<Tensor>> tensor21 = {I18, t2};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  task20->add_dep(task21);
  energy_->add_task(task21);


  vector<IndexRange> I21_index = {closed_, virt_, closed_, virt_};
  auto I21 = make_shared<Tensor>(I21_index, false);
  vector<shared_ptr<Tensor>> tensor22 = {I17, t2, I21};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  task19->add_dep(task22);
  energy_->add_task(task22);


  vector<shared_ptr<Tensor>> tensor23 = {I21, t2};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  task22->add_dep(task23);
  energy_->add_task(task23);


  vector<IndexRange> I23_index = {closed_, closed_};
  auto I23 = make_shared<Tensor>(I23_index, false);
  vector<shared_ptr<Tensor>> tensor24 = {I16, f1_, I23};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  task19->add_dep(task24);
  energy_->add_task(task24);


  vector<IndexRange> I24_index = {virt_, closed_, virt_, closed_};
  auto I24 = make_shared<Tensor>(I24_index, false);
  vector<shared_ptr<Tensor>> tensor25 = {I23, t2, I24};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  task24->add_dep(task25);
  energy_->add_task(task25);


  vector<shared_ptr<Tensor>> tensor26 = {I24, t2};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  task25->add_dep(task26);
  energy_->add_task(task26);


  vector<IndexRange> I27_index = {virt_, closed_, virt_, closed_};
  auto I27 = make_shared<Tensor>(I27_index, false);
  vector<shared_ptr<Tensor>> tensor27 = {I23, t2, I27};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  task24->add_dep(task27);
  energy_->add_task(task27);


  vector<shared_ptr<Tensor>> tensor28 = {I27, t2};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  task27->add_dep(task28);
  energy_->add_task(task28);


  vector<IndexRange> I29_index = {virt_, virt_};
  auto I29 = make_shared<Tensor>(I29_index, false);
  vector<shared_ptr<Tensor>> tensor29 = {I16, f1_, I29};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  task19->add_dep(task29);
  energy_->add_task(task29);


  vector<IndexRange> I30_index = {virt_, closed_, virt_, closed_};
  auto I30 = make_shared<Tensor>(I30_index, false);
  vector<shared_ptr<Tensor>> tensor30 = {I29, t2, I30};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  task29->add_dep(task30);
  energy_->add_task(task30);


  vector<shared_ptr<Tensor>> tensor31 = {I30, t2};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  task30->add_dep(task31);
  energy_->add_task(task31);


  vector<IndexRange> I33_index = {virt_, closed_, virt_, closed_};
  auto I33 = make_shared<Tensor>(I33_index, false);
  vector<shared_ptr<Tensor>> tensor32 = {I29, t2, I33};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  task29->add_dep(task32);
  energy_->add_task(task32);


  vector<shared_ptr<Tensor>> tensor33 = {I33, t2};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  task32->add_dep(task33);
  energy_->add_task(task33);


  vector<IndexRange> I35_index = {closed_, virt_, closed_, virt_};
  auto I35 = make_shared<Tensor>(I35_index, false);
  vector<shared_ptr<Tensor>> tensor34 = {I16, t2, I35};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  task19->add_dep(task34);
  energy_->add_task(task34);


  vector<shared_ptr<Tensor>> tensor35 = {I35, t2, v2_};
  auto task35 = make_shared<Task35>(tensor35, pindex, this->e0_);
  task34->add_dep(task35);
  energy_->add_task(task35);


  vector<IndexRange> I37_index = {virt_, closed_, virt_, closed_};
  auto I37 = make_shared<Tensor>(I37_index, false);
  vector<shared_ptr<Tensor>> tensor36 = {I16, t2, I37};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  task19->add_dep(task36);
  energy_->add_task(task36);


  vector<shared_ptr<Tensor>> tensor37 = {I37, t2};
  auto task37 = make_shared<Task37>(tensor37, pindex, this->e0_);
  task36->add_dep(task37);
  energy_->add_task(task37);


  vector<IndexRange> I39_index = {virt_, closed_, virt_, closed_};
  auto I39 = make_shared<Tensor>(I39_index, false);
  vector<shared_ptr<Tensor>> tensor38 = {I16, v2_, I39};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  task19->add_dep(task38);
  energy_->add_task(task38);


  vector<shared_ptr<Tensor>> tensor39 = {I39, t2};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  task38->add_dep(task39);
  energy_->add_task(task39);


  auto correction_ = make_shared<Queue>();
  vector<IndexRange> I42_index;
  auto I42 = make_shared<Tensor>(I42_index, false);
  vector<IndexRange> I43_index = {closed_, virt_, closed_, virt_};
  auto I43 = make_shared<Tensor>(I43_index, false);
  vector<shared_ptr<Tensor>> tensor40 = {I42, t2, I43};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  correction_->add_task(task40);


  vector<shared_ptr<Tensor>> tensor41 = {I43, t2};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  task40->add_dep(task41);
  correction_->add_task(task41);


  auto density_ = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor42 = {den2};
  auto task42 = make_shared<Task42>(tensor42);
  density_->add_task(task42);

  vector<IndexRange> I46_index = {active_, active_};
  auto I46 = make_shared<Tensor>(I46_index, false);
  vector<shared_ptr<Tensor>> tensor43 = {den2, I46};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  task43->add_dep(task42);
  density_->add_task(task43);


  vector<IndexRange> I47_index;
  auto I47 = make_shared<Tensor>(I47_index, false);
  vector<shared_ptr<Tensor>> tensor44 = {I46, Gamma4, I47};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  task43->add_dep(task44);
  task44->add_dep(task42);
  density_->add_task(task44);

  task44->add_dep(task2);

  vector<IndexRange> I48_index = {virt_, closed_, virt_, closed_};
  auto I48 = make_shared<Tensor>(I48_index, false);
  vector<shared_ptr<Tensor>> tensor45 = {I47, t2, I48};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  task44->add_dep(task45);
  task45->add_dep(task42);
  density_->add_task(task45);


  vector<shared_ptr<Tensor>> tensor46 = {I48, t2};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  task45->add_dep(task46);
  task46->add_dep(task42);
  density_->add_task(task46);


  vector<IndexRange> I51_index = {virt_, closed_, virt_, closed_};
  auto I51 = make_shared<Tensor>(I51_index, false);
  vector<shared_ptr<Tensor>> tensor47 = {I47, t2, I51};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  task44->add_dep(task47);
  task47->add_dep(task42);
  density_->add_task(task47);


  vector<shared_ptr<Tensor>> tensor48 = {I51, t2};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  task47->add_dep(task48);
  task48->add_dep(task42);
  density_->add_task(task48);


  vector<IndexRange> I52_index = {closed_, closed_};
  auto I52 = make_shared<Tensor>(I52_index, false);
  vector<shared_ptr<Tensor>> tensor49 = {den2, I52};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  task49->add_dep(task42);
  density_->add_task(task49);


  vector<IndexRange> I53_index = {virt_, closed_, virt_, closed_};
  auto I53 = make_shared<Tensor>(I53_index, false);
  vector<shared_ptr<Tensor>> tensor50 = {I52, t2, I53};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  task49->add_dep(task50);
  task50->add_dep(task42);
  density_->add_task(task50);


  vector<shared_ptr<Tensor>> tensor51 = {I53, t2};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  task50->add_dep(task51);
  task51->add_dep(task42);
  density_->add_task(task51);


  vector<IndexRange> I55_index = {closed_, virt_, closed_, virt_};
  auto I55 = make_shared<Tensor>(I55_index, false);
  vector<shared_ptr<Tensor>> tensor52 = {I52, t2, I55};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  task49->add_dep(task52);
  task52->add_dep(task42);
  density_->add_task(task52);


  vector<shared_ptr<Tensor>> tensor53 = {I55, t2};
  auto task53 = make_shared<Task53>(tensor53, pindex);
  task52->add_dep(task53);
  task53->add_dep(task42);
  density_->add_task(task53);


  vector<IndexRange> I56_index = {virt_, virt_};
  auto I56 = make_shared<Tensor>(I56_index, false);
  vector<shared_ptr<Tensor>> tensor54 = {den2, I56};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  task54->add_dep(task42);
  density_->add_task(task54);


  vector<IndexRange> I57_index = {virt_, closed_, virt_, closed_};
  auto I57 = make_shared<Tensor>(I57_index, false);
  vector<shared_ptr<Tensor>> tensor55 = {I56, t2, I57};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  task54->add_dep(task55);
  task55->add_dep(task42);
  density_->add_task(task55);


  vector<shared_ptr<Tensor>> tensor56 = {I57, t2};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  task55->add_dep(task56);
  task56->add_dep(task42);
  density_->add_task(task56);


  vector<IndexRange> I59_index = {closed_, virt_, closed_, virt_};
  auto I59 = make_shared<Tensor>(I59_index, false);
  vector<shared_ptr<Tensor>> tensor57 = {I56, t2, I59};
  auto task57 = make_shared<Task57>(tensor57, pindex);
  task54->add_dep(task57);
  task57->add_dep(task42);
  density_->add_task(task57);


  vector<shared_ptr<Tensor>> tensor58 = {I59, t2};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  task57->add_dep(task58);
  task58->add_dep(task42);
  density_->add_task(task58);


  auto density1_ = make_shared<Queue>();
  auto density2_ = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor59 = {Den1};
  auto task59 = make_shared<Task59>(tensor59);
  density2_->add_task(task59);

  vector<IndexRange> I60_index = {closed_, virt_, closed_, virt_};
  auto I60 = make_shared<Tensor>(I60_index, false);
  vector<shared_ptr<Tensor>> tensor60 = {Den1, I60};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  task60->add_dep(task59);
  density2_->add_task(task60);


  vector<shared_ptr<Tensor>> tensor61 = {I60, t2};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  task60->add_dep(task61);
  task61->add_dep(task59);
  density2_->add_task(task61);


  auto dedci_ = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor62 = {deci};
  auto task62 = make_shared<Task62>(tensor62);
  dedci_->add_task(task62);

  vector<IndexRange> I62_index = {ci_};
  auto I62 = make_shared<Tensor>(I62_index, false);
  vector<shared_ptr<Tensor>> tensor63 = {deci, I62};
  auto task63 = make_shared<Task63>(tensor63, cindex);
  task63->add_dep(task62);
  dedci_->add_task(task63);


  vector<IndexRange> I63_index;
  auto I63 = make_shared<Tensor>(I63_index, false);
  vector<shared_ptr<Tensor>> tensor64 = {I62, Gamma6, I63};
  auto task64 = make_shared<Task64>(tensor64, cindex);
  task63->add_dep(task64);
  task64->add_dep(task62);
  dedci_->add_task(task64);

  task64->add_dep(task3);

  vector<IndexRange> I64_index = {virt_, closed_, virt_, closed_};
  auto I64 = make_shared<Tensor>(I64_index, false);
  vector<shared_ptr<Tensor>> tensor65 = {I63, t2, I64};
  auto task65 = make_shared<Task65>(tensor65, cindex);
  task64->add_dep(task65);
  task65->add_dep(task62);
  dedci_->add_task(task65);


  vector<shared_ptr<Tensor>> tensor66 = {I64, t2};
  auto task66 = make_shared<Task66>(tensor66, cindex);
  task65->add_dep(task66);
  task66->add_dep(task62);
  dedci_->add_task(task66);


  vector<IndexRange> I67_index = {virt_, closed_, virt_, closed_};
  auto I67 = make_shared<Tensor>(I67_index, false);
  vector<shared_ptr<Tensor>> tensor67 = {I63, t2, I67};
  auto task67 = make_shared<Task67>(tensor67, cindex);
  task64->add_dep(task67);
  task67->add_dep(task62);
  dedci_->add_task(task67);


  vector<shared_ptr<Tensor>> tensor68 = {I67, t2};
  auto task68 = make_shared<Task68>(tensor68, cindex);
  task67->add_dep(task68);
  task68->add_dep(task62);
  dedci_->add_task(task68);


  return make_tuple(queue_, energy_, correction_, density_, density1_, density2_, dedci_);
}

CASPT2::CASPT2::CASPT2(shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {
  this->eig_ = f1_->diag();
  t2 = v2_->clone();
  e0_ = this->e0();
  this->update_amplitude(t2, v2_, true);
  t2->scale(2.0);
  r = t2->clone();
  den1 = h1_->clone();
  den2 = h1_->clone();
  Den1 = v2_->clone();
  deci = rdm0deriv_->clone();
}

void SMITH::CASPT2::CASPT2::solve() {
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
