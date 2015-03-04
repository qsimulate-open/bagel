//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_residualqq.cc
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/MRCI.h>
#include <src/smith/MRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_residualq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto residualq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor7 = {r};
  auto task7 = make_shared<Task7>(tensor7);
  residualq->add_task(task7);

  vector<IndexRange> I0_index = {active_, active_, virt_, virt_};
  auto I0 = make_shared<Tensor>(I0_index);
  vector<shared_ptr<Tensor>> tensor8 = {r, I0};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  task8->add_dep(task7);
  residualq->add_task(task8);

  vector<IndexRange> I1_index = {active_, active_, active_, active_};
  auto I1 = make_shared<Tensor>(I1_index);
  vector<shared_ptr<Tensor>> tensor9 = {I0, t2, I1};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  task8->add_dep(task9);
  task9->add_dep(task7);
  residualq->add_task(task9);

  vector<shared_ptr<Tensor>> tensor10 = {I1, Gamma0_()};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  task9->add_dep(task10);
  task10->add_dep(task7);
  residualq->add_task(task10);

  vector<IndexRange> I6_index = {active_, virt_, active_, virt_};
  auto I6 = make_shared<Tensor>(I6_index);
  vector<shared_ptr<Tensor>> tensor11 = {I0, Gamma2_(), I6};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  task8->add_dep(task11);
  task11->add_dep(task7);
  residualq->add_task(task11);

  vector<shared_ptr<Tensor>> tensor12 = {I6, t2};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  task11->add_dep(task12);
  task12->add_dep(task7);
  residualq->add_task(task12);

  vector<IndexRange> I20_index = {virt_, virt_, active_, active_};
  auto I20 = make_shared<Tensor>(I20_index);
  vector<shared_ptr<Tensor>> tensor13 = {I0, Gamma7_(), I20};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  task8->add_dep(task13);
  task13->add_dep(task7);
  residualq->add_task(task13);

  vector<IndexRange> I21_index = {virt_, virt_, virt_, virt_};
  auto I21 = make_shared<Tensor>(I21_index);
  vector<shared_ptr<Tensor>> tensor14 = {I20, t2, I21};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  task13->add_dep(task14);
  task14->add_dep(task7);
  residualq->add_task(task14);

  vector<shared_ptr<Tensor>> tensor15 = {I21, v2_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  task14->add_dep(task15);
  task15->add_dep(task7);
  residualq->add_task(task15);

  vector<IndexRange> I23_index = {active_, virt_, active_, virt_};
  auto I23 = make_shared<Tensor>(I23_index);
  vector<shared_ptr<Tensor>> tensor16 = {I0, Gamma8_(), I23};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  task8->add_dep(task16);
  task16->add_dep(task7);
  residualq->add_task(task16);

  vector<shared_ptr<Tensor>> tensor17 = {I23, t2};
  auto task17 = make_shared<Task17>(tensor17, pindex, this->e0_);
  task16->add_dep(task17);
  task17->add_dep(task7);
  residualq->add_task(task17);

  vector<IndexRange> I2_index = {virt_, virt_, active_, active_};
  auto I2 = make_shared<Tensor>(I2_index);
  vector<shared_ptr<Tensor>> tensor18 = {r, I2};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  task18->add_dep(task7);
  residualq->add_task(task18);

  vector<IndexRange> I3_index = {virt_, active_, virt_, active_};
  auto I3 = make_shared<Tensor>(I3_index);
  vector<shared_ptr<Tensor>> tensor19 = {I2, Gamma8_(), I3};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  task18->add_dep(task19);
  task19->add_dep(task7);
  residualq->add_task(task19);

  vector<IndexRange> I4_index = {virt_, virt_};
  auto I4 = make_shared<Tensor>(I4_index);
  vector<shared_ptr<Tensor>> tensor20 = {I3, t2, I4};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  task19->add_dep(task20);
  task20->add_dep(task7);
  residualq->add_task(task20);

  vector<shared_ptr<Tensor>> tensor21 = {I4, h1_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  task20->add_dep(task21);
  task21->add_dep(task7);
  residualq->add_task(task21);

  vector<IndexRange> I8_index = {virt_, virt_, active_, active_, active_, active_};
  auto I8 = make_shared<Tensor>(I8_index);
  vector<shared_ptr<Tensor>> tensor22 = {I2, v2_, I8};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  task18->add_dep(task22);
  task22->add_dep(task7);
  residualq->add_task(task22);

  vector<IndexRange> I9_index = {active_, virt_, active_, virt_};
  auto I9 = make_shared<Tensor>(I9_index);
  vector<shared_ptr<Tensor>> tensor23 = {I8, Gamma3_(), I9};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  task22->add_dep(task23);
  task23->add_dep(task7);
  residualq->add_task(task23);

  vector<shared_ptr<Tensor>> tensor24 = {I9, t2};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  task23->add_dep(task24);
  task24->add_dep(task7);
  residualq->add_task(task24);

  vector<IndexRange> I11_index = {active_, active_, virt_, active_, virt_, active_};
  auto I11 = make_shared<Tensor>(I11_index);
  vector<shared_ptr<Tensor>> tensor25 = {I2, Gamma4_(), I11};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  task18->add_dep(task25);
  task25->add_dep(task7);
  residualq->add_task(task25);

  vector<IndexRange> I12_index = {virt_, active_, active_, virt_};
  auto I12 = make_shared<Tensor>(I12_index);
  vector<shared_ptr<Tensor>> tensor26 = {I11, t2, I12};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  task25->add_dep(task26);
  task26->add_dep(task7);
  residualq->add_task(task26);

  vector<shared_ptr<Tensor>> tensor27 = {I12, v2_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  task26->add_dep(task27);
  task27->add_dep(task7);
  residualq->add_task(task27);

  vector<IndexRange> I14_index = {active_, virt_, active_, active_, virt_, active_};
  auto I14 = make_shared<Tensor>(I14_index);
  vector<shared_ptr<Tensor>> tensor28 = {I2, Gamma5_(), I14};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  task18->add_dep(task28);
  task28->add_dep(task7);
  residualq->add_task(task28);

  vector<IndexRange> I15_index = {active_, virt_, virt_, active_};
  auto I15 = make_shared<Tensor>(I15_index);
  vector<shared_ptr<Tensor>> tensor29 = {I14, t2, I15};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  task28->add_dep(task29);
  task29->add_dep(task7);
  residualq->add_task(task29);

  vector<shared_ptr<Tensor>> tensor30 = {I15, v2_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  task29->add_dep(task30);
  task30->add_dep(task7);
  residualq->add_task(task30);

  vector<IndexRange> I17_index = {active_, active_, virt_, active_, virt_, active_};
  auto I17 = make_shared<Tensor>(I17_index);
  vector<shared_ptr<Tensor>> tensor31 = {I2, Gamma3_(), I17};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  task18->add_dep(task31);
  task31->add_dep(task7);
  residualq->add_task(task31);

  vector<IndexRange> I18_index = {active_, active_, virt_, virt_};
  auto I18 = make_shared<Tensor>(I18_index);
  vector<shared_ptr<Tensor>> tensor32 = {I17, t2, I18};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  task31->add_dep(task32);
  task32->add_dep(task7);
  residualq->add_task(task32);

  vector<shared_ptr<Tensor>> tensor33 = {I18, v2_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  task32->add_dep(task33);
  task33->add_dep(task7);
  residualq->add_task(task33);

  return residualq;
}


#endif
