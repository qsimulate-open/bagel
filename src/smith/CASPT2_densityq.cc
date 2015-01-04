//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_densityqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_densityq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto densityq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor42 = {den2};
  auto task42 = make_shared<Task42>(tensor42);
  densityq->add_task(task42);

  vector<IndexRange> I46_index = {active_, active_};
  auto I46 = make_shared<Tensor>(I46_index);
  vector<shared_ptr<Tensor>> tensor43 = {den2, I46};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  task43->add_dep(task42);
  densityq->add_task(task43);

  vector<IndexRange> I47_index;
  auto I47 = make_shared<Tensor>(I47_index);
  vector<shared_ptr<Tensor>> tensor44 = {I46, Gamma4_(), I47};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  task43->add_dep(task44);
  task44->add_dep(task42);
  densityq->add_task(task44);

  vector<IndexRange> I48_index = {closed_, virt_, closed_, virt_};
  auto I48 = make_shared<Tensor>(I48_index);
  vector<shared_ptr<Tensor>> tensor45 = {I47, t2, I48};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  task44->add_dep(task45);
  task45->add_dep(task42);
  densityq->add_task(task45);

  vector<shared_ptr<Tensor>> tensor46 = {I48, t2};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  task45->add_dep(task46);
  task46->add_dep(task42);
  densityq->add_task(task46);

  vector<IndexRange> I51_index = {virt_, closed_, virt_, closed_};
  auto I51 = make_shared<Tensor>(I51_index);
  vector<shared_ptr<Tensor>> tensor47 = {I47, t2, I51};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  task44->add_dep(task47);
  task47->add_dep(task42);
  densityq->add_task(task47);

  vector<shared_ptr<Tensor>> tensor48 = {I51, t2};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  task47->add_dep(task48);
  task48->add_dep(task42);
  densityq->add_task(task48);

  vector<IndexRange> I52_index = {closed_, closed_};
  auto I52 = make_shared<Tensor>(I52_index);
  vector<shared_ptr<Tensor>> tensor49 = {den2, I52};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  task49->add_dep(task42);
  densityq->add_task(task49);

  vector<IndexRange> I53_index = {virt_, closed_, virt_, closed_};
  auto I53 = make_shared<Tensor>(I53_index);
  vector<shared_ptr<Tensor>> tensor50 = {I52, t2, I53};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  task49->add_dep(task50);
  task50->add_dep(task42);
  densityq->add_task(task50);

  vector<shared_ptr<Tensor>> tensor51 = {I53, t2};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  task50->add_dep(task51);
  task51->add_dep(task42);
  densityq->add_task(task51);

  vector<IndexRange> I55_index = {virt_, closed_, virt_, closed_};
  auto I55 = make_shared<Tensor>(I55_index);
  vector<shared_ptr<Tensor>> tensor52 = {I52, t2, I55};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  task49->add_dep(task52);
  task52->add_dep(task42);
  densityq->add_task(task52);

  vector<shared_ptr<Tensor>> tensor53 = {I55, t2};
  auto task53 = make_shared<Task53>(tensor53, pindex);
  task52->add_dep(task53);
  task53->add_dep(task42);
  densityq->add_task(task53);

  vector<IndexRange> I56_index = {virt_, virt_};
  auto I56 = make_shared<Tensor>(I56_index);
  vector<shared_ptr<Tensor>> tensor54 = {den2, I56};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  task54->add_dep(task42);
  densityq->add_task(task54);

  vector<IndexRange> I57_index = {virt_, closed_, virt_, closed_};
  auto I57 = make_shared<Tensor>(I57_index);
  vector<shared_ptr<Tensor>> tensor55 = {I56, t2, I57};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  task54->add_dep(task55);
  task55->add_dep(task42);
  densityq->add_task(task55);

  vector<shared_ptr<Tensor>> tensor56 = {I57, t2};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  task55->add_dep(task56);
  task56->add_dep(task42);
  densityq->add_task(task56);

  vector<IndexRange> I59_index = {virt_, closed_, virt_, closed_};
  auto I59 = make_shared<Tensor>(I59_index);
  vector<shared_ptr<Tensor>> tensor57 = {I56, t2, I59};
  auto task57 = make_shared<Task57>(tensor57, pindex);
  task54->add_dep(task57);
  task57->add_dep(task42);
  densityq->add_task(task57);

  vector<shared_ptr<Tensor>> tensor58 = {I59, t2};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  task57->add_dep(task58);
  task58->add_dep(task42);
  densityq->add_task(task58);

  return densityq;
}


