//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_corrqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_corrq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto corrq = make_shared<Queue>();
  vector<IndexRange> I832_index;
  auto I832 = make_shared<Tensor>(I832_index);
  vector<IndexRange> I833_index = {active_, active_, active_, active_};
  auto I833 = make_shared<Tensor>(I833_index);
  vector<shared_ptr<Tensor>> tensor706 = {I832, Gamma94_(), I833};
  auto task706 = make_shared<Task706>(tensor706, pindex);
  corrq->add_task(task706);

  vector<IndexRange> I834_index = {closed_, active_, closed_, active_};
  auto I834 = make_shared<Tensor>(I834_index);
  vector<shared_ptr<Tensor>> tensor707 = {I833, t2, I834};
  auto task707 = make_shared<Task707>(tensor707, pindex);
  task706->add_dep(task707);
  corrq->add_task(task707);

  vector<shared_ptr<Tensor>> tensor708 = {I834, t2};
  auto task708 = make_shared<Task708>(tensor708, pindex);
  task707->add_dep(task708);
  corrq->add_task(task708);

  vector<IndexRange> I836_index = {active_, active_, active_, closed_};
  auto I836 = make_shared<Tensor>(I836_index);
  vector<shared_ptr<Tensor>> tensor709 = {I832, t2, I836};
  auto task709 = make_shared<Task709>(tensor709, pindex);
  task706->add_dep(task709);
  corrq->add_task(task709);

  vector<IndexRange> I837_index = {active_, active_, active_, active_, active_, active_};
  auto I837 = make_shared<Tensor>(I837_index);
  vector<shared_ptr<Tensor>> tensor710 = {I836, t2, I837};
  auto task710 = make_shared<Task710>(tensor710, pindex);
  task709->add_dep(task710);
  corrq->add_task(task710);

  vector<shared_ptr<Tensor>> tensor711 = {I837, Gamma6_()};
  auto task711 = make_shared<Task711>(tensor711, pindex);
  task710->add_dep(task711);
  corrq->add_task(task711);

  vector<IndexRange> I839_index = {active_, active_};
  auto I839 = make_shared<Tensor>(I839_index);
  vector<shared_ptr<Tensor>> tensor712 = {I832, Gamma16_(), I839};
  auto task712 = make_shared<Task712>(tensor712, pindex);
  task706->add_dep(task712);
  corrq->add_task(task712);

  vector<IndexRange> I840_index = {closed_, virt_, closed_, active_};
  auto I840 = make_shared<Tensor>(I840_index);
  vector<shared_ptr<Tensor>> tensor713 = {I839, t2, I840};
  auto task713 = make_shared<Task713>(tensor713, pindex);
  task712->add_dep(task713);
  corrq->add_task(task713);

  vector<shared_ptr<Tensor>> tensor714 = {I840, t2};
  auto task714 = make_shared<Task714>(tensor714, pindex);
  task713->add_dep(task714);
  corrq->add_task(task714);

  vector<IndexRange> I845_index = {active_, active_, active_, active_};
  auto I845 = make_shared<Tensor>(I845_index);
  vector<shared_ptr<Tensor>> tensor715 = {I832, Gamma32_(), I845};
  auto task715 = make_shared<Task715>(tensor715, pindex);
  task706->add_dep(task715);
  corrq->add_task(task715);

  vector<IndexRange> I846_index = {active_, virt_, closed_, active_};
  auto I846 = make_shared<Tensor>(I846_index);
  vector<shared_ptr<Tensor>> tensor716 = {I845, t2, I846};
  auto task716 = make_shared<Task716>(tensor716, pindex);
  task715->add_dep(task716);
  corrq->add_task(task716);

  vector<shared_ptr<Tensor>> tensor717 = {I846, t2};
  auto task717 = make_shared<Task717>(tensor717, pindex);
  task716->add_dep(task717);
  corrq->add_task(task717);

  vector<IndexRange> I848_index = {active_, active_, active_, active_};
  auto I848 = make_shared<Tensor>(I848_index);
  vector<shared_ptr<Tensor>> tensor718 = {I832, Gamma35_(), I848};
  auto task718 = make_shared<Task718>(tensor718, pindex);
  task706->add_dep(task718);
  corrq->add_task(task718);

  vector<IndexRange> I849_index = {closed_, virt_, active_, active_};
  auto I849 = make_shared<Tensor>(I849_index);
  vector<shared_ptr<Tensor>> tensor719 = {I848, t2, I849};
  auto task719 = make_shared<Task719>(tensor719, pindex);
  task718->add_dep(task719);
  corrq->add_task(task719);

  vector<shared_ptr<Tensor>> tensor720 = {I849, t2};
  auto task720 = make_shared<Task720>(tensor720, pindex);
  task719->add_dep(task720);
  corrq->add_task(task720);

  vector<IndexRange> I852_index = {active_, active_, virt_, closed_};
  auto I852 = make_shared<Tensor>(I852_index);
  vector<shared_ptr<Tensor>> tensor721 = {I848, t2, I852};
  auto task721 = make_shared<Task721>(tensor721, pindex);
  task718->add_dep(task721);
  corrq->add_task(task721);

  vector<shared_ptr<Tensor>> tensor722 = {I852, t2};
  auto task722 = make_shared<Task722>(tensor722, pindex);
  task721->add_dep(task722);
  corrq->add_task(task722);

  vector<IndexRange> I855_index = {closed_, virt_, active_, active_};
  auto I855 = make_shared<Tensor>(I855_index);
  vector<shared_ptr<Tensor>> tensor723 = {I848, t2, I855};
  auto task723 = make_shared<Task723>(tensor723, pindex);
  task718->add_dep(task723);
  corrq->add_task(task723);

  vector<shared_ptr<Tensor>> tensor724 = {I855, t2};
  auto task724 = make_shared<Task724>(tensor724, pindex);
  task723->add_dep(task724);
  corrq->add_task(task724);

  vector<IndexRange> I857_index = {active_, active_, active_, active_, active_, active_};
  auto I857 = make_shared<Tensor>(I857_index);
  vector<shared_ptr<Tensor>> tensor725 = {I832, Gamma59_(), I857};
  auto task725 = make_shared<Task725>(tensor725, pindex);
  task706->add_dep(task725);
  corrq->add_task(task725);

  vector<IndexRange> I858_index = {active_, virt_, active_, active_};
  auto I858 = make_shared<Tensor>(I858_index);
  vector<shared_ptr<Tensor>> tensor726 = {I857, t2, I858};
  auto task726 = make_shared<Task726>(tensor726, pindex);
  task725->add_dep(task726);
  corrq->add_task(task726);

  vector<shared_ptr<Tensor>> tensor727 = {I858, t2};
  auto task727 = make_shared<Task727>(tensor727, pindex);
  task726->add_dep(task727);
  corrq->add_task(task727);

  vector<IndexRange> I860_index = {virt_, closed_, virt_, closed_};
  auto I860 = make_shared<Tensor>(I860_index);
  vector<shared_ptr<Tensor>> tensor728 = {I832, t2, I860};
  auto task728 = make_shared<Task728>(tensor728, pindex);
  task706->add_dep(task728);
  corrq->add_task(task728);

  vector<shared_ptr<Tensor>> tensor729 = {I860, t2};
  auto task729 = make_shared<Task729>(tensor729, pindex);
  task728->add_dep(task729);
  corrq->add_task(task729);

  vector<IndexRange> I862_index = {closed_, virt_, closed_, virt_};
  auto I862 = make_shared<Tensor>(I862_index);
  vector<shared_ptr<Tensor>> tensor730 = {I832, t2, I862};
  auto task730 = make_shared<Task730>(tensor730, pindex);
  task706->add_dep(task730);
  corrq->add_task(task730);

  vector<shared_ptr<Tensor>> tensor731 = {I862, t2};
  auto task731 = make_shared<Task731>(tensor731, pindex);
  task730->add_dep(task731);
  corrq->add_task(task731);

  vector<IndexRange> I864_index = {active_, active_};
  auto I864 = make_shared<Tensor>(I864_index);
  vector<shared_ptr<Tensor>> tensor732 = {I832, Gamma38_(), I864};
  auto task732 = make_shared<Task732>(tensor732, pindex);
  task706->add_dep(task732);
  corrq->add_task(task732);

  vector<IndexRange> I865_index = {active_, virt_, closed_, virt_};
  auto I865 = make_shared<Tensor>(I865_index);
  vector<shared_ptr<Tensor>> tensor733 = {I864, t2, I865};
  auto task733 = make_shared<Task733>(tensor733, pindex);
  task732->add_dep(task733);
  corrq->add_task(task733);

  vector<shared_ptr<Tensor>> tensor734 = {I865, t2};
  auto task734 = make_shared<Task734>(tensor734, pindex);
  task733->add_dep(task734);
  corrq->add_task(task734);

  vector<IndexRange> I868_index = {virt_, closed_, virt_, active_};
  auto I868 = make_shared<Tensor>(I868_index);
  vector<shared_ptr<Tensor>> tensor735 = {I864, t2, I868};
  auto task735 = make_shared<Task735>(tensor735, pindex);
  task732->add_dep(task735);
  corrq->add_task(task735);

  vector<shared_ptr<Tensor>> tensor736 = {I868, t2};
  auto task736 = make_shared<Task736>(tensor736, pindex);
  task735->add_dep(task736);
  corrq->add_task(task736);

  vector<IndexRange> I870_index = {active_, active_, active_, active_};
  auto I870 = make_shared<Tensor>(I870_index);
  vector<shared_ptr<Tensor>> tensor737 = {I832, Gamma60_(), I870};
  auto task737 = make_shared<Task737>(tensor737, pindex);
  task706->add_dep(task737);
  corrq->add_task(task737);

  vector<IndexRange> I871_index = {active_, virt_, active_, virt_};
  auto I871 = make_shared<Tensor>(I871_index);
  vector<shared_ptr<Tensor>> tensor738 = {I870, t2, I871};
  auto task738 = make_shared<Task738>(tensor738, pindex);
  task737->add_dep(task738);
  corrq->add_task(task738);

  vector<shared_ptr<Tensor>> tensor739 = {I871, t2};
  auto task739 = make_shared<Task739>(tensor739, pindex);
  task738->add_dep(task739);
  corrq->add_task(task739);

  return corrq;
}


