//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_density2qq.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density2q(const bool reset) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor740 = {Den1};
  auto task740 = make_shared<Task740>(tensor740, reset);
  density2q->add_task(task740);

  vector<IndexRange> I742_index = {closed_, closed_, active_, active_};
  auto I742 = make_shared<Tensor>(I742_index);
  vector<shared_ptr<Tensor>> tensor741 = {Den1, I742};
  auto task741 = make_shared<Task741>(tensor741, pindex);
  task741->add_dep(task740);
  density2q->add_task(task741);

  vector<IndexRange> I743_index = {closed_, active_, closed_, active_};
  auto I743 = make_shared<Tensor>(I743_index);
  vector<shared_ptr<Tensor>> tensor742 = {I742, Gamma92_(), I743};
  auto task742 = make_shared<Task742>(tensor742, pindex);
  task741->add_dep(task742);
  task742->add_dep(task740);
  density2q->add_task(task742);

  vector<shared_ptr<Tensor>> tensor743 = {I743, t2};
  auto task743 = make_shared<Task743>(tensor743, pindex);
  task742->add_dep(task743);
  task743->add_dep(task740);
  density2q->add_task(task743);

  vector<IndexRange> I744_index = {closed_, active_, active_, active_};
  auto I744 = make_shared<Tensor>(I744_index);
  vector<shared_ptr<Tensor>> tensor744 = {Den1, I744};
  auto task744 = make_shared<Task744>(tensor744, pindex);
  task744->add_dep(task740);
  density2q->add_task(task744);

  vector<IndexRange> I745_index = {active_, active_, closed_, active_};
  auto I745 = make_shared<Tensor>(I745_index);
  vector<shared_ptr<Tensor>> tensor745 = {I744, Gamma6_(), I745};
  auto task745 = make_shared<Task745>(tensor745, pindex);
  task744->add_dep(task745);
  task745->add_dep(task740);
  density2q->add_task(task745);

  vector<shared_ptr<Tensor>> tensor746 = {I745, t2};
  auto task746 = make_shared<Task746>(tensor746, pindex);
  task745->add_dep(task746);
  task746->add_dep(task740);
  density2q->add_task(task746);

  vector<IndexRange> I746_index = {closed_, virt_, closed_, active_};
  auto I746 = make_shared<Tensor>(I746_index);
  vector<shared_ptr<Tensor>> tensor747 = {Den1, I746};
  auto task747 = make_shared<Task747>(tensor747, pindex);
  task747->add_dep(task740);
  density2q->add_task(task747);

  vector<IndexRange> I747_index = {closed_, virt_, closed_, active_};
  auto I747 = make_shared<Tensor>(I747_index);
  vector<shared_ptr<Tensor>> tensor748 = {I746, Gamma16_(), I747};
  auto task748 = make_shared<Task748>(tensor748, pindex);
  task747->add_dep(task748);
  task748->add_dep(task740);
  density2q->add_task(task748);

  vector<shared_ptr<Tensor>> tensor749 = {I747, t2};
  auto task749 = make_shared<Task749>(tensor749, pindex);
  task748->add_dep(task749);
  task749->add_dep(task740);
  density2q->add_task(task749);

  vector<IndexRange> I750_index = {virt_, closed_, active_, active_};
  auto I750 = make_shared<Tensor>(I750_index);
  vector<shared_ptr<Tensor>> tensor750 = {Den1, I750};
  auto task750 = make_shared<Task750>(tensor750, pindex);
  task750->add_dep(task740);
  density2q->add_task(task750);

  vector<IndexRange> I751_index = {active_, virt_, closed_, active_};
  auto I751 = make_shared<Tensor>(I751_index);
  vector<shared_ptr<Tensor>> tensor751 = {I750, Gamma32_(), I751};
  auto task751 = make_shared<Task751>(tensor751, pindex);
  task750->add_dep(task751);
  task751->add_dep(task740);
  density2q->add_task(task751);

  vector<shared_ptr<Tensor>> tensor752 = {I751, t2};
  auto task752 = make_shared<Task752>(tensor752, pindex);
  task751->add_dep(task752);
  task752->add_dep(task740);
  density2q->add_task(task752);

  vector<IndexRange> I753_index = {closed_, virt_, active_, active_};
  auto I753 = make_shared<Tensor>(I753_index);
  vector<shared_ptr<Tensor>> tensor753 = {I750, Gamma35_(), I753};
  auto task753 = make_shared<Task753>(tensor753, pindex);
  task750->add_dep(task753);
  task753->add_dep(task740);
  density2q->add_task(task753);

  vector<shared_ptr<Tensor>> tensor754 = {I753, t2};
  auto task754 = make_shared<Task754>(tensor754, pindex);
  task753->add_dep(task754);
  task754->add_dep(task740);
  density2q->add_task(task754);

  vector<IndexRange> I754_index = {virt_, closed_, active_, active_};
  auto I754 = make_shared<Tensor>(I754_index);
  vector<shared_ptr<Tensor>> tensor755 = {Den1, I754};
  auto task755 = make_shared<Task755>(tensor755, pindex);
  task755->add_dep(task740);
  density2q->add_task(task755);

  vector<IndexRange> I755_index = {active_, virt_, closed_, active_};
  auto I755 = make_shared<Tensor>(I755_index);
  vector<shared_ptr<Tensor>> tensor756 = {I754, Gamma35_(), I755};
  auto task756 = make_shared<Task756>(tensor756, pindex);
  task755->add_dep(task756);
  task756->add_dep(task740);
  density2q->add_task(task756);

  vector<shared_ptr<Tensor>> tensor757 = {I755, t2};
  auto task757 = make_shared<Task757>(tensor757, pindex);
  task756->add_dep(task757);
  task757->add_dep(task740);
  density2q->add_task(task757);

  vector<IndexRange> I758_index = {virt_, active_, active_, active_};
  auto I758 = make_shared<Tensor>(I758_index);
  vector<shared_ptr<Tensor>> tensor758 = {Den1, I758};
  auto task758 = make_shared<Task758>(tensor758, pindex);
  task758->add_dep(task740);
  density2q->add_task(task758);

  vector<IndexRange> I759_index = {active_, virt_, active_, active_};
  auto I759 = make_shared<Tensor>(I759_index);
  vector<shared_ptr<Tensor>> tensor759 = {I758, Gamma59_(), I759};
  auto task759 = make_shared<Task759>(tensor759, pindex);
  task758->add_dep(task759);
  task759->add_dep(task740);
  density2q->add_task(task759);

  vector<shared_ptr<Tensor>> tensor760 = {I759, t2};
  auto task760 = make_shared<Task760>(tensor760, pindex);
  task759->add_dep(task760);
  task760->add_dep(task740);
  density2q->add_task(task760);

  vector<IndexRange> I760_index = {closed_, virt_, closed_, virt_};
  auto I760 = make_shared<Tensor>(I760_index);
  vector<shared_ptr<Tensor>> tensor761 = {Den1, I760};
  auto task761 = make_shared<Task761>(tensor761, pindex);
  task761->add_dep(task740);
  density2q->add_task(task761);

  vector<shared_ptr<Tensor>> tensor762 = {I760, t2};
  auto task762 = make_shared<Task762>(tensor762, pindex);
  task761->add_dep(task762);
  task762->add_dep(task740);
  density2q->add_task(task762);

  vector<IndexRange> I762_index = {virt_, closed_, virt_, active_};
  auto I762 = make_shared<Tensor>(I762_index);
  vector<shared_ptr<Tensor>> tensor763 = {Den1, I762};
  auto task763 = make_shared<Task763>(tensor763, pindex);
  task763->add_dep(task740);
  density2q->add_task(task763);

  vector<IndexRange> I763_index = {active_, virt_, closed_, virt_};
  auto I763 = make_shared<Tensor>(I763_index);
  vector<shared_ptr<Tensor>> tensor764 = {I762, Gamma38_(), I763};
  auto task764 = make_shared<Task764>(tensor764, pindex);
  task763->add_dep(task764);
  task764->add_dep(task740);
  density2q->add_task(task764);

  vector<shared_ptr<Tensor>> tensor765 = {I763, t2};
  auto task765 = make_shared<Task765>(tensor765, pindex);
  task764->add_dep(task765);
  task765->add_dep(task740);
  density2q->add_task(task765);

  vector<IndexRange> I766_index = {virt_, virt_, active_, active_};
  auto I766 = make_shared<Tensor>(I766_index);
  vector<shared_ptr<Tensor>> tensor766 = {Den1, I766};
  auto task766 = make_shared<Task766>(tensor766, pindex);
  task766->add_dep(task740);
  density2q->add_task(task766);

  vector<IndexRange> I767_index = {active_, virt_, active_, virt_};
  auto I767 = make_shared<Tensor>(I767_index);
  vector<shared_ptr<Tensor>> tensor767 = {I766, Gamma60_(), I767};
  auto task767 = make_shared<Task767>(tensor767, pindex);
  task766->add_dep(task767);
  task767->add_dep(task740);
  density2q->add_task(task767);

  vector<shared_ptr<Tensor>> tensor768 = {I767, t2};
  auto task768 = make_shared<Task768>(tensor768, pindex);
  task767->add_dep(task768);
  task768->add_dep(task740);
  density2q->add_task(task768);

  return density2q;
}


#endif
