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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density2q() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor753 = {Den1};
  auto task753 = make_shared<Task753>(tensor753);
  density2q->add_task(task753);

  vector<IndexRange> I756_index = {closed_, closed_, active_, active_};
  auto I756 = make_shared<Tensor>(I756_index);
  vector<shared_ptr<Tensor>> tensor754 = {Den1, I756};
  auto task754 = make_shared<Task754>(tensor754, pindex);
  task754->add_dep(task753);
  density2q->add_task(task754);

  vector<IndexRange> I757_index = {closed_, active_, closed_, active_};
  auto I757 = make_shared<Tensor>(I757_index);
  vector<shared_ptr<Tensor>> tensor755 = {I756, Gamma94_(), I757};
  auto task755 = make_shared<Task755>(tensor755, pindex);
  task754->add_dep(task755);
  task755->add_dep(task753);
  density2q->add_task(task755);

  vector<shared_ptr<Tensor>> tensor756 = {I757, t2};
  auto task756 = make_shared<Task756>(tensor756, pindex);
  task755->add_dep(task756);
  task756->add_dep(task753);
  density2q->add_task(task756);

  vector<IndexRange> I758_index = {closed_, active_, active_, active_};
  auto I758 = make_shared<Tensor>(I758_index);
  vector<shared_ptr<Tensor>> tensor757 = {Den1, I758};
  auto task757 = make_shared<Task757>(tensor757, pindex);
  task757->add_dep(task753);
  density2q->add_task(task757);

  vector<IndexRange> I759_index = {active_, active_, closed_, active_};
  auto I759 = make_shared<Tensor>(I759_index);
  vector<shared_ptr<Tensor>> tensor758 = {I758, Gamma6_(), I759};
  auto task758 = make_shared<Task758>(tensor758, pindex);
  task757->add_dep(task758);
  task758->add_dep(task753);
  density2q->add_task(task758);

  vector<shared_ptr<Tensor>> tensor759 = {I759, t2};
  auto task759 = make_shared<Task759>(tensor759, pindex);
  task758->add_dep(task759);
  task759->add_dep(task753);
  density2q->add_task(task759);

  vector<IndexRange> I760_index = {closed_, virt_, closed_, active_};
  auto I760 = make_shared<Tensor>(I760_index);
  vector<shared_ptr<Tensor>> tensor760 = {Den1, I760};
  auto task760 = make_shared<Task760>(tensor760, pindex);
  task760->add_dep(task753);
  density2q->add_task(task760);

  vector<IndexRange> I761_index = {closed_, virt_, closed_, active_};
  auto I761 = make_shared<Tensor>(I761_index);
  vector<shared_ptr<Tensor>> tensor761 = {I760, Gamma16_(), I761};
  auto task761 = make_shared<Task761>(tensor761, pindex);
  task760->add_dep(task761);
  task761->add_dep(task753);
  density2q->add_task(task761);

  vector<shared_ptr<Tensor>> tensor762 = {I761, t2};
  auto task762 = make_shared<Task762>(tensor762, pindex);
  task761->add_dep(task762);
  task762->add_dep(task753);
  density2q->add_task(task762);

  vector<IndexRange> I764_index = {virt_, closed_, active_, active_};
  auto I764 = make_shared<Tensor>(I764_index);
  vector<shared_ptr<Tensor>> tensor763 = {Den1, I764};
  auto task763 = make_shared<Task763>(tensor763, pindex);
  task763->add_dep(task753);
  density2q->add_task(task763);

  vector<IndexRange> I765_index = {active_, virt_, closed_, active_};
  auto I765 = make_shared<Tensor>(I765_index);
  vector<shared_ptr<Tensor>> tensor764 = {I764, Gamma32_(), I765};
  auto task764 = make_shared<Task764>(tensor764, pindex);
  task763->add_dep(task764);
  task764->add_dep(task753);
  density2q->add_task(task764);

  vector<shared_ptr<Tensor>> tensor765 = {I765, t2};
  auto task765 = make_shared<Task765>(tensor765, pindex);
  task764->add_dep(task765);
  task765->add_dep(task753);
  density2q->add_task(task765);

  vector<IndexRange> I767_index = {closed_, virt_, active_, active_};
  auto I767 = make_shared<Tensor>(I767_index);
  vector<shared_ptr<Tensor>> tensor766 = {I764, Gamma35_(), I767};
  auto task766 = make_shared<Task766>(tensor766, pindex);
  task763->add_dep(task766);
  task766->add_dep(task753);
  density2q->add_task(task766);

  vector<shared_ptr<Tensor>> tensor767 = {I767, t2};
  auto task767 = make_shared<Task767>(tensor767, pindex);
  task766->add_dep(task767);
  task767->add_dep(task753);
  density2q->add_task(task767);

  vector<IndexRange> I768_index = {virt_, closed_, active_, active_};
  auto I768 = make_shared<Tensor>(I768_index);
  vector<shared_ptr<Tensor>> tensor768 = {Den1, I768};
  auto task768 = make_shared<Task768>(tensor768, pindex);
  task768->add_dep(task753);
  density2q->add_task(task768);

  vector<IndexRange> I769_index = {active_, virt_, closed_, active_};
  auto I769 = make_shared<Tensor>(I769_index);
  vector<shared_ptr<Tensor>> tensor769 = {I768, Gamma35_(), I769};
  auto task769 = make_shared<Task769>(tensor769, pindex);
  task768->add_dep(task769);
  task769->add_dep(task753);
  density2q->add_task(task769);

  vector<shared_ptr<Tensor>> tensor770 = {I769, t2};
  auto task770 = make_shared<Task770>(tensor770, pindex);
  task769->add_dep(task770);
  task770->add_dep(task753);
  density2q->add_task(task770);

  vector<IndexRange> I772_index = {virt_, active_, active_, active_};
  auto I772 = make_shared<Tensor>(I772_index);
  vector<shared_ptr<Tensor>> tensor771 = {Den1, I772};
  auto task771 = make_shared<Task771>(tensor771, pindex);
  task771->add_dep(task753);
  density2q->add_task(task771);

  vector<IndexRange> I773_index = {active_, virt_, active_, active_};
  auto I773 = make_shared<Tensor>(I773_index);
  vector<shared_ptr<Tensor>> tensor772 = {I772, Gamma59_(), I773};
  auto task772 = make_shared<Task772>(tensor772, pindex);
  task771->add_dep(task772);
  task772->add_dep(task753);
  density2q->add_task(task772);

  vector<shared_ptr<Tensor>> tensor773 = {I773, t2};
  auto task773 = make_shared<Task773>(tensor773, pindex);
  task772->add_dep(task773);
  task773->add_dep(task753);
  density2q->add_task(task773);

  vector<IndexRange> I774_index = {closed_, virt_, closed_, virt_};
  auto I774 = make_shared<Tensor>(I774_index);
  vector<shared_ptr<Tensor>> tensor774 = {Den1, I774};
  auto task774 = make_shared<Task774>(tensor774, pindex);
  task774->add_dep(task753);
  density2q->add_task(task774);

  vector<shared_ptr<Tensor>> tensor775 = {I774, t2};
  auto task775 = make_shared<Task775>(tensor775, pindex);
  task774->add_dep(task775);
  task775->add_dep(task753);
  density2q->add_task(task775);

  vector<IndexRange> I776_index = {virt_, closed_, virt_, active_};
  auto I776 = make_shared<Tensor>(I776_index);
  vector<shared_ptr<Tensor>> tensor776 = {Den1, I776};
  auto task776 = make_shared<Task776>(tensor776, pindex);
  task776->add_dep(task753);
  density2q->add_task(task776);

  vector<IndexRange> I777_index = {active_, virt_, closed_, virt_};
  auto I777 = make_shared<Tensor>(I777_index);
  vector<shared_ptr<Tensor>> tensor777 = {I776, Gamma38_(), I777};
  auto task777 = make_shared<Task777>(tensor777, pindex);
  task776->add_dep(task777);
  task777->add_dep(task753);
  density2q->add_task(task777);

  vector<shared_ptr<Tensor>> tensor778 = {I777, t2};
  auto task778 = make_shared<Task778>(tensor778, pindex);
  task777->add_dep(task778);
  task778->add_dep(task753);
  density2q->add_task(task778);

  vector<IndexRange> I780_index = {virt_, virt_, active_, active_};
  auto I780 = make_shared<Tensor>(I780_index);
  vector<shared_ptr<Tensor>> tensor779 = {Den1, I780};
  auto task779 = make_shared<Task779>(tensor779, pindex);
  task779->add_dep(task753);
  density2q->add_task(task779);

  vector<IndexRange> I781_index = {active_, virt_, active_, virt_};
  auto I781 = make_shared<Tensor>(I781_index);
  vector<shared_ptr<Tensor>> tensor780 = {I780, Gamma60_(), I781};
  auto task780 = make_shared<Task780>(tensor780, pindex);
  task779->add_dep(task780);
  task780->add_dep(task753);
  density2q->add_task(task780);

  vector<shared_ptr<Tensor>> tensor781 = {I781, t2};
  auto task781 = make_shared<Task781>(tensor781, pindex);
  task780->add_dep(task781);
  task781->add_dep(task753);
  density2q->add_task(task781);

  return density2q;
}


