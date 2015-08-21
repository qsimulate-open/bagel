//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_normqq.cc
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


#include <src/smith/RelMRCI.h>
#include <src/smith/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor755 = vector<shared_ptr<Tensor>>{n};
  auto task755 = make_shared<Task755>(tensor755, reset);
  normq->add_task(task755);

  vector<IndexRange> I1348_index = {active_, active_, closed_, closed_};
  auto I1348 = make_shared<Tensor>(I1348_index);
  auto tensor756 = vector<shared_ptr<Tensor>>{n, I1348};
  auto task756 = make_shared<Task756>(tensor756, pindex);
  task756->add_dep(task755);
  normq->add_task(task756);

  auto tensor757 = vector<shared_ptr<Tensor>>{I1348, t2, Gamma0_()};
  auto task757 = make_shared<Task757>(tensor757, pindex);
  task756->add_dep(task757);
  task757->add_dep(task755);
  normq->add_task(task757);

  vector<IndexRange> I1350_index = {closed_, active_, active_, active_};
  auto I1350 = make_shared<Tensor>(I1350_index);
  auto tensor758 = vector<shared_ptr<Tensor>>{n, I1350};
  auto task758 = make_shared<Task758>(tensor758, pindex);
  task758->add_dep(task755);
  normq->add_task(task758);

  auto tensor759 = vector<shared_ptr<Tensor>>{I1350, Gamma4_(), t2};
  auto task759 = make_shared<Task759>(tensor759, pindex);
  task758->add_dep(task759);
  task759->add_dep(task755);
  normq->add_task(task759);

  vector<IndexRange> I1352_index = {closed_, virt_, closed_, active_};
  auto I1352 = make_shared<Tensor>(I1352_index);
  auto tensor760 = vector<shared_ptr<Tensor>>{n, I1352};
  auto task760 = make_shared<Task760>(tensor760, pindex);
  task760->add_dep(task755);
  normq->add_task(task760);

  auto tensor761 = vector<shared_ptr<Tensor>>{I1352, Gamma11_(), t2};
  auto task761 = make_shared<Task761>(tensor761, pindex);
  task760->add_dep(task761);
  task761->add_dep(task755);
  normq->add_task(task761);

  auto tensor762 = vector<shared_ptr<Tensor>>{I1352, t2, Gamma11_()};
  auto task762 = make_shared<Task762>(tensor762, pindex);
  task760->add_dep(task762);
  task762->add_dep(task755);
  normq->add_task(task762);

  vector<IndexRange> I1356_index = {active_, active_, closed_, virt_};
  auto I1356 = make_shared<Tensor>(I1356_index);
  auto tensor763 = vector<shared_ptr<Tensor>>{n, I1356};
  auto task763 = make_shared<Task763>(tensor763, pindex);
  task763->add_dep(task755);
  normq->add_task(task763);

  auto tensor764 = vector<shared_ptr<Tensor>>{I1356, t2, Gamma24_()};
  auto task764 = make_shared<Task764>(tensor764, pindex);
  task763->add_dep(task764);
  task764->add_dep(task755);
  normq->add_task(task764);

  vector<IndexRange> I1358_index = {virt_, active_, active_, active_};
  auto I1358 = make_shared<Tensor>(I1358_index);
  auto tensor765 = vector<shared_ptr<Tensor>>{n, I1358};
  auto task765 = make_shared<Task765>(tensor765, pindex);
  task765->add_dep(task755);
  normq->add_task(task765);

  auto tensor766 = vector<shared_ptr<Tensor>>{I1358, Gamma32_(), t2};
  auto task766 = make_shared<Task766>(tensor766, pindex);
  task765->add_dep(task766);
  task766->add_dep(task755);
  normq->add_task(task766);

  shared_ptr<Tensor> I1360;
  if (diagonal) {
    vector<IndexRange> I1360_index = {closed_, virt_, closed_, virt_};
    I1360 = make_shared<Tensor>(I1360_index);
  }
  shared_ptr<Task767> task767;
  if (diagonal) {
    auto tensor767 = vector<shared_ptr<Tensor>>{n, I1360};
    task767 = make_shared<Task767>(tensor767, pindex);
    task767->add_dep(task755);
    normq->add_task(task767);
  }

  shared_ptr<Task768> task768;
  if (diagonal) {
    auto tensor768 = vector<shared_ptr<Tensor>>{I1360, t2};
    task768 = make_shared<Task768>(tensor768, pindex);
    task767->add_dep(task768);
    task768->add_dep(task755);
    normq->add_task(task768);
  }

  vector<IndexRange> I1362_index = {active_, virt_, closed_, virt_};
  auto I1362 = make_shared<Tensor>(I1362_index);
  auto tensor769 = vector<shared_ptr<Tensor>>{n, I1362};
  auto task769 = make_shared<Task769>(tensor769, pindex);
  task769->add_dep(task755);
  normq->add_task(task769);

  auto tensor770 = vector<shared_ptr<Tensor>>{I1362, t2, Gamma27_()};
  auto task770 = make_shared<Task770>(tensor770, pindex);
  task769->add_dep(task770);
  task770->add_dep(task755);
  normq->add_task(task770);

  auto tensor771 = vector<shared_ptr<Tensor>>{I1362, t2, Gamma27_()};
  auto task771 = make_shared<Task771>(tensor771, pindex);
  task769->add_dep(task771);
  task771->add_dep(task755);
  normq->add_task(task771);

  vector<IndexRange> I1366_index = {active_, active_, virt_, virt_};
  auto I1366 = make_shared<Tensor>(I1366_index);
  auto tensor772 = vector<shared_ptr<Tensor>>{n, I1366};
  auto task772 = make_shared<Task772>(tensor772, pindex);
  task772->add_dep(task755);
  normq->add_task(task772);

  auto tensor773 = vector<shared_ptr<Tensor>>{I1366, t2, Gamma33_()};
  auto task773 = make_shared<Task773>(tensor773, pindex);
  task772->add_dep(task773);
  task773->add_dep(task755);
  normq->add_task(task773);

  return normq;
}


#endif
