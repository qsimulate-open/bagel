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
  auto tensor970 = vector<shared_ptr<Tensor>>{n};
  auto task970 = make_shared<Task970>(tensor970, reset);
  normq->add_task(task970);

  vector<IndexRange> I1778_index = {active_, active_, closed_, closed_};
  auto I1778 = make_shared<Tensor>(I1778_index);
  auto tensor971 = vector<shared_ptr<Tensor>>{n, I1778};
  auto task971 = make_shared<Task971>(tensor971, pindex);
  task971->add_dep(task970);
  normq->add_task(task971);

  auto tensor972 = vector<shared_ptr<Tensor>>{I1778, t2, Gamma0_()};
  auto task972 = make_shared<Task972>(tensor972, pindex);
  task971->add_dep(task972);
  task972->add_dep(task970);
  normq->add_task(task972);

  vector<IndexRange> I1780_index = {closed_, active_, active_, active_};
  auto I1780 = make_shared<Tensor>(I1780_index);
  auto tensor973 = vector<shared_ptr<Tensor>>{n, I1780};
  auto task973 = make_shared<Task973>(tensor973, pindex);
  task973->add_dep(task970);
  normq->add_task(task973);

  auto tensor974 = vector<shared_ptr<Tensor>>{I1780, Gamma4_(), t2};
  auto task974 = make_shared<Task974>(tensor974, pindex);
  task973->add_dep(task974);
  task974->add_dep(task970);
  normq->add_task(task974);

  vector<IndexRange> I1782_index = {active_, closed_, virt_, closed_};
  auto I1782 = make_shared<Tensor>(I1782_index);
  auto tensor975 = vector<shared_ptr<Tensor>>{n, I1782};
  auto task975 = make_shared<Task975>(tensor975, pindex);
  task975->add_dep(task970);
  normq->add_task(task975);

  auto tensor976 = vector<shared_ptr<Tensor>>{I1782, t2, Gamma12_()};
  auto task976 = make_shared<Task976>(tensor976, pindex);
  task975->add_dep(task976);
  task976->add_dep(task970);
  normq->add_task(task976);

  auto tensor977 = vector<shared_ptr<Tensor>>{I1782, t2, Gamma12_()};
  auto task977 = make_shared<Task977>(tensor977, pindex);
  task975->add_dep(task977);
  task977->add_dep(task970);
  normq->add_task(task977);

  vector<IndexRange> I1786_index = {virt_, closed_, active_, active_};
  auto I1786 = make_shared<Tensor>(I1786_index);
  auto tensor978 = vector<shared_ptr<Tensor>>{n, I1786};
  auto task978 = make_shared<Task978>(tensor978, pindex);
  task978->add_dep(task970);
  normq->add_task(task978);

  auto tensor979 = vector<shared_ptr<Tensor>>{I1786, Gamma28_(), t2};
  auto task979 = make_shared<Task979>(tensor979, pindex);
  task978->add_dep(task979);
  task979->add_dep(task970);
  normq->add_task(task979);

  auto tensor980 = vector<shared_ptr<Tensor>>{I1786, Gamma31_(), t2};
  auto task980 = make_shared<Task980>(tensor980, pindex);
  task978->add_dep(task980);
  task980->add_dep(task970);
  normq->add_task(task980);

  vector<IndexRange> I1790_index = {active_, active_, virt_, closed_};
  auto I1790 = make_shared<Tensor>(I1790_index);
  auto tensor981 = vector<shared_ptr<Tensor>>{n, I1790};
  auto task981 = make_shared<Task981>(tensor981, pindex);
  task981->add_dep(task970);
  normq->add_task(task981);

  auto tensor982 = vector<shared_ptr<Tensor>>{I1790, t2, Gamma31_()};
  auto task982 = make_shared<Task982>(tensor982, pindex);
  task981->add_dep(task982);
  task982->add_dep(task970);
  normq->add_task(task982);

  auto tensor983 = vector<shared_ptr<Tensor>>{I1790, t2, Gamma31_()};
  auto task983 = make_shared<Task983>(tensor983, pindex);
  task981->add_dep(task983);
  task983->add_dep(task970);
  normq->add_task(task983);

  vector<IndexRange> I1794_index = {active_, active_, active_, virt_};
  auto I1794 = make_shared<Tensor>(I1794_index);
  auto tensor984 = vector<shared_ptr<Tensor>>{n, I1794};
  auto task984 = make_shared<Task984>(tensor984, pindex);
  task984->add_dep(task970);
  normq->add_task(task984);

  auto tensor985 = vector<shared_ptr<Tensor>>{I1794, t2, Gamma54_()};
  auto task985 = make_shared<Task985>(tensor985, pindex);
  task984->add_dep(task985);
  task985->add_dep(task970);
  normq->add_task(task985);

  shared_ptr<Tensor> I1796;
  if (diagonal) {
    vector<IndexRange> I1796_index = {closed_, virt_, closed_, virt_};
    I1796 = make_shared<Tensor>(I1796_index);
  }
  shared_ptr<Task986> task986;
  if (diagonal) {
    auto tensor986 = vector<shared_ptr<Tensor>>{n, I1796};
    task986 = make_shared<Task986>(tensor986, pindex);
    task986->add_dep(task970);
    normq->add_task(task986);
  }

  shared_ptr<Task987> task987;
  if (diagonal) {
    auto tensor987 = vector<shared_ptr<Tensor>>{I1796, t2};
    task987 = make_shared<Task987>(tensor987, pindex);
    task986->add_dep(task987);
    task987->add_dep(task970);
    normq->add_task(task987);
  }

  vector<IndexRange> I1798_index = {active_, virt_, closed_, virt_};
  auto I1798 = make_shared<Tensor>(I1798_index);
  auto tensor988 = vector<shared_ptr<Tensor>>{n, I1798};
  auto task988 = make_shared<Task988>(tensor988, pindex);
  task988->add_dep(task970);
  normq->add_task(task988);

  auto tensor989 = vector<shared_ptr<Tensor>>{I1798, t2, Gamma34_()};
  auto task989 = make_shared<Task989>(tensor989, pindex);
  task988->add_dep(task989);
  task989->add_dep(task970);
  normq->add_task(task989);

  auto tensor990 = vector<shared_ptr<Tensor>>{I1798, Gamma34_(), t2};
  auto task990 = make_shared<Task990>(tensor990, pindex);
  task988->add_dep(task990);
  task990->add_dep(task970);
  normq->add_task(task990);

  vector<IndexRange> I1802_index = {virt_, virt_, active_, active_};
  auto I1802 = make_shared<Tensor>(I1802_index);
  auto tensor991 = vector<shared_ptr<Tensor>>{n, I1802};
  auto task991 = make_shared<Task991>(tensor991, pindex);
  task991->add_dep(task970);
  normq->add_task(task991);

  auto tensor992 = vector<shared_ptr<Tensor>>{I1802, Gamma55_(), t2};
  auto task992 = make_shared<Task992>(tensor992, pindex);
  task991->add_dep(task992);
  task992->add_dep(task970);
  normq->add_task(task992);

  return normq;
}


#endif
