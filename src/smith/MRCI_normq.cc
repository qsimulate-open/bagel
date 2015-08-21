//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_normqq.cc
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

shared_ptr<Queue> MRCI::MRCI::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor969 = vector<shared_ptr<Tensor>>{n};
  auto task969 = make_shared<Task969>(tensor969, reset);
  normq->add_task(task969);

  vector<IndexRange> I1778_index = {active_, active_, closed_, closed_};
  auto I1778 = make_shared<Tensor>(I1778_index);
  auto tensor970 = vector<shared_ptr<Tensor>>{n, I1778};
  auto task970 = make_shared<Task970>(tensor970, pindex);
  task970->add_dep(task969);
  normq->add_task(task970);

  auto tensor971 = vector<shared_ptr<Tensor>>{I1778, t2, Gamma0_()};
  auto task971 = make_shared<Task971>(tensor971, pindex);
  task970->add_dep(task971);
  task971->add_dep(task969);
  normq->add_task(task971);

  vector<IndexRange> I1780_index = {closed_, active_, active_, active_};
  auto I1780 = make_shared<Tensor>(I1780_index);
  auto tensor972 = vector<shared_ptr<Tensor>>{n, I1780};
  auto task972 = make_shared<Task972>(tensor972, pindex);
  task972->add_dep(task969);
  normq->add_task(task972);

  auto tensor973 = vector<shared_ptr<Tensor>>{I1780, Gamma4_(), t2};
  auto task973 = make_shared<Task973>(tensor973, pindex);
  task972->add_dep(task973);
  task973->add_dep(task969);
  normq->add_task(task973);

  vector<IndexRange> I1782_index = {active_, closed_, virt_, closed_};
  auto I1782 = make_shared<Tensor>(I1782_index);
  auto tensor974 = vector<shared_ptr<Tensor>>{n, I1782};
  auto task974 = make_shared<Task974>(tensor974, pindex);
  task974->add_dep(task969);
  normq->add_task(task974);

  auto tensor975 = vector<shared_ptr<Tensor>>{I1782, t2, Gamma12_()};
  auto task975 = make_shared<Task975>(tensor975, pindex);
  task974->add_dep(task975);
  task975->add_dep(task969);
  normq->add_task(task975);

  auto tensor976 = vector<shared_ptr<Tensor>>{I1782, Gamma12_(), t2};
  auto task976 = make_shared<Task976>(tensor976, pindex);
  task974->add_dep(task976);
  task976->add_dep(task969);
  normq->add_task(task976);

  vector<IndexRange> I1786_index = {active_, active_, virt_, closed_};
  auto I1786 = make_shared<Tensor>(I1786_index);
  auto tensor977 = vector<shared_ptr<Tensor>>{n, I1786};
  auto task977 = make_shared<Task977>(tensor977, pindex);
  task977->add_dep(task969);
  normq->add_task(task977);

  auto tensor978 = vector<shared_ptr<Tensor>>{I1786, t2, Gamma27_()};
  auto task978 = make_shared<Task978>(tensor978, pindex);
  task977->add_dep(task978);
  task978->add_dep(task969);
  normq->add_task(task978);

  auto tensor979 = vector<shared_ptr<Tensor>>{I1786, t2, Gamma29_()};
  auto task979 = make_shared<Task979>(tensor979, pindex);
  task977->add_dep(task979);
  task979->add_dep(task969);
  normq->add_task(task979);

  vector<IndexRange> I1790_index = {active_, active_, virt_, closed_};
  auto I1790 = make_shared<Tensor>(I1790_index);
  auto tensor980 = vector<shared_ptr<Tensor>>{n, I1790};
  auto task980 = make_shared<Task980>(tensor980, pindex);
  task980->add_dep(task969);
  normq->add_task(task980);

  auto tensor981 = vector<shared_ptr<Tensor>>{I1790, t2, Gamma29_()};
  auto task981 = make_shared<Task981>(tensor981, pindex);
  task980->add_dep(task981);
  task981->add_dep(task969);
  normq->add_task(task981);

  auto tensor982 = vector<shared_ptr<Tensor>>{I1790, t2, Gamma29_()};
  auto task982 = make_shared<Task982>(tensor982, pindex);
  task980->add_dep(task982);
  task982->add_dep(task969);
  normq->add_task(task982);

  vector<IndexRange> I1794_index = {active_, active_, active_, virt_};
  auto I1794 = make_shared<Tensor>(I1794_index);
  auto tensor983 = vector<shared_ptr<Tensor>>{n, I1794};
  auto task983 = make_shared<Task983>(tensor983, pindex);
  task983->add_dep(task969);
  normq->add_task(task983);

  auto tensor984 = vector<shared_ptr<Tensor>>{I1794, t2, Gamma50_()};
  auto task984 = make_shared<Task984>(tensor984, pindex);
  task983->add_dep(task984);
  task984->add_dep(task969);
  normq->add_task(task984);

  shared_ptr<Tensor> I1796;
  if (diagonal) {
    vector<IndexRange> I1796_index = {closed_, virt_, closed_, virt_};
    I1796 = make_shared<Tensor>(I1796_index);
  }
  shared_ptr<Task985> task985;
  if (diagonal) {
    auto tensor985 = vector<shared_ptr<Tensor>>{n, I1796};
    task985 = make_shared<Task985>(tensor985, pindex);
    task985->add_dep(task969);
    normq->add_task(task985);
  }

  shared_ptr<Task986> task986;
  if (diagonal) {
    auto tensor986 = vector<shared_ptr<Tensor>>{I1796, t2};
    task986 = make_shared<Task986>(tensor986, pindex);
    task985->add_dep(task986);
    task986->add_dep(task969);
    normq->add_task(task986);
  }

  vector<IndexRange> I1798_index = {active_, virt_, closed_, virt_};
  auto I1798 = make_shared<Tensor>(I1798_index);
  auto tensor987 = vector<shared_ptr<Tensor>>{n, I1798};
  auto task987 = make_shared<Task987>(tensor987, pindex);
  task987->add_dep(task969);
  normq->add_task(task987);

  auto tensor988 = vector<shared_ptr<Tensor>>{I1798, t2, Gamma32_()};
  auto task988 = make_shared<Task988>(tensor988, pindex);
  task987->add_dep(task988);
  task988->add_dep(task969);
  normq->add_task(task988);

  auto tensor989 = vector<shared_ptr<Tensor>>{I1798, t2, Gamma32_()};
  auto task989 = make_shared<Task989>(tensor989, pindex);
  task987->add_dep(task989);
  task989->add_dep(task969);
  normq->add_task(task989);

  vector<IndexRange> I1802_index = {virt_, virt_, active_, active_};
  auto I1802 = make_shared<Tensor>(I1802_index);
  auto tensor990 = vector<shared_ptr<Tensor>>{n, I1802};
  auto task990 = make_shared<Task990>(tensor990, pindex);
  task990->add_dep(task969);
  normq->add_task(task990);

  auto tensor991 = vector<shared_ptr<Tensor>>{I1802, Gamma51_(), t2};
  auto task991 = make_shared<Task991>(tensor991, pindex);
  task990->add_dep(task991);
  task991->add_dep(task969);
  normq->add_task(task991);

  return normq;
}


#endif
