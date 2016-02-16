//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_normqq.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks20.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor967 = vector<shared_ptr<Tensor>>{n};
  auto task967 = make_shared<Task967>(tensor967, reset);
  normq->add_task(task967);

  vector<IndexRange> I1774_index = {active_, active_, closed_, closed_};
  auto I1774 = make_shared<Tensor>(I1774_index);
  auto tensor968 = vector<shared_ptr<Tensor>>{n, I1774};
  auto task968 = make_shared<Task968>(tensor968, pindex);
  task968->add_dep(task967);
  normq->add_task(task968);

  auto tensor969 = vector<shared_ptr<Tensor>>{I1774, t2, Gamma0_()};
  auto task969 = make_shared<Task969>(tensor969, pindex);
  task968->add_dep(task969);
  task969->add_dep(task967);
  normq->add_task(task969);

  vector<IndexRange> I1776_index = {active_, active_, active_, closed_};
  auto I1776 = make_shared<Tensor>(I1776_index);
  auto tensor970 = vector<shared_ptr<Tensor>>{n, I1776};
  auto task970 = make_shared<Task970>(tensor970, pindex);
  task970->add_dep(task967);
  normq->add_task(task970);

  auto tensor971 = vector<shared_ptr<Tensor>>{I1776, t2, Gamma4_()};
  auto task971 = make_shared<Task971>(tensor971, pindex);
  task970->add_dep(task971);
  task971->add_dep(task967);
  normq->add_task(task971);

  vector<IndexRange> I1778_index = {active_, closed_, virt_, closed_};
  auto I1778 = make_shared<Tensor>(I1778_index);
  auto tensor972 = vector<shared_ptr<Tensor>>{n, I1778};
  auto task972 = make_shared<Task972>(tensor972, pindex);
  task972->add_dep(task967);
  normq->add_task(task972);

  auto tensor973 = vector<shared_ptr<Tensor>>{I1778, t2, Gamma12_()};
  auto task973 = make_shared<Task973>(tensor973, pindex);
  task972->add_dep(task973);
  task973->add_dep(task967);
  normq->add_task(task973);

  auto tensor974 = vector<shared_ptr<Tensor>>{I1778, t2, Gamma12_()};
  auto task974 = make_shared<Task974>(tensor974, pindex);
  task972->add_dep(task974);
  task974->add_dep(task967);
  normq->add_task(task974);

  vector<IndexRange> I1782_index = {active_, active_, virt_, closed_};
  auto I1782 = make_shared<Tensor>(I1782_index);
  auto tensor975 = vector<shared_ptr<Tensor>>{n, I1782};
  auto task975 = make_shared<Task975>(tensor975, pindex);
  task975->add_dep(task967);
  normq->add_task(task975);

  auto tensor976 = vector<shared_ptr<Tensor>>{I1782, t2, Gamma27_()};
  auto task976 = make_shared<Task976>(tensor976, pindex);
  task975->add_dep(task976);
  task976->add_dep(task967);
  normq->add_task(task976);

  auto tensor977 = vector<shared_ptr<Tensor>>{I1782, t2, Gamma29_()};
  auto task977 = make_shared<Task977>(tensor977, pindex);
  task975->add_dep(task977);
  task977->add_dep(task967);
  normq->add_task(task977);

  vector<IndexRange> I1786_index = {active_, active_, virt_, closed_};
  auto I1786 = make_shared<Tensor>(I1786_index);
  auto tensor978 = vector<shared_ptr<Tensor>>{n, I1786};
  auto task978 = make_shared<Task978>(tensor978, pindex);
  task978->add_dep(task967);
  normq->add_task(task978);

  auto tensor979 = vector<shared_ptr<Tensor>>{I1786, t2, Gamma29_()};
  auto task979 = make_shared<Task979>(tensor979, pindex);
  task978->add_dep(task979);
  task979->add_dep(task967);
  normq->add_task(task979);

  auto tensor980 = vector<shared_ptr<Tensor>>{I1786, t2, Gamma29_()};
  auto task980 = make_shared<Task980>(tensor980, pindex);
  task978->add_dep(task980);
  task980->add_dep(task967);
  normq->add_task(task980);

  vector<IndexRange> I1790_index = {active_, active_, active_, virt_};
  auto I1790 = make_shared<Tensor>(I1790_index);
  auto tensor981 = vector<shared_ptr<Tensor>>{n, I1790};
  auto task981 = make_shared<Task981>(tensor981, pindex);
  task981->add_dep(task967);
  normq->add_task(task981);

  auto tensor982 = vector<shared_ptr<Tensor>>{I1790, t2, Gamma50_()};
  auto task982 = make_shared<Task982>(tensor982, pindex);
  task981->add_dep(task982);
  task982->add_dep(task967);
  normq->add_task(task982);

  shared_ptr<Tensor> I1792;
  if (diagonal) {
    vector<IndexRange> I1792_index = {closed_, virt_, closed_, virt_};
    I1792 = make_shared<Tensor>(I1792_index);
  }
  shared_ptr<Task983> task983;
  if (diagonal) {
    auto tensor983 = vector<shared_ptr<Tensor>>{n, I1792};
    task983 = make_shared<Task983>(tensor983, pindex);
    task983->add_dep(task967);
    normq->add_task(task983);
  }

  shared_ptr<Task984> task984;
  if (diagonal) {
    auto tensor984 = vector<shared_ptr<Tensor>>{I1792, t2};
    task984 = make_shared<Task984>(tensor984, pindex);
    task983->add_dep(task984);
    task984->add_dep(task967);
    normq->add_task(task984);
  }

  vector<IndexRange> I1794_index = {active_, virt_, closed_, virt_};
  auto I1794 = make_shared<Tensor>(I1794_index);
  auto tensor985 = vector<shared_ptr<Tensor>>{n, I1794};
  auto task985 = make_shared<Task985>(tensor985, pindex);
  task985->add_dep(task967);
  normq->add_task(task985);

  auto tensor986 = vector<shared_ptr<Tensor>>{I1794, t2, Gamma32_()};
  auto task986 = make_shared<Task986>(tensor986, pindex);
  task985->add_dep(task986);
  task986->add_dep(task967);
  normq->add_task(task986);

  auto tensor987 = vector<shared_ptr<Tensor>>{I1794, Gamma32_(), t2};
  auto task987 = make_shared<Task987>(tensor987, pindex);
  task985->add_dep(task987);
  task987->add_dep(task967);
  normq->add_task(task987);

  vector<IndexRange> I1798_index = {active_, active_, virt_, virt_};
  auto I1798 = make_shared<Tensor>(I1798_index);
  auto tensor988 = vector<shared_ptr<Tensor>>{n, I1798};
  auto task988 = make_shared<Task988>(tensor988, pindex);
  task988->add_dep(task967);
  normq->add_task(task988);

  auto tensor989 = vector<shared_ptr<Tensor>>{I1798, t2, Gamma51_()};
  auto task989 = make_shared<Task989>(tensor989, pindex);
  task988->add_dep(task989);
  task989->add_dep(task967);
  normq->add_task(task989);

  return normq;
}


#endif
