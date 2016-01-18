//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_normqq.cc
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


#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relmrci/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor758 = vector<shared_ptr<Tensor>>{n};
  auto task758 = make_shared<Task758>(tensor758, reset);
  normq->add_task(task758);

  vector<IndexRange> I1344_index = {active_, active_, closed_, closed_};
  auto I1344 = make_shared<Tensor>(I1344_index);
  auto tensor759 = vector<shared_ptr<Tensor>>{n, I1344};
  auto task759 = make_shared<Task759>(tensor759, pindex);
  task759->add_dep(task758);
  normq->add_task(task759);

  auto tensor760 = vector<shared_ptr<Tensor>>{I1344, t2, Gamma0_()};
  auto task760 = make_shared<Task760>(tensor760, pindex);
  task759->add_dep(task760);
  task760->add_dep(task758);
  normq->add_task(task760);

  vector<IndexRange> I1346_index = {active_, active_, active_, closed_};
  auto I1346 = make_shared<Tensor>(I1346_index);
  auto tensor761 = vector<shared_ptr<Tensor>>{n, I1346};
  auto task761 = make_shared<Task761>(tensor761, pindex);
  task761->add_dep(task758);
  normq->add_task(task761);

  auto tensor762 = vector<shared_ptr<Tensor>>{I1346, t2, Gamma4_()};
  auto task762 = make_shared<Task762>(tensor762, pindex);
  task761->add_dep(task762);
  task762->add_dep(task758);
  normq->add_task(task762);

  vector<IndexRange> I1348_index = {closed_, virt_, closed_, active_};
  auto I1348 = make_shared<Tensor>(I1348_index);
  auto tensor763 = vector<shared_ptr<Tensor>>{n, I1348};
  auto task763 = make_shared<Task763>(tensor763, pindex);
  task763->add_dep(task758);
  normq->add_task(task763);

  vector<IndexRange> I1349_index = {closed_, virt_, closed_, active_};
  auto I1349 = make_shared<Tensor>(I1349_index);
  auto tensor764 = vector<shared_ptr<Tensor>>{I1348, Gamma11_(), I1349};
  auto task764 = make_shared<Task764>(tensor764, pindex);
  task763->add_dep(task764);
  task764->add_dep(task758);
  normq->add_task(task764);

  auto tensor765 = vector<shared_ptr<Tensor>>{I1349, t2};
  auto task765 = make_shared<Task765>(tensor765, pindex);
  task764->add_dep(task765);
  task765->add_dep(task758);
  normq->add_task(task765);

  vector<IndexRange> I1352_index = {active_, active_, closed_, virt_};
  auto I1352 = make_shared<Tensor>(I1352_index);
  auto tensor766 = vector<shared_ptr<Tensor>>{n, I1352};
  auto task766 = make_shared<Task766>(tensor766, pindex);
  task766->add_dep(task758);
  normq->add_task(task766);

  auto tensor767 = vector<shared_ptr<Tensor>>{I1352, t2, Gamma24_()};
  auto task767 = make_shared<Task767>(tensor767, pindex);
  task766->add_dep(task767);
  task767->add_dep(task758);
  normq->add_task(task767);

  vector<IndexRange> I1354_index = {active_, active_, active_, virt_};
  auto I1354 = make_shared<Tensor>(I1354_index);
  auto tensor768 = vector<shared_ptr<Tensor>>{n, I1354};
  auto task768 = make_shared<Task768>(tensor768, pindex);
  task768->add_dep(task758);
  normq->add_task(task768);

  auto tensor769 = vector<shared_ptr<Tensor>>{I1354, t2, Gamma32_()};
  auto task769 = make_shared<Task769>(tensor769, pindex);
  task768->add_dep(task769);
  task769->add_dep(task758);
  normq->add_task(task769);

  shared_ptr<Tensor> I1356;
  if (diagonal) {
    vector<IndexRange> I1356_index = {closed_, virt_, closed_, virt_};
    I1356 = make_shared<Tensor>(I1356_index);
  }
  shared_ptr<Task770> task770;
  if (diagonal) {
    auto tensor770 = vector<shared_ptr<Tensor>>{n, I1356};
    task770 = make_shared<Task770>(tensor770, pindex);
    task770->add_dep(task758);
    normq->add_task(task770);
  }

  shared_ptr<Task771> task771;
  if (diagonal) {
    auto tensor771 = vector<shared_ptr<Tensor>>{I1356, t2};
    task771 = make_shared<Task771>(tensor771, pindex);
    task770->add_dep(task771);
    task771->add_dep(task758);
    normq->add_task(task771);
  }

  vector<IndexRange> I1358_index = {active_, virt_, closed_, virt_};
  auto I1358 = make_shared<Tensor>(I1358_index);
  auto tensor772 = vector<shared_ptr<Tensor>>{n, I1358};
  auto task772 = make_shared<Task772>(tensor772, pindex);
  task772->add_dep(task758);
  normq->add_task(task772);

  auto tensor773 = vector<shared_ptr<Tensor>>{I1358, t2, Gamma27_()};
  auto task773 = make_shared<Task773>(tensor773, pindex);
  task772->add_dep(task773);
  task773->add_dep(task758);
  normq->add_task(task773);

  auto tensor774 = vector<shared_ptr<Tensor>>{I1358, t2, Gamma27_()};
  auto task774 = make_shared<Task774>(tensor774, pindex);
  task772->add_dep(task774);
  task774->add_dep(task758);
  normq->add_task(task774);

  vector<IndexRange> I1362_index = {active_, active_, virt_, virt_};
  auto I1362 = make_shared<Tensor>(I1362_index);
  auto tensor775 = vector<shared_ptr<Tensor>>{n, I1362};
  auto task775 = make_shared<Task775>(tensor775, pindex);
  task775->add_dep(task758);
  normq->add_task(task775);

  auto tensor776 = vector<shared_ptr<Tensor>>{I1362, t2, Gamma33_()};
  auto task776 = make_shared<Task776>(tensor776, pindex);
  task775->add_dep(task776);
  task776->add_dep(task758);
  normq->add_task(task776);

  return normq;
}


#endif
