//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_residualqq.cc
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
#include <src/smith/mrci/MRCI_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_residualq(const bool reset, const bool diagonal) {
  auto out = make_shared<Queue>();
  auto tensor108 = vector<shared_ptr<Tensor>>{r};
  auto task108 = make_shared<Task108>(tensor108, reset);
  out->add_task(task108);

  make_residualq1(out, task108, diagonal);
  make_residualq2(out, task108, diagonal);
  make_residualq3(out, task108, diagonal);
  make_residualq4(out, task108, diagonal);
  make_residualq5(out, task108, diagonal);
  make_residualq6(out, task108, diagonal);
  make_residualq7(out, task108, diagonal);
  make_residualq8(out, task108, diagonal);
  make_residualq9(out, task108, diagonal);
  return out;
}


void MRCI::MRCI::make_residualq1(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I0_index = {closed_, closed_, active_, active_};
  auto I0 = make_shared<Tensor>(I0_index);
  auto tensor109 = vector<shared_ptr<Tensor>>{r, I0};
  auto task109 = make_shared<Task109>(tensor109, pindex);
  task109->add_dep(task108);
  residualq->add_task(task109);

  vector<IndexRange> I1_index = {closed_, closed_, active_, active_};
  auto I1 = make_shared<Tensor>(I1_index);
  auto tensor110 = vector<shared_ptr<Tensor>>{I0, Gamma0_(), I1};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task109->add_dep(task110);
  task110->add_dep(task108);
  residualq->add_task(task110);

  auto tensor111 = vector<shared_ptr<Tensor>>{I1, t2, h1_};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task110->add_dep(task111);
  task111->add_dep(task108);
  residualq->add_task(task111);

  vector<IndexRange> I280_index = {virt_, active_, closed_, closed_};
  auto I280 = make_shared<Tensor>(I280_index);
  auto tensor112 = vector<shared_ptr<Tensor>>{I1, t2, I280};
  auto task112 = make_shared<Task112>(tensor112, pindex);
  task110->add_dep(task112);
  task112->add_dep(task108);
  residualq->add_task(task112);

  auto tensor113 = vector<shared_ptr<Tensor>>{I280, v2_};
  auto task113 = make_shared<Task113>(tensor113, pindex);
  task112->add_dep(task113);
  task113->add_dep(task108);
  residualq->add_task(task113);

  auto tensor114 = vector<shared_ptr<Tensor>>{I1, t2, v2_};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task110->add_dep(task114);
  task114->add_dep(task108);
  residualq->add_task(task114);

  vector<IndexRange> I4_index = {closed_, active_, active_, active_};
  auto I4 = make_shared<Tensor>(I4_index);
  auto tensor115 = vector<shared_ptr<Tensor>>{I0, h1_, I4};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task109->add_dep(task115);
  task115->add_dep(task108);
  residualq->add_task(task115);

  auto tensor116 = vector<shared_ptr<Tensor>>{I4, Gamma1_(), t2};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task115->add_dep(task116);
  task116->add_dep(task108);
  residualq->add_task(task116);

  vector<IndexRange> I7_index = {active_, closed_, closed_, active_};
  auto I7 = make_shared<Tensor>(I7_index);
  auto tensor117 = vector<shared_ptr<Tensor>>{I0, Gamma2_(), I7};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task109->add_dep(task117);
  task117->add_dep(task108);
  residualq->add_task(task117);

  auto tensor118 = vector<shared_ptr<Tensor>>{I7, t2, h1_};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task117->add_dep(task118);
  task118->add_dep(task108);
  residualq->add_task(task118);

  auto tensor119 = vector<shared_ptr<Tensor>>{I7, t2, v2_};
  auto task119 = make_shared<Task119>(tensor119, pindex);
  task117->add_dep(task119);
  task119->add_dep(task108);
  residualq->add_task(task119);

  vector<IndexRange> I249_index = {closed_, active_, active_, closed_, active_, active_};
  auto I249 = make_shared<Tensor>(I249_index);
  auto tensor120 = vector<shared_ptr<Tensor>>{I0, Gamma80_(), I249};
  auto task120 = make_shared<Task120>(tensor120, pindex);
  task109->add_dep(task120);
  task120->add_dep(task108);
  residualq->add_task(task120);

  vector<IndexRange> I250_index = {closed_, closed_, active_, active_};
  auto I250 = make_shared<Tensor>(I250_index);
  auto tensor121 = vector<shared_ptr<Tensor>>{I249, t2, I250};
  auto task121 = make_shared<Task121>(tensor121, pindex);
  task120->add_dep(task121);
  task121->add_dep(task108);
  residualq->add_task(task121);

  auto tensor122 = vector<shared_ptr<Tensor>>{I250, v2_};
  auto task122 = make_shared<Task122>(tensor122, pindex);
  task121->add_dep(task122);
  task122->add_dep(task108);
  residualq->add_task(task122);

  vector<IndexRange> I252_index = {closed_, active_, active_, closed_, active_, active_};
  auto I252 = make_shared<Tensor>(I252_index);
  auto tensor123 = vector<shared_ptr<Tensor>>{I0, Gamma81_(), I252};
  auto task123 = make_shared<Task123>(tensor123, pindex);
  task109->add_dep(task123);
  task123->add_dep(task108);
  residualq->add_task(task123);

  auto tensor124 = vector<shared_ptr<Tensor>>{I252, t2, v2_};
  auto task124 = make_shared<Task124>(tensor124, pindex);
  task123->add_dep(task124);
  task124->add_dep(task108);
  residualq->add_task(task124);

  vector<IndexRange> I255_index = {active_, closed_, active_, closed_, active_, active_};
  auto I255 = make_shared<Tensor>(I255_index);
  auto tensor125 = vector<shared_ptr<Tensor>>{I0, Gamma82_(), I255};
  auto task125 = make_shared<Task125>(tensor125, pindex);
  task109->add_dep(task125);
  task125->add_dep(task108);
  residualq->add_task(task125);

  auto tensor126 = vector<shared_ptr<Tensor>>{I255, t2, v2_};
  auto task126 = make_shared<Task126>(tensor126, pindex);
  task125->add_dep(task126);
  task126->add_dep(task108);
  residualq->add_task(task126);

  vector<IndexRange> I264_index = {closed_, active_, active_, active_, active_, active_};
  auto I264 = make_shared<Tensor>(I264_index);
  auto tensor127 = vector<shared_ptr<Tensor>>{I0, t2, I264};
  auto task127 = make_shared<Task127>(tensor127, pindex);
  task109->add_dep(task127);
  task127->add_dep(task108);
  residualq->add_task(task127);

  auto tensor128 = vector<shared_ptr<Tensor>>{I264, Gamma85_(), v2_};
  auto task128 = make_shared<Task128>(tensor128, pindex);
  task127->add_dep(task128);
  task128->add_dep(task108);
  residualq->add_task(task128);

  auto tensor129 = vector<shared_ptr<Tensor>>{I264, Gamma86_(), v2_};
  auto task129 = make_shared<Task129>(tensor129, pindex);
  task127->add_dep(task129);
  task129->add_dep(task108);
  residualq->add_task(task129);

  vector<IndexRange> I270_index = {closed_, active_, active_, active_};
  auto I270 = make_shared<Tensor>(I270_index);
  auto tensor130 = vector<shared_ptr<Tensor>>{I0, v2_, I270};
  auto task130 = make_shared<Task130>(tensor130, pindex);
  task109->add_dep(task130);
  task130->add_dep(task108);
  residualq->add_task(task130);

  auto tensor131 = vector<shared_ptr<Tensor>>{I270, Gamma87_(), t2};
  auto task131 = make_shared<Task131>(tensor131, pindex);
  task130->add_dep(task131);
  task131->add_dep(task108);
  residualq->add_task(task131);

  vector<IndexRange> I273_index = {virt_, active_, active_, active_};
  auto I273 = make_shared<Tensor>(I273_index);
  auto tensor132 = vector<shared_ptr<Tensor>>{I0, t2, I273};
  auto task132 = make_shared<Task132>(tensor132, pindex);
  task109->add_dep(task132);
  task132->add_dep(task108);
  residualq->add_task(task132);

  auto tensor133 = vector<shared_ptr<Tensor>>{I273, Gamma88_(), v2_};
  auto task133 = make_shared<Task133>(tensor133, pindex);
  task132->add_dep(task133);
  task133->add_dep(task108);
  residualq->add_task(task133);

  auto tensor134 = vector<shared_ptr<Tensor>>{I273, Gamma89_(), v2_};
  auto task134 = make_shared<Task134>(tensor134, pindex);
  task132->add_dep(task134);
  task134->add_dep(task108);
  residualq->add_task(task134);

  vector<IndexRange> I291_index = {closed_, active_, active_, active_, closed_, active_};
  auto I291 = make_shared<Tensor>(I291_index);
  auto tensor135 = vector<shared_ptr<Tensor>>{I0, Gamma94_(), I291};
  auto task135 = make_shared<Task135>(tensor135, pindex);
  task109->add_dep(task135);
  task135->add_dep(task108);
  residualq->add_task(task135);

  auto tensor136 = vector<shared_ptr<Tensor>>{I291, t2, v2_};
  auto task136 = make_shared<Task136>(tensor136, pindex);
  task135->add_dep(task136);
  task136->add_dep(task108);
  residualq->add_task(task136);

  vector<IndexRange> I294_index = {closed_, active_, active_, closed_, active_, active_};
  auto I294 = make_shared<Tensor>(I294_index);
  auto tensor137 = vector<shared_ptr<Tensor>>{I0, Gamma87_(), I294};
  auto task137 = make_shared<Task137>(tensor137, pindex);
  task109->add_dep(task137);
  task137->add_dep(task108);
  residualq->add_task(task137);

  auto tensor138 = vector<shared_ptr<Tensor>>{I294, t2, v2_};
  auto task138 = make_shared<Task138>(tensor138, pindex);
  task137->add_dep(task138);
  task138->add_dep(task108);
  residualq->add_task(task138);
}

#endif
