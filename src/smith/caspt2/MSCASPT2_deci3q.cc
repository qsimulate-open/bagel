//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci3q.cc
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci3q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto deci3q = make_shared<Queue>();
  auto tensor480 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task480 = make_shared<Task480>(tensor480, reset);
  deci3q->add_task(task480);

  vector<IndexRange> I723_index = {active_, active_, active_, active_};
  auto I723 = make_shared<Tensor>(I723_index);
  auto tensor482 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I723};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task482->add_dep(task480);
  deci3q->add_task(task482);

  auto tensor483 = vector<shared_ptr<Tensor>>{I723, v2_, l2};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task482->add_dep(task483);
  task483->add_dep(task480);
  deci3q->add_task(task483);

  vector<IndexRange> I726_index = {active_, active_, active_, active_, active_, active_};
  auto I726 = make_shared<Tensor>(I726_index);
  auto tensor484 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I726};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task484->add_dep(task480);
  deci3q->add_task(task484);

  auto tensor485 = vector<shared_ptr<Tensor>>{I726, v2_, l2};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task484->add_dep(task485);
  task485->add_dep(task480);
  deci3q->add_task(task485);

  vector<IndexRange> I729_index = {active_, active_, active_, active_, active_, active_};
  auto I729 = make_shared<Tensor>(I729_index);
  auto tensor486 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I729};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task486->add_dep(task480);
  deci3q->add_task(task486);

  auto tensor487 = vector<shared_ptr<Tensor>>{I729, v2_, l2};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task486->add_dep(task487);
  task487->add_dep(task480);
  deci3q->add_task(task487);

  vector<IndexRange> I732_index = {active_, active_};
  auto I732 = make_shared<Tensor>(I732_index);
  auto tensor488 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I732};
  auto task488 = make_shared<Task488>(tensor488, pindex);
  task488->add_dep(task480);
  deci3q->add_task(task488);

  auto tensor489 = vector<shared_ptr<Tensor>>{I732, v2_, l2};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task488->add_dep(task489);
  task489->add_dep(task480);
  deci3q->add_task(task489);

  auto tensor490 = vector<shared_ptr<Tensor>>{I732, v2_, l2};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task488->add_dep(task490);
  task490->add_dep(task480);
  deci3q->add_task(task490);

  vector<IndexRange> I738_index = {active_, active_, active_, active_};
  auto I738 = make_shared<Tensor>(I738_index);
  auto tensor491 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I738};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task491->add_dep(task480);
  deci3q->add_task(task491);

  auto tensor492 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task491->add_dep(task492);
  task492->add_dep(task480);
  deci3q->add_task(task492);

  auto tensor493 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task491->add_dep(task493);
  task493->add_dep(task480);
  deci3q->add_task(task493);

  auto tensor494 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task491->add_dep(task494);
  task494->add_dep(task480);
  deci3q->add_task(task494);

  auto tensor495 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task491->add_dep(task495);
  task495->add_dep(task480);
  deci3q->add_task(task495);

  auto tensor496 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task496 = make_shared<Task496>(tensor496, pindex);
  task491->add_dep(task496);
  task496->add_dep(task480);
  deci3q->add_task(task496);

  vector<IndexRange> I741_index = {active_, active_, active_, active_};
  auto I741 = make_shared<Tensor>(I741_index);
  auto tensor497 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I741};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task497->add_dep(task480);
  deci3q->add_task(task497);

  auto tensor498 = vector<shared_ptr<Tensor>>{I741, v2_, l2};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task497->add_dep(task498);
  task498->add_dep(task480);
  deci3q->add_task(task498);

  vector<IndexRange> I744_index = {active_, active_, active_, active_};
  auto I744 = make_shared<Tensor>(I744_index);
  auto tensor499 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I744};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task499->add_dep(task480);
  deci3q->add_task(task499);

  auto tensor500 = vector<shared_ptr<Tensor>>{I744, v2_, l2};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task499->add_dep(task500);
  task500->add_dep(task480);
  deci3q->add_task(task500);

  vector<IndexRange> I753_index = {active_, active_, active_, active_};
  auto I753 = make_shared<Tensor>(I753_index);
  auto tensor501 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I753};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task501->add_dep(task480);
  deci3q->add_task(task501);

  auto tensor502 = vector<shared_ptr<Tensor>>{I753, v2_, l2};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task501->add_dep(task502);
  task502->add_dep(task480);
  deci3q->add_task(task502);

  auto tensor503 = vector<shared_ptr<Tensor>>{I753, h1_, l2};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task501->add_dep(task503);
  task503->add_dep(task480);
  deci3q->add_task(task503);

  vector<IndexRange> I762_index = {active_, active_, active_, active_, active_, active_};
  auto I762 = make_shared<Tensor>(I762_index);
  auto tensor504 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I762};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task504->add_dep(task480);
  deci3q->add_task(task504);

  auto tensor505 = vector<shared_ptr<Tensor>>{I762, v2_, l2};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task504->add_dep(task505);
  task505->add_dep(task480);
  deci3q->add_task(task505);

  vector<IndexRange> I765_index = {active_, active_, active_, active_, active_, active_};
  auto I765 = make_shared<Tensor>(I765_index);
  auto tensor506 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I765};
  auto task506 = make_shared<Task506>(tensor506, pindex);
  task506->add_dep(task480);
  deci3q->add_task(task506);

  auto tensor507 = vector<shared_ptr<Tensor>>{I765, v2_, l2};
  auto task507 = make_shared<Task507>(tensor507, pindex);
  task506->add_dep(task507);
  task507->add_dep(task480);
  deci3q->add_task(task507);

  vector<IndexRange> I768_index = {active_, active_};
  auto I768 = make_shared<Tensor>(I768_index);
  auto tensor508 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I768};
  auto task508 = make_shared<Task508>(tensor508, pindex);
  task508->add_dep(task480);
  deci3q->add_task(task508);

  auto tensor509 = vector<shared_ptr<Tensor>>{I768, v2_, l2};
  auto task509 = make_shared<Task509>(tensor509, pindex);
  task508->add_dep(task509);
  task509->add_dep(task480);
  deci3q->add_task(task509);

  auto tensor510 = vector<shared_ptr<Tensor>>{I768, v2_, l2};
  auto task510 = make_shared<Task510>(tensor510, pindex);
  task508->add_dep(task510);
  task510->add_dep(task480);
  deci3q->add_task(task510);

  auto tensor511 = vector<shared_ptr<Tensor>>{I768, h1_, l2};
  auto task511 = make_shared<Task511>(tensor511, pindex);
  task508->add_dep(task511);
  task511->add_dep(task480);
  deci3q->add_task(task511);

  auto tensor512 = vector<shared_ptr<Tensor>>{I768, h1_, l2};
  auto task512 = make_shared<Task512>(tensor512, pindex);
  task508->add_dep(task512);
  task512->add_dep(task480);
  deci3q->add_task(task512);

  vector<IndexRange> I774_index = {active_, active_, active_, active_};
  auto I774 = make_shared<Tensor>(I774_index);
  auto tensor513 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I774};
  auto task513 = make_shared<Task513>(tensor513, pindex);
  task513->add_dep(task480);
  deci3q->add_task(task513);

  auto tensor514 = vector<shared_ptr<Tensor>>{I774, v2_, l2};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task513->add_dep(task514);
  task514->add_dep(task480);
  deci3q->add_task(task514);

  auto tensor515 = vector<shared_ptr<Tensor>>{I774, h1_, l2};
  auto task515 = make_shared<Task515>(tensor515, pindex);
  task513->add_dep(task515);
  task515->add_dep(task480);
  deci3q->add_task(task515);

  return deci3q;
}


#endif
