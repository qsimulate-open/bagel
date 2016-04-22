//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci2qq.cc
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

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deci2q = make_shared<Queue>();
  auto tensor502 = vector<shared_ptr<Tensor>>{deci};
  auto task502 = make_shared<Task502>(tensor502, reset);
  deci2q->add_task(task502);

  vector<IndexRange> I722_index = {ci_};
  auto I722 = make_shared<Tensor>(I722_index);
  auto tensor503 = vector<shared_ptr<Tensor>>{deci, I722};
  auto task503 = make_shared<Task503>(tensor503, cindex);
  task503->add_dep(task502);
  deci2q->add_task(task503);

  vector<IndexRange> I723_index = {active_, active_, active_, active_};
  auto I723 = make_shared<Tensor>(I723_index);
  auto tensor504 = vector<shared_ptr<Tensor>>{I722, Gamma111_(), I723};
  auto task504 = make_shared<Task504>(tensor504, cindex);
  task503->add_dep(task504);
  task504->add_dep(task502);
  deci2q->add_task(task504);

  auto tensor505 = vector<shared_ptr<Tensor>>{I723, v2_, l2};
  auto task505 = make_shared<Task505>(tensor505, cindex);
  task504->add_dep(task505);
  task505->add_dep(task502);
  deci2q->add_task(task505);

  vector<IndexRange> I726_index = {active_, active_, active_, active_, active_, active_};
  auto I726 = make_shared<Tensor>(I726_index);
  auto tensor506 = vector<shared_ptr<Tensor>>{I722, Gamma217_(), I726};
  auto task506 = make_shared<Task506>(tensor506, cindex);
  task503->add_dep(task506);
  task506->add_dep(task502);
  deci2q->add_task(task506);

  auto tensor507 = vector<shared_ptr<Tensor>>{I726, v2_, l2};
  auto task507 = make_shared<Task507>(tensor507, cindex);
  task506->add_dep(task507);
  task507->add_dep(task502);
  deci2q->add_task(task507);

  vector<IndexRange> I729_index = {active_, active_, active_, active_, active_, active_};
  auto I729 = make_shared<Tensor>(I729_index);
  auto tensor508 = vector<shared_ptr<Tensor>>{I722, Gamma116_(), I729};
  auto task508 = make_shared<Task508>(tensor508, cindex);
  task503->add_dep(task508);
  task508->add_dep(task502);
  deci2q->add_task(task508);

  auto tensor509 = vector<shared_ptr<Tensor>>{I729, v2_, l2};
  auto task509 = make_shared<Task509>(tensor509, cindex);
  task508->add_dep(task509);
  task509->add_dep(task502);
  deci2q->add_task(task509);

  vector<IndexRange> I732_index = {active_, active_};
  auto I732 = make_shared<Tensor>(I732_index);
  auto tensor510 = vector<shared_ptr<Tensor>>{I722, Gamma126_(), I732};
  auto task510 = make_shared<Task510>(tensor510, cindex);
  task503->add_dep(task510);
  task510->add_dep(task502);
  deci2q->add_task(task510);

  auto tensor511 = vector<shared_ptr<Tensor>>{I732, v2_, l2};
  auto task511 = make_shared<Task511>(tensor511, cindex);
  task510->add_dep(task511);
  task511->add_dep(task502);
  deci2q->add_task(task511);

  auto tensor512 = vector<shared_ptr<Tensor>>{I732, v2_, l2};
  auto task512 = make_shared<Task512>(tensor512, cindex);
  task510->add_dep(task512);
  task512->add_dep(task502);
  deci2q->add_task(task512);

  vector<IndexRange> I738_index = {active_, active_, active_, active_};
  auto I738 = make_shared<Tensor>(I738_index);
  auto tensor513 = vector<shared_ptr<Tensor>>{I722, Gamma145_(), I738};
  auto task513 = make_shared<Task513>(tensor513, cindex);
  task503->add_dep(task513);
  task513->add_dep(task502);
  deci2q->add_task(task513);

  auto tensor514 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task514 = make_shared<Task514>(tensor514, cindex);
  task513->add_dep(task514);
  task514->add_dep(task502);
  deci2q->add_task(task514);

  auto tensor515 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task515 = make_shared<Task515>(tensor515, cindex);
  task513->add_dep(task515);
  task515->add_dep(task502);
  deci2q->add_task(task515);

  auto tensor516 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task516 = make_shared<Task516>(tensor516, cindex);
  task513->add_dep(task516);
  task516->add_dep(task502);
  deci2q->add_task(task516);

  auto tensor517 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task517 = make_shared<Task517>(tensor517, cindex);
  task513->add_dep(task517);
  task517->add_dep(task502);
  deci2q->add_task(task517);

  auto tensor518 = vector<shared_ptr<Tensor>>{I738, v2_, l2};
  auto task518 = make_shared<Task518>(tensor518, cindex);
  task513->add_dep(task518);
  task518->add_dep(task502);
  deci2q->add_task(task518);

  vector<IndexRange> I741_index = {active_, active_, active_, active_};
  auto I741 = make_shared<Tensor>(I741_index);
  auto tensor519 = vector<shared_ptr<Tensor>>{I722, Gamma139_(), I741};
  auto task519 = make_shared<Task519>(tensor519, cindex);
  task503->add_dep(task519);
  task519->add_dep(task502);
  deci2q->add_task(task519);

  auto tensor520 = vector<shared_ptr<Tensor>>{I741, v2_, l2};
  auto task520 = make_shared<Task520>(tensor520, cindex);
  task519->add_dep(task520);
  task520->add_dep(task502);
  deci2q->add_task(task520);

  vector<IndexRange> I744_index = {active_, active_, active_, active_};
  auto I744 = make_shared<Tensor>(I744_index);
  auto tensor521 = vector<shared_ptr<Tensor>>{I722, Gamma142_(), I744};
  auto task521 = make_shared<Task521>(tensor521, cindex);
  task503->add_dep(task521);
  task521->add_dep(task502);
  deci2q->add_task(task521);

  auto tensor522 = vector<shared_ptr<Tensor>>{I744, v2_, l2};
  auto task522 = make_shared<Task522>(tensor522, cindex);
  task521->add_dep(task522);
  task522->add_dep(task502);
  deci2q->add_task(task522);

  vector<IndexRange> I753_index = {active_, active_, active_, active_};
  auto I753 = make_shared<Tensor>(I753_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{I722, Gamma117_(), I753};
  auto task523 = make_shared<Task523>(tensor523, cindex);
  task503->add_dep(task523);
  task523->add_dep(task502);
  deci2q->add_task(task523);

  auto tensor524 = vector<shared_ptr<Tensor>>{I753, v2_, l2};
  auto task524 = make_shared<Task524>(tensor524, cindex);
  task523->add_dep(task524);
  task524->add_dep(task502);
  deci2q->add_task(task524);

  auto tensor525 = vector<shared_ptr<Tensor>>{I753, h1_, l2};
  auto task525 = make_shared<Task525>(tensor525, cindex);
  task523->add_dep(task525);
  task525->add_dep(task502);
  deci2q->add_task(task525);

  vector<IndexRange> I762_index = {active_, active_, active_, active_, active_, active_};
  auto I762 = make_shared<Tensor>(I762_index);
  auto tensor526 = vector<shared_ptr<Tensor>>{I722, Gamma169_(), I762};
  auto task526 = make_shared<Task526>(tensor526, cindex);
  task503->add_dep(task526);
  task526->add_dep(task502);
  deci2q->add_task(task526);

  auto tensor527 = vector<shared_ptr<Tensor>>{I762, v2_, l2};
  auto task527 = make_shared<Task527>(tensor527, cindex);
  task526->add_dep(task527);
  task527->add_dep(task502);
  deci2q->add_task(task527);

  vector<IndexRange> I765_index = {active_, active_, active_, active_, active_, active_};
  auto I765 = make_shared<Tensor>(I765_index);
  auto tensor528 = vector<shared_ptr<Tensor>>{I722, Gamma167_(), I765};
  auto task528 = make_shared<Task528>(tensor528, cindex);
  task503->add_dep(task528);
  task528->add_dep(task502);
  deci2q->add_task(task528);

  auto tensor529 = vector<shared_ptr<Tensor>>{I765, v2_, l2};
  auto task529 = make_shared<Task529>(tensor529, cindex);
  task528->add_dep(task529);
  task529->add_dep(task502);
  deci2q->add_task(task529);

  vector<IndexRange> I768_index = {active_, active_};
  auto I768 = make_shared<Tensor>(I768_index);
  auto tensor530 = vector<shared_ptr<Tensor>>{I722, Gamma148_(), I768};
  auto task530 = make_shared<Task530>(tensor530, cindex);
  task503->add_dep(task530);
  task530->add_dep(task502);
  deci2q->add_task(task530);

  auto tensor531 = vector<shared_ptr<Tensor>>{I768, v2_, l2};
  auto task531 = make_shared<Task531>(tensor531, cindex);
  task530->add_dep(task531);
  task531->add_dep(task502);
  deci2q->add_task(task531);

  auto tensor532 = vector<shared_ptr<Tensor>>{I768, v2_, l2};
  auto task532 = make_shared<Task532>(tensor532, cindex);
  task530->add_dep(task532);
  task532->add_dep(task502);
  deci2q->add_task(task532);

  auto tensor533 = vector<shared_ptr<Tensor>>{I768, h1_, l2};
  auto task533 = make_shared<Task533>(tensor533, cindex);
  task530->add_dep(task533);
  task533->add_dep(task502);
  deci2q->add_task(task533);

  auto tensor534 = vector<shared_ptr<Tensor>>{I768, h1_, l2};
  auto task534 = make_shared<Task534>(tensor534, cindex);
  task530->add_dep(task534);
  task534->add_dep(task502);
  deci2q->add_task(task534);

  vector<IndexRange> I774_index = {active_, active_, active_, active_};
  auto I774 = make_shared<Tensor>(I774_index);
  auto tensor535 = vector<shared_ptr<Tensor>>{I722, Gamma170_(), I774};
  auto task535 = make_shared<Task535>(tensor535, cindex);
  task503->add_dep(task535);
  task535->add_dep(task502);
  deci2q->add_task(task535);

  auto tensor536 = vector<shared_ptr<Tensor>>{I774, v2_, l2};
  auto task536 = make_shared<Task536>(tensor536, cindex);
  task535->add_dep(task536);
  task536->add_dep(task502);
  deci2q->add_task(task536);

  auto tensor537 = vector<shared_ptr<Tensor>>{I774, h1_, l2};
  auto task537 = make_shared<Task537>(tensor537, cindex);
  task535->add_dep(task537);
  task537->add_dep(task502);
  deci2q->add_task(task537);

  return deci2q;
}


#endif
