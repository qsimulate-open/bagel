//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_sourceqq.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_sourceq(const bool reset, const bool diagonal) {

  auto sourceq = make_shared<Queue>();
  auto task229 = make_shared<Task229>(s, reset);
  sourceq->add_task(task229);

  auto I288 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task230 = make_shared<Task230>(s, I288);
  task230->add_dep(task229);
  sourceq->add_task(task230);

  auto task231 = make_shared<Task231>(I288, Gamma92_(), v2_);
  task230->add_dep(task231);
  task231->add_dep(task229);
  sourceq->add_task(task231);

  auto I290 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task232 = make_shared<Task232>(s, I290);
  task232->add_dep(task229);
  sourceq->add_task(task232);

  auto task233 = make_shared<Task233>(I290, Gamma105_(), v2_);
  task232->add_dep(task233);
  task233->add_dep(task229);
  sourceq->add_task(task233);

  auto task234 = make_shared<Task234>(I290, Gamma6_(), v2_);
  task232->add_dep(task234);
  task234->add_dep(task229);
  sourceq->add_task(task234);

  auto task235 = make_shared<Task235>(I290, Gamma7_(), h1_);
  task232->add_dep(task235);
  task235->add_dep(task229);
  sourceq->add_task(task235);

  auto I294 = make_shared<TATensor<double,4>>({closed_, closed_, virt_, active_});
  auto task236 = make_shared<Task236>(s, I294);
  task236->add_dep(task229);
  sourceq->add_task(task236);

  auto I295 = make_shared<TATensor<double,4>>({closed_, active_, closed_, virt_});
  auto task237 = make_shared<Task237>(I294, Gamma16_(), I295);
  task236->add_dep(task237);
  task237->add_dep(task229);
  sourceq->add_task(task237);

  auto task238 = make_shared<Task238>(I295, v2_);
  task237->add_dep(task238);
  task238->add_dep(task229);
  sourceq->add_task(task238);

  auto I298 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task239 = make_shared<Task239>(s, I298);
  task239->add_dep(task229);
  sourceq->add_task(task239);

  auto I299 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task240 = make_shared<Task240>(I298, Gamma35_(), I299);
  task239->add_dep(task240);
  task240->add_dep(task229);
  sourceq->add_task(task240);

  auto task241 = make_shared<Task241>(I299, v2_);
  task240->add_dep(task241);
  task241->add_dep(task229);
  sourceq->add_task(task241);

  auto task242 = make_shared<Task242>(I298, Gamma29_(), v2_);
  task239->add_dep(task242);
  task242->add_dep(task229);
  sourceq->add_task(task242);

  auto task243 = make_shared<Task243>(I298, Gamma32_(), v2_);
  task239->add_dep(task243);
  task243->add_dep(task229);
  sourceq->add_task(task243);

  auto task244 = make_shared<Task244>(I298, Gamma38_(), h1_);
  task239->add_dep(task244);
  task244->add_dep(task229);
  sourceq->add_task(task244);

  auto I306 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task245 = make_shared<Task245>(s, I306);
  task245->add_dep(task229);
  sourceq->add_task(task245);

  auto I307 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task246 = make_shared<Task246>(I306, Gamma35_(), I307);
  task245->add_dep(task246);
  task246->add_dep(task229);
  sourceq->add_task(task246);

  auto task247 = make_shared<Task247>(I307, v2_);
  task246->add_dep(task247);
  task247->add_dep(task229);
  sourceq->add_task(task247);

  auto task248 = make_shared<Task248>(I306, Gamma7_(), v2_);
  task245->add_dep(task248);
  task248->add_dep(task229);
  sourceq->add_task(task248);

  auto task249 = make_shared<Task249>(I306, Gamma38_(), h1_);
  task245->add_dep(task249);
  task249->add_dep(task229);
  sourceq->add_task(task249);

  auto I314 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task250 = make_shared<Task250>(s, I314);
  task250->add_dep(task229);
  sourceq->add_task(task250);

  auto task251 = make_shared<Task251>(I314, Gamma59_(), v2_);
  task250->add_dep(task251);
  task251->add_dep(task229);
  sourceq->add_task(task251);

  auto task252 = make_shared<Task252>(I314, Gamma57_(), v2_);
  task250->add_dep(task252);
  task252->add_dep(task229);
  sourceq->add_task(task252);

  auto task253 = make_shared<Task253>(I314, Gamma60_(), h1_);
  task250->add_dep(task253);
  task253->add_dep(task229);
  sourceq->add_task(task253);

  shared_ptr<TATensor<double,4>> I318;
  if (diagonal) {
    I318 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task254> task254;
  if (diagonal) {
    task254 = make_shared<Task254>(s, I318);
    task254->add_dep(task229);
    sourceq->add_task(task254);
  }

  shared_ptr<Task255> task255;
  if (diagonal) {
    task255 = make_shared<Task255>(I318, v2_);
    task254->add_dep(task255);
    task255->add_dep(task229);
    sourceq->add_task(task255);
  }

  auto I320 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task256 = make_shared<Task256>(s, I320);
  task256->add_dep(task229);
  sourceq->add_task(task256);

  auto I321 = make_shared<TATensor<double,4>>({active_, virt_, closed_, virt_});
  auto task257 = make_shared<Task257>(I320, Gamma38_(), I321);
  task256->add_dep(task257);
  task257->add_dep(task229);
  sourceq->add_task(task257);

  auto task258 = make_shared<Task258>(I321, v2_);
  task257->add_dep(task258);
  task258->add_dep(task229);
  sourceq->add_task(task258);

  auto I324 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task259 = make_shared<Task259>(s, I324);
  task259->add_dep(task229);
  sourceq->add_task(task259);

  auto task260 = make_shared<Task260>(I324, Gamma60_(), v2_);
  task259->add_dep(task260);
  task260->add_dep(task229);
  sourceq->add_task(task260);

  return sourceq;
}


#endif
