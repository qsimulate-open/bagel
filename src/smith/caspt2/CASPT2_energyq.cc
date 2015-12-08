//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_energyqq.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_energyq(const bool reset, const bool diagonal) {

  auto energyq = make_shared<Queue>();
  auto I335 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task241 = make_shared<Task241>(Gamma92_(), I335);
  energyq->add_task(task241);

  auto task242 = make_shared<Task242>(I335, v2_, t2);
  task241->add_dep(task242);
  energyq->add_task(task242);

  auto I338 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task243 = make_shared<Task243>(t2, I338);
  task241->add_dep(task243);
  energyq->add_task(task243);

  auto task244 = make_shared<Task244>(I338, Gamma105_(), v2_);
  task243->add_dep(task244);
  energyq->add_task(task244);

  auto task245 = make_shared<Task245>(I338, Gamma6_(), v2_);
  task243->add_dep(task245);
  energyq->add_task(task245);

  auto I344 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task246 = make_shared<Task246>(Gamma16_(), I344);
  task241->add_dep(task246);
  energyq->add_task(task246);

  auto task247 = make_shared<Task247>(I344, v2_, t2);
  task246->add_dep(task247);
  energyq->add_task(task247);

  auto task248 = make_shared<Task248>(I344, v2_, t2);
  task246->add_dep(task248);
  energyq->add_task(task248);

  auto I350 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task249 = make_shared<Task249>(Gamma35_(), I350);
  task241->add_dep(task249);
  energyq->add_task(task249);

  auto task250 = make_shared<Task250>(I350, v2_, t2);
  task249->add_dep(task250);
  energyq->add_task(task250);

  auto task251 = make_shared<Task251>(I350, v2_, t2);
  task249->add_dep(task251);
  energyq->add_task(task251);

  auto task252 = make_shared<Task252>(I350, v2_, t2);
  task249->add_dep(task252);
  energyq->add_task(task252);

  auto task253 = make_shared<Task253>(I350, v2_, t2);
  task249->add_dep(task253);
  energyq->add_task(task253);

  auto task254 = make_shared<Task254>(I350, v2_, t2);
  task249->add_dep(task254);
  energyq->add_task(task254);

  auto I353 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task255 = make_shared<Task255>(Gamma29_(), I353);
  task241->add_dep(task255);
  energyq->add_task(task255);

  auto task256 = make_shared<Task256>(I353, v2_, t2);
  task255->add_dep(task256);
  energyq->add_task(task256);

  auto I356 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task257 = make_shared<Task257>(Gamma32_(), I356);
  task241->add_dep(task257);
  energyq->add_task(task257);

  auto task258 = make_shared<Task258>(I356, v2_, t2);
  task257->add_dep(task258);
  energyq->add_task(task258);

  auto I365 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task259 = make_shared<Task259>(Gamma7_(), I365);
  task241->add_dep(task259);
  energyq->add_task(task259);

  auto task260 = make_shared<Task260>(I365, v2_, t2);
  task259->add_dep(task260);
  energyq->add_task(task260);

  auto I374 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task261 = make_shared<Task261>(Gamma59_(), I374);
  task241->add_dep(task261);
  energyq->add_task(task261);

  auto task262 = make_shared<Task262>(I374, v2_, t2);
  task261->add_dep(task262);
  energyq->add_task(task262);

  auto I377 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task263 = make_shared<Task263>(Gamma57_(), I377);
  task241->add_dep(task263);
  energyq->add_task(task263);

  auto task264 = make_shared<Task264>(I377, v2_, t2);
  task263->add_dep(task264);
  energyq->add_task(task264);

  shared_ptr<Task265> task265;
  if (diagonal) {
    task265 = make_shared<Task265>(v2_, t2);
    task241->add_dep(task265);
    energyq->add_task(task265);
  }

  shared_ptr<Task266> task266;
  if (diagonal) {
    task266 = make_shared<Task266>(v2_, t2);
    task241->add_dep(task266);
    energyq->add_task(task266);
  }

  auto I384 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task267 = make_shared<Task267>(Gamma38_(), I384);
  task241->add_dep(task267);
  energyq->add_task(task267);

  auto task268 = make_shared<Task268>(I384, v2_, t2);
  task267->add_dep(task268);
  energyq->add_task(task268);

  auto task269 = make_shared<Task269>(I384, v2_, t2);
  task267->add_dep(task269);
  energyq->add_task(task269);

  auto task270 = make_shared<Task270>(I384, h1_, t2);
  task267->add_dep(task270);
  energyq->add_task(task270);

  auto task271 = make_shared<Task271>(I384, h1_, t2);
  task267->add_dep(task271);
  energyq->add_task(task271);

  auto I390 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task272 = make_shared<Task272>(Gamma60_(), I390);
  task241->add_dep(task272);
  energyq->add_task(task272);

  auto task273 = make_shared<Task273>(I390, v2_, t2);
  task272->add_dep(task273);
  energyq->add_task(task273);

  auto task274 = make_shared<Task274>(I390, h1_, t2);
  task272->add_dep(task274);
  energyq->add_task(task274);

  auto I393 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task275 = make_shared<Task275>(h1_, I393);
  task241->add_dep(task275);
  energyq->add_task(task275);

  auto task276 = make_shared<Task276>(I393, Gamma7_(), t2);
  task275->add_dep(task276);
  energyq->add_task(task276);

  return energyq;
}


#endif
