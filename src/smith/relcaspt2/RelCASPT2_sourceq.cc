//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_sourceqq.cc
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


#include <src/smith/relcaspt2/RelCASPT2.h>
#include <src/smith/relcaspt2/RelCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASPT2::RelCASPT2::make_sourceq(const bool reset, const bool diagonal) {

  auto sourceq = make_shared<Queue>();
  auto task209 = make_shared<Task209>(s, reset);
  sourceq->add_task(task209);

  auto I348 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task210 = make_shared<Task210>(s, I348);
  task210->add_dep(task209);
  sourceq->add_task(task210);

  auto task211 = make_shared<Task211>(I348, Gamma7_(), h1_);
  task210->add_dep(task211);
  task211->add_dep(task209);
  sourceq->add_task(task211);

  auto task212 = make_shared<Task212>(I348, Gamma107_(), v2_);
  task210->add_dep(task212);
  task212->add_dep(task209);
  sourceq->add_task(task212);

  auto task213 = make_shared<Task213>(I348, Gamma6_(), v2_);
  task210->add_dep(task213);
  task213->add_dep(task209);
  sourceq->add_task(task213);

  auto I350 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task214 = make_shared<Task214>(s, I350);
  task214->add_dep(task209);
  sourceq->add_task(task214);

  auto task215 = make_shared<Task215>(I350, Gamma38_(), h1_);
  task214->add_dep(task215);
  task215->add_dep(task209);
  sourceq->add_task(task215);

  auto I367 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task216 = make_shared<Task216>(I350, Gamma35_(), I367);
  task214->add_dep(task216);
  task216->add_dep(task209);
  sourceq->add_task(task216);

  auto task217 = make_shared<Task217>(I367, v2_);
  task216->add_dep(task217);
  task217->add_dep(task209);
  sourceq->add_task(task217);

  auto task218 = make_shared<Task218>(I350, Gamma29_(), v2_);
  task214->add_dep(task218);
  task218->add_dep(task209);
  sourceq->add_task(task218);

  auto task219 = make_shared<Task219>(I350, Gamma32_(), v2_);
  task214->add_dep(task219);
  task219->add_dep(task209);
  sourceq->add_task(task219);

  auto I352 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task220 = make_shared<Task220>(s, I352);
  task220->add_dep(task209);
  sourceq->add_task(task220);

  auto task221 = make_shared<Task221>(I352, Gamma38_(), h1_);
  task220->add_dep(task221);
  task221->add_dep(task209);
  sourceq->add_task(task221);

  auto I375 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task222 = make_shared<Task222>(I352, Gamma35_(), I375);
  task220->add_dep(task222);
  task222->add_dep(task209);
  sourceq->add_task(task222);

  auto task223 = make_shared<Task223>(I375, v2_);
  task222->add_dep(task223);
  task223->add_dep(task209);
  sourceq->add_task(task223);

  auto task224 = make_shared<Task224>(I352, Gamma7_(), v2_);
  task220->add_dep(task224);
  task224->add_dep(task209);
  sourceq->add_task(task224);

  auto I354 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task225 = make_shared<Task225>(s, I354);
  task225->add_dep(task209);
  sourceq->add_task(task225);

  auto task226 = make_shared<Task226>(I354, Gamma60_(), h1_);
  task225->add_dep(task226);
  task226->add_dep(task209);
  sourceq->add_task(task226);

  auto task227 = make_shared<Task227>(I354, Gamma59_(), v2_);
  task225->add_dep(task227);
  task227->add_dep(task209);
  sourceq->add_task(task227);

  auto task228 = make_shared<Task228>(I354, Gamma57_(), v2_);
  task225->add_dep(task228);
  task228->add_dep(task209);
  sourceq->add_task(task228);

  auto I356 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task229 = make_shared<Task229>(s, I356);
  task229->add_dep(task209);
  sourceq->add_task(task229);

  auto task230 = make_shared<Task230>(I356, Gamma94_(), v2_);
  task229->add_dep(task230);
  task230->add_dep(task209);
  sourceq->add_task(task230);

  auto I362 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, virt_, active_});
  auto task231 = make_shared<Task231>(s, I362);
  task231->add_dep(task209);
  sourceq->add_task(task231);

  auto I363 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, closed_, virt_});
  auto task232 = make_shared<Task232>(I362, Gamma16_(), I363);
  task231->add_dep(task232);
  task232->add_dep(task209);
  sourceq->add_task(task232);

  auto task233 = make_shared<Task233>(I363, v2_);
  task232->add_dep(task233);
  task233->add_dep(task209);
  sourceq->add_task(task233);

  shared_ptr<TATensor<std::complex<double>,4>> I386;
  if (diagonal) {
    I386 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task234> task234;
  if (diagonal) {
    task234 = make_shared<Task234>(s, I386);
    task234->add_dep(task209);
    sourceq->add_task(task234);
  }

  shared_ptr<Task235> task235;
  if (diagonal) {
    task235 = make_shared<Task235>(I386, v2_);
    task234->add_dep(task235);
    task235->add_dep(task209);
    sourceq->add_task(task235);
  }

  auto I388 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task236 = make_shared<Task236>(s, I388);
  task236->add_dep(task209);
  sourceq->add_task(task236);

  auto I389 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task237 = make_shared<Task237>(I388, Gamma38_(), I389);
  task236->add_dep(task237);
  task237->add_dep(task209);
  sourceq->add_task(task237);

  auto task238 = make_shared<Task238>(I389, v2_);
  task237->add_dep(task238);
  task238->add_dep(task209);
  sourceq->add_task(task238);

  auto I392 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task239 = make_shared<Task239>(s, I392);
  task239->add_dep(task209);
  sourceq->add_task(task239);

  auto task240 = make_shared<Task240>(I392, Gamma60_(), v2_);
  task239->add_dep(task240);
  task240->add_dep(task209);
  sourceq->add_task(task240);

  return sourceq;
}


#endif
