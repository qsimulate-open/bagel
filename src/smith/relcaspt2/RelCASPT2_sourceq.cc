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
  auto task202 = make_shared<Task202>(s, reset);
  sourceq->add_task(task202);

  auto I334 = make_shared<TATensor<std::complex<double>,4>>({closed_, active_, active_, active_});
  auto task203 = make_shared<Task203>(s, I334);
  task203->add_dep(task202);
  sourceq->add_task(task203);

  auto task204 = make_shared<Task204>(I334, Gamma7_(), h1_);
  task203->add_dep(task204);
  task204->add_dep(task202);
  sourceq->add_task(task204);

  auto task205 = make_shared<Task205>(I334, Gamma105_(), v2_);
  task203->add_dep(task205);
  task205->add_dep(task202);
  sourceq->add_task(task205);

  auto task206 = make_shared<Task206>(I334, Gamma6_(), v2_);
  task203->add_dep(task206);
  task206->add_dep(task202);
  sourceq->add_task(task206);

  auto I336 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task207 = make_shared<Task207>(s, I336);
  task207->add_dep(task202);
  sourceq->add_task(task207);

  auto task208 = make_shared<Task208>(I336, Gamma38_(), h1_);
  task207->add_dep(task208);
  task208->add_dep(task202);
  sourceq->add_task(task208);

  auto I353 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task209 = make_shared<Task209>(I336, Gamma35_(), I353);
  task207->add_dep(task209);
  task209->add_dep(task202);
  sourceq->add_task(task209);

  auto task210 = make_shared<Task210>(I353, v2_);
  task209->add_dep(task210);
  task210->add_dep(task202);
  sourceq->add_task(task210);

  auto task211 = make_shared<Task211>(I336, Gamma29_(), v2_);
  task207->add_dep(task211);
  task211->add_dep(task202);
  sourceq->add_task(task211);

  auto task212 = make_shared<Task212>(I336, Gamma32_(), v2_);
  task207->add_dep(task212);
  task212->add_dep(task202);
  sourceq->add_task(task212);

  auto I338 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task213 = make_shared<Task213>(s, I338);
  task213->add_dep(task202);
  sourceq->add_task(task213);

  auto task214 = make_shared<Task214>(I338, Gamma38_(), h1_);
  task213->add_dep(task214);
  task214->add_dep(task202);
  sourceq->add_task(task214);

  auto I361 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task215 = make_shared<Task215>(I338, Gamma35_(), I361);
  task213->add_dep(task215);
  task215->add_dep(task202);
  sourceq->add_task(task215);

  auto task216 = make_shared<Task216>(I361, v2_);
  task215->add_dep(task216);
  task216->add_dep(task202);
  sourceq->add_task(task216);

  auto task217 = make_shared<Task217>(I338, Gamma7_(), v2_);
  task213->add_dep(task217);
  task217->add_dep(task202);
  sourceq->add_task(task217);

  auto I340 = make_shared<TATensor<std::complex<double>,4>>({virt_, active_, active_, active_});
  auto task218 = make_shared<Task218>(s, I340);
  task218->add_dep(task202);
  sourceq->add_task(task218);

  auto task219 = make_shared<Task219>(I340, Gamma60_(), h1_);
  task218->add_dep(task219);
  task219->add_dep(task202);
  sourceq->add_task(task219);

  auto task220 = make_shared<Task220>(I340, Gamma59_(), v2_);
  task218->add_dep(task220);
  task220->add_dep(task202);
  sourceq->add_task(task220);

  auto task221 = make_shared<Task221>(I340, Gamma57_(), v2_);
  task218->add_dep(task221);
  task221->add_dep(task202);
  sourceq->add_task(task221);

  auto I342 = make_shared<TATensor<std::complex<double>,4>>({closed_, closed_, active_, active_});
  auto task222 = make_shared<Task222>(s, I342);
  task222->add_dep(task202);
  sourceq->add_task(task222);

  auto task223 = make_shared<Task223>(I342, Gamma92_(), v2_);
  task222->add_dep(task223);
  task223->add_dep(task202);
  sourceq->add_task(task223);

  auto I348 = make_shared<TATensor<std::complex<double>,4>>({closed_, closed_, virt_, active_});
  auto task224 = make_shared<Task224>(s, I348);
  task224->add_dep(task202);
  sourceq->add_task(task224);

  auto I349 = make_shared<TATensor<std::complex<double>,4>>({closed_, active_, closed_, virt_});
  auto task225 = make_shared<Task225>(I348, Gamma16_(), I349);
  task224->add_dep(task225);
  task225->add_dep(task202);
  sourceq->add_task(task225);

  auto task226 = make_shared<Task226>(I349, v2_);
  task225->add_dep(task226);
  task226->add_dep(task202);
  sourceq->add_task(task226);

  shared_ptr<TATensor<std::complex<double>,4>> I372;
  if (diagonal) {
    I372 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task227> task227;
  if (diagonal) {
    task227 = make_shared<Task227>(s, I372);
    task227->add_dep(task202);
    sourceq->add_task(task227);
  }

  shared_ptr<Task228> task228;
  if (diagonal) {
    task228 = make_shared<Task228>(I372, v2_);
    task227->add_dep(task228);
    task228->add_dep(task202);
    sourceq->add_task(task228);
  }

  auto I374 = make_shared<TATensor<std::complex<double>,4>>({virt_, closed_, virt_, active_});
  auto task229 = make_shared<Task229>(s, I374);
  task229->add_dep(task202);
  sourceq->add_task(task229);

  auto I375 = make_shared<TATensor<std::complex<double>,4>>({active_, virt_, closed_, virt_});
  auto task230 = make_shared<Task230>(I374, Gamma38_(), I375);
  task229->add_dep(task230);
  task230->add_dep(task202);
  sourceq->add_task(task230);

  auto task231 = make_shared<Task231>(I375, v2_);
  task230->add_dep(task231);
  task231->add_dep(task202);
  sourceq->add_task(task231);

  auto I378 = make_shared<TATensor<std::complex<double>,4>>({virt_, virt_, active_, active_});
  auto task232 = make_shared<Task232>(s, I378);
  task232->add_dep(task202);
  sourceq->add_task(task232);

  auto task233 = make_shared<Task233>(I378, Gamma60_(), v2_);
  task232->add_dep(task233);
  task233->add_dep(task202);
  sourceq->add_task(task233);

  return sourceq;
}


#endif
