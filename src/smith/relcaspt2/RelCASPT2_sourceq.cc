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
  auto task190 = make_shared<Task190>(s, reset);
  sourceq->add_task(task190);

  auto I288 = make_shared<TATensor<std::complex<double>,4>>({closed_, active_, active_, active_});
  auto task191 = make_shared<Task191>(s, I288);
  task191->add_dep(task190);
  sourceq->add_task(task191);

  auto task192 = make_shared<Task192>(I288, Gamma7_(), h1_);
  task191->add_dep(task192);
  task192->add_dep(task190);
  sourceq->add_task(task192);

  auto task193 = make_shared<Task193>(I288, Gamma109_(), v2_);
  task191->add_dep(task193);
  task193->add_dep(task190);
  sourceq->add_task(task193);

  auto task194 = make_shared<Task194>(I288, Gamma6_(), v2_);
  task191->add_dep(task194);
  task194->add_dep(task190);
  sourceq->add_task(task194);

  auto I290 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task195 = make_shared<Task195>(s, I290);
  task195->add_dep(task190);
  sourceq->add_task(task195);

  auto task196 = make_shared<Task196>(I290, Gamma38_(), h1_);
  task195->add_dep(task196);
  task196->add_dep(task190);
  sourceq->add_task(task196);

  auto I307 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task197 = make_shared<Task197>(I290, Gamma35_(), I307);
  task195->add_dep(task197);
  task197->add_dep(task190);
  sourceq->add_task(task197);

  auto task198 = make_shared<Task198>(I307, v2_);
  task197->add_dep(task198);
  task198->add_dep(task190);
  sourceq->add_task(task198);

  auto task199 = make_shared<Task199>(I290, Gamma29_(), v2_);
  task195->add_dep(task199);
  task199->add_dep(task190);
  sourceq->add_task(task199);

  auto task200 = make_shared<Task200>(I290, Gamma32_(), v2_);
  task195->add_dep(task200);
  task200->add_dep(task190);
  sourceq->add_task(task200);

  auto I292 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task201 = make_shared<Task201>(s, I292);
  task201->add_dep(task190);
  sourceq->add_task(task201);

  auto task202 = make_shared<Task202>(I292, Gamma38_(), h1_);
  task201->add_dep(task202);
  task202->add_dep(task190);
  sourceq->add_task(task202);

  auto I315 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task203 = make_shared<Task203>(I292, Gamma35_(), I315);
  task201->add_dep(task203);
  task203->add_dep(task190);
  sourceq->add_task(task203);

  auto task204 = make_shared<Task204>(I315, v2_);
  task203->add_dep(task204);
  task204->add_dep(task190);
  sourceq->add_task(task204);

  auto task205 = make_shared<Task205>(I292, Gamma7_(), v2_);
  task201->add_dep(task205);
  task205->add_dep(task190);
  sourceq->add_task(task205);

  auto I294 = make_shared<TATensor<std::complex<double>,4>>({virt_, active_, active_, active_});
  auto task206 = make_shared<Task206>(s, I294);
  task206->add_dep(task190);
  sourceq->add_task(task206);

  auto task207 = make_shared<Task207>(I294, Gamma60_(), h1_);
  task206->add_dep(task207);
  task207->add_dep(task190);
  sourceq->add_task(task207);

  auto task208 = make_shared<Task208>(I294, Gamma59_(), v2_);
  task206->add_dep(task208);
  task208->add_dep(task190);
  sourceq->add_task(task208);

  auto task209 = make_shared<Task209>(I294, Gamma57_(), v2_);
  task206->add_dep(task209);
  task209->add_dep(task190);
  sourceq->add_task(task209);

  auto I296 = make_shared<TATensor<std::complex<double>,4>>({closed_, closed_, active_, active_});
  auto task210 = make_shared<Task210>(s, I296);
  task210->add_dep(task190);
  sourceq->add_task(task210);

  auto task211 = make_shared<Task211>(I296, Gamma92_(), v2_);
  task210->add_dep(task211);
  task211->add_dep(task190);
  sourceq->add_task(task211);

  auto I302 = make_shared<TATensor<std::complex<double>,4>>({closed_, closed_, virt_, active_});
  auto task212 = make_shared<Task212>(s, I302);
  task212->add_dep(task190);
  sourceq->add_task(task212);

  auto I303 = make_shared<TATensor<std::complex<double>,4>>({closed_, active_, closed_, virt_});
  auto task213 = make_shared<Task213>(I302, Gamma16_(), I303);
  task212->add_dep(task213);
  task213->add_dep(task190);
  sourceq->add_task(task213);

  auto task214 = make_shared<Task214>(I303, v2_);
  task213->add_dep(task214);
  task214->add_dep(task190);
  sourceq->add_task(task214);

  shared_ptr<TATensor<std::complex<double>,4>> I326;
  if (diagonal) {
    I326 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task215> task215;
  if (diagonal) {
    task215 = make_shared<Task215>(s, I326);
    task215->add_dep(task190);
    sourceq->add_task(task215);
  }

  shared_ptr<Task216> task216;
  if (diagonal) {
    task216 = make_shared<Task216>(I326, v2_);
    task215->add_dep(task216);
    task216->add_dep(task190);
    sourceq->add_task(task216);
  }

  auto I328 = make_shared<TATensor<std::complex<double>,4>>({virt_, closed_, virt_, active_});
  auto task217 = make_shared<Task217>(s, I328);
  task217->add_dep(task190);
  sourceq->add_task(task217);

  auto I329 = make_shared<TATensor<std::complex<double>,4>>({active_, virt_, closed_, virt_});
  auto task218 = make_shared<Task218>(I328, Gamma38_(), I329);
  task217->add_dep(task218);
  task218->add_dep(task190);
  sourceq->add_task(task218);

  auto task219 = make_shared<Task219>(I329, v2_);
  task218->add_dep(task219);
  task219->add_dep(task190);
  sourceq->add_task(task219);

  auto I332 = make_shared<TATensor<std::complex<double>,4>>({virt_, virt_, active_, active_});
  auto task220 = make_shared<Task220>(s, I332);
  task220->add_dep(task190);
  sourceq->add_task(task220);

  auto task221 = make_shared<Task221>(I332, Gamma60_(), v2_);
  task220->add_dep(task221);
  task221->add_dep(task190);
  sourceq->add_task(task221);

  return sourceq;
}


#endif
