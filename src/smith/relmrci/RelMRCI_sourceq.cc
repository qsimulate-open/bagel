//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_sourceqq.cc
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


#include <src/smith/RelMRCI.h>
#include <src/smith/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_sourceq(const bool reset, const bool diagonal) {

  auto sourceq = make_shared<Queue>();
  auto task722 = make_shared<Task722>(s, reset);
  sourceq->add_task(task722);

  auto I1299 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task723 = make_shared<Task723>(s, I1299);
  task723->add_dep(task722);
  sourceq->add_task(task723);

  auto task724 = make_shared<Task724>(I1299, Gamma5_(), h1_);
  task723->add_dep(task724);
  task724->add_dep(task722);
  sourceq->add_task(task724);

  auto task725 = make_shared<Task725>(I1299, Gamma81_(), v2_);
  task723->add_dep(task725);
  task725->add_dep(task722);
  sourceq->add_task(task725);

  auto task726 = make_shared<Task726>(I1299, Gamma4_(), v2_);
  task723->add_dep(task726);
  task726->add_dep(task722);
  sourceq->add_task(task726);

  auto I1301 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task727 = make_shared<Task727>(s, I1301);
  task727->add_dep(task722);
  sourceq->add_task(task727);

  auto task728 = make_shared<Task728>(I1301, Gamma27_(), h1_);
  task727->add_dep(task728);
  task728->add_dep(task722);
  sourceq->add_task(task728);

  auto I1316 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task729 = make_shared<Task729>(I1301, Gamma24_(), I1316);
  task727->add_dep(task729);
  task729->add_dep(task722);
  sourceq->add_task(task729);

  auto task730 = make_shared<Task730>(I1316, v2_);
  task729->add_dep(task730);
  task730->add_dep(task722);
  sourceq->add_task(task730);

  auto task731 = make_shared<Task731>(I1301, Gamma5_(), v2_);
  task727->add_dep(task731);
  task731->add_dep(task722);
  sourceq->add_task(task731);

  auto I1303 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task732 = make_shared<Task732>(s, I1303);
  task732->add_dep(task722);
  sourceq->add_task(task732);

  auto task733 = make_shared<Task733>(I1303, Gamma33_(), h1_);
  task732->add_dep(task733);
  task733->add_dep(task722);
  sourceq->add_task(task733);

  auto task734 = make_shared<Task734>(I1303, Gamma32_(), v2_);
  task732->add_dep(task734);
  task734->add_dep(task722);
  sourceq->add_task(task734);

  auto task735 = make_shared<Task735>(I1303, Gamma31_(), v2_);
  task732->add_dep(task735);
  task735->add_dep(task722);
  sourceq->add_task(task735);

  auto I1305 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task736 = make_shared<Task736>(s, I1305);
  task736->add_dep(task722);
  sourceq->add_task(task736);

  auto task737 = make_shared<Task737>(I1305, Gamma0_(), v2_);
  task736->add_dep(task737);
  task737->add_dep(task722);
  sourceq->add_task(task737);

  auto I1311 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, virt_, active_});
  auto task738 = make_shared<Task738>(s, I1311);
  task738->add_dep(task722);
  sourceq->add_task(task738);

  auto I1312 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, closed_, virt_});
  auto task739 = make_shared<Task739>(I1311, Gamma11_(), I1312);
  task738->add_dep(task739);
  task739->add_dep(task722);
  sourceq->add_task(task739);

  auto task740 = make_shared<Task740>(I1312, v2_);
  task739->add_dep(task740);
  task740->add_dep(task722);
  sourceq->add_task(task740);

  shared_ptr<TATensor<std::complex<double>,4>> I1327;
  if (diagonal) {
    I1327 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task741> task741;
  if (diagonal) {
    task741 = make_shared<Task741>(s, I1327);
    task741->add_dep(task722);
    sourceq->add_task(task741);
  }

  shared_ptr<Task742> task742;
  if (diagonal) {
    task742 = make_shared<Task742>(I1327, v2_);
    task741->add_dep(task742);
    task742->add_dep(task722);
    sourceq->add_task(task742);
  }

  auto I1329 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task743 = make_shared<Task743>(s, I1329);
  task743->add_dep(task722);
  sourceq->add_task(task743);

  auto I1330 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task744 = make_shared<Task744>(I1329, Gamma27_(), I1330);
  task743->add_dep(task744);
  task744->add_dep(task722);
  sourceq->add_task(task744);

  auto task745 = make_shared<Task745>(I1330, v2_);
  task744->add_dep(task745);
  task745->add_dep(task722);
  sourceq->add_task(task745);

  auto I1333 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task746 = make_shared<Task746>(s, I1333);
  task746->add_dep(task722);
  sourceq->add_task(task746);

  auto task747 = make_shared<Task747>(I1333, Gamma33_(), v2_);
  task746->add_dep(task747);
  task747->add_dep(task722);
  sourceq->add_task(task747);

  return sourceq;
}


#endif
