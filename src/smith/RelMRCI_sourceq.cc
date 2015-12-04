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
  auto task728 = make_shared<Task728>(s, reset);
  sourceq->add_task(task728);

  auto I1312 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task729 = make_shared<Task729>(s, I1312);
  task729->add_dep(task728);
  sourceq->add_task(task729);

  auto task730 = make_shared<Task730>(I1312, Gamma5_(), h1_);
  task729->add_dep(task730);
  task730->add_dep(task728);
  sourceq->add_task(task730);

  auto task731 = make_shared<Task731>(I1312, Gamma81_(), v2_);
  task729->add_dep(task731);
  task731->add_dep(task728);
  sourceq->add_task(task731);

  auto task732 = make_shared<Task732>(I1312, Gamma4_(), v2_);
  task729->add_dep(task732);
  task732->add_dep(task728);
  sourceq->add_task(task732);

  auto I1314 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task733 = make_shared<Task733>(s, I1314);
  task733->add_dep(task728);
  sourceq->add_task(task733);

  auto task734 = make_shared<Task734>(I1314, Gamma27_(), h1_);
  task733->add_dep(task734);
  task734->add_dep(task728);
  sourceq->add_task(task734);

  auto I1329 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task735 = make_shared<Task735>(I1314, Gamma24_(), I1329);
  task733->add_dep(task735);
  task735->add_dep(task728);
  sourceq->add_task(task735);

  auto task736 = make_shared<Task736>(I1329, v2_);
  task735->add_dep(task736);
  task736->add_dep(task728);
  sourceq->add_task(task736);

  auto task737 = make_shared<Task737>(I1314, Gamma5_(), v2_);
  task733->add_dep(task737);
  task737->add_dep(task728);
  sourceq->add_task(task737);

  auto I1316 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task738 = make_shared<Task738>(s, I1316);
  task738->add_dep(task728);
  sourceq->add_task(task738);

  auto task739 = make_shared<Task739>(I1316, Gamma33_(), h1_);
  task738->add_dep(task739);
  task739->add_dep(task728);
  sourceq->add_task(task739);

  auto task740 = make_shared<Task740>(I1316, Gamma32_(), v2_);
  task738->add_dep(task740);
  task740->add_dep(task728);
  sourceq->add_task(task740);

  auto task741 = make_shared<Task741>(I1316, Gamma31_(), v2_);
  task738->add_dep(task741);
  task741->add_dep(task728);
  sourceq->add_task(task741);

  auto I1318 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task742 = make_shared<Task742>(s, I1318);
  task742->add_dep(task728);
  sourceq->add_task(task742);

  auto task743 = make_shared<Task743>(I1318, Gamma0_(), v2_);
  task742->add_dep(task743);
  task743->add_dep(task728);
  sourceq->add_task(task743);

  auto I1324 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, virt_, active_});
  auto task744 = make_shared<Task744>(s, I1324);
  task744->add_dep(task728);
  sourceq->add_task(task744);

  auto I1325 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, closed_, virt_});
  auto task745 = make_shared<Task745>(I1324, Gamma11_(), I1325);
  task744->add_dep(task745);
  task745->add_dep(task728);
  sourceq->add_task(task745);

  auto task746 = make_shared<Task746>(I1325, v2_);
  task745->add_dep(task746);
  task746->add_dep(task728);
  sourceq->add_task(task746);

  shared_ptr<TATensor<std::complex<double>,4>> I1340;
  if (diagonal) {
    I1340 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task747> task747;
  if (diagonal) {
    task747 = make_shared<Task747>(s, I1340);
    task747->add_dep(task728);
    sourceq->add_task(task747);
  }

  shared_ptr<Task748> task748;
  if (diagonal) {
    task748 = make_shared<Task748>(I1340, v2_);
    task747->add_dep(task748);
    task748->add_dep(task728);
    sourceq->add_task(task748);
  }

  auto I1342 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task749 = make_shared<Task749>(s, I1342);
  task749->add_dep(task728);
  sourceq->add_task(task749);

  auto I1343 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task750 = make_shared<Task750>(I1342, Gamma27_(), I1343);
  task749->add_dep(task750);
  task750->add_dep(task728);
  sourceq->add_task(task750);

  auto task751 = make_shared<Task751>(I1343, v2_);
  task750->add_dep(task751);
  task751->add_dep(task728);
  sourceq->add_task(task751);

  auto I1346 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task752 = make_shared<Task752>(s, I1346);
  task752->add_dep(task728);
  sourceq->add_task(task752);

  auto task753 = make_shared<Task753>(I1346, Gamma33_(), v2_);
  task752->add_dep(task753);
  task753->add_dep(task728);
  sourceq->add_task(task753);

  return sourceq;
}


#endif
