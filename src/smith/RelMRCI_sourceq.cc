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
  auto task729 = make_shared<Task729>(s, reset);
  sourceq->add_task(task729);

  auto I1312 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task730 = make_shared<Task730>(s, I1312);
  task730->add_dep(task729);
  sourceq->add_task(task730);

  auto task731 = make_shared<Task731>(I1312, Gamma5_(), h1_);
  task730->add_dep(task731);
  task731->add_dep(task729);
  sourceq->add_task(task731);

  auto task732 = make_shared<Task732>(I1312, Gamma81_(), v2_);
  task730->add_dep(task732);
  task732->add_dep(task729);
  sourceq->add_task(task732);

  auto task733 = make_shared<Task733>(I1312, Gamma4_(), v2_);
  task730->add_dep(task733);
  task733->add_dep(task729);
  sourceq->add_task(task733);

  auto I1314 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task734 = make_shared<Task734>(s, I1314);
  task734->add_dep(task729);
  sourceq->add_task(task734);

  auto task735 = make_shared<Task735>(I1314, Gamma27_(), h1_);
  task734->add_dep(task735);
  task735->add_dep(task729);
  sourceq->add_task(task735);

  auto I1329 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task736 = make_shared<Task736>(I1314, Gamma24_(), I1329);
  task734->add_dep(task736);
  task736->add_dep(task729);
  sourceq->add_task(task736);

  auto task737 = make_shared<Task737>(I1329, v2_);
  task736->add_dep(task737);
  task737->add_dep(task729);
  sourceq->add_task(task737);

  auto task738 = make_shared<Task738>(I1314, Gamma5_(), v2_);
  task734->add_dep(task738);
  task738->add_dep(task729);
  sourceq->add_task(task738);

  auto I1316 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task739 = make_shared<Task739>(s, I1316);
  task739->add_dep(task729);
  sourceq->add_task(task739);

  auto task740 = make_shared<Task740>(I1316, Gamma33_(), h1_);
  task739->add_dep(task740);
  task740->add_dep(task729);
  sourceq->add_task(task740);

  auto task741 = make_shared<Task741>(I1316, Gamma32_(), v2_);
  task739->add_dep(task741);
  task741->add_dep(task729);
  sourceq->add_task(task741);

  auto task742 = make_shared<Task742>(I1316, Gamma31_(), v2_);
  task739->add_dep(task742);
  task742->add_dep(task729);
  sourceq->add_task(task742);

  auto I1318 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task743 = make_shared<Task743>(s, I1318);
  task743->add_dep(task729);
  sourceq->add_task(task743);

  auto task744 = make_shared<Task744>(I1318, Gamma0_(), v2_);
  task743->add_dep(task744);
  task744->add_dep(task729);
  sourceq->add_task(task744);

  auto I1324 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, virt_, active_});
  auto task745 = make_shared<Task745>(s, I1324);
  task745->add_dep(task729);
  sourceq->add_task(task745);

  auto I1325 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, closed_, virt_});
  auto task746 = make_shared<Task746>(I1324, Gamma11_(), I1325);
  task745->add_dep(task746);
  task746->add_dep(task729);
  sourceq->add_task(task746);

  auto task747 = make_shared<Task747>(I1325, v2_);
  task746->add_dep(task747);
  task747->add_dep(task729);
  sourceq->add_task(task747);

  shared_ptr<TATensor<std::complex<double>,4>> I1340;
  if (diagonal) {
    I1340 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task748> task748;
  if (diagonal) {
    task748 = make_shared<Task748>(s, I1340);
    task748->add_dep(task729);
    sourceq->add_task(task748);
  }

  shared_ptr<Task749> task749;
  if (diagonal) {
    task749 = make_shared<Task749>(I1340, v2_);
    task748->add_dep(task749);
    task749->add_dep(task729);
    sourceq->add_task(task749);
  }

  auto I1342 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task750 = make_shared<Task750>(s, I1342);
  task750->add_dep(task729);
  sourceq->add_task(task750);

  auto I1343 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task751 = make_shared<Task751>(I1342, Gamma27_(), I1343);
  task750->add_dep(task751);
  task751->add_dep(task729);
  sourceq->add_task(task751);

  auto task752 = make_shared<Task752>(I1343, v2_);
  task751->add_dep(task752);
  task752->add_dep(task729);
  sourceq->add_task(task752);

  auto I1346 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task753 = make_shared<Task753>(s, I1346);
  task753->add_dep(task729);
  sourceq->add_task(task753);

  auto task754 = make_shared<Task754>(I1346, Gamma33_(), v2_);
  task753->add_dep(task754);
  task754->add_dep(task729);
  sourceq->add_task(task754);

  return sourceq;
}


#endif
