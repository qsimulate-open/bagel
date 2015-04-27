//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_sourceqq.cc
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


#include <src/smith/MRCI.h>
#include <src/smith/MRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor936 = vector<shared_ptr<Tensor>>{s};
  auto task936 = make_shared<Task936>(tensor936, reset);
  sourceq->add_task(task936);

  vector<IndexRange> I1732_index = {closed_, active_, active_, active_};
  auto I1732 = make_shared<Tensor>(I1732_index);
  auto tensor937 = vector<shared_ptr<Tensor>>{s, I1732};
  auto task937 = make_shared<Task937>(tensor937, pindex);
  task937->add_dep(task936);
  sourceq->add_task(task937);

  auto tensor938 = vector<shared_ptr<Tensor>>{I1732, Gamma5_(), h1_};
  auto task938 = make_shared<Task938>(tensor938, pindex);
  task937->add_dep(task938);
  task938->add_dep(task936);
  sourceq->add_task(task938);

  auto tensor939 = vector<shared_ptr<Tensor>>{I1732, v2_, Gamma104_()};
  auto task939 = make_shared<Task939>(tensor939, pindex);
  task937->add_dep(task939);
  task939->add_dep(task936);
  sourceq->add_task(task939);

  auto tensor940 = vector<shared_ptr<Tensor>>{I1732, v2_, Gamma4_()};
  auto task940 = make_shared<Task940>(tensor940, pindex);
  task937->add_dep(task940);
  task940->add_dep(task936);
  sourceq->add_task(task940);

  vector<IndexRange> I1734_index = {active_, active_, closed_, virt_};
  auto I1734 = make_shared<Tensor>(I1734_index);
  auto tensor941 = vector<shared_ptr<Tensor>>{s, I1734};
  auto task941 = make_shared<Task941>(tensor941, pindex);
  task941->add_dep(task936);
  sourceq->add_task(task941);

  auto tensor942 = vector<shared_ptr<Tensor>>{I1734, h1_, Gamma32_()};
  auto task942 = make_shared<Task942>(tensor942, pindex);
  task941->add_dep(task942);
  task942->add_dep(task936);
  sourceq->add_task(task942);

  auto tensor943 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma29_()};
  auto task943 = make_shared<Task943>(tensor943, pindex);
  task941->add_dep(task943);
  task943->add_dep(task936);
  sourceq->add_task(task943);

  auto tensor944 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma25_()};
  auto task944 = make_shared<Task944>(tensor944, pindex);
  task941->add_dep(task944);
  task944->add_dep(task936);
  sourceq->add_task(task944);

  auto tensor945 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma27_()};
  auto task945 = make_shared<Task945>(tensor945, pindex);
  task941->add_dep(task945);
  task945->add_dep(task936);
  sourceq->add_task(task945);

  auto tensor946 = vector<shared_ptr<Tensor>>{I1734, Gamma29_(), v2_};
  auto task946 = make_shared<Task946>(tensor946, pindex);
  task941->add_dep(task946);
  task946->add_dep(task936);
  sourceq->add_task(task946);

  vector<IndexRange> I1736_index = {active_, active_, closed_, virt_};
  auto I1736 = make_shared<Tensor>(I1736_index);
  auto tensor947 = vector<shared_ptr<Tensor>>{s, I1736};
  auto task947 = make_shared<Task947>(tensor947, pindex);
  task947->add_dep(task936);
  sourceq->add_task(task947);

  auto tensor948 = vector<shared_ptr<Tensor>>{I1736, h1_, Gamma32_()};
  auto task948 = make_shared<Task948>(tensor948, pindex);
  task947->add_dep(task948);
  task948->add_dep(task936);
  sourceq->add_task(task948);

  auto tensor949 = vector<shared_ptr<Tensor>>{I1736, v2_, Gamma29_()};
  auto task949 = make_shared<Task949>(tensor949, pindex);
  task947->add_dep(task949);
  task949->add_dep(task936);
  sourceq->add_task(task949);

  auto tensor950 = vector<shared_ptr<Tensor>>{I1736, Gamma5_(), v2_};
  auto task950 = make_shared<Task950>(tensor950, pindex);
  task947->add_dep(task950);
  task950->add_dep(task936);
  sourceq->add_task(task950);

  auto tensor951 = vector<shared_ptr<Tensor>>{I1736, v2_, Gamma29_()};
  auto task951 = make_shared<Task951>(tensor951, pindex);
  task947->add_dep(task951);
  task951->add_dep(task936);
  sourceq->add_task(task951);

  auto tensor952 = vector<shared_ptr<Tensor>>{I1736, Gamma29_(), v2_};
  auto task952 = make_shared<Task952>(tensor952, pindex);
  task947->add_dep(task952);
  task952->add_dep(task936);
  sourceq->add_task(task952);

  vector<IndexRange> I1738_index = {virt_, active_, active_, active_};
  auto I1738 = make_shared<Tensor>(I1738_index);
  auto tensor953 = vector<shared_ptr<Tensor>>{s, I1738};
  auto task953 = make_shared<Task953>(tensor953, pindex);
  task953->add_dep(task936);
  sourceq->add_task(task953);

  auto tensor954 = vector<shared_ptr<Tensor>>{I1738, Gamma51_(), h1_};
  auto task954 = make_shared<Task954>(tensor954, pindex);
  task953->add_dep(task954);
  task954->add_dep(task936);
  sourceq->add_task(task954);

  auto tensor955 = vector<shared_ptr<Tensor>>{I1738, Gamma50_(), v2_};
  auto task955 = make_shared<Task955>(tensor955, pindex);
  task953->add_dep(task955);
  task955->add_dep(task936);
  sourceq->add_task(task955);

  auto tensor956 = vector<shared_ptr<Tensor>>{I1738, v2_, Gamma49_()};
  auto task956 = make_shared<Task956>(tensor956, pindex);
  task953->add_dep(task956);
  task956->add_dep(task936);
  sourceq->add_task(task956);

  vector<IndexRange> I1740_index = {active_, active_, closed_, closed_};
  auto I1740 = make_shared<Tensor>(I1740_index);
  auto tensor957 = vector<shared_ptr<Tensor>>{s, I1740};
  auto task957 = make_shared<Task957>(tensor957, pindex);
  task957->add_dep(task936);
  sourceq->add_task(task957);

  auto tensor958 = vector<shared_ptr<Tensor>>{I1740, v2_, Gamma0_()};
  auto task958 = make_shared<Task958>(tensor958, pindex);
  task957->add_dep(task958);
  task958->add_dep(task936);
  sourceq->add_task(task958);

  vector<IndexRange> I1746_index = {active_, closed_, closed_, virt_};
  auto I1746 = make_shared<Tensor>(I1746_index);
  auto tensor959 = vector<shared_ptr<Tensor>>{s, I1746};
  auto task959 = make_shared<Task959>(tensor959, pindex);
  task959->add_dep(task936);
  sourceq->add_task(task959);

  auto tensor960 = vector<shared_ptr<Tensor>>{I1746, v2_, Gamma12_()};
  auto task960 = make_shared<Task960>(tensor960, pindex);
  task959->add_dep(task960);
  task960->add_dep(task936);
  sourceq->add_task(task960);

  auto tensor961 = vector<shared_ptr<Tensor>>{I1746, v2_, Gamma12_()};
  auto task961 = make_shared<Task961>(tensor961, pindex);
  task959->add_dep(task961);
  task961->add_dep(task936);
  sourceq->add_task(task961);

  shared_ptr<Tensor> I1770;
  if (diagonal) {
    vector<IndexRange> I1770_index = {closed_, virt_, closed_, virt_};
    I1770 = make_shared<Tensor>(I1770_index);
  }
  shared_ptr<Task962> task962;
  if (diagonal) {
    auto tensor962 = vector<shared_ptr<Tensor>>{s, I1770};
    task962 = make_shared<Task962>(tensor962, pindex);
    task962->add_dep(task936);
    sourceq->add_task(task962);
  }

  shared_ptr<Task963> task963;
  if (diagonal) {
    auto tensor963 = vector<shared_ptr<Tensor>>{I1770, v2_};
    task963 = make_shared<Task963>(tensor963, pindex);
    task962->add_dep(task963);
    task963->add_dep(task936);
    sourceq->add_task(task963);
  }

  vector<IndexRange> I1772_index = {active_, virt_, closed_, virt_};
  auto I1772 = make_shared<Tensor>(I1772_index);
  auto tensor964 = vector<shared_ptr<Tensor>>{s, I1772};
  auto task964 = make_shared<Task964>(tensor964, pindex);
  task964->add_dep(task936);
  sourceq->add_task(task964);

  auto tensor965 = vector<shared_ptr<Tensor>>{I1772, v2_, Gamma32_()};
  auto task965 = make_shared<Task965>(tensor965, pindex);
  task964->add_dep(task965);
  task965->add_dep(task936);
  sourceq->add_task(task965);

  auto tensor966 = vector<shared_ptr<Tensor>>{I1772, v2_, Gamma32_()};
  auto task966 = make_shared<Task966>(tensor966, pindex);
  task964->add_dep(task966);
  task966->add_dep(task936);
  sourceq->add_task(task966);

  vector<IndexRange> I1776_index = {virt_, virt_, active_, active_};
  auto I1776 = make_shared<Tensor>(I1776_index);
  auto tensor967 = vector<shared_ptr<Tensor>>{s, I1776};
  auto task967 = make_shared<Task967>(tensor967, pindex);
  task967->add_dep(task936);
  sourceq->add_task(task967);

  auto tensor968 = vector<shared_ptr<Tensor>>{I1776, Gamma51_(), v2_};
  auto task968 = make_shared<Task968>(tensor968, pindex);
  task967->add_dep(task968);
  task968->add_dep(task936);
  sourceq->add_task(task968);

  return sourceq;
}


#endif
