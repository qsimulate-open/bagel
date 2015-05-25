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

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor937 = vector<shared_ptr<Tensor>>{s};
  auto task937 = make_shared<Task937>(tensor937, reset);
  sourceq->add_task(task937);

  vector<IndexRange> I1732_index = {active_, active_, active_, closed_};
  auto I1732 = make_shared<Tensor>(I1732_index);
  auto tensor938 = vector<shared_ptr<Tensor>>{s, I1732};
  auto task938 = make_shared<Task938>(tensor938, pindex);
  task938->add_dep(task937);
  sourceq->add_task(task938);

  auto tensor939 = vector<shared_ptr<Tensor>>{I1732, h1_, Gamma5_()};
  auto task939 = make_shared<Task939>(tensor939, pindex);
  task938->add_dep(task939);
  task939->add_dep(task937);
  sourceq->add_task(task939);

  auto tensor940 = vector<shared_ptr<Tensor>>{I1732, v2_, Gamma108_()};
  auto task940 = make_shared<Task940>(tensor940, pindex);
  task938->add_dep(task940);
  task940->add_dep(task937);
  sourceq->add_task(task940);

  auto tensor941 = vector<shared_ptr<Tensor>>{I1732, v2_, Gamma4_()};
  auto task941 = make_shared<Task941>(tensor941, pindex);
  task938->add_dep(task941);
  task941->add_dep(task937);
  sourceq->add_task(task941);

  vector<IndexRange> I1734_index = {active_, active_, closed_, virt_};
  auto I1734 = make_shared<Tensor>(I1734_index);
  auto tensor942 = vector<shared_ptr<Tensor>>{s, I1734};
  auto task942 = make_shared<Task942>(tensor942, pindex);
  task942->add_dep(task937);
  sourceq->add_task(task942);

  auto tensor943 = vector<shared_ptr<Tensor>>{I1734, h1_, Gamma34_()};
  auto task943 = make_shared<Task943>(tensor943, pindex);
  task942->add_dep(task943);
  task943->add_dep(task937);
  sourceq->add_task(task943);

  auto tensor944 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma31_()};
  auto task944 = make_shared<Task944>(tensor944, pindex);
  task942->add_dep(task944);
  task944->add_dep(task937);
  sourceq->add_task(task944);

  auto tensor945 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma25_()};
  auto task945 = make_shared<Task945>(tensor945, pindex);
  task942->add_dep(task945);
  task945->add_dep(task937);
  sourceq->add_task(task945);

  auto tensor946 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma28_()};
  auto task946 = make_shared<Task946>(tensor946, pindex);
  task942->add_dep(task946);
  task946->add_dep(task937);
  sourceq->add_task(task946);

  auto tensor947 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma31_()};
  auto task947 = make_shared<Task947>(tensor947, pindex);
  task942->add_dep(task947);
  task947->add_dep(task937);
  sourceq->add_task(task947);

  vector<IndexRange> I1736_index = {active_, active_, closed_, virt_};
  auto I1736 = make_shared<Tensor>(I1736_index);
  auto tensor948 = vector<shared_ptr<Tensor>>{s, I1736};
  auto task948 = make_shared<Task948>(tensor948, pindex);
  task948->add_dep(task937);
  sourceq->add_task(task948);

  auto tensor949 = vector<shared_ptr<Tensor>>{I1736, h1_, Gamma34_()};
  auto task949 = make_shared<Task949>(tensor949, pindex);
  task948->add_dep(task949);
  task949->add_dep(task937);
  sourceq->add_task(task949);

  auto tensor950 = vector<shared_ptr<Tensor>>{I1736, v2_, Gamma31_()};
  auto task950 = make_shared<Task950>(tensor950, pindex);
  task948->add_dep(task950);
  task950->add_dep(task937);
  sourceq->add_task(task950);

  auto tensor951 = vector<shared_ptr<Tensor>>{I1736, v2_, Gamma5_()};
  auto task951 = make_shared<Task951>(tensor951, pindex);
  task948->add_dep(task951);
  task951->add_dep(task937);
  sourceq->add_task(task951);

  auto tensor952 = vector<shared_ptr<Tensor>>{I1736, v2_, Gamma31_()};
  auto task952 = make_shared<Task952>(tensor952, pindex);
  task948->add_dep(task952);
  task952->add_dep(task937);
  sourceq->add_task(task952);

  auto tensor953 = vector<shared_ptr<Tensor>>{I1736, v2_, Gamma31_()};
  auto task953 = make_shared<Task953>(tensor953, pindex);
  task948->add_dep(task953);
  task953->add_dep(task937);
  sourceq->add_task(task953);

  vector<IndexRange> I1738_index = {virt_, active_, active_, active_};
  auto I1738 = make_shared<Tensor>(I1738_index);
  auto tensor954 = vector<shared_ptr<Tensor>>{s, I1738};
  auto task954 = make_shared<Task954>(tensor954, pindex);
  task954->add_dep(task937);
  sourceq->add_task(task954);

  auto tensor955 = vector<shared_ptr<Tensor>>{I1738, Gamma55_(), h1_};
  auto task955 = make_shared<Task955>(tensor955, pindex);
  task954->add_dep(task955);
  task955->add_dep(task937);
  sourceq->add_task(task955);

  auto tensor956 = vector<shared_ptr<Tensor>>{I1738, v2_, Gamma54_()};
  auto task956 = make_shared<Task956>(tensor956, pindex);
  task954->add_dep(task956);
  task956->add_dep(task937);
  sourceq->add_task(task956);

  auto tensor957 = vector<shared_ptr<Tensor>>{I1738, Gamma53_(), v2_};
  auto task957 = make_shared<Task957>(tensor957, pindex);
  task954->add_dep(task957);
  task957->add_dep(task937);
  sourceq->add_task(task957);

  vector<IndexRange> I1740_index = {closed_, closed_, active_, active_};
  auto I1740 = make_shared<Tensor>(I1740_index);
  auto tensor958 = vector<shared_ptr<Tensor>>{s, I1740};
  auto task958 = make_shared<Task958>(tensor958, pindex);
  task958->add_dep(task937);
  sourceq->add_task(task958);

  auto tensor959 = vector<shared_ptr<Tensor>>{I1740, Gamma0_(), v2_};
  auto task959 = make_shared<Task959>(tensor959, pindex);
  task958->add_dep(task959);
  task959->add_dep(task937);
  sourceq->add_task(task959);

  vector<IndexRange> I1746_index = {active_, closed_, closed_, virt_};
  auto I1746 = make_shared<Tensor>(I1746_index);
  auto tensor960 = vector<shared_ptr<Tensor>>{s, I1746};
  auto task960 = make_shared<Task960>(tensor960, pindex);
  task960->add_dep(task937);
  sourceq->add_task(task960);

  auto tensor961 = vector<shared_ptr<Tensor>>{I1746, v2_, Gamma12_()};
  auto task961 = make_shared<Task961>(tensor961, pindex);
  task960->add_dep(task961);
  task961->add_dep(task937);
  sourceq->add_task(task961);

  auto tensor962 = vector<shared_ptr<Tensor>>{I1746, v2_, Gamma12_()};
  auto task962 = make_shared<Task962>(tensor962, pindex);
  task960->add_dep(task962);
  task962->add_dep(task937);
  sourceq->add_task(task962);

  shared_ptr<Tensor> I1770;
  if (diagonal) {
    vector<IndexRange> I1770_index = {closed_, virt_, closed_, virt_};
    I1770 = make_shared<Tensor>(I1770_index);
  }
  shared_ptr<Task963> task963;
  if (diagonal) {
    auto tensor963 = vector<shared_ptr<Tensor>>{s, I1770};
    task963 = make_shared<Task963>(tensor963, pindex);
    task963->add_dep(task937);
    sourceq->add_task(task963);
  }

  shared_ptr<Task964> task964;
  if (diagonal) {
    auto tensor964 = vector<shared_ptr<Tensor>>{I1770, v2_};
    task964 = make_shared<Task964>(tensor964, pindex);
    task963->add_dep(task964);
    task964->add_dep(task937);
    sourceq->add_task(task964);
  }

  vector<IndexRange> I1772_index = {active_, virt_, closed_, virt_};
  auto I1772 = make_shared<Tensor>(I1772_index);
  auto tensor965 = vector<shared_ptr<Tensor>>{s, I1772};
  auto task965 = make_shared<Task965>(tensor965, pindex);
  task965->add_dep(task937);
  sourceq->add_task(task965);

  auto tensor966 = vector<shared_ptr<Tensor>>{I1772, v2_, Gamma34_()};
  auto task966 = make_shared<Task966>(tensor966, pindex);
  task965->add_dep(task966);
  task966->add_dep(task937);
  sourceq->add_task(task966);

  auto tensor967 = vector<shared_ptr<Tensor>>{I1772, v2_, Gamma34_()};
  auto task967 = make_shared<Task967>(tensor967, pindex);
  task965->add_dep(task967);
  task967->add_dep(task937);
  sourceq->add_task(task967);

  vector<IndexRange> I1776_index = {active_, active_, virt_, virt_};
  auto I1776 = make_shared<Tensor>(I1776_index);
  auto tensor968 = vector<shared_ptr<Tensor>>{s, I1776};
  auto task968 = make_shared<Task968>(tensor968, pindex);
  task968->add_dep(task937);
  sourceq->add_task(task968);

  auto tensor969 = vector<shared_ptr<Tensor>>{I1776, v2_, Gamma55_()};
  auto task969 = make_shared<Task969>(tensor969, pindex);
  task968->add_dep(task969);
  task969->add_dep(task937);
  sourceq->add_task(task969);

  return sourceq;
}


#endif
