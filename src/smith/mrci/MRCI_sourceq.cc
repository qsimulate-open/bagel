//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_sourceqq.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks19.h>
#include <src/smith/mrci/MRCI_tasks20.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor934 = vector<shared_ptr<Tensor>>{s};
  auto task934 = make_shared<Task934>(tensor934, reset);
  sourceq->add_task(task934);

  vector<IndexRange> I1728_index = {active_, active_, active_, closed_};
  auto I1728 = make_shared<Tensor>(I1728_index);
  auto tensor935 = vector<shared_ptr<Tensor>>{s, I1728};
  auto task935 = make_shared<Task935>(tensor935, pindex);
  task935->add_dep(task934);
  sourceq->add_task(task935);

  auto tensor936 = vector<shared_ptr<Tensor>>{I1728, h1_, Gamma5_()};
  auto task936 = make_shared<Task936>(tensor936, pindex);
  task935->add_dep(task936);
  task936->add_dep(task934);
  sourceq->add_task(task936);

  auto tensor937 = vector<shared_ptr<Tensor>>{I1728, v2_, Gamma104_()};
  auto task937 = make_shared<Task937>(tensor937, pindex);
  task935->add_dep(task937);
  task937->add_dep(task934);
  sourceq->add_task(task937);

  auto tensor938 = vector<shared_ptr<Tensor>>{I1728, v2_, Gamma4_()};
  auto task938 = make_shared<Task938>(tensor938, pindex);
  task935->add_dep(task938);
  task938->add_dep(task934);
  sourceq->add_task(task938);

  vector<IndexRange> I1730_index = {active_, active_, closed_, virt_};
  auto I1730 = make_shared<Tensor>(I1730_index);
  auto tensor939 = vector<shared_ptr<Tensor>>{s, I1730};
  auto task939 = make_shared<Task939>(tensor939, pindex);
  task939->add_dep(task934);
  sourceq->add_task(task939);

  auto tensor940 = vector<shared_ptr<Tensor>>{I1730, h1_, Gamma32_()};
  auto task940 = make_shared<Task940>(tensor940, pindex);
  task939->add_dep(task940);
  task940->add_dep(task934);
  sourceq->add_task(task940);

  auto tensor941 = vector<shared_ptr<Tensor>>{I1730, v2_, Gamma29_()};
  auto task941 = make_shared<Task941>(tensor941, pindex);
  task939->add_dep(task941);
  task941->add_dep(task934);
  sourceq->add_task(task941);

  auto tensor942 = vector<shared_ptr<Tensor>>{I1730, v2_, Gamma25_()};
  auto task942 = make_shared<Task942>(tensor942, pindex);
  task939->add_dep(task942);
  task942->add_dep(task934);
  sourceq->add_task(task942);

  auto tensor943 = vector<shared_ptr<Tensor>>{I1730, v2_, Gamma27_()};
  auto task943 = make_shared<Task943>(tensor943, pindex);
  task939->add_dep(task943);
  task943->add_dep(task934);
  sourceq->add_task(task943);

  auto tensor944 = vector<shared_ptr<Tensor>>{I1730, v2_, Gamma29_()};
  auto task944 = make_shared<Task944>(tensor944, pindex);
  task939->add_dep(task944);
  task944->add_dep(task934);
  sourceq->add_task(task944);

  vector<IndexRange> I1732_index = {active_, active_, closed_, virt_};
  auto I1732 = make_shared<Tensor>(I1732_index);
  auto tensor945 = vector<shared_ptr<Tensor>>{s, I1732};
  auto task945 = make_shared<Task945>(tensor945, pindex);
  task945->add_dep(task934);
  sourceq->add_task(task945);

  auto tensor946 = vector<shared_ptr<Tensor>>{I1732, h1_, Gamma32_()};
  auto task946 = make_shared<Task946>(tensor946, pindex);
  task945->add_dep(task946);
  task946->add_dep(task934);
  sourceq->add_task(task946);

  auto tensor947 = vector<shared_ptr<Tensor>>{I1732, v2_, Gamma29_()};
  auto task947 = make_shared<Task947>(tensor947, pindex);
  task945->add_dep(task947);
  task947->add_dep(task934);
  sourceq->add_task(task947);

  auto tensor948 = vector<shared_ptr<Tensor>>{I1732, v2_, Gamma5_()};
  auto task948 = make_shared<Task948>(tensor948, pindex);
  task945->add_dep(task948);
  task948->add_dep(task934);
  sourceq->add_task(task948);

  auto tensor949 = vector<shared_ptr<Tensor>>{I1732, Gamma29_(), v2_};
  auto task949 = make_shared<Task949>(tensor949, pindex);
  task945->add_dep(task949);
  task949->add_dep(task934);
  sourceq->add_task(task949);

  auto tensor950 = vector<shared_ptr<Tensor>>{I1732, v2_, Gamma29_()};
  auto task950 = make_shared<Task950>(tensor950, pindex);
  task945->add_dep(task950);
  task950->add_dep(task934);
  sourceq->add_task(task950);

  vector<IndexRange> I1734_index = {active_, active_, active_, virt_};
  auto I1734 = make_shared<Tensor>(I1734_index);
  auto tensor951 = vector<shared_ptr<Tensor>>{s, I1734};
  auto task951 = make_shared<Task951>(tensor951, pindex);
  task951->add_dep(task934);
  sourceq->add_task(task951);

  auto tensor952 = vector<shared_ptr<Tensor>>{I1734, h1_, Gamma51_()};
  auto task952 = make_shared<Task952>(tensor952, pindex);
  task951->add_dep(task952);
  task952->add_dep(task934);
  sourceq->add_task(task952);

  auto tensor953 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma50_()};
  auto task953 = make_shared<Task953>(tensor953, pindex);
  task951->add_dep(task953);
  task953->add_dep(task934);
  sourceq->add_task(task953);

  auto tensor954 = vector<shared_ptr<Tensor>>{I1734, v2_, Gamma49_()};
  auto task954 = make_shared<Task954>(tensor954, pindex);
  task951->add_dep(task954);
  task954->add_dep(task934);
  sourceq->add_task(task954);

  vector<IndexRange> I1736_index = {active_, active_, closed_, closed_};
  auto I1736 = make_shared<Tensor>(I1736_index);
  auto tensor955 = vector<shared_ptr<Tensor>>{s, I1736};
  auto task955 = make_shared<Task955>(tensor955, pindex);
  task955->add_dep(task934);
  sourceq->add_task(task955);

  auto tensor956 = vector<shared_ptr<Tensor>>{I1736, v2_, Gamma0_()};
  auto task956 = make_shared<Task956>(tensor956, pindex);
  task955->add_dep(task956);
  task956->add_dep(task934);
  sourceq->add_task(task956);

  vector<IndexRange> I1742_index = {active_, closed_, closed_, virt_};
  auto I1742 = make_shared<Tensor>(I1742_index);
  auto tensor957 = vector<shared_ptr<Tensor>>{s, I1742};
  auto task957 = make_shared<Task957>(tensor957, pindex);
  task957->add_dep(task934);
  sourceq->add_task(task957);

  auto tensor958 = vector<shared_ptr<Tensor>>{I1742, v2_, Gamma12_()};
  auto task958 = make_shared<Task958>(tensor958, pindex);
  task957->add_dep(task958);
  task958->add_dep(task934);
  sourceq->add_task(task958);

  auto tensor959 = vector<shared_ptr<Tensor>>{I1742, v2_, Gamma12_()};
  auto task959 = make_shared<Task959>(tensor959, pindex);
  task957->add_dep(task959);
  task959->add_dep(task934);
  sourceq->add_task(task959);

  shared_ptr<Tensor> I1766;
  if (diagonal) {
    vector<IndexRange> I1766_index = {closed_, virt_, closed_, virt_};
    I1766 = make_shared<Tensor>(I1766_index);
  }
  shared_ptr<Task960> task960;
  if (diagonal) {
    auto tensor960 = vector<shared_ptr<Tensor>>{s, I1766};
    task960 = make_shared<Task960>(tensor960, pindex);
    task960->add_dep(task934);
    sourceq->add_task(task960);
  }

  shared_ptr<Task961> task961;
  if (diagonal) {
    auto tensor961 = vector<shared_ptr<Tensor>>{I1766, v2_};
    task961 = make_shared<Task961>(tensor961, pindex);
    task960->add_dep(task961);
    task961->add_dep(task934);
    sourceq->add_task(task961);
  }

  vector<IndexRange> I1768_index = {virt_, closed_, virt_, active_};
  auto I1768 = make_shared<Tensor>(I1768_index);
  auto tensor962 = vector<shared_ptr<Tensor>>{s, I1768};
  auto task962 = make_shared<Task962>(tensor962, pindex);
  task962->add_dep(task934);
  sourceq->add_task(task962);

  auto tensor963 = vector<shared_ptr<Tensor>>{I1768, Gamma32_(), v2_};
  auto task963 = make_shared<Task963>(tensor963, pindex);
  task962->add_dep(task963);
  task963->add_dep(task934);
  sourceq->add_task(task963);

  auto tensor964 = vector<shared_ptr<Tensor>>{I1768, v2_, Gamma32_()};
  auto task964 = make_shared<Task964>(tensor964, pindex);
  task962->add_dep(task964);
  task964->add_dep(task934);
  sourceq->add_task(task964);

  vector<IndexRange> I1772_index = {active_, active_, virt_, virt_};
  auto I1772 = make_shared<Tensor>(I1772_index);
  auto tensor965 = vector<shared_ptr<Tensor>>{s, I1772};
  auto task965 = make_shared<Task965>(tensor965, pindex);
  task965->add_dep(task934);
  sourceq->add_task(task965);

  auto tensor966 = vector<shared_ptr<Tensor>>{I1772, v2_, Gamma51_()};
  auto task966 = make_shared<Task966>(tensor966, pindex);
  task965->add_dep(task966);
  task966->add_dep(task934);
  sourceq->add_task(task966);

  return sourceq;
}


#endif
