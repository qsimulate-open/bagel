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

  auto sourceq = make_shared<Queue>();
  auto task937 = make_shared<Task937>(s, reset);
  sourceq->add_task(task937);

  auto I1732 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task938 = make_shared<Task938>(s, I1732);
  task938->add_dep(task937);
  sourceq->add_task(task938);

  auto task939 = make_shared<Task939>(I1732, Gamma5_(), h1_);
  task938->add_dep(task939);
  task939->add_dep(task937);
  sourceq->add_task(task939);

  auto task940 = make_shared<Task940>(I1732, Gamma104_(), v2_);
  task938->add_dep(task940);
  task940->add_dep(task937);
  sourceq->add_task(task940);

  auto task941 = make_shared<Task941>(I1732, Gamma4_(), v2_);
  task938->add_dep(task941);
  task941->add_dep(task937);
  sourceq->add_task(task941);

  auto I1734 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task942 = make_shared<Task942>(s, I1734);
  task942->add_dep(task937);
  sourceq->add_task(task942);

  auto task943 = make_shared<Task943>(I1734, Gamma32_(), h1_);
  task942->add_dep(task943);
  task943->add_dep(task937);
  sourceq->add_task(task943);

  auto I1751 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task944 = make_shared<Task944>(I1734, Gamma29_(), I1751);
  task942->add_dep(task944);
  task944->add_dep(task937);
  sourceq->add_task(task944);

  auto task945 = make_shared<Task945>(I1751, v2_);
  task944->add_dep(task945);
  task945->add_dep(task937);
  sourceq->add_task(task945);

  auto task946 = make_shared<Task946>(I1734, Gamma25_(), v2_);
  task942->add_dep(task946);
  task946->add_dep(task937);
  sourceq->add_task(task946);

  auto task947 = make_shared<Task947>(I1734, Gamma27_(), v2_);
  task942->add_dep(task947);
  task947->add_dep(task937);
  sourceq->add_task(task947);

  auto I1736 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task948 = make_shared<Task948>(s, I1736);
  task948->add_dep(task937);
  sourceq->add_task(task948);

  auto task949 = make_shared<Task949>(I1736, Gamma32_(), h1_);
  task948->add_dep(task949);
  task949->add_dep(task937);
  sourceq->add_task(task949);

  auto I1759 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task950 = make_shared<Task950>(I1736, Gamma29_(), I1759);
  task948->add_dep(task950);
  task950->add_dep(task937);
  sourceq->add_task(task950);

  auto task951 = make_shared<Task951>(I1759, v2_);
  task950->add_dep(task951);
  task951->add_dep(task937);
  sourceq->add_task(task951);

  auto task952 = make_shared<Task952>(I1736, Gamma5_(), v2_);
  task948->add_dep(task952);
  task952->add_dep(task937);
  sourceq->add_task(task952);

  auto I1738 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task953 = make_shared<Task953>(s, I1738);
  task953->add_dep(task937);
  sourceq->add_task(task953);

  auto task954 = make_shared<Task954>(I1738, Gamma51_(), h1_);
  task953->add_dep(task954);
  task954->add_dep(task937);
  sourceq->add_task(task954);

  auto task955 = make_shared<Task955>(I1738, Gamma50_(), v2_);
  task953->add_dep(task955);
  task955->add_dep(task937);
  sourceq->add_task(task955);

  auto task956 = make_shared<Task956>(I1738, Gamma49_(), v2_);
  task953->add_dep(task956);
  task956->add_dep(task937);
  sourceq->add_task(task956);

  auto I1740 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task957 = make_shared<Task957>(s, I1740);
  task957->add_dep(task937);
  sourceq->add_task(task957);

  auto task958 = make_shared<Task958>(I1740, Gamma0_(), v2_);
  task957->add_dep(task958);
  task958->add_dep(task937);
  sourceq->add_task(task958);

  auto I1746 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, closed_, virt_, active_});
  auto task959 = make_shared<Task959>(s, I1746);
  task959->add_dep(task937);
  sourceq->add_task(task959);

  auto I1747 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, closed_, virt_});
  auto task960 = make_shared<Task960>(I1746, Gamma12_(), I1747);
  task959->add_dep(task960);
  task960->add_dep(task937);
  sourceq->add_task(task960);

  auto task961 = make_shared<Task961>(I1747, v2_);
  task960->add_dep(task961);
  task961->add_dep(task937);
  sourceq->add_task(task961);

  shared_ptr<TATensor<double,4>> I1770;
  if (diagonal) {
    I1770 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task962> task962;
  if (diagonal) {
    task962 = make_shared<Task962>(s, I1770);
    task962->add_dep(task937);
    sourceq->add_task(task962);
  }

  shared_ptr<Task963> task963;
  if (diagonal) {
    task963 = make_shared<Task963>(I1770, v2_);
    task962->add_dep(task963);
    task963->add_dep(task937);
    sourceq->add_task(task963);
  }

  auto I1772 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task964 = make_shared<Task964>(s, I1772);
  task964->add_dep(task937);
  sourceq->add_task(task964);

  auto I1773 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task965 = make_shared<Task965>(I1772, Gamma32_(), I1773);
  task964->add_dep(task965);
  task965->add_dep(task937);
  sourceq->add_task(task965);

  auto task966 = make_shared<Task966>(I1773, v2_);
  task965->add_dep(task966);
  task966->add_dep(task937);
  sourceq->add_task(task966);

  auto I1776 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task967 = make_shared<Task967>(s, I1776);
  task967->add_dep(task937);
  sourceq->add_task(task967);

  auto task968 = make_shared<Task968>(I1776, Gamma51_(), v2_);
  task967->add_dep(task968);
  task968->add_dep(task937);
  sourceq->add_task(task968);

  return sourceq;
}


#endif
