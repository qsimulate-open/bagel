//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_sourceqq.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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
#include <src/smith/mrci/MRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_sourceq(const bool reset, const bool diagonal) {

  auto sourceq = make_shared<Queue>();
  auto task930 = make_shared<Task930>(s, reset);
  sourceq->add_task(task930);

  auto I1719 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task931 = make_shared<Task931>(s, I1719);
  task931->add_dep(task930);
  sourceq->add_task(task931);

  auto task932 = make_shared<Task932>(I1719, Gamma5_(), h1_);
  task931->add_dep(task932);
  task932->add_dep(task930);
  sourceq->add_task(task932);

  auto task933 = make_shared<Task933>(I1719, Gamma104_(), v2_);
  task931->add_dep(task933);
  task933->add_dep(task930);
  sourceq->add_task(task933);

  auto task934 = make_shared<Task934>(I1719, Gamma4_(), v2_);
  task931->add_dep(task934);
  task934->add_dep(task930);
  sourceq->add_task(task934);

  auto I1721 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task935 = make_shared<Task935>(s, I1721);
  task935->add_dep(task930);
  sourceq->add_task(task935);

  auto task936 = make_shared<Task936>(I1721, Gamma32_(), h1_);
  task935->add_dep(task936);
  task936->add_dep(task930);
  sourceq->add_task(task936);

  auto I1738 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task937 = make_shared<Task937>(I1721, Gamma29_(), I1738);
  task935->add_dep(task937);
  task937->add_dep(task930);
  sourceq->add_task(task937);

  auto task938 = make_shared<Task938>(I1738, v2_);
  task937->add_dep(task938);
  task938->add_dep(task930);
  sourceq->add_task(task938);

  auto task939 = make_shared<Task939>(I1721, Gamma25_(), v2_);
  task935->add_dep(task939);
  task939->add_dep(task930);
  sourceq->add_task(task939);

  auto task940 = make_shared<Task940>(I1721, Gamma27_(), v2_);
  task935->add_dep(task940);
  task940->add_dep(task930);
  sourceq->add_task(task940);

  auto I1723 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task941 = make_shared<Task941>(s, I1723);
  task941->add_dep(task930);
  sourceq->add_task(task941);

  auto task942 = make_shared<Task942>(I1723, Gamma32_(), h1_);
  task941->add_dep(task942);
  task942->add_dep(task930);
  sourceq->add_task(task942);

  auto I1746 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task943 = make_shared<Task943>(I1723, Gamma29_(), I1746);
  task941->add_dep(task943);
  task943->add_dep(task930);
  sourceq->add_task(task943);

  auto task944 = make_shared<Task944>(I1746, v2_);
  task943->add_dep(task944);
  task944->add_dep(task930);
  sourceq->add_task(task944);

  auto task945 = make_shared<Task945>(I1723, Gamma5_(), v2_);
  task941->add_dep(task945);
  task945->add_dep(task930);
  sourceq->add_task(task945);

  auto I1725 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task946 = make_shared<Task946>(s, I1725);
  task946->add_dep(task930);
  sourceq->add_task(task946);

  auto task947 = make_shared<Task947>(I1725, Gamma51_(), h1_);
  task946->add_dep(task947);
  task947->add_dep(task930);
  sourceq->add_task(task947);

  auto task948 = make_shared<Task948>(I1725, Gamma50_(), v2_);
  task946->add_dep(task948);
  task948->add_dep(task930);
  sourceq->add_task(task948);

  auto task949 = make_shared<Task949>(I1725, Gamma49_(), v2_);
  task946->add_dep(task949);
  task949->add_dep(task930);
  sourceq->add_task(task949);

  auto I1727 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task950 = make_shared<Task950>(s, I1727);
  task950->add_dep(task930);
  sourceq->add_task(task950);

  auto task951 = make_shared<Task951>(I1727, Gamma0_(), v2_);
  task950->add_dep(task951);
  task951->add_dep(task930);
  sourceq->add_task(task951);

  auto I1733 = make_shared<TATensor<double,4>>({closed_, closed_, virt_, active_});
  auto task952 = make_shared<Task952>(s, I1733);
  task952->add_dep(task930);
  sourceq->add_task(task952);

  auto I1734 = make_shared<TATensor<double,4>>({closed_, active_, closed_, virt_});
  auto task953 = make_shared<Task953>(I1733, Gamma12_(), I1734);
  task952->add_dep(task953);
  task953->add_dep(task930);
  sourceq->add_task(task953);

  auto task954 = make_shared<Task954>(I1734, v2_);
  task953->add_dep(task954);
  task954->add_dep(task930);
  sourceq->add_task(task954);

  shared_ptr<TATensor<double,4>> I1757;
  if (diagonal) {
    I1757 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task955> task955;
  if (diagonal) {
    task955 = make_shared<Task955>(s, I1757);
    task955->add_dep(task930);
    sourceq->add_task(task955);
  }

  shared_ptr<Task956> task956;
  if (diagonal) {
    task956 = make_shared<Task956>(I1757, v2_);
    task955->add_dep(task956);
    task956->add_dep(task930);
    sourceq->add_task(task956);
  }

  auto I1759 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task957 = make_shared<Task957>(s, I1759);
  task957->add_dep(task930);
  sourceq->add_task(task957);

  auto I1760 = make_shared<TATensor<double,4>>({active_, virt_, closed_, virt_});
  auto task958 = make_shared<Task958>(I1759, Gamma32_(), I1760);
  task957->add_dep(task958);
  task958->add_dep(task930);
  sourceq->add_task(task958);

  auto task959 = make_shared<Task959>(I1760, v2_);
  task958->add_dep(task959);
  task959->add_dep(task930);
  sourceq->add_task(task959);

  auto I1763 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task960 = make_shared<Task960>(s, I1763);
  task960->add_dep(task930);
  sourceq->add_task(task960);

  auto task961 = make_shared<Task961>(I1763, Gamma51_(), v2_);
  task960->add_dep(task961);
  task961->add_dep(task930);
  sourceq->add_task(task961);

  return sourceq;
}


#endif
