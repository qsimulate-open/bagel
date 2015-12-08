//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_normqq.cc
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


#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_normq(const bool reset, const bool diagonal) {

  auto normq = make_shared<Queue>();
  auto task962 = make_shared<Task962>(n, reset);
  normq->add_task(task962);

  auto I1765 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task963 = make_shared<Task963>(n, I1765);
  task963->add_dep(task962);
  normq->add_task(task963);

  auto task964 = make_shared<Task964>(I1765, Gamma0_(), t2);
  task963->add_dep(task964);
  task964->add_dep(task962);
  normq->add_task(task964);

  auto I1767 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task965 = make_shared<Task965>(n, I1767);
  task965->add_dep(task962);
  normq->add_task(task965);

  auto task966 = make_shared<Task966>(I1767, Gamma4_(), t2);
  task965->add_dep(task966);
  task966->add_dep(task962);
  normq->add_task(task966);

  auto I1769 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task967 = make_shared<Task967>(n, I1769);
  task967->add_dep(task962);
  normq->add_task(task967);

  auto I1770 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task968 = make_shared<Task968>(I1769, Gamma12_(), I1770);
  task967->add_dep(task968);
  task968->add_dep(task962);
  normq->add_task(task968);

  auto task969 = make_shared<Task969>(I1770, t2);
  task968->add_dep(task969);
  task969->add_dep(task962);
  normq->add_task(task969);

  auto I1773 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task970 = make_shared<Task970>(n, I1773);
  task970->add_dep(task962);
  normq->add_task(task970);

  auto task971 = make_shared<Task971>(I1773, Gamma27_(), t2);
  task970->add_dep(task971);
  task971->add_dep(task962);
  normq->add_task(task971);

  auto task972 = make_shared<Task972>(I1773, Gamma29_(), t2);
  task970->add_dep(task972);
  task972->add_dep(task962);
  normq->add_task(task972);

  auto I1777 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task973 = make_shared<Task973>(n, I1777);
  task973->add_dep(task962);
  normq->add_task(task973);

  auto I1778 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task974 = make_shared<Task974>(I1777, Gamma29_(), I1778);
  task973->add_dep(task974);
  task974->add_dep(task962);
  normq->add_task(task974);

  auto task975 = make_shared<Task975>(I1778, t2);
  task974->add_dep(task975);
  task975->add_dep(task962);
  normq->add_task(task975);

  auto I1781 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task976 = make_shared<Task976>(n, I1781);
  task976->add_dep(task962);
  normq->add_task(task976);

  auto task977 = make_shared<Task977>(I1781, Gamma50_(), t2);
  task976->add_dep(task977);
  task977->add_dep(task962);
  normq->add_task(task977);

  shared_ptr<TATensor<double,4>> I1783;
  if (diagonal) {
    I1783 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task978> task978;
  if (diagonal) {
    task978 = make_shared<Task978>(n, I1783);
    task978->add_dep(task962);
    normq->add_task(task978);
  }

  shared_ptr<Task979> task979;
  if (diagonal) {
    task979 = make_shared<Task979>(I1783, t2);
    task978->add_dep(task979);
    task979->add_dep(task962);
    normq->add_task(task979);
  }

  auto I1785 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task980 = make_shared<Task980>(n, I1785);
  task980->add_dep(task962);
  normq->add_task(task980);

  auto I1786 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task981 = make_shared<Task981>(I1785, Gamma32_(), I1786);
  task980->add_dep(task981);
  task981->add_dep(task962);
  normq->add_task(task981);

  auto task982 = make_shared<Task982>(I1786, t2);
  task981->add_dep(task982);
  task982->add_dep(task962);
  normq->add_task(task982);

  auto I1789 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task983 = make_shared<Task983>(n, I1789);
  task983->add_dep(task962);
  normq->add_task(task983);

  auto task984 = make_shared<Task984>(I1789, Gamma51_(), t2);
  task983->add_dep(task984);
  task984->add_dep(task962);
  normq->add_task(task984);

  return normq;
}


#endif
