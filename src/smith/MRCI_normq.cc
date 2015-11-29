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


#include <src/smith/MRCI.h>
#include <src/smith/MRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_normq(const bool reset, const bool diagonal) {

  auto normq = make_shared<Queue>();
  auto task969 = make_shared<Task969>(n, reset);
  normq->add_task(task969);

  auto I1778 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task970 = make_shared<Task970>(n, I1778);
  task970->add_dep(task969);
  normq->add_task(task970);

  auto task971 = make_shared<Task971>(I1778, Gamma0_(), t2);
  task970->add_dep(task971);
  task971->add_dep(task969);
  normq->add_task(task971);

  auto I1780 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task972 = make_shared<Task972>(n, I1780);
  task972->add_dep(task969);
  normq->add_task(task972);

  auto task973 = make_shared<Task973>(I1780, Gamma4_(), t2);
  task972->add_dep(task973);
  task973->add_dep(task969);
  normq->add_task(task973);

  auto I1782 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task974 = make_shared<Task974>(n, I1782);
  task974->add_dep(task969);
  normq->add_task(task974);

  auto I1783 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task975 = make_shared<Task975>(I1782, Gamma12_(), I1783);
  task974->add_dep(task975);
  task975->add_dep(task969);
  normq->add_task(task975);

  auto task976 = make_shared<Task976>(I1783, t2);
  task975->add_dep(task976);
  task976->add_dep(task969);
  normq->add_task(task976);

  auto I1786 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task977 = make_shared<Task977>(n, I1786);
  task977->add_dep(task969);
  normq->add_task(task977);

  auto task978 = make_shared<Task978>(I1786, Gamma27_(), t2);
  task977->add_dep(task978);
  task978->add_dep(task969);
  normq->add_task(task978);

  auto task979 = make_shared<Task979>(I1786, Gamma29_(), t2);
  task977->add_dep(task979);
  task979->add_dep(task969);
  normq->add_task(task979);

  auto I1790 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task980 = make_shared<Task980>(n, I1790);
  task980->add_dep(task969);
  normq->add_task(task980);

  auto I1791 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task981 = make_shared<Task981>(I1790, Gamma29_(), I1791);
  task980->add_dep(task981);
  task981->add_dep(task969);
  normq->add_task(task981);

  auto task982 = make_shared<Task982>(I1791, t2);
  task981->add_dep(task982);
  task982->add_dep(task969);
  normq->add_task(task982);

  auto I1794 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task983 = make_shared<Task983>(n, I1794);
  task983->add_dep(task969);
  normq->add_task(task983);

  auto task984 = make_shared<Task984>(I1794, Gamma50_(), t2);
  task983->add_dep(task984);
  task984->add_dep(task969);
  normq->add_task(task984);

  shared_ptr<TATensor<double,4>> I1796;
  if (diagonal) {
    I1796 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task985> task985;
  if (diagonal) {
    task985 = make_shared<Task985>(n, I1796);
    task985->add_dep(task969);
    normq->add_task(task985);
  }

  shared_ptr<Task986> task986;
  if (diagonal) {
    task986 = make_shared<Task986>(I1796, t2);
    task985->add_dep(task986);
    task986->add_dep(task969);
    normq->add_task(task986);
  }

  auto I1798 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task987 = make_shared<Task987>(n, I1798);
  task987->add_dep(task969);
  normq->add_task(task987);

  auto I1799 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task988 = make_shared<Task988>(I1798, Gamma32_(), I1799);
  task987->add_dep(task988);
  task988->add_dep(task969);
  normq->add_task(task988);

  auto task989 = make_shared<Task989>(I1799, t2);
  task988->add_dep(task989);
  task989->add_dep(task969);
  normq->add_task(task989);

  auto I1802 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task990 = make_shared<Task990>(n, I1802);
  task990->add_dep(task969);
  normq->add_task(task990);

  auto task991 = make_shared<Task991>(I1802, Gamma51_(), t2);
  task990->add_dep(task991);
  task991->add_dep(task969);
  normq->add_task(task991);

  return normq;
}


#endif
