//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_normqq.cc
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

shared_ptr<Queue> RelMRCI::RelMRCI::make_normq(const bool reset, const bool diagonal) {

  auto normq = make_shared<Queue>();
  auto task754 = make_shared<Task754>(n, reset);
  normq->add_task(task754);

  auto I1348 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task755 = make_shared<Task755>(n, I1348);
  task755->add_dep(task754);
  normq->add_task(task755);

  auto task756 = make_shared<Task756>(I1348, Gamma0_(), t2);
  task755->add_dep(task756);
  task756->add_dep(task754);
  normq->add_task(task756);

  auto I1350 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task757 = make_shared<Task757>(n, I1350);
  task757->add_dep(task754);
  normq->add_task(task757);

  auto task758 = make_shared<Task758>(I1350, Gamma4_(), t2);
  task757->add_dep(task758);
  task758->add_dep(task754);
  normq->add_task(task758);

  auto I1352 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task759 = make_shared<Task759>(n, I1352);
  task759->add_dep(task754);
  normq->add_task(task759);

  auto I1353 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task760 = make_shared<Task760>(I1352, Gamma11_(), I1353);
  task759->add_dep(task760);
  task760->add_dep(task754);
  normq->add_task(task760);

  auto task761 = make_shared<Task761>(I1353, t2);
  task760->add_dep(task761);
  task761->add_dep(task754);
  normq->add_task(task761);

  auto I1356 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task762 = make_shared<Task762>(n, I1356);
  task762->add_dep(task754);
  normq->add_task(task762);

  auto task763 = make_shared<Task763>(I1356, Gamma24_(), t2);
  task762->add_dep(task763);
  task763->add_dep(task754);
  normq->add_task(task763);

  auto I1358 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task764 = make_shared<Task764>(n, I1358);
  task764->add_dep(task754);
  normq->add_task(task764);

  auto task765 = make_shared<Task765>(I1358, Gamma32_(), t2);
  task764->add_dep(task765);
  task765->add_dep(task754);
  normq->add_task(task765);

  shared_ptr<TATensor<std::complex<double>,4>> I1360;
  if (diagonal) {
    I1360 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task766> task766;
  if (diagonal) {
    task766 = make_shared<Task766>(n, I1360);
    task766->add_dep(task754);
    normq->add_task(task766);
  }

  shared_ptr<Task767> task767;
  if (diagonal) {
    task767 = make_shared<Task767>(I1360, t2);
    task766->add_dep(task767);
    task767->add_dep(task754);
    normq->add_task(task767);
  }

  auto I1362 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task768 = make_shared<Task768>(n, I1362);
  task768->add_dep(task754);
  normq->add_task(task768);

  auto I1363 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task769 = make_shared<Task769>(I1362, Gamma27_(), I1363);
  task768->add_dep(task769);
  task769->add_dep(task754);
  normq->add_task(task769);

  auto task770 = make_shared<Task770>(I1363, t2);
  task769->add_dep(task770);
  task770->add_dep(task754);
  normq->add_task(task770);

  auto I1366 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task771 = make_shared<Task771>(n, I1366);
  task771->add_dep(task754);
  normq->add_task(task771);

  auto task772 = make_shared<Task772>(I1366, Gamma33_(), t2);
  task771->add_dep(task772);
  task772->add_dep(task754);
  normq->add_task(task772);

  return normq;
}


#endif
