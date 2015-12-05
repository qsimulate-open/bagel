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
  auto task755 = make_shared<Task755>(n, reset);
  normq->add_task(task755);

  auto I1348 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task756 = make_shared<Task756>(n, I1348);
  task756->add_dep(task755);
  normq->add_task(task756);

  auto task757 = make_shared<Task757>(I1348, Gamma0_(), t2);
  task756->add_dep(task757);
  task757->add_dep(task755);
  normq->add_task(task757);

  auto I1350 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task758 = make_shared<Task758>(n, I1350);
  task758->add_dep(task755);
  normq->add_task(task758);

  auto task759 = make_shared<Task759>(I1350, Gamma4_(), t2);
  task758->add_dep(task759);
  task759->add_dep(task755);
  normq->add_task(task759);

  auto I1352 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task760 = make_shared<Task760>(n, I1352);
  task760->add_dep(task755);
  normq->add_task(task760);

  auto I1353 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task761 = make_shared<Task761>(I1352, Gamma11_(), I1353);
  task760->add_dep(task761);
  task761->add_dep(task755);
  normq->add_task(task761);

  auto task762 = make_shared<Task762>(I1353, t2);
  task761->add_dep(task762);
  task762->add_dep(task755);
  normq->add_task(task762);

  auto I1356 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task763 = make_shared<Task763>(n, I1356);
  task763->add_dep(task755);
  normq->add_task(task763);

  auto task764 = make_shared<Task764>(I1356, Gamma24_(), t2);
  task763->add_dep(task764);
  task764->add_dep(task755);
  normq->add_task(task764);

  auto I1358 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task765 = make_shared<Task765>(n, I1358);
  task765->add_dep(task755);
  normq->add_task(task765);

  auto task766 = make_shared<Task766>(I1358, Gamma32_(), t2);
  task765->add_dep(task766);
  task766->add_dep(task755);
  normq->add_task(task766);

  shared_ptr<TATensor<std::complex<double>,4>> I1360;
  if (diagonal) {
    I1360 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task767> task767;
  if (diagonal) {
    task767 = make_shared<Task767>(n, I1360);
    task767->add_dep(task755);
    normq->add_task(task767);
  }

  shared_ptr<Task768> task768;
  if (diagonal) {
    task768 = make_shared<Task768>(I1360, t2);
    task767->add_dep(task768);
    task768->add_dep(task755);
    normq->add_task(task768);
  }

  auto I1362 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task769 = make_shared<Task769>(n, I1362);
  task769->add_dep(task755);
  normq->add_task(task769);

  auto I1363 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task770 = make_shared<Task770>(I1362, Gamma27_(), I1363);
  task769->add_dep(task770);
  task770->add_dep(task755);
  normq->add_task(task770);

  auto task771 = make_shared<Task771>(I1363, t2);
  task770->add_dep(task771);
  task771->add_dep(task755);
  normq->add_task(task771);

  auto I1366 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task772 = make_shared<Task772>(n, I1366);
  task772->add_dep(task755);
  normq->add_task(task772);

  auto task773 = make_shared<Task773>(I1366, Gamma33_(), t2);
  task772->add_dep(task773);
  task773->add_dep(task755);
  normq->add_task(task773);

  return normq;
}


#endif
