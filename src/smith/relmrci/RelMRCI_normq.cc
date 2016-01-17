//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_normqq.cc
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


#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relmrci/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_normq(const bool reset, const bool diagonal) {

  auto normq = make_shared<Queue>();
  auto task747 = make_shared<Task747>(n, reset);
  normq->add_task(task747);

  auto I1335 = make_shared<TATensor<std::complex<double>,4>>({closed_, closed_, active_, active_});
  auto task748 = make_shared<Task748>(n, I1335);
  task748->add_dep(task747);
  normq->add_task(task748);

  auto task749 = make_shared<Task749>(I1335, Gamma0_(), t2);
  task748->add_dep(task749);
  task749->add_dep(task747);
  normq->add_task(task749);

  auto I1337 = make_shared<TATensor<std::complex<double>,4>>({closed_, active_, active_, active_});
  auto task750 = make_shared<Task750>(n, I1337);
  task750->add_dep(task747);
  normq->add_task(task750);

  auto task751 = make_shared<Task751>(I1337, Gamma4_(), t2);
  task750->add_dep(task751);
  task751->add_dep(task747);
  normq->add_task(task751);

  auto I1339 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, closed_, active_});
  auto task752 = make_shared<Task752>(n, I1339);
  task752->add_dep(task747);
  normq->add_task(task752);

  auto I1340 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, closed_, active_});
  auto task753 = make_shared<Task753>(I1339, Gamma11_(), I1340);
  task752->add_dep(task753);
  task753->add_dep(task747);
  normq->add_task(task753);

  auto task754 = make_shared<Task754>(I1340, t2);
  task753->add_dep(task754);
  task754->add_dep(task747);
  normq->add_task(task754);

  auto I1343 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, active_, active_});
  auto task755 = make_shared<Task755>(n, I1343);
  task755->add_dep(task747);
  normq->add_task(task755);

  auto task756 = make_shared<Task756>(I1343, Gamma24_(), t2);
  task755->add_dep(task756);
  task756->add_dep(task747);
  normq->add_task(task756);

  auto I1345 = make_shared<TATensor<std::complex<double>,4>>({virt_, active_, active_, active_});
  auto task757 = make_shared<Task757>(n, I1345);
  task757->add_dep(task747);
  normq->add_task(task757);

  auto task758 = make_shared<Task758>(I1345, Gamma32_(), t2);
  task757->add_dep(task758);
  task758->add_dep(task747);
  normq->add_task(task758);

  shared_ptr<TATensor<std::complex<double>,4>> I1347;
  if (diagonal) {
    I1347 = make_shared<TATensor<std::complex<double>,4>>({closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task759> task759;
  if (diagonal) {
    task759 = make_shared<Task759>(n, I1347);
    task759->add_dep(task747);
    normq->add_task(task759);
  }

  shared_ptr<Task760> task760;
  if (diagonal) {
    task760 = make_shared<Task760>(I1347, t2);
    task759->add_dep(task760);
    task760->add_dep(task747);
    normq->add_task(task760);
  }

  auto I1349 = make_shared<TATensor<std::complex<double>,4>>({virt_, closed_, virt_, active_});
  auto task761 = make_shared<Task761>(n, I1349);
  task761->add_dep(task747);
  normq->add_task(task761);

  auto I1350 = make_shared<TATensor<std::complex<double>,4>>({active_, virt_, closed_, virt_});
  auto task762 = make_shared<Task762>(I1349, Gamma27_(), I1350);
  task761->add_dep(task762);
  task762->add_dep(task747);
  normq->add_task(task762);

  auto task763 = make_shared<Task763>(I1350, t2);
  task762->add_dep(task763);
  task763->add_dep(task747);
  normq->add_task(task763);

  auto I1353 = make_shared<TATensor<std::complex<double>,4>>({virt_, virt_, active_, active_});
  auto task764 = make_shared<Task764>(n, I1353);
  task764->add_dep(task747);
  normq->add_task(task764);

  auto task765 = make_shared<Task765>(I1353, Gamma33_(), t2);
  task764->add_dep(task765);
  task765->add_dep(task747);
  normq->add_task(task765);

  return normq;
}


#endif
