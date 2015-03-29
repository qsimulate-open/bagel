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

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor1410 = {n};
  auto task1410 = make_shared<Task1410>(tensor1410, reset);
  normq->add_task(task1410);

  vector<IndexRange> I1778_index = {closed_, closed_, active_, active_};
  auto I1778 = make_shared<Tensor>(I1778_index);
  vector<shared_ptr<Tensor>> tensor1411 = {n, I1778};
  auto task1411 = make_shared<Task1411>(tensor1411, pindex);
  task1411->add_dep(task1410);
  normq->add_task(task1411);

  vector<IndexRange> I1779_index = {closed_, active_, closed_, active_};
  auto I1779 = make_shared<Tensor>(I1779_index);
  vector<shared_ptr<Tensor>> tensor1412 = {I1778, Gamma0_(), I1779};
  auto task1412 = make_shared<Task1412>(tensor1412, pindex);
  task1411->add_dep(task1412);
  task1412->add_dep(task1410);
  normq->add_task(task1412);

  vector<shared_ptr<Tensor>> tensor1413 = {I1779, t2};
  auto task1413 = make_shared<Task1413>(tensor1413, pindex);
  task1412->add_dep(task1413);
  task1413->add_dep(task1410);
  normq->add_task(task1413);

  vector<IndexRange> I1780_index = {active_, active_, active_, closed_};
  auto I1780 = make_shared<Tensor>(I1780_index);
  vector<shared_ptr<Tensor>> tensor1414 = {n, I1780};
  auto task1414 = make_shared<Task1414>(tensor1414, pindex);
  task1414->add_dep(task1410);
  normq->add_task(task1414);

  vector<IndexRange> I1781_index = {active_, active_, active_, active_, active_, active_};
  auto I1781 = make_shared<Tensor>(I1781_index);
  vector<shared_ptr<Tensor>> tensor1415 = {I1780, t2, I1781};
  auto task1415 = make_shared<Task1415>(tensor1415, pindex);
  task1414->add_dep(task1415);
  task1415->add_dep(task1410);
  normq->add_task(task1415);

  vector<shared_ptr<Tensor>> tensor1416 = {I1781, Gamma4_()};
  auto task1416 = make_shared<Task1416>(tensor1416, pindex);
  task1415->add_dep(task1416);
  task1416->add_dep(task1410);
  normq->add_task(task1416);

  vector<IndexRange> I1782_index = {active_, closed_, virt_, closed_};
  auto I1782 = make_shared<Tensor>(I1782_index);
  vector<shared_ptr<Tensor>> tensor1417 = {n, I1782};
  auto task1417 = make_shared<Task1417>(tensor1417, pindex);
  task1417->add_dep(task1410);
  normq->add_task(task1417);

  vector<IndexRange> I1783_index = {active_, active_};
  auto I1783 = make_shared<Tensor>(I1783_index);
  vector<shared_ptr<Tensor>> tensor1418 = {I1782, t2, I1783};
  auto task1418 = make_shared<Task1418>(tensor1418, pindex);
  task1417->add_dep(task1418);
  task1418->add_dep(task1410);
  normq->add_task(task1418);

  vector<shared_ptr<Tensor>> tensor1419 = {I1783, Gamma12_()};
  auto task1419 = make_shared<Task1419>(tensor1419, pindex);
  task1418->add_dep(task1419);
  task1419->add_dep(task1410);
  normq->add_task(task1419);

  vector<IndexRange> I1785_index = {active_, active_};
  auto I1785 = make_shared<Tensor>(I1785_index);
  vector<shared_ptr<Tensor>> tensor1420 = {I1782, t2, I1785};
  auto task1420 = make_shared<Task1420>(tensor1420, pindex);
  task1417->add_dep(task1420);
  task1420->add_dep(task1410);
  normq->add_task(task1420);

  vector<shared_ptr<Tensor>> tensor1421 = {I1785, Gamma12_()};
  auto task1421 = make_shared<Task1421>(tensor1421, pindex);
  task1420->add_dep(task1421);
  task1421->add_dep(task1410);
  normq->add_task(task1421);

  vector<IndexRange> I1786_index = {virt_, closed_, active_, active_};
  auto I1786 = make_shared<Tensor>(I1786_index);
  vector<shared_ptr<Tensor>> tensor1422 = {n, I1786};
  auto task1422 = make_shared<Task1422>(tensor1422, pindex);
  task1422->add_dep(task1410);
  normq->add_task(task1422);

  vector<IndexRange> I1787_index = {active_, virt_, closed_, active_};
  auto I1787 = make_shared<Tensor>(I1787_index);
  vector<shared_ptr<Tensor>> tensor1423 = {I1786, Gamma27_(), I1787};
  auto task1423 = make_shared<Task1423>(tensor1423, pindex);
  task1422->add_dep(task1423);
  task1423->add_dep(task1410);
  normq->add_task(task1423);

  vector<shared_ptr<Tensor>> tensor1424 = {I1787, t2};
  auto task1424 = make_shared<Task1424>(tensor1424, pindex);
  task1423->add_dep(task1424);
  task1424->add_dep(task1410);
  normq->add_task(task1424);

  vector<IndexRange> I1789_index = {active_, active_, active_, active_};
  auto I1789 = make_shared<Tensor>(I1789_index);
  vector<shared_ptr<Tensor>> tensor1425 = {I1786, t2, I1789};
  auto task1425 = make_shared<Task1425>(tensor1425, pindex);
  task1422->add_dep(task1425);
  task1425->add_dep(task1410);
  normq->add_task(task1425);

  vector<shared_ptr<Tensor>> tensor1426 = {I1789, Gamma29_()};
  auto task1426 = make_shared<Task1426>(tensor1426, pindex);
  task1425->add_dep(task1426);
  task1426->add_dep(task1410);
  normq->add_task(task1426);

  vector<IndexRange> I1790_index = {virt_, closed_, active_, active_};
  auto I1790 = make_shared<Tensor>(I1790_index);
  vector<shared_ptr<Tensor>> tensor1427 = {n, I1790};
  auto task1427 = make_shared<Task1427>(tensor1427, pindex);
  task1427->add_dep(task1410);
  normq->add_task(task1427);

  vector<IndexRange> I1791_index = {active_, virt_, closed_, active_};
  auto I1791 = make_shared<Tensor>(I1791_index);
  vector<shared_ptr<Tensor>> tensor1428 = {I1790, Gamma29_(), I1791};
  auto task1428 = make_shared<Task1428>(tensor1428, pindex);
  task1427->add_dep(task1428);
  task1428->add_dep(task1410);
  normq->add_task(task1428);

  vector<shared_ptr<Tensor>> tensor1429 = {I1791, t2};
  auto task1429 = make_shared<Task1429>(tensor1429, pindex);
  task1428->add_dep(task1429);
  task1429->add_dep(task1410);
  normq->add_task(task1429);

  vector<IndexRange> I1793_index = {active_, active_, active_, active_};
  auto I1793 = make_shared<Tensor>(I1793_index);
  vector<shared_ptr<Tensor>> tensor1430 = {I1790, t2, I1793};
  auto task1430 = make_shared<Task1430>(tensor1430, pindex);
  task1427->add_dep(task1430);
  task1430->add_dep(task1410);
  normq->add_task(task1430);

  vector<shared_ptr<Tensor>> tensor1431 = {I1793, Gamma29_()};
  auto task1431 = make_shared<Task1431>(tensor1431, pindex);
  task1430->add_dep(task1431);
  task1431->add_dep(task1410);
  normq->add_task(task1431);

  vector<IndexRange> I1794_index = {active_, active_, active_, virt_};
  auto I1794 = make_shared<Tensor>(I1794_index);
  vector<shared_ptr<Tensor>> tensor1432 = {n, I1794};
  auto task1432 = make_shared<Task1432>(tensor1432, pindex);
  task1432->add_dep(task1410);
  normq->add_task(task1432);

  vector<IndexRange> I1795_index = {active_, active_, active_, active_, active_, active_};
  auto I1795 = make_shared<Tensor>(I1795_index);
  vector<shared_ptr<Tensor>> tensor1433 = {I1794, t2, I1795};
  auto task1433 = make_shared<Task1433>(tensor1433, pindex);
  task1432->add_dep(task1433);
  task1433->add_dep(task1410);
  normq->add_task(task1433);

  vector<shared_ptr<Tensor>> tensor1434 = {I1795, Gamma50_()};
  auto task1434 = make_shared<Task1434>(tensor1434, pindex);
  task1433->add_dep(task1434);
  task1434->add_dep(task1410);
  normq->add_task(task1434);

  shared_ptr<Tensor> I1796;
  if (diagonal) {
    vector<IndexRange> I1796_index = {closed_, virt_, closed_, virt_};
    I1796 = make_shared<Tensor>(I1796_index);
  }
  shared_ptr<Task1435> task1435;
  if (diagonal) {
    vector<shared_ptr<Tensor>> tensor1435 = {n, I1796};
    task1435 = make_shared<Task1435>(tensor1435, pindex);
    task1435->add_dep(task1410);
    normq->add_task(task1435);
  }

  shared_ptr<Task1436> task1436;
  if (diagonal) {
    vector<shared_ptr<Tensor>> tensor1436 = {I1796, t2};
    task1436 = make_shared<Task1436>(tensor1436, pindex);
    task1435->add_dep(task1436);
    task1436->add_dep(task1410);
    normq->add_task(task1436);
  }

  vector<IndexRange> I1798_index = {active_, virt_, closed_, virt_};
  auto I1798 = make_shared<Tensor>(I1798_index);
  vector<shared_ptr<Tensor>> tensor1437 = {n, I1798};
  auto task1437 = make_shared<Task1437>(tensor1437, pindex);
  task1437->add_dep(task1410);
  normq->add_task(task1437);

  vector<IndexRange> I1799_index = {active_, active_};
  auto I1799 = make_shared<Tensor>(I1799_index);
  vector<shared_ptr<Tensor>> tensor1438 = {I1798, t2, I1799};
  auto task1438 = make_shared<Task1438>(tensor1438, pindex);
  task1437->add_dep(task1438);
  task1438->add_dep(task1410);
  normq->add_task(task1438);

  vector<shared_ptr<Tensor>> tensor1439 = {I1799, Gamma32_()};
  auto task1439 = make_shared<Task1439>(tensor1439, pindex);
  task1438->add_dep(task1439);
  task1439->add_dep(task1410);
  normq->add_task(task1439);

  vector<IndexRange> I1801_index = {active_, active_};
  auto I1801 = make_shared<Tensor>(I1801_index);
  vector<shared_ptr<Tensor>> tensor1440 = {I1798, t2, I1801};
  auto task1440 = make_shared<Task1440>(tensor1440, pindex);
  task1437->add_dep(task1440);
  task1440->add_dep(task1410);
  normq->add_task(task1440);

  vector<shared_ptr<Tensor>> tensor1441 = {I1801, Gamma32_()};
  auto task1441 = make_shared<Task1441>(tensor1441, pindex);
  task1440->add_dep(task1441);
  task1441->add_dep(task1410);
  normq->add_task(task1441);

  vector<IndexRange> I1802_index = {virt_, virt_, active_, active_};
  auto I1802 = make_shared<Tensor>(I1802_index);
  vector<shared_ptr<Tensor>> tensor1442 = {n, I1802};
  auto task1442 = make_shared<Task1442>(tensor1442, pindex);
  task1442->add_dep(task1410);
  normq->add_task(task1442);

  vector<IndexRange> I1803_index = {active_, virt_, active_, virt_};
  auto I1803 = make_shared<Tensor>(I1803_index);
  vector<shared_ptr<Tensor>> tensor1443 = {I1802, Gamma51_(), I1803};
  auto task1443 = make_shared<Task1443>(tensor1443, pindex);
  task1442->add_dep(task1443);
  task1443->add_dep(task1410);
  normq->add_task(task1443);

  vector<shared_ptr<Tensor>> tensor1444 = {I1803, t2};
  auto task1444 = make_shared<Task1444>(tensor1444, pindex);
  task1443->add_dep(task1444);
  task1444->add_dep(task1410);
  normq->add_task(task1444);

  return normq;
}


#endif
