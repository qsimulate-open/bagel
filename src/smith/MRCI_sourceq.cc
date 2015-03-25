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

shared_ptr<Queue> MRCI::MRCI::make_sourceq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor1358 = {s};
  auto task1358 = make_shared<Task1358>(tensor1358);
  sourceq->add_task(task1358);

  vector<IndexRange> I1732_index = {closed_, active_, active_, active_};
  auto I1732 = make_shared<Tensor>(I1732_index);
  vector<shared_ptr<Tensor>> tensor1359 = {s, I1732};
  auto task1359 = make_shared<Task1359>(tensor1359, pindex);
  task1359->add_dep(task1358);
  sourceq->add_task(task1359);

  vector<IndexRange> I1733_index = {closed_, active_};
  auto I1733 = make_shared<Tensor>(I1733_index);
  vector<shared_ptr<Tensor>> tensor1360 = {I1732, Gamma5_(), I1733};
  auto task1360 = make_shared<Task1360>(tensor1360, pindex);
  task1359->add_dep(task1360);
  task1360->add_dep(task1358);
  sourceq->add_task(task1360);

  vector<shared_ptr<Tensor>> tensor1361 = {I1733, h1_};
  auto task1361 = make_shared<Task1361>(tensor1361, pindex);
  task1360->add_dep(task1361);
  task1361->add_dep(task1358);
  sourceq->add_task(task1361);

  vector<IndexRange> I1743_index = {active_, active_, active_, active_, active_, active_};
  auto I1743 = make_shared<Tensor>(I1743_index);
  vector<shared_ptr<Tensor>> tensor1362 = {I1732, v2_, I1743};
  auto task1362 = make_shared<Task1362>(tensor1362, pindex);
  task1359->add_dep(task1362);
  task1362->add_dep(task1358);
  sourceq->add_task(task1362);

  vector<shared_ptr<Tensor>> tensor1363 = {I1743, Gamma104_()};
  auto task1363 = make_shared<Task1363>(tensor1363, pindex);
  task1362->add_dep(task1363);
  task1363->add_dep(task1358);
  sourceq->add_task(task1363);

  vector<IndexRange> I1745_index = {active_, active_, active_, active_, active_, active_};
  auto I1745 = make_shared<Tensor>(I1745_index);
  vector<shared_ptr<Tensor>> tensor1364 = {I1732, v2_, I1745};
  auto task1364 = make_shared<Task1364>(tensor1364, pindex);
  task1359->add_dep(task1364);
  task1364->add_dep(task1358);
  sourceq->add_task(task1364);

  vector<shared_ptr<Tensor>> tensor1365 = {I1745, Gamma4_()};
  auto task1365 = make_shared<Task1365>(tensor1365, pindex);
  task1364->add_dep(task1365);
  task1365->add_dep(task1358);
  sourceq->add_task(task1365);

  vector<IndexRange> I1734_index = {active_, active_, closed_, virt_};
  auto I1734 = make_shared<Tensor>(I1734_index);
  vector<shared_ptr<Tensor>> tensor1366 = {s, I1734};
  auto task1366 = make_shared<Task1366>(tensor1366, pindex);
  task1366->add_dep(task1358);
  sourceq->add_task(task1366);

  vector<IndexRange> I1735_index = {active_, active_};
  auto I1735 = make_shared<Tensor>(I1735_index);
  vector<shared_ptr<Tensor>> tensor1367 = {I1734, h1_, I1735};
  auto task1367 = make_shared<Task1367>(tensor1367, pindex);
  task1366->add_dep(task1367);
  task1367->add_dep(task1358);
  sourceq->add_task(task1367);

  vector<shared_ptr<Tensor>> tensor1368 = {I1735, Gamma32_()};
  auto task1368 = make_shared<Task1368>(tensor1368, pindex);
  task1367->add_dep(task1368);
  task1368->add_dep(task1358);
  sourceq->add_task(task1368);

  vector<IndexRange> I1751_index = {active_, active_, active_, active_};
  auto I1751 = make_shared<Tensor>(I1751_index);
  vector<shared_ptr<Tensor>> tensor1369 = {I1734, v2_, I1751};
  auto task1369 = make_shared<Task1369>(tensor1369, pindex);
  task1366->add_dep(task1369);
  task1369->add_dep(task1358);
  sourceq->add_task(task1369);

  vector<shared_ptr<Tensor>> tensor1370 = {I1751, Gamma29_()};
  auto task1370 = make_shared<Task1370>(tensor1370, pindex);
  task1369->add_dep(task1370);
  task1370->add_dep(task1358);
  sourceq->add_task(task1370);

  vector<IndexRange> I1753_index = {active_, active_, active_, active_};
  auto I1753 = make_shared<Tensor>(I1753_index);
  vector<shared_ptr<Tensor>> tensor1371 = {I1734, v2_, I1753};
  auto task1371 = make_shared<Task1371>(tensor1371, pindex);
  task1366->add_dep(task1371);
  task1371->add_dep(task1358);
  sourceq->add_task(task1371);

  vector<shared_ptr<Tensor>> tensor1372 = {I1753, Gamma25_()};
  auto task1372 = make_shared<Task1372>(tensor1372, pindex);
  task1371->add_dep(task1372);
  task1372->add_dep(task1358);
  sourceq->add_task(task1372);

  vector<IndexRange> I1755_index = {active_, active_, active_, active_};
  auto I1755 = make_shared<Tensor>(I1755_index);
  vector<shared_ptr<Tensor>> tensor1373 = {I1734, v2_, I1755};
  auto task1373 = make_shared<Task1373>(tensor1373, pindex);
  task1366->add_dep(task1373);
  task1373->add_dep(task1358);
  sourceq->add_task(task1373);

  vector<shared_ptr<Tensor>> tensor1374 = {I1755, Gamma27_()};
  auto task1374 = make_shared<Task1374>(tensor1374, pindex);
  task1373->add_dep(task1374);
  task1374->add_dep(task1358);
  sourceq->add_task(task1374);

  vector<IndexRange> I1757_index = {active_, active_, active_, active_};
  auto I1757 = make_shared<Tensor>(I1757_index);
  vector<shared_ptr<Tensor>> tensor1375 = {I1734, v2_, I1757};
  auto task1375 = make_shared<Task1375>(tensor1375, pindex);
  task1366->add_dep(task1375);
  task1375->add_dep(task1358);
  sourceq->add_task(task1375);

  vector<shared_ptr<Tensor>> tensor1376 = {I1757, Gamma29_()};
  auto task1376 = make_shared<Task1376>(tensor1376, pindex);
  task1375->add_dep(task1376);
  task1376->add_dep(task1358);
  sourceq->add_task(task1376);

  vector<IndexRange> I1736_index = {active_, active_, closed_, virt_};
  auto I1736 = make_shared<Tensor>(I1736_index);
  vector<shared_ptr<Tensor>> tensor1377 = {s, I1736};
  auto task1377 = make_shared<Task1377>(tensor1377, pindex);
  task1377->add_dep(task1358);
  sourceq->add_task(task1377);

  vector<IndexRange> I1737_index = {active_, active_};
  auto I1737 = make_shared<Tensor>(I1737_index);
  vector<shared_ptr<Tensor>> tensor1378 = {I1736, h1_, I1737};
  auto task1378 = make_shared<Task1378>(tensor1378, pindex);
  task1377->add_dep(task1378);
  task1378->add_dep(task1358);
  sourceq->add_task(task1378);

  vector<shared_ptr<Tensor>> tensor1379 = {I1737, Gamma32_()};
  auto task1379 = make_shared<Task1379>(tensor1379, pindex);
  task1378->add_dep(task1379);
  task1379->add_dep(task1358);
  sourceq->add_task(task1379);

  vector<IndexRange> I1759_index = {active_, active_, active_, active_};
  auto I1759 = make_shared<Tensor>(I1759_index);
  vector<shared_ptr<Tensor>> tensor1380 = {I1736, v2_, I1759};
  auto task1380 = make_shared<Task1380>(tensor1380, pindex);
  task1377->add_dep(task1380);
  task1380->add_dep(task1358);
  sourceq->add_task(task1380);

  vector<shared_ptr<Tensor>> tensor1381 = {I1759, Gamma29_()};
  auto task1381 = make_shared<Task1381>(tensor1381, pindex);
  task1380->add_dep(task1381);
  task1381->add_dep(task1358);
  sourceq->add_task(task1381);

  vector<IndexRange> I1761_index = {active_, active_, active_, active_};
  auto I1761 = make_shared<Tensor>(I1761_index);
  vector<shared_ptr<Tensor>> tensor1382 = {I1736, v2_, I1761};
  auto task1382 = make_shared<Task1382>(tensor1382, pindex);
  task1377->add_dep(task1382);
  task1382->add_dep(task1358);
  sourceq->add_task(task1382);

  vector<shared_ptr<Tensor>> tensor1383 = {I1761, Gamma5_()};
  auto task1383 = make_shared<Task1383>(tensor1383, pindex);
  task1382->add_dep(task1383);
  task1383->add_dep(task1358);
  sourceq->add_task(task1383);

  vector<IndexRange> I1763_index = {active_, virt_, closed_, active_};
  auto I1763 = make_shared<Tensor>(I1763_index);
  vector<shared_ptr<Tensor>> tensor1384 = {I1736, Gamma29_(), I1763};
  auto task1384 = make_shared<Task1384>(tensor1384, pindex);
  task1377->add_dep(task1384);
  task1384->add_dep(task1358);
  sourceq->add_task(task1384);

  vector<shared_ptr<Tensor>> tensor1385 = {I1763, v2_};
  auto task1385 = make_shared<Task1385>(tensor1385, pindex);
  task1384->add_dep(task1385);
  task1385->add_dep(task1358);
  sourceq->add_task(task1385);

  vector<IndexRange> I1765_index = {active_, active_, active_, active_};
  auto I1765 = make_shared<Tensor>(I1765_index);
  vector<shared_ptr<Tensor>> tensor1386 = {I1736, v2_, I1765};
  auto task1386 = make_shared<Task1386>(tensor1386, pindex);
  task1377->add_dep(task1386);
  task1386->add_dep(task1358);
  sourceq->add_task(task1386);

  vector<shared_ptr<Tensor>> tensor1387 = {I1765, Gamma29_()};
  auto task1387 = make_shared<Task1387>(tensor1387, pindex);
  task1386->add_dep(task1387);
  task1387->add_dep(task1358);
  sourceq->add_task(task1387);

  vector<IndexRange> I1738_index = {active_, active_, active_, virt_};
  auto I1738 = make_shared<Tensor>(I1738_index);
  vector<shared_ptr<Tensor>> tensor1388 = {s, I1738};
  auto task1388 = make_shared<Task1388>(tensor1388, pindex);
  task1388->add_dep(task1358);
  sourceq->add_task(task1388);

  vector<IndexRange> I1739_index = {active_, active_, active_, active_};
  auto I1739 = make_shared<Tensor>(I1739_index);
  vector<shared_ptr<Tensor>> tensor1389 = {I1738, h1_, I1739};
  auto task1389 = make_shared<Task1389>(tensor1389, pindex);
  task1388->add_dep(task1389);
  task1389->add_dep(task1358);
  sourceq->add_task(task1389);

  vector<shared_ptr<Tensor>> tensor1390 = {I1739, Gamma51_()};
  auto task1390 = make_shared<Task1390>(tensor1390, pindex);
  task1389->add_dep(task1390);
  task1390->add_dep(task1358);
  sourceq->add_task(task1390);

  vector<IndexRange> I1767_index = {active_, virt_, active_, active_};
  auto I1767 = make_shared<Tensor>(I1767_index);
  vector<shared_ptr<Tensor>> tensor1391 = {I1738, Gamma50_(), I1767};
  auto task1391 = make_shared<Task1391>(tensor1391, pindex);
  task1388->add_dep(task1391);
  task1391->add_dep(task1358);
  sourceq->add_task(task1391);

  vector<shared_ptr<Tensor>> tensor1392 = {I1767, v2_};
  auto task1392 = make_shared<Task1392>(tensor1392, pindex);
  task1391->add_dep(task1392);
  task1392->add_dep(task1358);
  sourceq->add_task(task1392);

  vector<IndexRange> I1769_index = {active_, active_, active_, active_, active_, active_};
  auto I1769 = make_shared<Tensor>(I1769_index);
  vector<shared_ptr<Tensor>> tensor1393 = {I1738, v2_, I1769};
  auto task1393 = make_shared<Task1393>(tensor1393, pindex);
  task1388->add_dep(task1393);
  task1393->add_dep(task1358);
  sourceq->add_task(task1393);

  vector<shared_ptr<Tensor>> tensor1394 = {I1769, Gamma49_()};
  auto task1394 = make_shared<Task1394>(tensor1394, pindex);
  task1393->add_dep(task1394);
  task1394->add_dep(task1358);
  sourceq->add_task(task1394);

  vector<IndexRange> I1740_index = {active_, active_, closed_, closed_};
  auto I1740 = make_shared<Tensor>(I1740_index);
  vector<shared_ptr<Tensor>> tensor1395 = {s, I1740};
  auto task1395 = make_shared<Task1395>(tensor1395, pindex);
  task1395->add_dep(task1358);
  sourceq->add_task(task1395);

  vector<IndexRange> I1741_index = {active_, active_, active_, active_};
  auto I1741 = make_shared<Tensor>(I1741_index);
  vector<shared_ptr<Tensor>> tensor1396 = {I1740, v2_, I1741};
  auto task1396 = make_shared<Task1396>(tensor1396, pindex);
  task1395->add_dep(task1396);
  task1396->add_dep(task1358);
  sourceq->add_task(task1396);

  vector<shared_ptr<Tensor>> tensor1397 = {I1741, Gamma0_()};
  auto task1397 = make_shared<Task1397>(tensor1397, pindex);
  task1396->add_dep(task1397);
  task1397->add_dep(task1358);
  sourceq->add_task(task1397);

  vector<IndexRange> I1746_index = {active_, closed_, closed_, virt_};
  auto I1746 = make_shared<Tensor>(I1746_index);
  vector<shared_ptr<Tensor>> tensor1398 = {s, I1746};
  auto task1398 = make_shared<Task1398>(tensor1398, pindex);
  task1398->add_dep(task1358);
  sourceq->add_task(task1398);

  vector<IndexRange> I1747_index = {active_, active_};
  auto I1747 = make_shared<Tensor>(I1747_index);
  vector<shared_ptr<Tensor>> tensor1399 = {I1746, v2_, I1747};
  auto task1399 = make_shared<Task1399>(tensor1399, pindex);
  task1398->add_dep(task1399);
  task1399->add_dep(task1358);
  sourceq->add_task(task1399);

  vector<shared_ptr<Tensor>> tensor1400 = {I1747, Gamma12_()};
  auto task1400 = make_shared<Task1400>(tensor1400, pindex);
  task1399->add_dep(task1400);
  task1400->add_dep(task1358);
  sourceq->add_task(task1400);

  vector<IndexRange> I1749_index = {active_, active_};
  auto I1749 = make_shared<Tensor>(I1749_index);
  vector<shared_ptr<Tensor>> tensor1401 = {I1746, v2_, I1749};
  auto task1401 = make_shared<Task1401>(tensor1401, pindex);
  task1398->add_dep(task1401);
  task1401->add_dep(task1358);
  sourceq->add_task(task1401);

  vector<shared_ptr<Tensor>> tensor1402 = {I1749, Gamma12_()};
  auto task1402 = make_shared<Task1402>(tensor1402, pindex);
  task1401->add_dep(task1402);
  task1402->add_dep(task1358);
  sourceq->add_task(task1402);

  vector<IndexRange> I1770_index = {closed_, virt_, closed_, virt_};
  auto I1770 = make_shared<Tensor>(I1770_index);
  vector<shared_ptr<Tensor>> tensor1403 = {s, I1770};
  auto task1403 = make_shared<Task1403>(tensor1403, pindex);
  task1403->add_dep(task1358);
  sourceq->add_task(task1403);

  vector<shared_ptr<Tensor>> tensor1404 = {I1770, v2_};
  auto task1404 = make_shared<Task1404>(tensor1404, pindex);
  task1403->add_dep(task1404);
  task1404->add_dep(task1358);
  sourceq->add_task(task1404);

  vector<IndexRange> I1772_index = {active_, virt_, closed_, virt_};
  auto I1772 = make_shared<Tensor>(I1772_index);
  vector<shared_ptr<Tensor>> tensor1405 = {s, I1772};
  auto task1405 = make_shared<Task1405>(tensor1405, pindex);
  task1405->add_dep(task1358);
  sourceq->add_task(task1405);

  vector<IndexRange> I1773_index = {active_, active_};
  auto I1773 = make_shared<Tensor>(I1773_index);
  vector<shared_ptr<Tensor>> tensor1406 = {I1772, v2_, I1773};
  auto task1406 = make_shared<Task1406>(tensor1406, pindex);
  task1405->add_dep(task1406);
  task1406->add_dep(task1358);
  sourceq->add_task(task1406);

  vector<shared_ptr<Tensor>> tensor1407 = {I1773, Gamma32_()};
  auto task1407 = make_shared<Task1407>(tensor1407, pindex);
  task1406->add_dep(task1407);
  task1407->add_dep(task1358);
  sourceq->add_task(task1407);

  vector<IndexRange> I1775_index = {active_, active_};
  auto I1775 = make_shared<Tensor>(I1775_index);
  vector<shared_ptr<Tensor>> tensor1408 = {I1772, v2_, I1775};
  auto task1408 = make_shared<Task1408>(tensor1408, pindex);
  task1405->add_dep(task1408);
  task1408->add_dep(task1358);
  sourceq->add_task(task1408);

  vector<shared_ptr<Tensor>> tensor1409 = {I1775, Gamma32_()};
  auto task1409 = make_shared<Task1409>(tensor1409, pindex);
  task1408->add_dep(task1409);
  task1409->add_dep(task1358);
  sourceq->add_task(task1409);

  vector<IndexRange> I1776_index = {virt_, virt_, active_, active_};
  auto I1776 = make_shared<Tensor>(I1776_index);
  vector<shared_ptr<Tensor>> tensor1410 = {s, I1776};
  auto task1410 = make_shared<Task1410>(tensor1410, pindex);
  task1410->add_dep(task1358);
  sourceq->add_task(task1410);

  vector<IndexRange> I1777_index = {active_, virt_, active_, virt_};
  auto I1777 = make_shared<Tensor>(I1777_index);
  vector<shared_ptr<Tensor>> tensor1411 = {I1776, Gamma51_(), I1777};
  auto task1411 = make_shared<Task1411>(tensor1411, pindex);
  task1410->add_dep(task1411);
  task1411->add_dep(task1358);
  sourceq->add_task(task1411);

  vector<shared_ptr<Tensor>> tensor1412 = {I1777, v2_};
  auto task1412 = make_shared<Task1412>(tensor1412, pindex);
  task1411->add_dep(task1412);
  task1412->add_dep(task1358);
  sourceq->add_task(task1412);

  return sourceq;
}


#endif
