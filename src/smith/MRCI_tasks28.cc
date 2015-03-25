//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks28.cc
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

#include <src/smith/MRCI_tasks28.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1350::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I1642
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, x2, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), c4.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, x3, c4, x2, c3);
}

void Task1351::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1640
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1674
      std::unique_ptr<double[]> i1data = in(1)->get_block(a1, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), a2.size());
  out()->put_block(odata, x1, x0, a1, a2);
}

void Task1352::Task_local::compute() {
  const Index a1 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1674
  std::unique_ptr<double[]> odata = out()->move_block(a1, a2, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a2, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a2, x3, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, x2, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), x2.size(), a4.size());
      // tensor label: I1675
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, a3, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), a1.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), a1.size(), a2.size());
  out()->put_block(odata, a1, a2, x3, x2);
}

void Task1353::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1675
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1, a3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a4.size(), a1.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, a4, a1, a3, a2);
}

void Task1354::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1640
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma566
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x1.size());
      // tensor label: I1713
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, x1, x0, a1, a2);
}

void Task1355::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index a2 = b(3);
  // tensor label: I1713
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x5.size(), a1.size(), x4.size(), a2.size());
  }
  out()->put_block(odata, x5, a1, x4, a2);
}

void Task1356::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1640
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma567
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x1.size());
      // tensor label: I1715
      std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a1, x6, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a1, x6, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, x7.size(), a1.size(), x6.size(), a2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*a2.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, x1, x0, a1, a2);
}

void Task1357::Task_local::compute() {
  const Index x7 = b(0);
  const Index a1 = b(1);
  const Index x6 = b(2);
  const Index a2 = b(3);
  // tensor label: I1715
  std::unique_ptr<double[]> odata = out()->move_block(x7, a1, x6, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x7.size(), a1.size(), x6.size(), a2.size());
  }
  out()->put_block(odata, x7, a1, x6, a2);
}

void Task1359::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x0, x1);
  {
    // tensor label: I1732
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, x1, x0);
    sort_indices<0,1,3,2,1,1,1,1>(i0data, odata, c1.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, c1, x2, x0, x1);
}

void Task1360::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1732
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma5
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    // tensor label: I1733
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size());
    dgemm_("T", "N", x2.size()*x1.size()*x0.size(), c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x2.size()*x1.size()*x0.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task1361::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  // tensor label: I1733
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c1.size(), x3.size());
  }
  out()->put_block(odata, c1, x3);
}

void Task1362::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1732
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        // tensor label: I1743
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task1363::Task_local::compute() {
  const Index x2 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1743
  std::unique_ptr<double[]> odata = out()->move_block(x2, x5, x4, x3, x1, x0);
  {
    // tensor label: Gamma104
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x5, x4, x3, x1, x0);
}

void Task1364::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1732
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        // tensor label: I1745
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x2, x1, x0);
}

void Task1365::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I1745
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma4
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task1366::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  {
    // tensor label: I1734
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, c2, a1);
    sort_indices<2,0,1,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), c2.size(), a1.size());
  }
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task1367::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1734
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I1735
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task1368::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1735
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task1369::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1734
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I1751
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task1370::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1751
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma29
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1371::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1734
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3, x2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3, x2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size(), x2.size(), a1.size());
      // tensor label: I1753
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x3, x2, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x3, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task1372::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I1753
  std::unique_ptr<double[]> odata = out()->move_block(x1, x3, x2, x0);
  {
    // tensor label: Gamma25
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x1, x3, x2, x0);
}

void Task1373::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1734
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I1755
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task1374::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I1755
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x1, x2);
  {
    // tensor label: Gamma27
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x3, x0, x1, x2);
}

void Task1375::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1734
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c2.size(), a1.size());
      // tensor label: I1757
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task1376::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1757
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma29
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1377::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, a2);
  {
    // tensor label: I1736
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, c1, a2);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), c1.size(), a2.size());
  }
  out()->put_block(odata, x0, x1, c1, a2);
}

void Task1378::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1736
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: I1737
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task1379::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1737
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task1380::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1736
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I1759
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task1381::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1759
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma29
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1382::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1736
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, x2, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), x2.size(), a2.size());
      // tensor label: I1761
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task1383::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1761
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma5
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task1384::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1736
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1763
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task1385::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I1763
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, c1, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, c1, x2);
}

void Task1386::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1736
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), a2.size());
      // tensor label: I1765
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, a2);
}

void Task1387::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1765
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma29
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1388::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, a1);
  {
    // tensor label: I1738
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x2, x1, a1);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, x0.size(), x2.size(), x1.size(), a1.size());
  }
  out()->put_block(odata, x1, x2, x0, a1);
}

void Task1389::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1738
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());
    // tensor label: I1739
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}

void Task1390::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1739
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma51
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task1391::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1738
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1767
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}

void Task1392::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1767
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

void Task1393::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1738
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, a1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), a1.size());
        // tensor label: I1769
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a1);
}

void Task1394::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1769
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma49
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x5, x4, x3, x0, x2, x1);
}

void Task1395::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, c1, x0);
  {
    // tensor label: I1740
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, c2);
    sort_indices<3,1,2,0,1,1,1,1>(i0data, odata, x0.size(), x1.size(), c1.size(), c2.size());
  }
  out()->put_block(odata, c2, x1, c1, x0);
}

void Task1396::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c2 = b(3);
  // tensor label: I1740
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c2, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      // tensor label: I1741
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3, x1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size(), x1.size(), x2.size());
      dgemm_("T", "N", c1.size()*c2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, c1, c2);
}

void Task1397::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I1741
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x1, x2);
  {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x0, x3, x1, x2);
}

void Task1398::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0, c1, a2);
  {
    // tensor label: I1746
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c3, c1, a2);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, x0.size(), c3.size(), c1.size(), a2.size());
  }
  out()->put_block(odata, c3, x0, c1, a2);
}

void Task1399::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1746
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), c1.size(), a2.size());
    // tensor label: I1747
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*c1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, c3, c1, a2);
}

#endif
