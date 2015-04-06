//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks25.cc
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

#include <src/smith/MRCI_tasks25.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task1200::Task_local::compute() {
  const Index x5 = b(0);
  const Index a1 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1456
  std::unique_ptr<double[]> odata = out()->move_block(x5, a1, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a1, x4, x3);
}

void Task1201::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, a4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), a4.size(), a1.size());
    // tensor label: I1458
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a3.size()*a1.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c2.size()*a3.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1202::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I1458
  std::unique_ptr<double[]> odata = out()->move_block(a4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1459
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->put_block(odata, a4, x0);
}

void Task1203::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1459
  std::unique_ptr<double[]> odata = out()->move_block(x3, a4, x2, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), a4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a4, x2, x1);
}

void Task1204::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, a4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), a4.size(), a3.size());
    // tensor label: I1461
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a1.size()*a3.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, c2.size()*a1.size()*a3.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1205::Task_local::compute() {
  const Index a4 = b(0);
  const Index x0 = b(1);
  // tensor label: I1461
  std::unique_ptr<double[]> odata = out()->move_block(a4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1462
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a4, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a4, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a4.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a4.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a4.size());
  out()->put_block(odata, a4, x0);
}

void Task1206::Task_local::compute() {
  const Index x3 = b(0);
  const Index a4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1462
  std::unique_ptr<double[]> odata = out()->move_block(x3, a4, x2, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a4, x2, x1);
}

void Task1207::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, c5, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, c5, a1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), c5.size(), a1.size());
      // tensor label: I1476
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c2, c4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c2, c4, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c5.size(), c2.size(), c4.size(), x0.size());
      dgemm_("T", "N", a3.size()*a1.size(), c2.size()*x0.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a3.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1208::Task_local::compute() {
  const Index c5 = b(0);
  const Index c2 = b(1);
  const Index c4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1476
  std::unique_ptr<double[]> odata = out()->move_block(c5, c2, c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c2, c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c2, c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1477
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c5, c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c5, c2, c4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c5.size(), c2.size(), c4.size());
    dgemm_("T", "N", x0.size(), c5.size()*c2.size()*c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c5.size(), c2.size(), c4.size());
  out()->put_block(odata, c5, c2, c4, x0);
}

void Task1209::Task_local::compute() {
  const Index x1 = b(0);
  const Index c5 = b(1);
  const Index c2 = b(2);
  const Index c4 = b(3);
  // tensor label: I1477
  std::unique_ptr<double[]> odata = out()->move_block(x1, c5, c2, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c5, c2, c4);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, x1.size(), c5.size(), c2.size(), c4.size());
  }
  out()->put_block(odata, x1, c5, c2, c4);
}

void Task1210::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c5 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a1, c5, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a1, c5, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a1.size(), c5.size(), a3.size());
      // tensor label: I1479
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c2, c4, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c2, c4, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c5.size(), c2.size(), c4.size(), x0.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), c5.size()*c4.size(),
             1.0, i0data_sorted, c5.size()*c4.size(), i1data_sorted, c5.size()*c4.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1211::Task_local::compute() {
  const Index c5 = b(0);
  const Index c2 = b(1);
  const Index c4 = b(2);
  const Index x0 = b(3);
  // tensor label: I1479
  std::unique_ptr<double[]> odata = out()->move_block(c5, c2, c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c2, c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c2, c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I1480
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c5, c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c5, c2, c4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c5.size(), c2.size(), c4.size());
    dgemm_("T", "N", x0.size(), c5.size()*c2.size()*c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c5.size(), c2.size(), c4.size());
  out()->put_block(odata, c5, c2, c4, x0);
}

void Task1212::Task_local::compute() {
  const Index x1 = b(0);
  const Index c5 = b(1);
  const Index c2 = b(2);
  const Index c4 = b(3);
  // tensor label: I1480
  std::unique_ptr<double[]> odata = out()->move_block(x1, c5, c2, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c5, c2, c4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x1.size(), c5.size(), c2.size(), c4.size());
  }
  out()->put_block(odata, x1, c5, c2, c4);
}

void Task1213::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c4, a1)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c4.size(), a1.size());
      // tensor label: I1506
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a3.size()*a1.size(), c2.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, a3.size()*a1.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1214::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1506
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1507
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c2.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c2.size(), c4.size());
  out()->put_block(odata, c2, c4, x3, x0);
}

void Task1215::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1507
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c2.size(), c4.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, c2, c4);
    sort_indices<2,3,0,1,1,1,1,2>(i1data, odata, x2.size(), x1.size(), c2.size(), c4.size());
  }
  out()->put_block(odata, c2, c4, x2, x1);
}

void Task1216::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1506
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1525
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x2, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x2, x1, c4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x2.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), c2.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c2.size(), c4.size());
  out()->put_block(odata, c2, c4, x3, x0);
}

void Task1217::Task_local::compute() {
  const Index c2 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1525
  std::unique_ptr<double[]> odata = out()->move_block(c2, x2, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, x1, c4);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c2.size(), x2.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, c2, x2, x1, c4);
}

void Task1218::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1506
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1543
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c4.size()*c2.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c4.size(), c2.size());
  out()->put_block(odata, c2, c4, x3, x0);
}

void Task1219::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index c2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1543
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, c2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, c2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x2.size(), c4.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x2, c4, c2, x1);
}

void Task1220::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c4, a3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c4.size(), a3.size());
      // tensor label: I1509
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a1.size()*a3.size(), c2.size()*x0.size(), c4.size()*x3.size(),
             1.0, i0data_sorted, c4.size()*x3.size(), i1data_sorted, c4.size()*x3.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1221::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1509
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1510
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c2.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c2.size(), c4.size());
  out()->put_block(odata, c2, c4, x3, x0);
}

void Task1222::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1510
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), c4.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, c4, c2, x1);
    sort_indices<2,1,0,3,1,1,1,2>(i1data, odata, x2.size(), c4.size(), c2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x2, x1, c2, c4);
    sort_indices<2,3,0,1,1,1,-1,1>(i2data, odata, x2.size(), x1.size(), c2.size(), c4.size());
  }
  out()->put_block(odata, c2, c4, x2, x1);
}

void Task1223::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1509
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I1528
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x2, x1, c4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x2, x1, c4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x2.size(), x1.size(), c4.size());
      dgemm_("T", "N", x3.size()*x0.size(), c2.size()*c4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c2.size(), c4.size());
  out()->put_block(odata, c2, c4, x3, x0);
}

void Task1224::Task_local::compute() {
  const Index c2 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c4 = b(3);
  // tensor label: I1528
  std::unique_ptr<double[]> odata = out()->move_block(c2, x2, x1, c4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2, x1, c4);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c2.size(), x2.size(), x1.size(), c4.size());
  }
  out()->put_block(odata, c2, x2, x1, c4);
}

void Task1225::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c2, a3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a3.size());
      // tensor label: I1512
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), x3.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1226::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1512
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1513
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a1.size());
  out()->put_block(odata, a4, a1, x3, x0);
}

void Task1227::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1513
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), a1.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, a4, a1);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x2.size(), x1.size(), a4.size(), a1.size());
  }
  out()->put_block(odata, a4, a1, x2, x1);
}

void Task1228::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1512
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1531
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a1.size());
  out()->put_block(odata, a4, a1, x3, x0);
}

void Task1229::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1531
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, x1, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, x1, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), x2.size(), x1.size(), a1.size());
  }
  out()->put_block(odata, a4, x2, x1, a1);
}

void Task1230::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1512
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1549
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a1.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a1.size(), a4.size());
  out()->put_block(odata, a4, a1, x3, x0);
}

void Task1231::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1549
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, a4, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, a4, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), a1.size(), a4.size(), x1.size());
  }
  out()->put_block(odata, x2, a1, a4, x1);
}

void Task1232::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a4.size());
      // tensor label: I1515
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), x3.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), a1.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1233::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1515
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1516
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a1.size());
  out()->put_block(odata, a4, a1, x3, x0);
}

void Task1234::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1516
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), a1.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, a4, a1);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x2.size(), x1.size(), a4.size(), a1.size());
  }
  out()->put_block(odata, a4, a1, x2, x1);
}

void Task1235::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1515
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1534
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a1.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a1.size());
  out()->put_block(odata, a4, a1, x3, x0);
}

void Task1236::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  // tensor label: I1534
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, x1, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, x1, a1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), x2.size(), x1.size(), a1.size());
  }
  out()->put_block(odata, a4, x2, x1, a1);
}

void Task1237::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1515
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a1, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a1, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1552
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a1.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a1.size(), a4.size());
  out()->put_block(odata, a4, a1, x3, x0);
}

void Task1238::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1552
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, a4, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, a4, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x2.size(), a1.size(), a4.size(), x1.size());
  }
  out()->put_block(odata, x2, a1, a4, x1);
}

void Task1239::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c2, a1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a1.size());
      // tensor label: I1518
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x3.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1240::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1518
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1519
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->put_block(odata, a4, a3, x3, x0);
}

void Task1241::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1519
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), a3.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, x1, a4, a3);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x2.size(), x1.size(), a4.size(), a3.size());
  }
  out()->put_block(odata, a4, a3, x2, x1);
}

void Task1242::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1518
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I1537
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a3.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->put_block(odata, a4, a3, x3, x0);
}

void Task1243::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  // tensor label: I1537
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, x1, a3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, x1, a3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a4.size(), x2.size(), x1.size(), a3.size());
  }
  out()->put_block(odata, a4, x2, x1, a3);
}

void Task1244::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1518
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma503
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x2.size(), x0.size());
      // tensor label: I1555
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a3, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a3, a4, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a3.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a3.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a3.size(), a4.size());
  out()->put_block(odata, a4, a3, x3, x0);
}

void Task1245::Task_local::compute() {
  const Index x2 = b(0);
  const Index a3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1555
  std::unique_ptr<double[]> odata = out()->move_block(x2, a3, a4, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a3, a4, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x2.size(), a3.size(), a4.size(), x1.size());
  }
  out()->put_block(odata, x2, a3, a4, x1);
}

void Task1246::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a4.size());
      // tensor label: I1521
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3, x3, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x3.size(), x0.size());
      dgemm_("T", "N", a1.size()*c2.size(), a3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task1247::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1521
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma51
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I1522
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3, x2, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size(), x2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->put_block(odata, a4, a3, x3, x0);
}

void Task1248::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1522
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3, x2, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), a3.size(), x2.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x2, a3, a4, x1);
    sort_indices<2,1,0,3,1,1,-1,2>(i1data, odata, x2.size(), a3.size(), a4.size(), x1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x2, x1, a4, a3);
    sort_indices<2,3,0,1,1,1,1,1>(i2data, odata, x2.size(), x1.size(), a4.size(), a3.size());
  }
  out()->put_block(odata, a4, a3, x2, x1);
}

void Task1249::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1521
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a3, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a3, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I1540
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, x1, a3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, x1, a3)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), x1.size(), a3.size());
      dgemm_("T", "N", x3.size()*x0.size(), a4.size()*a3.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), a4.size(), a3.size());
  out()->put_block(odata, a4, a3, x3, x0);
}

#endif
