//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks29.cc
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


#include <src/smith/CASPT2_tasks29.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1400::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1533
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
    // tensor label: I1534
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c2, a1);
}

void Task1401::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1534
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1402::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1533
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());
    // tensor label: I1538
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c2, a1);
}

void Task1403::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1538
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1404::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x1.size());
    // tensor label: I1541
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a3.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x0.size()*a1.size()*a3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,4,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x0.size(), a1.size(), a3.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1405::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);
  // tensor label: I1541
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
      // tensor label: I1542
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a3);
}

void Task1406::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1542
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1407::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1949
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1408::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1949
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1409::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1952
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1410::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1952
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1411::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I2039
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1412::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I2039
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1413::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I2042
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1414::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I2042
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1415::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x1 : *range_[1]) {
        for (auto& a2 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
          // tensor label: I1544
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a2)]);
          sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a2.size());
          dgemm_("T", "N", 1, ci0.size(), x0.size()*x1.size()*a1.size()*a2.size(),
                 1.0, i0data_sorted, x0.size()*x1.size()*a1.size()*a2.size(), i1data_sorted, x0.size()*x1.size()*a1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1416::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1544
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a2.size());
    // tensor label: I1545
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x0.size()*x1.size()*a1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}

void Task1417::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  // tensor label: I1545
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        // tensor label: I1546
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<4,3,1,0,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1);
}

void Task1418::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1546
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma437
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}

void Task1419::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1544
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
      // tensor label: I1549
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x2, x1, a1, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x2, x1, a1, c3, a2)]);
      sort_indices<5,2,0,1,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, ci0.size()*x0.size()*x1.size()*a1.size()*a2.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), x1.size(), a1.size(), a2.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}

void Task1420::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);
  // tensor label: I1549
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x2, x1, a1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x2, x1, a1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x2, x1, a1, c3, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());
    // tensor label: I1550
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size()*c3.size()*a2.size(), ci0.size()*x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c3.size()*a2.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), ci0.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x0, x2, x1, a1, c3, a2);
}

void Task1421::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1550
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1422::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1544
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
      // tensor label: I1553
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}

void Task1423::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  // tensor label: I1553
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  {
    // tensor label: Gamma470
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}

void Task1424::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1544
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I1556
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, a3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), a3.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x0.size()*x1.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}

void Task1425::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);
  // tensor label: I1556
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a3.size());
      // tensor label: I1557
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a3);
}

void Task1426::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1557
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1427::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1544
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      // tensor label: I1955
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}

void Task1428::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1955
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    dscal_(ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1429::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1544
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      // tensor label: I2045
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, a2);
}

void Task1430::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I2045
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1431::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& c2 : *range_[0]) {
        for (auto& x4 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c2, x4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
          // tensor label: I1559
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, c1, c2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, c1, c2)]);
          sort_indices<3,1,4,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), c1.size(), c2.size());
          dgemm_("T", "N", 1, ci0.size(), x5.size()*x4.size()*c1.size()*c2.size(),
                 1.0, i0data_sorted, x5.size()*x4.size()*c1.size()*c2.size(), i1data_sorted, x5.size()*x4.size()*c1.size()*c2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1432::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1559
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c2, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
      // tensor label: I1560
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x5, x1, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x5, x1, x4)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x5.size(), x1.size(), x4.size());
      dgemm_("T", "N", c1.size()*c2.size(), ci0.size()*x5.size()*x4.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), ci0.size(), x5.size(), x4.size());
  out()->put_block(odata, ci0, x5, x4, c1, c2);
}

void Task1433::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x5 = b(2);
  const Index x1 = b(3);
  const Index x4 = b(4);
  // tensor label: I1560
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x5, x1, x4);
  {
    // tensor label: Gamma378
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x5, x1, x4);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x0.size(), x5.size(), x1.size(), x4.size());
  }
  out()->put_block(odata, ci0, x0, x5, x1, x4);
}

void Task1434::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1559
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size());
    // tensor label: I1567
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x3, x4, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x3, x4, c1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x3.size(), x4.size(), c1.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x5.size()*x4.size()*c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x5.size(), x4.size(), c1.size());
  out()->put_block(odata, ci0, x5, x4, c1, c2);
}

void Task1435::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index c1 = b(4);
  // tensor label: I1567
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x3, x4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x3, x4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x3, x4, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        // tensor label: I1568
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x2, x5, x3, x4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x2, x5, x3, x4, x1, x0)]);
        sort_indices<1,5,6,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), ci0.size()*x5.size()*x3.size()*x4.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x5.size(), x3.size(), x4.size());
  out()->put_block(odata, ci0, x5, x3, x4, c1);
}

void Task1436::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x5 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: I1568
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x5, x3, x4, x1, x0);
  {
    // tensor label: Gamma382
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x5, x3, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x2, x5, x3, x4, x1, x0);
}

void Task1437::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      for (auto& c2 : *range_[0]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c2, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c2, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
          // tensor label: I1562
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, c1, c2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, c1, c2)]);
          sort_indices<3,1,4,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), c1.size(), c2.size());
          dgemm_("T", "N", 1, ci0.size(), x3.size()*x2.size()*c1.size()*c2.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*c1.size()*c2.size(), i1data_sorted, x3.size()*x2.size()*c1.size()*c2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1438::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1562
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, c1, c2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c3.size());
    // tensor label: I1563
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, c1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, c1, c3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), c1.size(), c3.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x3.size()*x2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x3.size(), x2.size(), c1.size());
  out()->put_block(odata, ci0, x3, x2, c1, c2);
}

void Task1439::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index c3 = b(4);
  // tensor label: I1563
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, c1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, c1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, c1, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c3, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c3.size(), x1.size());
      // tensor label: I1564
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x3, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x3, x1, x2)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x3.size(), x1.size(), x2.size());
      dgemm_("T", "N", c1.size()*c3.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, c1, c3);
}

void Task1440::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1564
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x3, x1, x2);
  {
    // tensor label: Gamma379
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x3, x1, x2);
    sort_indices<0,1,2,3,4,1,1,-1,1>(i0data, odata, ci0.size(), x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x0, x3, x1, x2);
}

void Task1441::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1562
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, c1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x1)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), x1.size());
      // tensor label: I1571
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x3, x2, c1, a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x3, x2, c1, a3, c2)]);
      sort_indices<1,5,0,2,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x3.size(), x2.size(), c1.size(), a3.size(), c2.size());
      dgemm_("T", "N", 1, ci0.size()*x3.size()*x2.size()*c1.size()*c2.size(), x1.size()*a3.size(),
             1.0, i0data_sorted, x1.size()*a3.size(), i1data_sorted, x1.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x3.size(), x2.size(), c1.size(), c2.size());
  out()->put_block(odata, ci0, x3, x2, c1, c2);
}

void Task1442::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index c1 = b(4);
  const Index a3 = b(5);
  const Index c2 = b(6);
  // tensor label: I1571
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x3, x2, c1, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x3, x2, c1, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x3, x2, c1, a3, c2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x0.size());
    // tensor label: I1572
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x3, x0, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x3, x0, x2)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x3.size(), x0.size(), x2.size());
    dgemm_("T", "N", c1.size()*a3.size()*c2.size(), ci0.size()*x1.size()*x3.size()*x2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*a3.size()*c2.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a3.size(), c2.size(), ci0.size(), x1.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x1, x3, x2, c1, a3, c2);
}

void Task1443::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  // tensor label: I1572
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x3, x0, x2);
  {
    // tensor label: Gamma381
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x3, x0, x2);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x3.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, ci0, x1, x3, x0, x2);
}

void Task1444::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1562
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, c1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c2, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
      // tensor label: I1958
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x3, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x3, x1, x2)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x3.size(), x1.size(), x2.size());
      dgemm_("T", "N", c1.size()*c2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, c1, c2);
}

void Task1445::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1958
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x3, x1, x2);
  {
    // tensor label: Gamma379
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x3, x1, x2);
    dscal_(ci0.size()*x0.size()*x3.size()*x1.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x0, x3, x1, x2);
}

void Task1446::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1562
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, c1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c2, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c2.size(), x1.size());
      // tensor label: I2048
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x3, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x3, x1, x2)]);
      sort_indices<3,1,0,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x3.size(), x1.size(), x2.size());
      dgemm_("T", "N", c1.size()*c2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, c1, c2);
}

void Task1447::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I2048
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x3, x1, x2);
  {
    // tensor label: Gamma379
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x3, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x0, x3, x1, x2);
}

void Task1448::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& c1 : *range_[0]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
          // tensor label: I1574
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, c1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, c1)]);
          sort_indices<1,2,4,3,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), c1.size());
          dgemm_("T", "N", 1, ci0.size(), x5.size()*x4.size()*x3.size()*c1.size(),
                 1.0, i0data_sorted, x5.size()*x4.size()*x3.size()*c1.size(), i1data_sorted, x5.size()*x4.size()*x3.size()*c1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1449::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1574
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, c1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c2)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x2.size(), c2.size());
      // tensor label: I1575
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, x2, c1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, x2, c1, c2)]);
      sort_indices<6,4,0,1,2,3,5,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), x2.size(), c1.size(), c2.size());
      dgemm_("T", "N", 1, ci0.size()*x5.size()*x4.size()*x3.size()*c1.size(), x2.size()*c2.size(),
             1.0, i0data_sorted, x2.size()*c2.size(), i1data_sorted, x2.size()*c2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x5.size(), x4.size(), x3.size(), c1.size());
  out()->put_block(odata, ci0, x5, x4, x3, c1);
}

