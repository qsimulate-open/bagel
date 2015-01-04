//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks35.cc
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


#include <src/smith/CASPT2_tasks35.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task1700::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1542
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1701::Task_local::compute() {
  const Index a3 = b(0);
  const Index a1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I1541
  std::unique_ptr<double[]> odata = out()->move_block(a3, a1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a1, x0, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I1557
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1, a1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size(), a1.size(), x0.size());
    dgemm_("T", "N", a3.size(), x1.size()*a1.size()*x0.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), x1.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, a1, x0, x1);
}

void Task1702::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1557
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x0.size(), a1.size(), x1.size(), a2.size());
  }
  out()->put_block(odata, a2, x1, a1, x0);
}

void Task1703::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: I1549
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a1, a2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a1.size(), a2.size());
      dgemm_("T", "N", x1.size()*x0.size(), x2.size()*x3.size(), a1.size()*a2.size(),
             1.0, i0data_sorted, a1.size()*a2.size(), i1data_sorted, a1.size()*a2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1704::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1549
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());
    // tensor label: I1550
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", x3.size()*a1.size()*a2.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), a2.size(), x2.size());
  out()->put_block(odata, x2, x3, a1, a2);
}

void Task1705::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I1550
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x2.size(), c3.size());
  }
  out()->put_block(odata, x2, c3);
}

void Task1706::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
    // tensor label: I1791
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a1.size());
    dgemm_("T", "N", x1.size()*x2.size()*x3.size(), x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x1.size()*x2.size()*x3.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1707::Task_local::compute() {
  const Index x0 = b(0);
  const Index a1 = b(1);
  // tensor label: I1791
  std::unique_ptr<double[]> odata = out()->move_block(x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I1792
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", x0.size()*a1.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x0.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size());
  out()->put_block(odata, x0, a1);
}

void Task1708::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I1792
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2);
}

void Task1709::Task_local::compute() {
  const Index x0 = b(0);
  const Index a1 = b(1);
  // tensor label: I1791
  std::unique_ptr<double[]> odata = out()->move_block(x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I1796
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", x0.size()*a1.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x0.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size());
  out()->put_block(odata, x0, a1);
}

void Task1710::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I1796
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,1,2>(i0data, odata, a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2);
}

void Task1711::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), x1.size(), x2.size());
    // tensor label: I1857
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1712::Task_local::compute() {
  const Index a3 = b(0);
  const Index x3 = b(1);
  // tensor label: I1857
  std::unique_ptr<double[]> odata = out()->move_block(a3, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
      // tensor label: I1858
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x3.size());
      dgemm_("T", "N", 1, a3.size()*x3.size(), c2.size()*a1.size(),
             1.0, i0data_sorted, c2.size()*a1.size(), i1data_sorted, c2.size()*a1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a3.size(), x3.size());
  out()->put_block(odata, a3, x3);
}

void Task1713::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1858
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x3);
}

void Task1714::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
    // tensor label: I1861
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x3.size());
    dgemm_("T", "N", x0.size()*x1.size()*x2.size(), x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size()*x1.size()*x2.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1715::Task_local::compute() {
  const Index a1 = b(0);
  const Index x3 = b(1);
  // tensor label: I1861
  std::unique_ptr<double[]> odata = out()->move_block(a1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
      // tensor label: I1862
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x3)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x3.size());
      dgemm_("T", "N", 1, a1.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size());
  out()->put_block(odata, a1, x3);
}

void Task1716::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1862
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x3);
}

void Task1717::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a3.size());
      // tensor label: I1903
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a1, x3, x2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), x3.size()*x2.size(), a3.size()*a1.size(),
             1.0, i0data_sorted, a3.size()*a1.size(), i1data_sorted, a3.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x3.size(), x2.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1718::Task_local::compute() {
  const Index a3 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1903
  std::unique_ptr<double[]> odata = out()->move_block(a3, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a1, x3, x2), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I1904
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x3.size());
    dgemm_("T", "N", x2.size(), a3.size()*a1.size()*x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x2.size(), a3.size(), a1.size(), x3.size());
  out()->put_block(odata, a3, a1, x3, x2);
}

void Task1719::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1904
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x3);
}

void Task1720::Task_local::compute() {
  const Index a3 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1903
  std::unique_ptr<double[]> odata = out()->move_block(a3, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, a1, x3, x2), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I1919
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x2, a1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x2, a1, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), a1.size(), x3.size());
    dgemm_("T", "N", a3.size(), x2.size()*a1.size()*x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), x2.size(), a1.size(), x3.size());
  out()->put_block(odata, a3, a1, x3, x2);
}

void Task1721::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1919
  std::unique_ptr<double[]> odata = out()->move_block(a2, x2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x3.size(), a1.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, a2, x2, a1, x3);
}

void Task1722::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& a1 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      // tensor label: I1911
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a1, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a1, a2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a1.size(), a2.size());
      dgemm_("T", "N", x2.size()*x3.size(), x1.size()*x0.size(), a1.size()*a2.size(),
             1.0, i0data_sorted, a1.size()*a2.size(), i1data_sorted, a1.size()*a2.size(),
             1.0, odata_sorted, x2.size()*x3.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1723::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  const Index a2 = b(3);
  // tensor label: I1911
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c3, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c3.size(), a2.size());
    // tensor label: I1912
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), c3.size());
    dgemm_("T", "N", x0.size()*a1.size()*a2.size(), x1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x0.size()*a1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size(), a2.size(), x1.size());
  out()->put_block(odata, x1, x0, a1, a2);
}

void Task1724::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  // tensor label: I1912
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x1.size(), c3.size());
  }
  out()->put_block(odata, x1, c3);
}

void Task1725::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      // tensor label: I1955
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1, a1, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1726::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1955
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
    dscal_(a2.size()*x1.size()*a1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), x1.size(), a2.size());
  }
  out()->put_block(odata, a2, x1, a1, x0);
}

void Task1727::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: I1991
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x2, a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x2, a1, x3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), a1.size(), x3.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1728::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1991
  std::unique_ptr<double[]> odata = out()->move_block(a2, x2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
    dscal_(a2.size()*x2.size()*a1.size()*x3.size(), e0_, i0data.get(), 1);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x3.size(), a1.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, a2, x2, a1, x3);
}

void Task1729::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
      // tensor label: I2045
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1, a1, x0)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), x1.size()*x0.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1730::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2045
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x0.size(), a1.size(), x1.size(), a2.size());
  }
  out()->put_block(odata, a2, x1, a1, x0);
}

void Task1731::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: I2099
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x2, a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x2, a1, x3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), a1.size(), x3.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*a1.size(),
             1.0, i0data_sorted, a2.size()*a1.size(), i1data_sorted, a2.size()*a1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1732::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I2099
  std::unique_ptr<double[]> odata = out()->move_block(a2, x2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x3.size(), a1.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, a2, x2, a1, x3);
}

void Task1733::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());
    // tensor label: I2111
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, a1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), a1.size(), x0.size());
    dgemm_("T", "N", x3.size(), x2.size()*x1.size()*x0.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x3.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1734::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I2111
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x0.size(), a1.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x2, x1, a1, x0);
}

void Task1735::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1428
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0), 0.0);
  for (auto& a1 : *range_[2]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size());
    // tensor label: I2123
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x2, a1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x2, a1, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x2.size(), a1.size(), x3.size());
    dgemm_("T", "N", x0.size(), x1.size()*x2.size()*x3.size(), a1.size(),
           1.0, i0data_sorted, a1.size(), i1data_sorted, a1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task1736::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I2123
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x3.size(), a1.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x1, x2, a1, x3);
}

void Task1737::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  // tensor label: Gamma447
  std::unique_ptr<double[]> i0data = in(0)->get_block(ci0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0)]);
  sort_indices<0,0,1,1,1>(i0data, i0data_sorted, ci0.size());
  // tensor label: I1464
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);
  dgemm_("T", "N", ci0.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, ci0.size());
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1738::Task_local::compute() {
  // tensor label: I1464
  std::unique_ptr<double[]> odata = out()->move_block();
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size()]);
  std::fill_n(odata_sorted.get(), out()->get_size(), 0.0);
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I1465
  std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
  sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
  odata_sorted[0] += ddot_(a4.size()*c3.size()*a2.size()*c1.size(), i0data_sorted, 1, i1data_sorted, 1);
  sort_indices<1,1,1,1>(odata_sorted, odata);
  out()->put_block(odata);
}

void Task1739::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1465
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task1740::Task_local::compute() {
  // tensor label: I1464
  std::unique_ptr<double[]> odata = out()->move_block();
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size()]);
  std::fill_n(odata_sorted.get(), out()->get_size(), 0.0);
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I1468
  std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
  sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
  odata_sorted[0] += ddot_(a4.size()*c3.size()*a2.size()*c1.size(), i0data_sorted, 1, i1data_sorted, 1);
  sort_indices<1,1,1,1>(odata_sorted, odata);
  out()->put_block(odata);
}

void Task1741::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1468
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,2,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task1742::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma459
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x0)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x0.size());
      // tensor label: I1510
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size());
      dgemm_("T", "N", ci0.size(), 1, x0.size()*x3.size(),
             1.0, i0data_sorted, x0.size()*x3.size(), i1data_sorted, x0.size()*x3.size(),
             1.0, odata_sorted, ci0.size());
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1743::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I1510
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I1511
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x3.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size());
  out()->put_block(odata, x0, x3);
}

void Task1744::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1511
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1745::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I1510
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I1514
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x3.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size());
  out()->put_block(odata, x0, x3);
}

void Task1746::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1514
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task1747::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I1510
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I1873
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x3.size());
        dgemm_("T", "N", x0.size(), x3.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->put_block(odata, x0, x3);
}

void Task1748::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1873
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x3);
}

void Task1749::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I1510
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I1876
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x3)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x3.size());
        dgemm_("T", "N", x0.size(), x3.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x3.size());
  out()->put_block(odata, x0, x3);
}

