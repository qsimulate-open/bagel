//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks31.cc
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


#include <src/smith/CAS_test_tasks31.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task1500::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x2 = b(1);
  const Index x3 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1657
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x2, x3, x1, x0);
  {
    // tensor label: Gamma385
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x2, x3, x1, x0);
}

void Task1501::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
          // tensor label: I1619
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a2, c3, c1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a2, c3, c1)]);
          sort_indices<4,2,3,1,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a2.size(), c3.size(), c1.size());
          dgemm_("T", "N", 1, ci0.size(), x1.size()*a2.size()*c3.size()*c1.size(),
                 1.0, i0data_sorted, x1.size()*a2.size()*c3.size()*c1.size(), i1data_sorted, x1.size()*a2.size()*c3.size()*c1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1502::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1619
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a2, c3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a2, c3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a2, c3, c1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I1620
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, c4, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, c4, a2, c3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), c4.size(), a2.size(), c3.size());
    dgemm_("T", "N", c1.size(), ci0.size()*x1.size()*a2.size()*c3.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x1.size(), a2.size(), c3.size());
  out()->put_block(odata, ci0, x1, a2, c3, c1);
}

void Task1503::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c4 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1620
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c4, a2, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x0.size());
    // tensor label: I1621
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c4.size()*a2.size()*c3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c4.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c4.size(), a2.size(), c3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c4, a2, c3);
}

void Task1504::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1621
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1505::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c4 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1620
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c4, a2, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x0.size());
    // tensor label: I1625
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c4.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c4.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c4, a2, c3);
}

void Task1506::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1625
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1507::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1619
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a2, c3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a2, c3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a2, c3, c1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I1628
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, c4, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, c4, a2, c1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), c4.size(), a2.size(), c1.size());
    dgemm_("T", "N", c3.size(), ci0.size()*x1.size()*a2.size()*c1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), x1.size(), a2.size(), c1.size());
  out()->put_block(odata, ci0, x1, a2, c3, c1);
}

void Task1508::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c4 = b(2);
  const Index a2 = b(3);
  const Index c1 = b(4);
  // tensor label: I1628
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c4, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c4, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c4, a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x0.size());
    // tensor label: I1629
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c4.size()*a2.size()*c1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c4.size()*a2.size()*c1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c4.size(), a2.size(), c1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c4, a2, c1);
}

void Task1509::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1629
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1510::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c4 = b(2);
  const Index a2 = b(3);
  const Index c1 = b(4);
  // tensor label: I1628
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c4, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c4, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c4, a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x0.size());
    // tensor label: I1637
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c4.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c4.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c4, a2, c1);
}

void Task1511::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1637
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1512::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1619
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a2, c3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a2, c3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a2, c3, c1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I1632
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, c3, a4, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, c3, a4, c1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), c3.size(), a4.size(), c1.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x1.size()*c3.size()*c1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,0,3,4,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x1.size(), c3.size(), c1.size());
  out()->put_block(odata, ci0, x1, a2, c3, c1);
}

void Task1513::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  const Index c1 = b(4);
  // tensor label: I1632
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c3, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c3, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c3, a4, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x0.size());
    // tensor label: I1633
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a4.size()*c1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size()*a4.size()*c1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c3, a4, c1);
}

void Task1514::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1633
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1515::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  const Index c1 = b(4);
  // tensor label: I1632
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c3, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c3, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c3, a4, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x0.size());
    // tensor label: I1641
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c3, a4, c1);
}

void Task1516::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1641
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1517::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1619
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a2, c3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a2, c3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a2, c3, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x0)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), x0.size());
      // tensor label: I1660
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, c1, a4, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, c1, a4, c3, a2)]);
      sort_indices<1,4,0,2,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), c1.size(), a4.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*c1.size()*c3.size()*a2.size(), x0.size()*a4.size(),
             1.0, i0data_sorted, x0.size()*a4.size(), i1data_sorted, x0.size()*a4.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,4,3,2,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, ci0, x1, a2, c3, c1);
}

void Task1518::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c1 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);
  // tensor label: I1660
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I1661
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), ci0.size()*x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, c1, a4, c3, a2);
}

void Task1519::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1661
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1520::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index c1 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);
  // tensor label: I1660
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I1665
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), ci0.size()*x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, c1, a4, c3, a2);
}

void Task1521::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1665
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1522::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1619
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a2, c3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a2, c3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a2, c3, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x0.size());
    // tensor label: I1964
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a2, c3, c1);
}

void Task1523::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1964
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    dscal_(ci0.size()*x0.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1524::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1619
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a2, c3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a2, c3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a2, c3, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    // tensor label: I1967
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a2, c3, c1);
}

void Task1525::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1967
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    dscal_(ci0.size()*x0.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1526::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1619
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a2, c3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a2, c3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a2, c3, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size(), c1.size(), a2.size());
    // tensor label: I2057
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*c1.size()*a2.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size()*c1.size()*a2.size());
  }
  sort_indices<3,4,2,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a2, c3, c1);
}

void Task1527::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I2057
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1528::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1619
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a2, c3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a2, c3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a2, c3, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), c3.size(), a2.size());
    // tensor label: I2060
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a2, c3, c1);
}

void Task1529::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I2060
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1530::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& x4 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
          // tensor label: I1667
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, c2, a1)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), c2.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x5.size()*x4.size()*c2.size()*a1.size(),
                 1.0, i0data_sorted, x5.size()*x4.size()*c2.size()*a1.size(), i1data_sorted, x5.size()*x4.size()*c2.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1531::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1667
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size());
    // tensor label: I1668
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x3, x4, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x3, x4, c2)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x3.size(), x4.size(), c2.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x5.size()*x4.size()*c2.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x5.size(), x4.size(), c2.size());
  out()->put_block(odata, ci0, x5, x4, c2, a1);
}

void Task1532::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index c2 = b(4);
  // tensor label: I1668
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x3, x4, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x3, x4, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x3, x4, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c2, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c2, x2)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c2.size(), x2.size());
        // tensor label: I1669
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x3, x2, x4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x3, x2, x4, x1, x0)]);
        sort_indices<3,5,6,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
        dgemm_("T", "N", c2.size(), ci0.size()*x5.size()*x3.size()*x4.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x5.size(), x3.size(), x4.size());
  out()->put_block(odata, ci0, x5, x3, x4, c2);
}

void Task1533::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  const Index x4 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: I1669
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x3, x2, x4, x1, x0);
  {
    // tensor label: Gamma387
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x3, x2, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x3, x2, x4, x1, x0);
}

void Task1534::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1667
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I1680
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x1, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x1, x4)]);
      sort_indices<3,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x1.size(), x4.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x5.size()*x4.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x5.size(), x4.size());
  out()->put_block(odata, ci0, x5, x4, c2, a1);
}

void Task1535::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x4 = b(4);
  // tensor label: I1680
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x1, x4);
  {
    // tensor label: Gamma409
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x1, x4);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x1.size(), x4.size());
  }
  out()->put_block(odata, ci0, x5, x0, x1, x4);
}

void Task1536::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1667
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      // tensor label: I1691
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x5.size()*x4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x5.size(), x4.size());
  out()->put_block(odata, ci0, x5, x4, c2, a1);
}

void Task1537::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1691
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x0);
  {
    // tensor label: Gamma412
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x0);
}

void Task1538::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1667
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x3.size());
    // tensor label: I1702
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x3, x4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x3, x4, a1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x3.size(), x4.size(), a1.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x5.size()*x4.size()*a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x5.size(), x4.size(), a1.size());
  out()->put_block(odata, ci0, x5, x4, c2, a1);
}

void Task1539::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index a1 = b(4);
  // tensor label: I1702
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x3, x4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x3, x4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x3, x4, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        // tensor label: I1703
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x3, x4, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x3, x4, x2, x1)]);
        sort_indices<5,6,2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x5.size()*x3.size()*x4.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x5.size(), x3.size(), x4.size());
  out()->put_block(odata, ci0, x5, x3, x4, a1);
}

void Task1540::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x4 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1703
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x3, x4, x2, x1);
  {
    // tensor label: Gamma434
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x3, x4, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x3, x4, x2, x1);
}

void Task1541::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
          // tensor label: I1671
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, c2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, c2)]);
          sort_indices<1,3,4,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), c2.size());
          dgemm_("T", "N", 1, ci0.size(), x3.size()*x2.size()*a1.size()*c2.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*a1.size()*c2.size(), i1data_sorted, x3.size()*x2.size()*a1.size()*c2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1542::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c3.size());
      // tensor label: I1672
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x1, x2, c3, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x1, x2, c3, a1, c2)]);
      sort_indices<4,2,0,1,3,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x1.size(), x2.size(), c3.size(), a1.size(), c2.size());
      dgemm_("T", "N", 1, ci0.size()*x3.size()*x2.size()*a1.size()*c2.size(), x1.size()*c3.size(),
             1.0, i0data_sorted, x1.size()*c3.size(), i1data_sorted, x1.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x3.size(), x2.size(), a1.size(), c2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1543::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  const Index c3 = b(4);
  const Index a1 = b(5);
  const Index c2 = b(6);
  // tensor label: I1672
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x1, x2, c3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x1, x2, c3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x1, x2, c3, a1, c2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x0.size());
    // tensor label: I1673
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x1, x0, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x1, x0, x2)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
    dgemm_("T", "N", c3.size()*a1.size()*c2.size(), ci0.size()*x3.size()*x1.size()*x2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size()*a1.size()*c2.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a1.size(), c2.size(), ci0.size(), x3.size(), x1.size(), x2.size());
  out()->put_block(odata, ci0, x3, x1, x2, c3, a1, c2);
}

void Task1544::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  // tensor label: I1673
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x1, x0, x2);
  {
    // tensor label: Gamma400
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x1, x0, x2);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x1, x0, x2);
}

void Task1545::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  const Index c3 = b(4);
  const Index a1 = b(5);
  const Index c2 = b(6);
  // tensor label: I1672
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x1, x2, c3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x1, x2, c3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x1, x2, c3, a1, c2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x0.size());
    // tensor label: I1677
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x0, x1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
    dgemm_("T", "N", c2.size()*a1.size()*c3.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c2.size()*a1.size()*c3.size());
  }
  sort_indices<3,4,6,5,2,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c3.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x1, x2, c3, a1, c2);
}

void Task1546::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  // tensor label: I1677
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x0, x1);
  {
    // tensor label: Gamma390
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x2, x0, x1);
}

void Task1547::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c3.size());
    // tensor label: I1683
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, c3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), c3.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x3.size()*x2.size()*a1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x3.size(), x2.size(), a1.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1548::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c3 = b(4);
  // tensor label: I1683
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c3, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c3.size(), x1.size());
      // tensor label: I1684
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<3,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c3.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*c3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c3);
}

void Task1549::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1684
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

