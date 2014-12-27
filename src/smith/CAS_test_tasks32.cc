//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks32.cc
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


#include <src/smith/CAS_test_tasks32.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task1550::Task_local::compute() {
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
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), x0.size(), x1.size());
      // tensor label: I1695
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a1.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c3.size()*a1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a1.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c3);
}

void Task1551::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1695
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1552::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a1.size());
    // tensor label: I1687
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a3, c2)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a3.size(), c2.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x3.size()*x2.size()*c2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x3.size(), x2.size(), c2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1553::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  // tensor label: I1687
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a3, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), x1.size());
      // tensor label: I1688
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<3,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a3.size()*c2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a3, c2);
}

void Task1554::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1688
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1555::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  // tensor label: I1687
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a3, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x0.size(), x1.size());
      // tensor label: I1699
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a3, c2);
}

void Task1556::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1699
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1557::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x1)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a3.size(), x1.size());
      // tensor label: I1714
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, a3, c2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, a3, c2, a1)]);
      sort_indices<3,4,0,1,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), a3.size(), c2.size(), a1.size());
      dgemm_("T", "N", 1, ci0.size()*x3.size()*x2.size()*c2.size()*a1.size(), x1.size()*a3.size(),
             1.0, i0data_sorted, x1.size()*a3.size(), i1data_sorted, x1.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,4,3,1,1,1,1>(odata_sorted, odata, ci0.size(), x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1558::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  const Index a1 = b(6);
  // tensor label: I1714
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a3, c2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1715
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, a3, c2, a1);
}

void Task1559::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1715
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1560::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  const Index a1 = b(6);
  // tensor label: I1714
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a3, c2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1719
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x3.size()*x1.size()*x2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,6,5,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x3.size(), x1.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, x1, a3, c2, a1);
}

void Task1561::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1719
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1562::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I1970
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<3,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1563::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1970
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    dscal_(ci0.size()*x3.size()*x0.size()*x1.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1564::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      // tensor label: I1973
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1565::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1973
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    dscal_(ci0.size()*x3.size()*x2.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1566::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x0.size(), x1.size());
      // tensor label: I2063
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1567::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I2063
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1568::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x0, x1, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x0, x1, a1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), x0.size(), x1.size(), a1.size());
      // tensor label: I2066
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x1, x0, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x1, x0, x2)]);
      sort_indices<2,3,0,1,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1569::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  // tensor label: I2066
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x1, x0, x2);
  {
    // tensor label: Gamma400
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x1, x0, x2);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x1.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x1, x0, x2);
}

void Task1570::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
      // tensor label: I2069
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<3,2,0,1,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1571::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I2069
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1572::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1671
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, c2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c2, a1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c2.size(), a1.size());
      // tensor label: I2072
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x3.size()*x2.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, c2);
}

void Task1573::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I2072
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1574::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& x0 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, x0)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), x0.size());
          // tensor label: I1705
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a1)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x1.size()*x0.size()*c2.size()*a1.size(),
                 1.0, i0data_sorted, x1.size()*x0.size()*c2.size()*a1.size(), i1data_sorted, x1.size()*x0.size()*c2.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1575::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1705
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I1706
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a4, c3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a4, c3, a1)]);
      sort_indices<5,4,0,1,2,3,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a4.size(), c3.size(), a1.size());
      dgemm_("T", "N", 1, ci0.size()*x1.size()*x0.size()*c2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1576::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a1 = b(6);
  // tensor label: I1706
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a4, c3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a4, c3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a4, c3, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
  // tensor label: I1707
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a4.size()*c3.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a4.size()*c3.size()*a1.size());
  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a4.size(), c3.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a4, c3, a1);
}

void Task1577::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1707
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1578::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a4 = b(4);
  const Index c3 = b(5);
  const Index a1 = b(6);
  // tensor label: I1706
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a4, c3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a4, c3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a4, c3, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), a4.size());
  // tensor label: I1711
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size()*c3.size()*a4.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c3.size()*a4.size());
  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c3.size(), a4.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a4, c3, a1);
}

void Task1579::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1711
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1580::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1705
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a1), 0.0);
  // tensor label: h1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I2117
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a1);
}

void Task1581::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I2117
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1582::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& x5 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
          // tensor label: I1721
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, c1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, c1, a2)]);
          sort_indices<3,4,1,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), c1.size(), a2.size());
          dgemm_("T", "N", 1, ci0.size(), x5.size()*x4.size()*c1.size()*a2.size(),
                 1.0, i0data_sorted, x5.size()*x4.size()*c1.size()*a2.size(), i1data_sorted, x5.size()*x4.size()*c1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1583::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1721
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size());
    // tensor label: I1722
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, c1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), c1.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x5.size()*x4.size()*c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x5.size(), x4.size(), c1.size());
  out()->put_block(odata, ci0, x5, x4, c1, a2);
}

void Task1584::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index c1 = b(4);
  // tensor label: I1722
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
        // tensor label: I1723
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x2, x3, x1, x0)]);
        sort_indices<3,5,6,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), ci0.size()*x5.size()*x4.size()*x3.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, c1);
}

void Task1585::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x2 = b(3);
  const Index x3 = b(4);
  const Index x1 = b(5);
  const Index x0 = b(6);
  // tensor label: I1723
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma384
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,6,1,1,1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x2, x3, x1, x0);
}

void Task1586::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1721
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, x1)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      // tensor label: I1734
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), ci0.size()*x5.size()*x4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), ci0.size(), x5.size(), x4.size());
  out()->put_block(odata, ci0, x5, x4, c1, a2);
}

void Task1587::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1734
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x0);
  {
    // tensor label: Gamma412
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x0);
}

void Task1588::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1721
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x0, x1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x0.size(), x1.size());
      // tensor label: I1745
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x1, x0)]);
      sort_indices<3,4,0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), ci0.size()*x5.size()*x4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), ci0.size(), x5.size(), x4.size());
  out()->put_block(odata, ci0, x5, x4, c1, a2);
}

void Task1589::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1745
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x1, x0);
  {
    // tensor label: Gamma412
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x5, x4, x1, x0);
}

void Task1590::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1721
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, c1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size());
    // tensor label: I1756
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, a2)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), a2.size());
    dgemm_("T", "N", c1.size(), ci0.size()*x5.size()*x4.size()*a2.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x5.size(), x4.size(), a2.size());
  out()->put_block(odata, ci0, x5, x4, c1, a2);
}

void Task1591::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1756
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), x1.size(), x2.size());
        // tensor label: I1757
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, x0, x2, x1)]);
        sort_indices<5,6,4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a2.size(), ci0.size()*x5.size()*x4.size()*x3.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, a2);
}

void Task1592::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x0 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1757
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, x0, x2, x1);
  {
    // tensor label: Gamma435
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x4, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,-1,4>(i0data, odata, ci0.size(), x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x4, x3, x0, x2, x1);
}

void Task1593::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
          // tensor label: I1725
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a2, c1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a2, c1)]);
          sort_indices<4,3,1,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a2.size(), c1.size());
          dgemm_("T", "N", 1, ci0.size(), x3.size()*x2.size()*a2.size()*c1.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*a2.size()*c1.size(), i1data_sorted, x3.size()*x2.size()*a2.size()*c1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1594::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  const Index c1 = b(4);
  // tensor label: I1725
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a2, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c3.size());
      // tensor label: I1726
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, c3, a2, c1)]);
      sort_indices<4,3,0,1,2,5,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, ci0.size()*x3.size()*x2.size()*a2.size()*c1.size(), x1.size()*c3.size(),
             1.0, i0data_sorted, x1.size()*c3.size(), i1data_sorted, x1.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x3.size(), x2.size(), a2.size(), c1.size());
  out()->put_block(odata, ci0, x3, x2, a2, c1);
}

void Task1595::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  const Index c1 = b(6);
  // tensor label: I1726
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, c3, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, c3, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, c3, a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x0.size());
    // tensor label: I1727
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x0, x1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, c3, a2, c1);
}

void Task1596::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  // tensor label: I1727
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x0, x1);
  {
    // tensor label: Gamma390
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x2, x0, x1);
}

void Task1597::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  const Index c1 = b(6);
  // tensor label: I1726
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, c3, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, c3, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, c3, a2, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    // tensor label: I1731
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x0, x1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,5,6,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, c3, a2, c1);
}

void Task1598::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  const Index x1 = b(4);
  // tensor label: I1731
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x0, x1);
  {
    // tensor label: Gamma390
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x2, x0, x1);
}

void Task1599::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  const Index c1 = b(4);
  // tensor label: I1725
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a2, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c3.size());
    // tensor label: I1737
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a2, c3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a2.size(), c3.size());
    dgemm_("T", "N", c1.size(), ci0.size()*x3.size()*x2.size()*a2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, c1.size(), ci0.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, ci0, x3, x2, a2, c1);
}

