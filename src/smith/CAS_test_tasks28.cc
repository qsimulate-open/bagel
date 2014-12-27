//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks28.cc
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


#include <src/smith/CAS_test_tasks28.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task1350::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1471
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I1476
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c1, a2);
}

void Task1351::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1476
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1352::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size());
    // tensor label: I1803
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, c1, a4, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, c1, a4, c3)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), c1.size(), a4.size(), c3.size());
    dgemm_("T", "N", a2.size(), ci0.size()*c1.size()*a4.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), c1.size(), a4.size(), c3.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1353::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  // tensor label: I1803
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x0.size());
    // tensor label: I1804
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c1, a4, c3);
}

void Task1354::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1804
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1355::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size());
    // tensor label: I1807
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, c1, a2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, c1, a2, c3)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), c1.size(), a2.size(), c3.size());
    dgemm_("T", "N", a4.size(), ci0.size()*c1.size()*a2.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,0,4,3,1,1,1,1>(odata_sorted, odata, a4.size(), ci0.size(), c1.size(), a2.size(), c3.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1356::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  const Index c3 = b(4);
  // tensor label: I1807
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c1, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c1, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c1, a2, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
    // tensor label: I1808
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c1, a2, c3);
}

void Task1357::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  // tensor label: I1808
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1);
  {
    // tensor label: Gamma394
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
    sort_indices<0,1,2,1,1,1,1>(i0data, odata, ci0.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, ci0, x0, x1);
}

void Task1358::Task_local::compute() {
  const Index ci0 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  // tensor label: I1440
  std::unique_ptr<double[]> odata = out()->move_block(ci0, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, c1, a4, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I1833
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c1, a2)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c1.size(), a2.size());
    dgemm_("T", "N", c3.size(), ci0.size()*a4.size()*c1.size()*a2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,3,2,0,4,1,1,1,1>(odata_sorted, odata, c3.size(), ci0.size(), a4.size(), c1.size(), a2.size());
  out()->put_block(odata, ci0, c1, a4, c3, a2);
}

void Task1359::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1833
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c1, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I1834
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c1, a2);
}

void Task1360::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1834
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1361::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1833
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c1, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I1838
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c1, a2);
}

void Task1362::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1838
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1363::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
          // tensor label: I1478
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a3, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a3, c2, a1)]);
          sort_indices<1,4,3,2,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a3.size(), c2.size(), a1.size());
          dgemm_("T", "N", 1, ci0.size(), x0.size()*a3.size()*c2.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*a3.size()*c2.size()*a1.size(), i1data_sorted, x0.size()*a3.size()*c2.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1364::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size());
    // tensor label: I1479
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, a3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, a3, c2)]);
    sort_indices<1,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), a3.size(), c2.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*a3.size()*c2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), a3.size(), c2.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1365::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  // tensor label: I1479
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, a3, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      // tensor label: I1480
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, a3, c2);
}

void Task1366::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1480
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1367::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  // tensor label: I1479
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, a3, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I1488
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<2,3,4,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, a3, c2);
}

void Task1368::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1488
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1369::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: I1483
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, x1, a1, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, x1, a1, c2)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), x1.size(), a1.size(), c2.size());
    dgemm_("T", "N", a3.size(), ci0.size()*x0.size()*a1.size()*c2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,2,0,4,3,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x0.size(), a1.size(), c2.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1370::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1483
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I1484
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x1, x2)]);
      sort_indices<4,1,0,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), ci0.size()*x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), ci0.size(), x0.size(), x1.size());
  out()->put_block(odata, ci0, x0, x1, a1, c2);
}

void Task1371::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  const Index x2 = b(4);
  // tensor label: I1484
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x1, x2);
  {
    // tensor label: Gamma410
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x1, x2);
    sort_indices<0,1,2,3,4,1,1,1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, ci0, x3, x0, x1, x2);
}

void Task1372::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index a1 = b(3);
  const Index c2 = b(4);
  // tensor label: I1483
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, x1, a1, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I1492
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, x0)]);
      sort_indices<2,1,0,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,4,3,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x0, x1, a1, c2);
}

void Task1373::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index x0 = b(4);
  // tensor label: I1492
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, x0);
  {
    // tensor label: Gamma413
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x1, x0);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x2, x1, x0);
}

void Task1374::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I1495
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a3);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a3)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a3.size());
  dgemm_("T", "N", c2.size()*a1.size(), ci0.size()*x0.size()*a3.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), ci0.size(), x0.size(), a3.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1375::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  // tensor label: I1495
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, x2, x1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        // tensor label: I1496
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
        sort_indices<4,3,1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a3.size(), ci0.size()*x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,2,0,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3);
}

void Task1376::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1496
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1377::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
  // tensor label: I1499
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a1)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a1.size());
  dgemm_("T", "N", c2.size()*a3.size(), ci0.size()*x0.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size());
  sort_indices<2,3,1,0,4,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), ci0.size(), x0.size(), a1.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1378::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a1 = b(2);
  // tensor label: I1499
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, x1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        // tensor label: I1500
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
        sort_indices<4,3,1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a1);
}

void Task1379::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1500
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1380::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());
      // tensor label: I1503
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0, c2, a3, c4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
      sort_indices<5,1,0,2,3,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());
      dgemm_("T", "N", 1, ci0.size()*x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
             1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,3,2,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1381::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c4 = b(5);
  const Index a1 = b(6);
  // tensor label: I1503
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
  // tensor label: I1504
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());
  sort_indices<4,5,6,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
}

void Task1382::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1504
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1383::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c4 = b(5);
  const Index a1 = b(6);
  // tensor label: I1503
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, x0, c2, a3, c4, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
  // tensor label: I1508
  std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
  sort_indices<0,1,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), ci0.size()*x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());
  sort_indices<4,5,6,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), ci0.size(), x1.size(), x0.size());
  out()->put_block(odata, ci0, x1, x0, c2, a3, c4, a1);
}

void Task1384::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1508
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1385::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1511
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1386::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  // tensor label: I1511
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  {
    // tensor label: Gamma459
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x0);
}

void Task1387::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1514
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1388::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  // tensor label: I1514
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0);
  {
    // tensor label: Gamma459
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x3.size(), x0.size());
  }
  out()->put_block(odata, ci0, x3, x0);
}

void Task1389::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c4.size());
    // tensor label: I1517
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a3, c4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a3, c4, a1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a3.size(), c4.size(), a1.size());
    dgemm_("T", "N", c2.size(), ci0.size()*x0.size()*a3.size()*a1.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, c2.size(), ci0.size(), x0.size(), a3.size(), a1.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1390::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c4 = b(3);
  const Index a1 = b(4);
  // tensor label: I1517
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c4, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I1518
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c4.size()*a1.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c4.size()*a1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c4.size(), a1.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c4, a1);
}

void Task1391::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1518
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1392::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c4 = b(3);
  const Index a1 = b(4);
  // tensor label: I1517
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c4, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I1522
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c4.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size()*c4.size()*a3.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c4.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a3, c4, a1);
}

void Task1393::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1522
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1394::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a1.size());
    // tensor label: I1525
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c2, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c2, a3)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c2.size(), a3.size());
    dgemm_("T", "N", a1.size(), ci0.size()*x0.size()*c2.size()*a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,4,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x0.size(), c2.size(), a3.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

void Task1395::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  // tensor label: I1525
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
    // tensor label: I1526
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c2.size()*a3.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a3.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c2, a3);
}

void Task1396::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1526
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1397::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  // tensor label: I1525
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a4, c2, a3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());
    // tensor label: I1530
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<1,0,2,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a4.size(), ci0.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), ci0.size(), x0.size());
  out()->put_block(odata, ci0, x0, a4, c2, a3);
}

void Task1398::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1530
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1399::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1478
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x0, a3, c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());
    // tensor label: I1533
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x0, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x0, a4, c2, a1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x0.size(), a4.size(), c2.size(), a1.size());
    dgemm_("T", "N", a3.size(), ci0.size()*x0.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,2,0,3,4,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, ci0, x0, a3, c2, a1);
}

