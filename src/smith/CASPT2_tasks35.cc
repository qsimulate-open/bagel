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
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1888
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1701::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  // tensor label: I1887
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a4.size());
    // tensor label: I1892
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a4.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c2, a3);
}

void Task1702::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1892
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1703::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1864
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), a3.size());
    // tensor label: I1895
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, a4, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, a4, c2, a1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), a4.size(), c2.size(), a1.size());
    dgemm_("T", "N", a3.size(), ci0.size()*x1.size()*c2.size()*a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,2,3,0,4,1,1,1,1>(odata_sorted, odata, a3.size(), ci0.size(), x1.size(), c2.size(), a1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}

void Task1704::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1895
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c2.size(), a1.size());
    // tensor label: I1896
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c2.size(), a1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c2, a1);
}

void Task1705::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1896
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1706::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index a4 = b(2);
  const Index c2 = b(3);
  const Index a1 = b(4);
  // tensor label: I1895
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, a4, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, a4, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, a4, c2, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a4.size());
    // tensor label: I1900
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a4.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a4.size());
  }
  sort_indices<3,4,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a4.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, a4, c2, a1);
}

void Task1707::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1900
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1708::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1864
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1985
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}

void Task1709::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1985
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,1,4>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1710::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1864
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1988
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}

void Task1711::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I1988
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    dscal_(ci0.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,1,1,-1,2>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1712::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1864
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I2093
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,1,0,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}

void Task1713::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I2093
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,-1,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1714::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1864
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, c2, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x1, c2, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x1, c2, a3, a1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I2096
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x1, x0)]);
    sort_indices<2,0,1,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), ci0.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,1,2,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), ci0.size(), x1.size());
  out()->put_block(odata, ci0, x1, c2, a3, a1);
}

void Task1715::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  // tensor label: I2096
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x1, x0);
  {
    // tensor label: Gamma416
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x1, x0);
    sort_indices<0,1,2,1,1,2,1>(i0data, odata, ci0.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, ci0, x1, x0);
}

void Task1716::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x4 : *range_[1]) {
        for (auto& a2 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
          // tensor label: I1906
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, a1, a2)]);
          sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), a1.size(), a2.size());
          dgemm_("T", "N", 1, ci0.size(), x5.size()*x4.size()*a1.size()*a2.size(),
                 1.0, i0data_sorted, x5.size()*x4.size()*a1.size()*a2.size(), i1data_sorted, x5.size()*x4.size()*a1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1717::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1906
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, a1, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size());
    // tensor label: I1907
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x4, x3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x4, x3, a1)]);
    sort_indices<3,0,1,2,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x4.size(), x3.size(), a1.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x5.size()*x4.size()*a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x5.size(), x4.size(), a1.size());
  out()->put_block(odata, ci0, x5, x4, a1, a2);
}

void Task1718::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index a1 = b(4);
  // tensor label: I1907
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, x3, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, x2)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), x2.size());
        // tensor label: I1908
        std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x3, x2, x1)]);
        sort_indices<5,6,2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), ci0.size()*x5.size()*x4.size()*x3.size(), x0.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x0.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a1.size(), ci0.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, ci0, x5, x4, x3, a1);
}

void Task1719::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x3 = b(4);
  const Index x2 = b(5);
  const Index x1 = b(6);
  // tensor label: I1908
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x3, x2, x1);
  {
    // tensor label: Gamma437
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x3, x2, x1);
    sort_indices<0,1,2,3,4,5,6,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x3, x2, x1);
}

void Task1720::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1906
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x4, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x5, x4, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x5, x4, a1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: I1915
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x5, x0, x4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x5, x0, x4, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x5.size()*x4.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x5.size(), x4.size());
  out()->put_block(odata, ci0, x5, x4, a1, a2);
}

void Task1721::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x5 = b(1);
  const Index x0 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  // tensor label: I1915
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x5, x0, x4, x1);
  {
    // tensor label: Gamma470
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x5, x0, x4, x1);
    sort_indices<0,1,2,3,4,1,1,1,2>(i0data, odata, ci0.size(), x5.size(), x0.size(), x4.size(), x1.size());
  }
  out()->put_block(odata, ci0, x5, x0, x4, x1);
}

void Task1722::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I1196
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x2 : *range_[1]) {
        for (auto& a2 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, a2);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, a2)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), a2.size());
          // tensor label: I1910
          std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a2)]);
          sort_indices<1,3,2,4,0,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a2.size());
          dgemm_("T", "N", 1, ci0.size(), x3.size()*x2.size()*a1.size()*a2.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*a1.size()*a2.size(), i1data_sorted, x3.size()*x2.size()*a1.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task1723::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1910
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c3)]);
      sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), c3.size());
      // tensor label: I1911
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, x1, a1, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, x1, a1, c3, a2)]);
      sort_indices<5,3,0,1,2,4,6,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), x1.size(), a1.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, ci0.size()*x3.size()*x2.size()*a1.size()*a2.size(), x1.size()*c3.size(),
             1.0, i0data_sorted, x1.size()*c3.size(), i1data_sorted, x1.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,4,1,1,1,1>(odata_sorted, odata, ci0.size(), x3.size(), x2.size(), a1.size(), a2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a2);
}

void Task1724::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  const Index a1 = b(4);
  const Index c3 = b(5);
  const Index a2 = b(6);
  // tensor label: I1911
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, x1, a1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, x1, a1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, x1, a1, c3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c3.size(), a2.size());
    // tensor label: I1912
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
    sort_indices<2,0,1,3,4,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size()*c3.size()*a2.size(), ci0.size()*x3.size()*x2.size()*x1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a1.size()*c3.size()*a2.size());
  }
  sort_indices<3,4,5,6,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), a2.size(), ci0.size(), x3.size(), x2.size(), x1.size());
  out()->put_block(odata, ci0, x3, x2, x1, a1, c3, a2);
}

void Task1725::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1912
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1726::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1910
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a2.size());
    // tensor label: I1918
    std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x2, a1, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x2, a1, a3)]);
    sort_indices<4,0,1,2,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x2.size(), a1.size(), a3.size());
    dgemm_("T", "N", a2.size(), ci0.size()*x3.size()*x2.size()*a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,4,0,1,1,1,1>(odata_sorted, odata, a2.size(), ci0.size(), x3.size(), x2.size(), a1.size());
  out()->put_block(odata, ci0, x3, x2, a1, a2);
}

void Task1727::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a3 = b(4);
  // tensor label: I1918
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a3)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a3.size());
      // tensor label: I1919
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a3.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a3.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a3);
}

void Task1728::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1919
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1729::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1910
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: I1991
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a2);
}

void Task1730::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I1991
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    dscal_(ci0.size()*x3.size()*x0.size()*x2.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,1,1,-1,2>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

void Task1731::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  const Index a2 = b(4);
  // tensor label: I1910
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x2, a1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0, x3, x2, a1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0, x3, x2, a1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, x1, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), x1.size(), a2.size());
      // tensor label: I2099
      std::unique_ptr<double[]> i1data = in(1)->get_block(ci0, x3, x0, x2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(ci0, x3, x0, x2, x1)]);
      sort_indices<4,2,0,1,3,0,1,1,1>(i1data, i1data_sorted, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
      dgemm_("T", "N", a1.size()*a2.size(), ci0.size()*x3.size()*x2.size(), x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, a1.size()*a2.size());
    }
  }
  sort_indices<2,3,4,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), a2.size(), ci0.size(), x3.size(), x2.size());
  out()->put_block(odata, ci0, x3, x2, a1, a2);
}

void Task1732::Task_local::compute() {
  const Index ci0 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  const Index x1 = b(4);
  // tensor label: I2099
  std::unique_ptr<double[]> odata = out()->move_block(ci0, x3, x0, x2, x1);
  {
    // tensor label: Gamma438
    std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x0, x2, x1);
    sort_indices<0,1,2,3,4,1,1,1,1>(i0data, odata, ci0.size(), x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, ci0, x3, x0, x2, x1);
}

