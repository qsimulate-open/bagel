//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks4.cc
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


#include <src/smith/CASPT2_tasks4.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task150::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I66
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I67
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", a2.size()*c1.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a2, c1);
}

void Task151::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I67
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task152::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I66
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I73
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a2, c1);
}

void Task153::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I73
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task154::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I31
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size());
      // tensor label: I75
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, a4, c3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, a4, c3, a2)]);
      sort_indices<3,1,0,2,4,5,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), a4.size(), c3.size(), a2.size());
      dgemm_("T", "N", 1, x0.size()*c1.size()*c3.size()*a2.size(), x1.size()*a4.size(),
             1.0, i0data_sorted, x1.size()*a4.size(), i1data_sorted, x1.size()*a4.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task155::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I75
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I76
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, c1, a4, c3, a2);
}

void Task156::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I76
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task157::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I75
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I79
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, c1, a4, c3, a2);
}

void Task158::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I79
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,4,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task159::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I31
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
    // tensor label: I281
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task160::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I281
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    dscal_(x0.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task161::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I31
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I283
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task162::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I283
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    dscal_(x0.size()*x1.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task163::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I31
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size(), c1.size(), a2.size());
    // tensor label: I309
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*c1.size()*a2.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task164::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I309
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,4,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task165::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I31
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size(), c3.size(), a2.size());
    // tensor label: I311
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*c3.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*c3.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task166::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I311
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task167::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  {
    // tensor label: I80
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0, c2, a1);
    sort_indices<2,0,1,3,1,1,1,1>(i0data, odata, x1.size(), x0.size(), c2.size(), a1.size());
  }
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task168::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: I81
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x2, x0, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x2, x0, c2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x2.size(), x0.size(), c2.size());
    dgemm_("T", "N", a1.size(), x1.size()*x0.size()*c2.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x1.size(), x0.size(), c2.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task169::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  // tensor label: I81
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x2, x0, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x2, x0, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        // tensor label: I82
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x1, x3, x2, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x1, x3, x2, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        dgemm_("T", "N", c2.size(), x1.size()*x2.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c2.size(), x1.size(), x2.size(), x0.size());
  out()->put_block(odata, x1, x2, x0, c2);
}

void Task170::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x3 = b(3);
  const Index x2 = b(4);
  const Index x0 = b(5);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x1, x3, x2, x0);
  {
    // tensor label: Gamma28
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x1, x3, x2, x0);
}

void Task171::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), c3.size());
      // tensor label: I84
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x2, x0, c3, a1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x2, x0, c3, a1, c2)]);
      sort_indices<1,3,0,2,4,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x2.size(), x0.size(), c3.size(), a1.size(), c2.size());
      dgemm_("T", "N", 1, x1.size()*x0.size()*a1.size()*c2.size(), x2.size()*c3.size(),
             1.0, i0data_sorted, x2.size()*c3.size(), i1data_sorted, x2.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task172::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index c3 = b(3);
  const Index a1 = b(4);
  const Index c2 = b(5);
  // tensor label: I84
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, c3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x2, x0, c3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x2, x0, c3, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x3.size());
    // tensor label: I85
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x3, x2, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x3, x2, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
    dgemm_("T", "N", c3.size()*a1.size()*c2.size(), x1.size()*x2.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c3.size()*a1.size()*c2.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a1.size(), c2.size(), x1.size(), x2.size(), x0.size());
  out()->put_block(odata, x1, x2, x0, c3, a1, c2);
}

void Task173::Task_local::compute() {
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I85
  std::unique_ptr<double[]> odata = out()->move_block(x1, x3, x2, x0);
  {
    // tensor label: Gamma29
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), x3.size(), x2.size(), x0.size());
  }
  out()->put_block(odata, x1, x3, x2, x0);
}

void Task174::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index c3 = b(3);
  const Index a1 = b(4);
  const Index c2 = b(5);
  // tensor label: I84
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, c3, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x2, x0, c3, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x2, x0, c3, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x3.size());
    // tensor label: I88
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c2.size()*a1.size()*c3.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c2.size()*a1.size()*c3.size());
  }
  sort_indices<4,3,5,2,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x2, x0, c3, a1, c2);
}

void Task175::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I88
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task176::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
      // tensor label: I90
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x1, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x1.size(), x4.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task177::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x4 = b(3);
  // tensor label: I90
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x1, x4);
  {
    // tensor label: Gamma31
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x1, x4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x1.size(), x4.size());
  }
  out()->put_block(odata, x5, x0, x1, x4);
}

void Task178::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), c3.size());
    // tensor label: I92
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, c3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), c3.size());
    dgemm_("T", "N", c2.size(), x0.size()*x1.size()*a1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task179::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index c3 = b(3);
  // tensor label: I92
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), x2.size());
      // tensor label: I93
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c3.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c3.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, c3);
}

void Task180::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I93
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x1, x2);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x3, x0, x1, x2);
}

void Task181::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index c3 = b(3);
  // tensor label: I92
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I101
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c3.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c3.size()*a1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x0, x1, a1, c3);
}

void Task182::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I101
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task183::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), a1.size());
    // tensor label: I95
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a3, c2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a3.size(), c2.size());
    dgemm_("T", "N", a1.size(), x0.size()*x1.size()*c2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x1.size(), c2.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task184::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I95
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a3, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      // tensor label: I96
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a3.size()*c2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a3, c2);
}

void Task185::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I96
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x1, x2);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x3, x0, x1, x2);
}

void Task186::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I95
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a3, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I104
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x0, x1, a3, c2);
}

void Task187::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I104
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task188::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
      // tensor label: I98
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task189::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I98
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x1, x0);
  {
    // tensor label: Gamma34
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x1, x0);
}

void Task190::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I106
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), a1.size());
    dgemm_("T", "N", c2.size(), x0.size()*x1.size()*a1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), x1.size(), a1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task191::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index a1 = b(3);
  // tensor label: I106
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, x2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x2, a1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), x3.size());
        // tensor label: I107
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x0, x4, x3, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x0, x4, x3, x1, x2)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
        dgemm_("T", "N", a1.size(), x0.size()*x1.size()*x2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x1.size(), x2.size());
  out()->put_block(odata, x0, x1, x2, a1);
}

void Task192::Task_local::compute() {
  const Index x5 = b(0);
  const Index x0 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x2 = b(5);
  // tensor label: I107
  std::unique_ptr<double[]> odata = out()->move_block(x5, x0, x4, x3, x1, x2);
  {
    // tensor label: Gamma37
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,1>(i0data, odata, x5.size(), x0.size(), x4.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x5, x0, x4, x3, x1, x2);
}

void Task193::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), c3.size());
      // tensor label: I109
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a4, c3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a4, c3, a1)]);
      sort_indices<3,4,0,1,2,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a4.size(), c3.size(), a1.size());
      dgemm_("T", "N", 1, x1.size()*x0.size()*c2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task194::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a1 = b(5);
  // tensor label: I109
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a4, c3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a4, c3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a4, c3, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
  // tensor label: I110
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a4.size()*c3.size()*a1.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a4.size()*c3.size()*a1.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a4.size(), c3.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a4, c3, a1);
}

void Task195::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I110
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task196::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a1 = b(5);
  // tensor label: I109
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a4, c3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a4, c3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a4, c3, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), a4.size());
  // tensor label: I113
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size()*c3.size()*a4.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c3.size()*a4.size());
  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c3.size(), a4.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a4, c3, a1);
}

void Task197::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I113
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task198::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
      // tensor label: I115
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a3, c2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a3, c2, a1)]);
      sort_indices<3,0,1,2,4,5,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a3.size(), c2.size(), a1.size());
      dgemm_("T", "N", 1, x1.size()*x0.size()*c2.size()*a1.size(), x2.size()*a3.size(),
             1.0, i0data_sorted, x2.size()*a3.size(), i1data_sorted, x2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, x1, x0, c2, a1);
}

void Task199::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  const Index a1 = b(5);
  // tensor label: I115
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I116
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, a3, c2, a1);
}

