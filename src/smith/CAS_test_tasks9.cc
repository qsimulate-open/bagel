//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks9.cc
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


#include <src/smith/CAS_test_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task400::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I354
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x1, x2);
  {
    // tensor label: Gamma94
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x0, x3, x1, x2);
}

void Task401::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c2 = b(3);
  // tensor label: I349
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, c2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c2.size(), x2.size());
    // tensor label: I357
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, x2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, x2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), x2.size(), c1.size());
    dgemm_("T", "N", c2.size(), x0.size()*x1.size()*c1.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c2.size(), x0.size(), x1.size(), c1.size());
  out()->put_block(odata, x0, x1, c1, c2);
}

void Task402::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  // tensor label: I357
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, x2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, x2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, x2, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        // tensor label: I358
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x0, x3, x1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x0, x3, x1, x2)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x0.size(), x3.size(), x1.size(), x2.size());
        dgemm_("T", "N", c1.size(), x0.size()*x1.size()*x2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size(), x1.size(), x2.size());
  out()->put_block(odata, x0, x1, x2, c1);
}

void Task403::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x0 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x2 = b(5);
  // tensor label: I358
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x0, x3, x1, x2);
  {
    // tensor label: Gamma2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x0, x3, x1, x2);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x4.size(), x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x5, x4, x0, x3, x1, x2);
}

void Task404::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c2 = b(3);
  // tensor label: I349
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, c2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size());
      // tensor label: I361
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, x2, c1, a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, x2, c1, a3, c2)]);
      sort_indices<4,2,0,1,3,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), x2.size(), c1.size(), a3.size(), c2.size());
      dgemm_("T", "N", 1, x1.size()*x0.size()*c1.size()*c2.size(), x2.size()*a3.size(),
             1.0, i0data_sorted, x2.size()*a3.size(), i1data_sorted, x2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), c2.size());
  out()->put_block(odata, x0, x1, c1, c2);
}

void Task405::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index c1 = b(3);
  const Index a3 = b(4);
  const Index c2 = b(5);
  // tensor label: I361
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, x2, c1, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, x2, c1, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, x2, c1, a3, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
    // tensor label: I362
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x3, x0, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x3, x0, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
    dgemm_("T", "N", c1.size()*a3.size()*c2.size(), x1.size()*x0.size()*x2.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a3.size()*c2.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a3.size(), c2.size(), x1.size(), x0.size(), x2.size());
  out()->put_block(odata, x1, x0, x2, c1, a3, c2);
}

void Task406::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I362
  std::unique_ptr<double[]> odata = out()->move_block(x1, x3, x0, x2);
  {
    // tensor label: Gamma3
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x1.size(), x3.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, x1, x3, x0, x2);
}

void Task407::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c2 = b(3);
  // tensor label: I349
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c2, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      // tensor label: I724
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

void Task408::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I724
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x1, x2);
  {
    // tensor label: Gamma94
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
    dscal_(x0.size()*x3.size()*x1.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x0, x3, x1, x2);
}

void Task409::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index c2 = b(3);
  // tensor label: I349
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c2, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c2.size(), x2.size());
      // tensor label: I764
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

void Task410::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I764
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x1, x2);
  {
    // tensor label: Gamma94
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x1, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), x3.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x0, x3, x1, x2);
}

void Task411::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1, c1, x2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1, c1, x2)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size(), c1.size(), x2.size());
  // tensor label: I364
  std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, c1)]);
  sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), c1.size());
  energy_ += ddot_(x2.size()*x1.size()*x0.size()*c1.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task412::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, c2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), c2.size());
      // tensor label: I365
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0, c1, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0, c1, c2)]);
      sort_indices<1,5,0,2,3,4,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size(), c1.size(), c2.size());
      dgemm_("T", "N", 1, x2.size()*x1.size()*x0.size()*c1.size(), x3.size()*c2.size(),
             1.0, i0data_sorted, x3.size()*c2.size(), i1data_sorted, x3.size()*c2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task413::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index c1 = b(4);
  const Index c2 = b(5);
  // tensor label: I365
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0, c1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x3, x1, x0, c1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x3, x1, x0, c1, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, c2, x4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), c2.size(), x4.size());
      // tensor label: I366
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x5, x3, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x5, x3, x4, x1, x0)]);
      sort_indices<1,3,0,2,4,5,0,1,1,1>(i1data, i1data_sorted, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*c2.size(), x2.size()*x3.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*c2.size());
    }
  }
  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c2.size(), x2.size(), x3.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x3, x1, x0, c1, c2);
}

void Task414::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x5 = b(1);
  const Index x3 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I366
  std::unique_ptr<double[]> odata = out()->move_block(x2, x5, x3, x4, x1, x0);
  {
    // tensor label: Gamma4
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x3, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x2.size(), x5.size(), x3.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x5, x3, x4, x1, x0);
}

void Task415::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c1, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, c1, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), c1.size(), x5.size());
        // tensor label: I369
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, x2, x5, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, x2, x5, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), x2.size(), x5.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task416::Task_local::compute() {
  energy_ = 0.0;
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I369
  std::unique_ptr<double[]> odata = out()->move_block(x7, x6, x2, x5, x1, x0);
  {
    // tensor label: Gamma5
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x2, x5, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x7.size(), x6.size(), x2.size(), x5.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x7, x6, x2, x5, x1, x0);
}

void Task417::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), c2.size());
    // tensor label: I372
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, c2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), c2.size());
    dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task418::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  // tensor label: I372
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        // tensor label: I373
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c2.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c2.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c2);
}

void Task419::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I373
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task420::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a3.size(), c2.size());
      // tensor label: I376
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, c2, a3, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, c2, a3, c1)]);
      sort_indices<4,3,0,1,2,5,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), c2.size(), a3.size(), c1.size());
      dgemm_("T", "N", 1, x2.size()*x1.size()*x0.size()*c1.size(), c2.size()*a3.size(),
             1.0, i0data_sorted, c2.size()*a3.size(), i1data_sorted, c2.size()*a3.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task421::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c1 = b(5);
  // tensor label: I376
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c2, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c2, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c2, a3, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c1.size(), x3.size());
    // tensor label: I377
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c2.size()*a3.size()*c1.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c2.size()*a3.size()*c1.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c2, a3, c1);
}

void Task422::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I377
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task423::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c2 = b(3);
  const Index a3 = b(4);
  const Index c1 = b(5);
  // tensor label: I376
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c2, a3, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c2, a3, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c2, a3, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c2, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c2.size(), x3.size());
    // tensor label: I381
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size()*a3.size()*c2.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a3.size()*c2.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a3.size(), c2.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c2, a3, c1);
}

void Task424::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I381
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task425::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a2, x3)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a2.size(), x3.size());
      // tensor label: I384
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0, a2, c1)]);
      sort_indices<4,0,1,2,3,5,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, x2.size()*x1.size()*x0.size()*c1.size(), x3.size()*a2.size(),
             1.0, i0data_sorted, x3.size()*a2.size(), i1data_sorted, x3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task426::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a2 = b(4);
  const Index c1 = b(5);
  // tensor label: I384
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0, a2, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      // tensor label: I385
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x3, x2, x4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x3, x2, x4, x1, x0)]);
      sort_indices<0,3,1,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), x3.size()*x2.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,4,5,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0, a2, c1);
}

void Task427::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I385
  std::unique_ptr<double[]> odata = out()->move_block(x5, x3, x2, x4, x1, x0);
  {
    // tensor label: Gamma9
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x2, x4, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x3.size(), x2.size(), x4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x3, x2, x4, x1, x0);
}

void Task428::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  const Index a2 = b(4);
  const Index c1 = b(5);
  // tensor label: I384
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, x1, x0, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, x1, x0, a2, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x5, x4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I389
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x2.size()*x3.size()*x1.size()*x0.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<3,2,4,5,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x2.size(), x3.size(), x1.size(), x0.size());
  out()->put_block(odata, x3, x2, x1, x0, a2, c1);
}

void Task429::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I389
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,4>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task430::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        // tensor label: I727
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task431::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I727
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
    dscal_(x5.size()*x4.size()*x2.size()*x3.size()*x1.size()*x0.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,4,5,1,1,-1,4>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task432::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x5, x4, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x5.size(), x4.size(), x3.size());
        // tensor label: I767
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x5, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task433::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I767
  std::unique_ptr<double[]> odata = out()->move_block(x2, x5, x4, x3, x1, x0);
  {
    // tensor label: Gamma107
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x5, x4, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x2.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x5, x4, x3, x1, x0);
}

void Task434::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, c1, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), c1.size(), x3.size());
        // tensor label: I770
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x2, x3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x2, x3, x1, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task435::Task_local::compute() {
  energy_ = 0.0;
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I770
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x2, x3, x1, x0);
  {
    // tensor label: Gamma6
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
    sort_indices<0,1,2,3,4,5,1,1,1,2>(i0data, odata, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x5, x4, x2, x3, x1, x0);
}

void Task436::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index c1 = b(3);
  // tensor label: I364
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1, x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1, x0, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size());
    // tensor label: I822
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size(), x2.size()*x1.size()*x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c1.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x2, x1, x0, c1);
}

void Task437::Task_local::compute() {
  energy_ = 0.0;
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I822
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, x1, x0);
  {
    // tensor label: Gamma7
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), x3.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, x3, x1, x0);
}

void Task438::Task_local::compute() {
  energy_ = 0.0;
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
  // tensor label: I391
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c1, c3, a2)]);
  sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), c1.size(), c3.size(), a2.size());
  energy_ += ddot_(x0.size()*c1.size()*c3.size()*a2.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task439::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size());
    // tensor label: I392
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c1, c3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c1.size(), c3.size());
    dgemm_("T", "N", a2.size(), x0.size()*c1.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), x0.size(), c1.size(), c3.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task440::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I392
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c1, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x3, c3, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x3.size(), c3.size(), x2.size());
      // tensor label: I393
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x3, x0, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x3, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), x3.size(), x0.size(), x2.size());
      dgemm_("T", "N", c1.size()*c3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*c3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), c3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c1, c3);
}

void Task441::Task_local::compute() {
  energy_ = 0.0;
  const Index x1 = b(0);
  const Index x3 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I393
  std::unique_ptr<double[]> odata = out()->move_block(x1, x3, x0, x2);
  {
    // tensor label: Gamma3
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x0, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x1.size(), x3.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, x1, x3, x0, x2);
}

void Task442::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
  // tensor label: I396
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size());
  dgemm_("T", "N", c1.size()*a2.size(), x0.size()*c3.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size());
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x0.size(), c3.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task443::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c3 = b(1);
  // tensor label: I396
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c3, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
        // tensor label: I397
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        dgemm_("T", "N", c3.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x0.size());
  out()->put_block(odata, x0, c3);
}

void Task444::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I397
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task445::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I400
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), c1.size());
  dgemm_("T", "N", c3.size()*a2.size(), x0.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x0.size(), c1.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task446::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  // tensor label: I400
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        // tensor label: I401
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        dgemm_("T", "N", c1.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x0.size());
  out()->put_block(odata, x0, c1);
}

void Task447::Task_local::compute() {
  energy_ = 0.0;
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I401
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task448::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index c1 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I391
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1, c3, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    // tensor label: I404
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, x0, c1, c3, a2);
}

void Task449::Task_local::compute() {
  energy_ = 0.0;
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I404
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  {
    // tensor label: Gamma14
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x0.size(), x3.size());
  }
  out()->put_block(odata, x0, x3);
}

