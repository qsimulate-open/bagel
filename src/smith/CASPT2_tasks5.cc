//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks5.cc
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


#include <src/smith/CASPT2_tasks5.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task200::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I132
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, a2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, a2, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), a2.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task201::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    dscal_(x3.size()*a2.size()*c1.size()*x2.size(), e0_, i0data.get(), 1);
    sort_indices<2,0,1,3,1,1,1,1>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x3, x2);
    dscal_(c1.size()*a2.size()*x3.size()*x2.size(), e0_, i1data.get(), 1);
    sort_indices<0,2,1,3,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(1)->get_block(c1, a2, x3, x2);
    sort_indices<0,2,1,3,1,1,2,1>(i2data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i3data = in(1)->get_block(x3, a2, c1, x2);
    sort_indices<2,0,1,3,1,1,-1,1>(i3data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i4data = in(1)->get_block(x3, x2, c1, a2);
    sort_indices<2,0,3,1,1,1,2,1>(i4data, odata, x3.size(), x2.size(), c1.size(), a2.size());
  }
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task202::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
    // tensor label: I133
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size());
    dgemm_("T", "N", x3.size()*a2.size()*x2.size(), c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*a2.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), x2.size(), c1.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task203::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  // tensor label: I133
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c1.size(), c3.size());
  }
  out()->put_block(odata, c1, c3);
}

void Task204::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
    // tensor label: I136
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size());
    dgemm_("T", "N", x3.size()*c1.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*c1.size()*x2.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x2.size(), a2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task205::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  // tensor label: I136
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2);
}

void Task206::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
    // tensor label: I141
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size());
    dgemm_("T", "N", a2.size()*x3.size()*x2.size(), c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a2.size()*x3.size()*x2.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a2.size(), x3.size(), x2.size(), c1.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task207::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  // tensor label: I141
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c1.size(), c3.size());
  }
  out()->put_block(odata, c1, c3);
}

void Task208::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x3, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
    // tensor label: I144
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), a2.size());
    dgemm_("T", "N", c1.size()*x3.size()*x2.size(), a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c1.size()*x3.size()*x2.size());
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task209::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2);
}

void Task210::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    // tensor label: I156
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*c1.size()*a2.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*c1.size()*a2.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), a2.size(), x2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task211::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I156
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task212::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I132
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a3.size());
    // tensor label: I159
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*a2.size()*c1.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a2.size()*c1.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), c1.size(), x2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task213::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I159
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task214::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I146
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x2, x1, x0)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), x2.size(), x1.size(), x0.size());
    dgemm_("T", "N", c1.size(), a2.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, c1.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task215::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I146
  std::unique_ptr<double[]> odata = out()->move_block(a2, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma51
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I147
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a2, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), a2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a2.size());
  out()->put_block(odata, a2, x2, x1, x0);
}

void Task216::Task_local::compute() {
  const Index x5 = b(0);
  const Index a2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I147
  std::unique_ptr<double[]> odata = out()->move_block(x5, a2, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), a2.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a2, x4, x3);
}

void Task217::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I120
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  // tensor label: Gamma38
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I149
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
  dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task218::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I149
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
    sort_indices<0,1,1,1,4,1>(i0data, odata, c1.size(), a2.size());
  }
  out()->put_block(odata, c1, a2);
}

void Task219::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I149
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I150
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, c1, a2);
}

void Task220::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I150
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task221::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  // tensor label: I149
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I153
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, c1, a2);
}

void Task222::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  // tensor label: I153
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3);
    sort_indices<0,1,1,1,8,1>(i0data, odata, a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3);
}

void Task223::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, a1);
  {
    // tensor label: I160
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, x0, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x1, x2, x0, a1);
}

void Task224::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I160
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma56
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x3, x4, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I161
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x5, a1, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x5, a1, x4)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x5.size(), a1.size(), x4.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, x3.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task225::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I161
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
    // tensor label: I162
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, x5, a1, x4);
}

void Task226::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x3.size(), c2.size());
  }
  out()->put_block(odata, x3, c2);
}

void Task227::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I160
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma57
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I164
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x5, x4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x5, x4)]);
        sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x5.size(), x4.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, x3.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task228::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I164
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x5, x4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, a1);
    sort_indices<2,3,0,1,1,1,1,1>(i0data, odata, x5.size(), x4.size(), x3.size(), a1.size());
  }
  out()->put_block(odata, x3, a1, x5, x4);
}

void Task229::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I164
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
    // tensor label: I165
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c2)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c2.size());
    dgemm_("T", "N", a1.size()*x5.size()*x4.size(), x3.size(), c2.size(),
           1.0, i0data_sorted, c2.size(), i1data_sorted, c2.size(),
           1.0, odata_sorted, a1.size()*x5.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, a1, x5, x4);
}

void Task230::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I165
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x3.size(), c2.size());
  }
  out()->put_block(odata, x3, c2);
}

void Task231::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I160
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma58
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x6, x5, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x6, x5, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x6.size(), x5.size(), x2.size(), x1.size());
        // tensor label: I167
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a1, x6, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a1, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x7.size(), a1.size(), x6.size(), x5.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task232::Task_local::compute() {
  const Index x7 = b(0);
  const Index a1 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  // tensor label: I167
  std::unique_ptr<double[]> odata = out()->move_block(x7, a1, x6, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, x6, x5);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x7.size(), a1.size(), x6.size(), x5.size());
  }
  out()->put_block(odata, x7, a1, x6, x5);
}

void Task233::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I160
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma59
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I169
        std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x5, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x5, x4, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, a1.size(), x5.size(), x4.size(), x3.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task234::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I169
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, x3);
    dscal_(x5.size()*a1.size()*x4.size()*x3.size(), e0_, i0data.get(), 1);
    sort_indices<1,0,2,3,1,1,-1,1>(i0data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a1, x4, x3);
    sort_indices<1,0,2,3,1,1,1,1>(i1data, odata, x5.size(), a1.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, a1, x5, x4, x3);
}

void Task235::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I169
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
    // tensor label: I170
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), a1.size());
    dgemm_("T", "N", x5.size()*x4.size()*x3.size(), a1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*x4.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->put_block(odata, a1, x5, x4, x3);
}

void Task236::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  // tensor label: I170
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a2.size(), a1.size());
  }
  out()->put_block(odata, a2, a1);
}

void Task237::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I169
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
    // tensor label: I179
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", x5.size()*a1.size()*x4.size(), x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*a1.size()*x4.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, a1, x5, x4, x3);
}

void Task238::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I179
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a2.size(), x3.size());
  }
  out()->put_block(odata, a2, x3);
}

void Task239::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I160
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    // tensor label: I172
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size());
    dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size()*x2.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task240::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I172
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x3.size(), a1.size());
  }
  out()->put_block(odata, x3, a1);
}

void Task241::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I172
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I173
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", x3.size()*a1.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->put_block(odata, x3, a1);
}

void Task242::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I173
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2);
}

void Task243::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I172
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I176
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size());
      dgemm_("T", "N", x3.size()*a1.size(), 1, a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x3.size()*a1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->put_block(odata, x3, a1);
}

void Task244::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I176
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2);
}

void Task245::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I180
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, c1, a4, c3);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a2.size(), c1.size(), a4.size(), c3.size());
  }
  {
    // tensor label: I180
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3, a2, c1);
    sort_indices<1,0,3,2,1,1,1,1>(i0data, odata, a4.size(), c3.size(), a2.size(), c1.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task246::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I181
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task247::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  // tensor label: I181
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I182
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size());
    dgemm_("T", "N", x1.size(), a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a2.size());
  out()->put_block(odata, a2, x1);
}

void Task248::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  // tensor label: I182
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x0.size(), a2.size());
  }
  out()->put_block(odata, x0, a2);
}

void Task249::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I184
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

