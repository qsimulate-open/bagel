//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks6.cc
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


#include <src/smith/CASPT2_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task250::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I184
  std::unique_ptr<double[]> odata = out()->move_block(a4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I185
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size());
    dgemm_("T", "N", x1.size(), a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), a4.size());
  out()->put_block(odata, a4, x1);
}

void Task251::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x0.size(), a4.size());
  }
  out()->put_block(odata, x0, a4);
}

void Task252::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
  // tensor label: I187
  std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
  dgemm_("T", "N", c3.size()*a2.size(), a4.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a2.size());
  sort_indices<1,3,2,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), a4.size(), c1.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task253::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma38
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I188
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a4, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a4.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->put_block(odata, a4, c1);
}

void Task254::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I188
  std::unique_ptr<double[]> odata = out()->move_block(x1, a4, c1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x1.size(), a4.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a4, x1, x0);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x1, a4, c1, x0);
}

void Task255::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size());
  // tensor label: I190
  std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size());
  dgemm_("T", "N", c3.size()*a4.size(), a2.size()*c1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c3.size()*a4.size());
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), a2.size(), c1.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task256::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I190
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma38
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I191
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, a2, c1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, a2, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), a2.size(), c1.size(), x0.size());
      dgemm_("T", "N", 1, a2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, a2, c1);
}

void Task257::Task_local::compute() {
  const Index x1 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x0 = b(3);
  // tensor label: I191
  std::unique_ptr<double[]> odata = out()->move_block(x1, a2, c1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x1.size(), a2.size(), c1.size(), x0.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x1, x0);
    sort_indices<2,1,0,3,1,1,4,1>(i1data, odata, c1.size(), a2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x1, a2, c1, x0);
}

void Task258::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I199
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task259::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I199
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I200
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->put_block(odata, c3, x1);
}

void Task260::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I200
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c3.size(), x0.size());
  }
  out()->put_block(odata, c3, x0);
}

void Task261::Task_local::compute() {
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
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I202
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task262::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I202
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I203
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x0.size());
    dgemm_("T", "N", x1.size(), c3.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->put_block(odata, c3, x1);
}

void Task263::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I203
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c3.size(), x0.size());
  }
  out()->put_block(odata, c3, x0);
}

void Task264::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3, x0, a1);
  {
    // tensor label: I204
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, x0, a1);
    sort_indices<1,0,2,3,1,1,1,1>(i0data, odata, a3.size(), c2.size(), x0.size(), a1.size());
  }
  out()->put_block(odata, c2, a3, x0, a1);
}

void Task265::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I204
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size());
    // tensor label: I205
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x1, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x1.size(), x0.size());
    dgemm_("T", "N", a1.size(), a3.size()*c2.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), a3.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task266::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I205
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I206
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a3.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a3.size(), c2.size());
  out()->put_block(odata, a3, c2, x1, x0);
}

void Task267::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I206
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, c2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c2, a3, x3, x2);
    sort_indices<2,1,0,3,1,1,2,1>(i1data, odata, c2.size(), a3.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x3, a3, c2, x2);
}

void Task268::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I204
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: I208
    std::unique_ptr<double[]> i1data = in(1)->get_block(a1, c2, x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, c2, x0, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a1.size(), c2.size(), x0.size(), x1.size());
    dgemm_("T", "N", a3.size(), a1.size()*c2.size()*x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size(), c2.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task269::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I208
  std::unique_ptr<double[]> odata = out()->move_block(a1, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I209
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a1.size()*c2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a1.size(), c2.size());
  out()->put_block(odata, a1, c2, x0, x1);
}

void Task270::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x2 = b(3);
  // tensor label: I209
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, c2, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), a1.size(), c2.size(), x2.size());
  }
  out()->put_block(odata, x3, a1, c2, x2);
}

void Task271::Task_local::compute() {
  const Index a1 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I208
  std::unique_ptr<double[]> odata = out()->move_block(a1, c2, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, c2, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, c2, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma35
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I215
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, a1, c2, x0, x1);
}

void Task272::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I215
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), a1.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task273::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I204
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I217
  std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size(), a3.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task274::Task_local::compute() {
  const Index a3 = b(0);
  const Index x0 = b(1);
  // tensor label: I217
  std::unique_ptr<double[]> odata = out()->move_block(a3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma60
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I218
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a3.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size());
  out()->put_block(odata, a3, x0);
}

void Task275::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I218
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, x2, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a3, x2, x1);
}

void Task276::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I204
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
  // tensor label: I220
  std::unique_ptr<double[]> i1data = in(1)->get_block(a1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a3.size(), a1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size());
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task277::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  // tensor label: I220
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma60
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I221
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        dgemm_("T", "N", x0.size(), a1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), a1.size());
  out()->put_block(odata, a1, x0);
}

void Task278::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I221
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x2, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), a1.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, a1, x2, x1);
}

void Task279::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I204
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I223
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a3.size()*a1.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a3.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), a1.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task280::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I223
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I224
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), c4.size());
    dgemm_("T", "N", x0.size(), c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task281::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);
  // tensor label: I224
  std::unique_ptr<double[]> odata = out()->move_block(x1, c4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, x1.size(), c4.size());
  }
  out()->put_block(odata, x1, c4);
}

void Task282::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I204
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I226
    std::unique_ptr<double[]> i1data = in(1)->get_block(c4, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c4, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c4.size(), x0.size());
    dgemm_("T", "N", c2.size()*a1.size()*a3.size(), x0.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, c2.size()*a1.size()*a3.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), a3.size(), x0.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task283::Task_local::compute() {
  const Index c4 = b(0);
  const Index x0 = b(1);
  // tensor label: I226
  std::unique_ptr<double[]> odata = out()->move_block(c4, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I227
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), c4.size());
    dgemm_("T", "N", x0.size(), c4.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c4.size());
  out()->put_block(odata, c4, x0);
}

void Task284::Task_local::compute() {
  const Index x1 = b(0);
  const Index c4 = b(1);
  // tensor label: I227
  std::unique_ptr<double[]> odata = out()->move_block(x1, c4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), c4.size());
  }
  out()->put_block(odata, x1, c4);
}

void Task285::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I204
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma79
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size());
    // tensor label: I229
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    dgemm_("T", "N", x0.size(), a3.size()*c2.size()*a1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task286::Task_local::compute() {
  const Index x3 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I229
  std::unique_ptr<double[]> odata = out()->move_block(x3, a3, c2, a1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, a1, c2, a3);
    sort_indices<0,3,2,1,1,1,2,1>(i1data, odata, x3.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, x3, a3, c2, a1);
}

void Task287::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I204
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, x0, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I233
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x1, a3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x1, a3, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), a3.size(), a1.size());
    dgemm_("T", "N", x0.size(), c2.size()*a3.size()*a1.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, a3, c2, x0, a1);
}

void Task288::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I233
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
    dscal_(x1.size()*a3.size()*c2.size()*a1.size(), e0_, i0data.get(), 1);
    sort_indices<2,0,1,3,1,1,1,1>(i0data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a1, c2, a3);
    dscal_(x1.size()*a1.size()*c2.size()*a3.size(), e0_, i1data.get(), 1);
    sort_indices<2,0,3,1,1,1,-2,1>(i1data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(1)->get_block(x1, a3, c2, a1);
    sort_indices<2,0,1,3,1,1,-2,1>(i2data, odata, x1.size(), a3.size(), c2.size(), a1.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i3data = in(1)->get_block(x1, a1, c2, a3);
    sort_indices<2,0,3,1,1,1,4,1>(i3data, odata, x1.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task289::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I233
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
    // tensor label: I234
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size());
    dgemm_("T", "N", x1.size()*a3.size()*a1.size(), c2.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x1.size()*a3.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), a1.size(), c2.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task290::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I234
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->put_block(odata, c2, c4);
}

void Task291::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I233
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& c4 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
    // tensor label: I237
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size());
    dgemm_("T", "N", x1.size()*a1.size()*a3.size(), c2.size(), c4.size(),
           1.0, i0data_sorted, c4.size(), i1data_sorted, c4.size(),
           1.0, odata_sorted, x1.size()*a1.size()*a3.size());
  }
  sort_indices<3,0,2,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), a3.size(), c2.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task292::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I237
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->put_block(odata, c2, c4);
}

void Task293::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I233
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a3.size());
    // tensor label: I240
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size());
    dgemm_("T", "N", x1.size()*c2.size()*a3.size(), a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*c2.size()*a3.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task294::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  // tensor label: I240
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a4.size(), a1.size());
  }
  out()->put_block(odata, a4, a1);
}

void Task295::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I233
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a4.size());
    // tensor label: I243
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size());
    dgemm_("T", "N", x1.size()*a3.size()*c2.size(), a1.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*a3.size()*c2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task296::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  // tensor label: I243
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), a1.size());
  }
  out()->put_block(odata, a4, a1);
}

void Task297::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I233
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c2, a1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c2.size(), a1.size());
    // tensor label: I246
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size());
    dgemm_("T", "N", x1.size()*c2.size()*a1.size(), a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*c2.size()*a1.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x1.size(), c2.size(), a1.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

void Task298::Task_local::compute() {
  const Index a4 = b(0);
  const Index a3 = b(1);
  // tensor label: I246
  std::unique_ptr<double[]> odata = out()->move_block(a4, a3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a4.size(), a3.size());
  }
  out()->put_block(odata, a4, a3);
}

void Task299::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I233
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, a3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, a3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, a3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a4)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a4.size());
    // tensor label: I249
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a3.size());
    dgemm_("T", "N", x1.size()*a1.size()*c2.size(), a3.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size()*a1.size()*c2.size());
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x1.size(), a1.size(), c2.size(), a3.size());
  out()->put_block(odata, c2, x1, a3, a1);
}

