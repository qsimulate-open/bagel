//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_tasks7.cc
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


#include <src/smith/CAS_test_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

void Task300::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task301::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);
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
  sort_indices<3,2,0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), a4.size(), c1.size());
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task302::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      // tensor label: I188
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", a4.size()*c1.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->put_block(odata, a4, c1);
}

void Task303::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I188
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task304::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), x1.size(), x0.size());
      // tensor label: I194
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a4.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size());
  out()->put_block(odata, a4, c1);
}

void Task305::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I194
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task306::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);
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
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task307::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I190
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), x0.size());
      // tensor label: I191
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size());
  out()->put_block(odata, a2, c1);
}

void Task308::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I191
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task309::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  // tensor label: I190
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x1.size(), x0.size());
      // tensor label: I197
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size());
  out()->put_block(odata, a2, c1);
}

void Task310::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I197
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,4,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task311::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  for (auto& c5 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c5)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), c5.size());
    // tensor label: I203
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());
    dgemm_("T", "N", c3.size(), c1.size()*a4.size()*a2.size(), c5.size(),
           1.0, i0data_sorted, c5.size(), i1data_sorted, c5.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), c1.size(), a4.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task312::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c5 = b(2);
  const Index a2 = b(3);
  // tensor label: I203
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c5, a4);
    sort_indices<0,3,2,1,1,1,-8,1>(i1data, odata, c1.size(), a2.size(), c5.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c5, a2);
}

void Task313::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a5.size(), a4.size());
    // tensor label: I207
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), a2.size());
    dgemm_("T", "N", a4.size(), c1.size()*c3.size()*a2.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task314::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I207
  std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, c1.size(), a5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a5);
    sort_indices<0,3,2,1,1,1,8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a5.size());
  }
  out()->put_block(odata, c1, a5, c3, a2);
}

void Task315::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I180
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), x0.size());
    // tensor label: I211
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4, c1, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size(), c1.size(), a2.size());
    dgemm_("T", "N", c3.size(), a4.size()*c1.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, c3.size());
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task316::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I211
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I212
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c1.size()*a2.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a4.size()*c1.size()*a2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), a2.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}

void Task317::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I212
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task318::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: I211
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4, c1, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4, c1, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4, c1, a2), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I215
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
    dgemm_("T", "N", a2.size()*c1.size()*a4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a2.size()*c1.size()*a4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), a4.size(), x0.size());
  out()->put_block(odata, x0, a4, c1, a2);
}

void Task319::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I215
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task320::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index c1 = b(2);
  const Index a2 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, c1, a2);
  {
    // tensor label: I198
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  out()->put_block(odata, c3, a4, c1, a2);
}

void Task321::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I198
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
    dscal_(c1.size()*a4.size()*c3.size()*a2.size(), e0_, i0data.get(), 1);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, a4);
    dscal_(c1.size()*a2.size()*c3.size()*a4.size(), e0_, i1data.get(), 1);
    sort_indices<0,3,2,1,1,1,-8,1>(i1data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(1)->get_block(c1, a4, c3, a2);
    sort_indices<0,1,2,3,1,1,-4,1>(i2data, odata, c1.size(), a4.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i3data = in(1)->get_block(c1, a2, c3, a4);
    sort_indices<0,3,2,1,1,1,8,1>(i3data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task322::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I198
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I199
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task323::Task_local::compute() {
  // tensor label: I199
  std::unique_ptr<double[]> odata = out()->move_block();
  {
    // scalar
    // tensor label: Gamma69
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    sort_indices<1,1,-4,1>(i0data, odata);
  }
  out()->put_block(odata);
}

void Task324::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I198
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I201
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);
  dgemm_("T", "N", c1.size()*a2.size()*c3.size()*a4.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a2.size()*c3.size()*a4.size());
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, c1, a4, c3, a2);
}

void Task325::Task_local::compute() {
  // tensor label: I201
  std::unique_ptr<double[]> odata = out()->move_block();
  {
    // scalar
    // tensor label: Gamma69
    std::unique_ptr<double[]> i0data = in(0)->get_block();
    sort_indices<1,1,8,1>(i0data, odata);
  }
  out()->put_block(odata);
}

void Task326::Task_local::compute() {
  const Index c2 = b(0);
  const Index a3 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, a3, x0, a1);
  {
    // tensor label: I216
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a3, c2, a1);
    sort_indices<2,1,0,3,1,1,1,1>(i0data, odata, x0.size(), a3.size(), c2.size(), a1.size());
  }
  out()->put_block(odata, c2, a3, x0, a1);
}

void Task327::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size());
    // tensor label: I217
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a3, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a3, c2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a3.size(), c2.size());
    dgemm_("T", "N", a1.size(), x0.size()*a3.size()*c2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), a3.size(), c2.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task328::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I217
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
      // tensor label: I218
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", a3.size()*c2.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a3.size()*c2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a3, c2);
}

void Task329::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I218
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task330::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I217
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a3, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a3, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a3, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I224
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a3.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a3, c2);
}

void Task331::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I224
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task332::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size());
    // tensor label: I220
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a1, c2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a1, c2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a1.size(), c2.size());
    dgemm_("T", "N", a3.size(), x0.size()*a1.size()*c2.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, a3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size(), a1.size(), c2.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task333::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I220
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), x2.size());
      // tensor label: I221
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      dgemm_("T", "N", a1.size()*c2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a1.size()*c2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a1, c2);
}

void Task334::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x1 = b(2);
  const Index x2 = b(3);
  // tensor label: I221
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x1, x2);
  {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), x0.size(), x1.size(), x2.size());
  }
  out()->put_block(odata, x3, x0, x1, x2);
}

void Task335::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index c2 = b(3);
  // tensor label: I220
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a1, c2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a1, c2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a1, c2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
      // tensor label: I227
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c2.size()*a1.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x0, x1, a1, c2);
}

void Task336::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I227
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x1, x0);
  {
    // tensor label: Gamma35
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x3, x2, x1, x0);
}

void Task337::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size());
  // tensor label: I229
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a3);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a3)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a3.size());
  dgemm_("T", "N", c2.size()*a1.size(), x0.size()*a3.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size());
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x0.size(), a3.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task338::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  // tensor label: I229
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), x2.size(), x1.size());
        // tensor label: I230
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a3.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), x0.size());
  out()->put_block(odata, x0, a3);
}

void Task339::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I230
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task340::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  // tensor label: f1
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size());
  // tensor label: I232
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), a1.size());
  dgemm_("T", "N", c2.size()*a3.size(), x0.size()*a1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size());
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), x0.size(), a1.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task341::Task_local::compute() {
  const Index x0 = b(0);
  const Index a1 = b(1);
  // tensor label: I232
  std::unique_ptr<double[]> odata = out()->move_block(x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), x2.size(), x1.size());
        // tensor label: I233
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size());
  out()->put_block(odata, x0, a1);
}

void Task342::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I233
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task343::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, c4)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), c4.size());
      // tensor label: I235
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, c2, a3, c4, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, c2, a3, c4, a1)]);
      sort_indices<0,4,1,2,3,5,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), c2.size(), a3.size(), c4.size(), a1.size());
      dgemm_("T", "N", 1, x0.size()*c2.size()*a3.size()*a1.size(), x1.size()*c4.size(),
             1.0, i0data_sorted, x1.size()*c4.size(), i1data_sorted, x1.size()*c4.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, x0.size(), c2.size(), a3.size(), a1.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task344::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);
  // tensor label: I235
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
  // tensor label: I236
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a3.size()*c4.size()*a1.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a3.size()*c4.size()*a1.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c2.size(), a3.size(), c4.size(), a1.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}

void Task345::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I236
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task346::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index c2 = b(2);
  const Index a3 = b(3);
  const Index c4 = b(4);
  const Index a1 = b(5);
  // tensor label: I235
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, c2, a3, c4, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, c2, a3, c4, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, c2, a3, c4, a1), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
  // tensor label: I239
  std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
  dgemm_("T", "N", c2.size()*a1.size()*c4.size()*a3.size(), x1.size()*x0.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c2.size()*a1.size()*c4.size()*a3.size());
  sort_indices<4,5,0,3,2,1,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), c4.size(), a3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, c2, a3, c4, a1);
}

void Task347::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I239
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,2,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task348::Task_local::compute() {
  const Index x0 = b(0);
  const Index a3 = b(1);
  const Index c2 = b(2);
  const Index a1 = b(3);
  // tensor label: I216
  std::unique_ptr<double[]> odata = out()->move_block(x0, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I241
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size());
  out()->put_block(odata, x0, a3, c2, a1);
}

void Task349::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  // tensor label: I241
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0);
  {
    // tensor label: Gamma81
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x3.size(), x0.size());
  }
  out()->put_block(odata, x3, x0);
}

