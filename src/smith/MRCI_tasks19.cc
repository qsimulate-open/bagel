//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks19.cc
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/MRCI_tasks19.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task900::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<0,1,1,1,-4,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->put_block(odata, a5, a4);
}

void Task901::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1249
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a5, a4);
}

void Task902::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1249
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4, x1, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a5.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a4, a5, x0);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, x1.size(), a4.size(), a5.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, a5, a4);
    sort_indices<2,3,0,1,1,1,-2,1>(i2data, odata, x1.size(), x0.size(), a5.size(), a4.size());
  }
  out()->put_block(odata, a5, a4, x1, x0);
}

void Task903::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I185
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: I1261
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, x0, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, x0, a4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), x0.size(), a4.size());
      dgemm_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a5, a4);
}

void Task904::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a4 = b(3);
  // tensor label: I1261
  std::unique_ptr<double[]> odata = out()->move_block(a5, x1, x0, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, x0, a4);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a5.size(), x1.size(), x0.size(), a4.size());
  }
  out()->put_block(odata, a5, x1, x0, a4);
}

void Task905::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& a5 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a5.size());
    // tensor label: I187
    std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), a5.size(),
           1.0, i0data_sorted, a5.size(), i1data_sorted, a5.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task906::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4);
    sort_indices<0,1,1,1,8,1>(i0data, odata, a5.size(), a4.size());
  }
  out()->put_block(odata, a5, a4);
}

void Task907::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
      // tensor label: I1252
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, a4, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, a4, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a5.size(), a4.size(), x1.size(), x0.size());
      dgemm_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a5, a4);
}

void Task908::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1252
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4, x1, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, a5.size(), a4.size(), x1.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x1, a4, a5, x0);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, x1.size(), a4.size(), a5.size(), x0.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x1, x0, a5, a4);
    sort_indices<2,3,0,1,1,1,4,1>(i2data, odata, x1.size(), x0.size(), a5.size(), a4.size());
  }
  out()->put_block(odata, a5, a4, x1, x0);
}

void Task909::Task_local::compute() {
  const Index a5 = b(0);
  const Index a4 = b(1);
  // tensor label: I187
  std::unique_ptr<double[]> odata = out()->move_block(a5, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a5, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a5, a4), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
      // tensor label: I1264
      std::unique_ptr<double[]> i1data = in(1)->get_block(a5, x1, x0, a4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a5, x1, x0, a4)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a5.size(), x1.size(), x0.size(), a4.size());
      dgemm_("T", "N", 1, a5.size()*a4.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a5, a4);
}

void Task910::Task_local::compute() {
  const Index a5 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a4 = b(3);
  // tensor label: I1264
  std::unique_ptr<double[]> odata = out()->move_block(a5, x1, x0, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a5, x1, x0, a4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a5.size(), x1.size(), x0.size(), a4.size());
  }
  out()->put_block(odata, a5, x1, x0, a4);
}

void Task911::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
    // tensor label: I189
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

void Task912::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I189
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I190
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

void Task913::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I190
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c3.size(), x0.size());
  }
  out()->put_block(odata, c3, x0);
}

void Task914::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
    // tensor label: I192
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

void Task915::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I192
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma32
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I193
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

void Task916::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I193
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c3.size(), x0.size());
  }
  out()->put_block(odata, c3, x0);
}

void Task917::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, c3, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), c3.size(), a2.size());
    // tensor label: I1116
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
    dgemm_("T", "N", a4.size()*c3.size()*a2.size(), c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a4.size()*c3.size()*a2.size());
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c3.size(), a2.size(), c1.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task918::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1116
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I1117
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->put_block(odata, c1, x0);
}

void Task919::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1117
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), x2.size(), c1.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c1, x1);
}

void Task920::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, c3, a4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), c3.size(), a4.size());
    // tensor label: I1119
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x0)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x0.size());
    dgemm_("T", "N", a2.size()*c3.size()*a4.size(), c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, a2.size()*c3.size()*a4.size());
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), a4.size(), c1.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task921::Task_local::compute() {
  const Index c1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1119
  std::unique_ptr<double[]> odata = out()->move_block(c1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma10
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
        // tensor label: I1120
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, c1, x1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
        dgemm_("T", "N", x0.size(), c1.size(), x3.size()*x2.size()*x1.size(),
               1.0, i0data_sorted, x3.size()*x2.size()*x1.size(), i1data_sorted, x3.size()*x2.size()*x1.size(),
               1.0, odata_sorted, x0.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->put_block(odata, c1, x0);
}

void Task922::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I1120
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, c1, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), x2.size(), c1.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, c1, x1);
}

void Task923::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
    // tensor label: I1122
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), a2.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task924::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I1122
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma5
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I1123
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), x1.size(), x0.size());
        dgemm_("T", "N", x3.size(), a2.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size());
  out()->put_block(odata, a2, x3);
}

void Task925::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1123
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, x1, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2, x1, x0);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x2.size(), a2.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, a2, x1, x0);
}

void Task926::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I1122
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma197
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1129
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a2)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a2.size());
        dgemm_("T", "N", x3.size(), a2.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size());
  out()->put_block(odata, a2, x3);
}

void Task927::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I1129
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, x0, a2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x2.size(), x1.size(), x0.size(), a2.size());
  }
  out()->put_block(odata, x2, x1, x0, a2);
}

void Task928::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I1125
    std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), a4.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task929::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  // tensor label: I1125
  std::unique_ptr<double[]> odata = out()->move_block(a4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: Gamma5
        std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I1126
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a4, x1, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a4.size(), x1.size(), x0.size());
        dgemm_("T", "N", x3.size(), a4.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a4.size());
  out()->put_block(odata, a4, x3);
}

void Task930::Task_local::compute() {
  const Index x2 = b(0);
  const Index a4 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I1126
  std::unique_ptr<double[]> odata = out()->move_block(x2, a4, x1, x0);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a4, x1, x0);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), a4.size(), x1.size(), x0.size());
  }
  out()->put_block(odata, x2, a4, x1, x0);
}

void Task931::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  // tensor label: I1125
  std::unique_ptr<double[]> odata = out()->move_block(a4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: Gamma197
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3, x2, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1132
        std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x1, x0, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x1, x0, a4)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), x1.size(), x0.size(), a4.size());
        dgemm_("T", "N", x3.size(), a4.size(), x2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, x2.size()*x1.size()*x0.size(), i1data_sorted, x2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), a4.size());
  out()->put_block(odata, a4, x3);
}

void Task932::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a4 = b(3);
  // tensor label: I1132
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1, x0, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1, x0, a4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), x1.size(), x0.size(), a4.size());
  }
  out()->put_block(odata, x2, x1, x0, a4);
}

void Task933::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a4, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a4, c1, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c5.size(), a4.size(), c1.size(), x1.size());
      // tensor label: I1134
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c3, a2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", a4.size()*c1.size(), c3.size()*a2.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task934::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1134
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3, a2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1135
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c5, c3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c5, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a2.size());
    dgemm_("T", "N", x1.size(), c5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a2.size());
  out()->put_block(odata, c5, c3, a2, x1);
}

void Task935::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1135
  std::unique_ptr<double[]> odata = out()->move_block(x0, c5, c3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c5, c3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a2, c3, c5);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, x0.size(), a2.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, x0, c5, c3, a2);
}

void Task936::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c5, a2, c1, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c5, a2, c1, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c5.size(), a2.size(), c1.size(), x1.size());
      // tensor label: I1137
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c3, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c3, a4, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a4.size(), x1.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*a4.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<0,1,3,2,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task937::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1137
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3, a4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1138
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c5, c3, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c5, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a4.size());
    dgemm_("T", "N", x1.size(), c5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a4.size());
  out()->put_block(odata, c5, c3, a4, x1);
}

void Task938::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I1138
  std::unique_ptr<double[]> odata = out()->move_block(x0, c5, c3, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c5, c3, a4);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a4, c3, c5);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, x0.size(), a4.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, x0, c5, c3, a4);
}

void Task939::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c5, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c5.size(), x1.size());
      // tensor label: I1146
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c3, a2, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", c1.size()*a4.size(), c3.size()*a2.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task940::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1146
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3, a2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1147
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c5, c3, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c5, c3, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a2.size());
    dgemm_("T", "N", x1.size(), c5.size()*c3.size()*a2.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a2.size());
  out()->put_block(odata, c5, c3, a2, x1);
}

void Task941::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1147
  std::unique_ptr<double[]> odata = out()->move_block(x0, c5, c3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c5, c3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a2, c3, c5);
    sort_indices<0,3,2,1,1,1,1,1>(i1data, odata, x0.size(), a2.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, x0, c5, c3, a2);
}

void Task942::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& c5 : *range_[0]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c5, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c5, x1)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), x1.size());
      // tensor label: I1149
      std::unique_ptr<double[]> i1data = in(1)->get_block(c5, c3, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c5, c3, a4, x1)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, c5.size(), c3.size(), a4.size(), x1.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*a4.size(), c5.size()*x1.size(),
             1.0, i0data_sorted, c5.size()*x1.size(), i1data_sorted, c5.size()*x1.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), a4.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task943::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I1149
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3, a4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3, a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3, a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1150
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c5, c3, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c5, c3, a4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c5.size(), c3.size(), a4.size());
    dgemm_("T", "N", x1.size(), c5.size()*c3.size()*a4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c5.size(), c3.size(), a4.size());
  out()->put_block(odata, c5, c3, a4, x1);
}

void Task944::Task_local::compute() {
  const Index x0 = b(0);
  const Index c5 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I1150
  std::unique_ptr<double[]> odata = out()->move_block(x0, c5, c3, a4);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, c5, c3, a4);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, x0.size(), c5.size(), c3.size(), a4.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x0, a4, c3, c5);
    sort_indices<0,3,2,1,1,1,-2,1>(i1data, odata, x0.size(), a4.size(), c3.size(), c5.size());
  }
  out()->put_block(odata, x0, c5, c3, a4);
}

void Task945::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a4, a5, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a4, a5, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a4.size(), a5.size(), a2.size());
      // tensor label: I1158
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, x0)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x0.size());
      dgemm_("T", "N", a4.size()*a2.size(), c1.size()*c3.size(), a5.size()*x0.size(),
             1.0, i0data_sorted, a5.size()*x0.size(), i1data_sorted, a5.size()*x0.size(),
             1.0, odata_sorted, a4.size()*a2.size());
    }
  }
  sort_indices<1,2,0,3,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size(), c1.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task946::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1158
  std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1159
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x1.size());
    dgemm_("T", "N", x0.size(), c1.size()*a5.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a5.size(), c3.size());
  out()->put_block(odata, c1, a5, c3, x0);
}

void Task947::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index x1 = b(3);
  // tensor label: I1159
  std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, x1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, x1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a5.size(), c3.size(), x1.size());
  }
  out()->put_block(odata, c1, a5, c3, x1);
}

void Task948::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I162
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, a4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, a4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, a4, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a5 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a2, a5, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a2, a5, a4)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a2.size(), a5.size(), a4.size());
      // tensor label: I1161
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, x0)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x0.size());
      dgemm_("T", "N", a2.size()*a4.size(), c1.size()*c3.size(), a5.size()*x0.size(),
             1.0, i0data_sorted, a5.size()*x0.size(), i1data_sorted, a5.size()*x0.size(),
             1.0, odata_sorted, a2.size()*a4.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size(), c1.size(), c3.size());
  out()->put_block(odata, a2, c1, a4, c3);
}

void Task949::Task_local::compute() {
  const Index c1 = b(0);
  const Index a5 = b(1);
  const Index c3 = b(2);
  const Index x0 = b(3);
  // tensor label: I1161
  std::unique_ptr<double[]> odata = out()->move_block(c1, a5, c3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a5, c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a5, c3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I1162
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a5, c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a5, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a5.size(), c3.size(), x1.size());
    dgemm_("T", "N", x0.size(), c1.size()*a5.size()*c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size(), a5.size(), c3.size());
  out()->put_block(odata, c1, a5, c3, x0);
}

#endif
