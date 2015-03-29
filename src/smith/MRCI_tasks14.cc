//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks14.cc
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

#include <src/smith/MRCI_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task650::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I877
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2, a4, c3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x2.size(), a2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, x2, a2, a4, c3);
}

void Task651::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c4, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c4.size(), x3.size());
      // tensor label: I880
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task652::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I880
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, a3, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x2.size(), c4.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, x2, c4, a3, a2);
}

void Task653::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x3.size());
      // tensor label: I883
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", c1.size()*x3.size(), x2.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c1.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task654::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I883
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a2, a4, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), a2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, x2, a2, a4, c3);
}

void Task655::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I964
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, x2, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task656::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I964
  std::unique_ptr<double[]> odata = out()->move_block(a4, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, a4, x3, x2, c3);
}

void Task657::Task_local::compute() {
  const Index x2 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I112
  std::unique_ptr<double[]> odata = out()->move_block(x2, a2, c1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a2, c1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a2, c1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I967
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, x2, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, a2, c1, x3);
}

void Task658::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I967
  std::unique_ptr<double[]> odata = out()->move_block(a4, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a4.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, a4, x3, x2, c3);
}

void Task659::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I118
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

void Task660::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
    // tensor label: I119
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

void Task661::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  // tensor label: I119
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c1.size(), c3.size());
  }
  out()->put_block(odata, c1, c3);
}

void Task662::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), x2.size());
    // tensor label: I122
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

void Task663::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  // tensor label: I122
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2);
}

void Task664::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
    // tensor label: I125
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

void Task665::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  // tensor label: I125
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3);
    sort_indices<0,1,1,1,-2,1>(i0data, odata, c1.size(), c3.size());
  }
  out()->put_block(odata, c1, c3);
}

void Task666::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x3, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x3.size(), x2.size());
    // tensor label: I128
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

void Task667::Task_local::compute() {
  const Index a3 = b(0);
  const Index a2 = b(1);
  // tensor label: I128
  std::unique_ptr<double[]> odata = out()->move_block(a3, a2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a3.size(), a2.size());
  }
  out()->put_block(odata, a3, a2);
}

void Task668::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a2.size());
    // tensor label: I140
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

void Task669::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I140
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task670::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a3.size());
    // tensor label: I143
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

void Task671::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I143
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task672::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c4, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c4, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c4.size(), x2.size());
      // tensor label: I910
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, a3, a2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task673::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I910
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4, a3, a2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c1.size(), c4.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, c1, c4, a3, a2);
}

void Task674::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), x2.size());
      // tensor label: I913
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, a4, c3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task675::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I913
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, a4, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), a2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c1, a2, a4, c3);
}

void Task676::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), x3.size(), x2.size());
      // tensor label: I940
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c4, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c4, a3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), c4.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task677::Task_local::compute() {
  const Index c1 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I940
  std::unique_ptr<double[]> odata = out()->move_block(c1, c4, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4, a3, a2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c1.size(), c4.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, c1, c4, a3, a2);
}

void Task678::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), x3.size(), x2.size());
      // tensor label: I943
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, a4, c3)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), a4.size(), c3.size());
      dgemm_("T", "N", x3.size()*x2.size(), c1.size()*a2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task679::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I943
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, a4, c3);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, c1.size(), a2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, c1, a2, a4, c3);
}

void Task680::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I958
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task681::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I958
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a4.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a4, c3);
    sort_indices<2,3,0,1,1,1,-2,1>(i1data, odata, x3.size(), x2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3, x3, x2);
}

void Task682::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I961
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, x3, x2)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task683::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I961
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,4,1>(i0data, odata, a4.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, x2, a4, c3);
    sort_indices<2,3,0,1,1,1,4,1>(i1data, odata, x3.size(), x2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a4, c3, x3, x2);
}

void Task684::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, c4, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, c4, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), c4.size(), a2.size());
      // tensor label: I970
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task685::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I970
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, a3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), c4.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x3, c4, a3, x2);
}

void Task686::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), a3.size());
      // tensor label: I973
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), x2.size());
      dgemm_("T", "N", c1.size()*a2.size(), x3.size()*x2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x3.size(), x2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task687::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I973
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, a3, x2);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x3.size(), c4.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x3, c4, a3, x2);
}

void Task688::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c3.size(), a2.size());
      // tensor label: I1006
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, c1, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c1.size(), c3.size());
      dgemm_("T", "N", x3.size()*a2.size(), x2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a2.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), x2.size(), c1.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task689::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I1006
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, c1, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, c1, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, a4.size(), x2.size(), c1.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,1,1>(i1data, odata, c1.size(), x2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a4, x2, c1, c3);
}

void Task690::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I1009
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, c1, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, c1, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), c1.size(), c3.size());
      dgemm_("T", "N", x3.size()*a2.size(), x2.size()*c1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, x3.size()*a2.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a2.size(), x2.size(), c1.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task691::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index c1 = b(2);
  const Index c3 = b(3);
  // tensor label: I1009
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, c1, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, c1, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), x2.size(), c1.size(), c3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, x2, a4, c3);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), x2.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, a4, x2, c1, c3);
}

void Task692::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c1, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c1.size(), a3.size());
      // tensor label: I1018
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, a3, a2)]);
      sort_indices<0,2,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*c1.size(), x2.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*c1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x2.size(), a2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task693::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1018
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, a3, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a4.size(), x2.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, a4, x2, a3, a2);
}

void Task694::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index x2 = b(3);
  // tensor label: I118
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a2, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, a2, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, a2, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c1, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c1, a4)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c1.size(), a4.size());
      // tensor label: I1021
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x2, a3, a2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x2, a3, a2)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a4.size(), x2.size(), a3.size(), a2.size());
      dgemm_("T", "N", x3.size()*c1.size(), x2.size()*a2.size(), a4.size()*a3.size(),
             1.0, i0data_sorted, a4.size()*a3.size(), i1data_sorted, a4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*c1.size());
    }
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), c1.size(), x2.size(), a2.size());
  out()->put_block(odata, c1, x3, a2, x2);
}

void Task695::Task_local::compute() {
  const Index a4 = b(0);
  const Index x2 = b(1);
  const Index a3 = b(2);
  const Index a2 = b(3);
  // tensor label: I1021
  std::unique_ptr<double[]> odata = out()->move_block(a4, x2, a3, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x2, a3, a2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), x2.size(), a3.size(), a2.size());
  }
  out()->put_block(odata, a4, x2, a3, a2);
}

void Task696::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size());
    // tensor label: I130
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

void Task697::Task_local::compute() {
  const Index a2 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I130
  std::unique_ptr<double[]> odata = out()->move_block(a2, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I131
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

void Task698::Task_local::compute() {
  const Index x5 = b(0);
  const Index a2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I131
  std::unique_ptr<double[]> odata = out()->move_block(x5, a2, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), a2.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a2, x4, x3);
}

void Task699::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  // tensor label: Gamma32
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I133
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size());
  dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

#endif
