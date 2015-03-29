//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks10.cc
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

#include <src/smith/MRCI_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task450::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I645
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a2.size()*c3.size(), c1.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,1,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), c1.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task451::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I645
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, a4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, a4, x3, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma10
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I646
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x2, a4, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), x2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c1.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c1.size(), a4.size());
  out()->put_block(odata, c1, a4, x3, x0);
}

void Task452::Task_local::compute() {
  const Index c1 = b(0);
  const Index x2 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I646
  std::unique_ptr<double[]> odata = out()->move_block(c1, x2, a4, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, a4, x1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), x2.size(), a4.size(), x1.size());
  }
  out()->put_block(odata, c1, x2, a4, x1);
}

void Task453::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c1.size(), a2.size());
      // tensor label: I648
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), c3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task454::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I648
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I649
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, a4, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c3.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c3.size(), a4.size());
  out()->put_block(odata, c3, a4, x3, x0);
}

void Task455::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I649
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2, a4, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x2, a4, x1);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, c3.size(), x2.size(), a4.size(), x1.size());
  }
  out()->put_block(odata, c3, x2, a4, x1);
}

void Task456::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), a4.size());
      // tensor label: I651
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a4, x3, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a4, x3, x0)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c3.size(), a4.size(), x3.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), c3.size()*x0.size(), a4.size()*x3.size(),
             1.0, i0data_sorted, a4.size()*x3.size(), i1data_sorted, a4.size()*x3.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), c3.size(), x0.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task457::Task_local::compute() {
  const Index c3 = b(0);
  const Index a4 = b(1);
  const Index x3 = b(2);
  const Index x0 = b(3);
  // tensor label: I651
  std::unique_ptr<double[]> odata = out()->move_block(c3, a4, x3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a4, x3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a4, x3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma18
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I652
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x2, a4, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x2, a4, x1)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), x2.size(), a4.size(), x1.size());
      dgemm_("T", "N", x3.size()*x0.size(), c3.size()*a4.size(), x2.size()*x1.size(),
             1.0, i0data_sorted, x2.size()*x1.size(), i1data_sorted, x2.size()*x1.size(),
             1.0, odata_sorted, x3.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size(), c3.size(), a4.size());
  out()->put_block(odata, c3, a4, x3, x0);
}

void Task458::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  const Index a4 = b(2);
  const Index x1 = b(3);
  // tensor label: I652
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2, a4, x1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x2, a4, x1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), x2.size(), a4.size(), x1.size());
  }
  out()->put_block(odata, c3, x2, a4, x1);
}

void Task459::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma552
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x3)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x3.size());
    // tensor label: I1685
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    dgemm_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task460::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x3 = b(3);
  // tensor label: I1685
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, c1, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c3.size(), a2.size(), c1.size(), x3.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x3);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x3.size());
  }
  out()->put_block(odata, c3, a2, c1, x3);
}

void Task461::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I27
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    // tensor label: Gamma554
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x5);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x5)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x0.size(), x5.size());
    // tensor label: I1689
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x5);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x5)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x5.size());
    dgemm_("T", "N", x0.size(), c3.size()*a2.size()*c1.size(), x5.size(),
           1.0, i0data_sorted, x5.size(), i1data_sorted, x5.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, c3, x0, a2);
}

void Task462::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1689
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, c1, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x5);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, c3.size(), a2.size(), c1.size(), x5.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, c3, x5);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), c3.size(), x5.size());
  }
  out()->put_block(odata, c3, a2, c1, x5);
}

void Task463::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  {
    // tensor label: I72
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, x1, x0, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, c2.size(), x1.size(), x0.size(), a1.size());
  }
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task464::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, a1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), a1.size());
    // tensor label: I73
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x1, x2, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x1, x2, x0)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), x2.size(), x0.size());
    dgemm_("T", "N", a1.size(), c2.size()*x1.size()*x0.size(), x2.size(),
           1.0, i0data_sorted, x2.size(), i1data_sorted, x2.size(),
           1.0, odata_sorted, a1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x1.size(), x0.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task465::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I73
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x2, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma24
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x2, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        // tensor label: I74
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, c2, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, c2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), c2.size(), x3.size());
        dgemm_("T", "N", x1.size()*x2.size()*x0.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x1.size()*x2.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x2.size(), x0.size(), c2.size());
  out()->put_block(odata, c2, x1, x2, x0);
}

void Task466::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c2, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), x4.size(), c2.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c2, x3);
}

void Task467::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma25
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x2, x0)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: I76
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, c2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, c2, x3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), c2.size(), x3.size());
      dgemm_("T", "N", x1.size()*x0.size(), a1.size()*c2.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a1.size(), c2.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task468::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I76
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, c2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x3.size());
    // tensor label: I77
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", a1.size()*c2.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, a1, c2, x3);
}

void Task469::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I77
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x2.size(), c3.size());
  }
  out()->put_block(odata, x2, c3);
}

void Task470::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I76
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, c2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a1, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a1, c3, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a1.size(), c3.size(), x3.size());
      // tensor label: I682
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c2, c3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c2.size(), c3.size());
      dgemm_("T", "N", a1.size()*x3.size(), x2.size()*c2.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a1.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size(), x2.size(), c2.size());
  out()->put_block(odata, x2, a1, c2, x3);
}

void Task471::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index c2 = b(2);
  const Index c3 = b(3);
  // tensor label: I682
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, c2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, c2, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), c4.size(), c2.size(), c3.size());
  }
  out()->put_block(odata, x2, c4, c2, c3);
}

void Task472::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I76
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, c2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a3, c2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a3, c2, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a3.size(), c2.size(), x3.size());
      // tensor label: I688
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a1.size());
      dgemm_("T", "N", c2.size()*x3.size(), x2.size()*a1.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c2.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size(), x2.size(), a1.size());
  out()->put_block(odata, x2, a1, c2, x3);
}

void Task473::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I688
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, a3, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, a3, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), c4.size(), a3.size(), a1.size());
  }
  out()->put_block(odata, x2, c4, a3, a1);
}

void Task474::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I76
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, c2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c2, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c2, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c2.size(), x3.size());
      // tensor label: I691
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, a4, c3)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), a4.size(), c3.size());
      dgemm_("T", "N", c2.size()*x3.size(), x2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c2.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size(), x2.size(), a1.size());
  out()->put_block(odata, x2, a1, c2, x3);
}

void Task475::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I691
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, a4, c3);
    sort_indices<0,1,2,3,1,1,-2,1>(i0data, odata, x2.size(), a1.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, x2, a1, a4, c3);
}

void Task476::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I76
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, c2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), x3.size());
      // tensor label: I697
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, a1, a4, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, a1, a4, c3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), a1.size(), a4.size(), c3.size());
      dgemm_("T", "N", c2.size()*x3.size(), x2.size()*a1.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c2.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size(), x2.size(), a1.size());
  out()->put_block(odata, x2, a1, c2, x3);
}

void Task477::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index a4 = b(2);
  const Index c3 = b(3);
  // tensor label: I697
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, a4, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1, a4, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), a1.size(), a4.size(), c3.size());
  }
  out()->put_block(odata, x2, a1, a4, c3);
}

void Task478::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I76
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1, c2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1, c2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1, c2, x3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a4, c3, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a4, c3, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a4.size(), c3.size(), a1.size());
      // tensor label: I778
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, x2, c3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", c2.size()*a1.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, a1, c2, x3);
}

void Task479::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I778
  std::unique_ptr<double[]> odata = out()->move_block(a4, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a4.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, a4, x3, x2, c3);
}

void Task480::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma5
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I79
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c2, a1, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c2, a1, x3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c2.size(), a1.size(), x3.size());
      dgemm_("T", "N", x1.size()*x0.size(), c2.size()*a1.size(), x2.size()*x3.size(),
             1.0, i0data_sorted, x2.size()*x3.size(), i1data_sorted, x2.size()*x3.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task481::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I79
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c2, a1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c2, a1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x3)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x3.size());
    // tensor label: I80
    std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x2.size(), c3.size());
    dgemm_("T", "N", c2.size()*a1.size()*x3.size(), x2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c2.size()*a1.size()*x3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task482::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x2.size(), c3.size());
  }
  out()->put_block(odata, x2, c3);
}

void Task483::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I79
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c2, a1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c2, a1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c4, x3)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c4.size(), x3.size());
      // tensor label: I685
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, c2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, c2, c3)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), c2.size(), c3.size());
      dgemm_("T", "N", a1.size()*x3.size(), x2.size()*c2.size(), c4.size()*c3.size(),
             1.0, i0data_sorted, c4.size()*c3.size(), i1data_sorted, c4.size()*c3.size(),
             1.0, odata_sorted, a1.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size(), x2.size(), c2.size());
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task484::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index c2 = b(2);
  const Index c3 = b(3);
  // tensor label: I685
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, c2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, c2, c3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x2.size(), c4.size(), c2.size(), c3.size());
  }
  out()->put_block(odata, x2, c4, c2, c3);
}

void Task485::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I79
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c2, a1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c2, a1, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), x3.size());
      // tensor label: I694
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, c4, a3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, c4, a3, a1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x2.size(), c4.size(), a3.size(), a1.size());
      dgemm_("T", "N", c2.size()*x3.size(), x2.size()*a1.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c2.size()*x3.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, c2.size(), x3.size(), x2.size(), a1.size());
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task486::Task_local::compute() {
  const Index x2 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I694
  std::unique_ptr<double[]> odata = out()->move_block(x2, c4, a3, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c4, a3, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x2.size(), c4.size(), a3.size(), a1.size());
  }
  out()->put_block(odata, x2, c4, a3, a1);
}

void Task487::Task_local::compute() {
  const Index x2 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x3 = b(3);
  // tensor label: I79
  std::unique_ptr<double[]> odata = out()->move_block(x2, c2, a1, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c2, a1, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c2, a1, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, a4)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), a4.size());
      // tensor label: I781
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, x3, x2, c3)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", c2.size()*a1.size(), x3.size()*x2.size(), a4.size()*c3.size(),
             1.0, i0data_sorted, a4.size()*c3.size(), i1data_sorted, a4.size()*c3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x2.size());
  out()->put_block(odata, x2, c2, a1, x3);
}

void Task488::Task_local::compute() {
  const Index a4 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I781
  std::unique_ptr<double[]> odata = out()->move_block(a4, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, a4, x3, x2, c3);
}

void Task489::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: I72
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x0, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x0, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x0, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma27
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I82
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, a1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, a1, x2)]);
      sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), a1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c2.size()*a1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,1,0,3,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x1, x0, a1);
}

void Task490::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, a1, x2), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), x2.size());
    // tensor label: I83
    std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c2.size(), c3.size());
    dgemm_("T", "N", x3.size()*a1.size()*x2.size(), c2.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*x2.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), x2.size(), c2.size());
  out()->put_block(odata, c2, x3, a1, x2);
}

void Task491::Task_local::compute() {
  const Index c2 = b(0);
  const Index c3 = b(1);
  // tensor label: I83
  std::unique_ptr<double[]> odata = out()->move_block(c2, c3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c3);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, c2.size(), c3.size());
  }
  out()->put_block(odata, c2, c3);
}

void Task492::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, a1, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
    // tensor label: I86
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, a1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, a1)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), a1.size());
    dgemm_("T", "N", x3.size()*c2.size()*x2.size(), a1.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*c2.size()*x2.size());
  }
  sort_indices<1,0,3,2,1,1,1,1>(odata_sorted, odata, x3.size(), c2.size(), x2.size(), a1.size());
  out()->put_block(odata, c2, x3, a1, x2);
}

void Task493::Task_local::compute() {
  const Index a3 = b(0);
  const Index a1 = b(1);
  // tensor label: I86
  std::unique_ptr<double[]> odata = out()->move_block(a3, a1);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, a1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a3.size(), a1.size());
  }
  out()->put_block(odata, a3, a1);
}

void Task494::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, a1, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I107
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size());
    dgemm_("T", "N", x3.size()*a1.size()*c2.size(), x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x3.size()*a1.size()*c2.size());
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size(), c2.size(), x2.size());
  out()->put_block(odata, c2, x3, a1, x2);
}

void Task495::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  // tensor label: I107
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a3.size(), x2.size());
  }
  out()->put_block(odata, a3, x2);
}

void Task496::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, a1, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c4, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c4, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c4.size(), x2.size());
      // tensor label: I724
      std::unique_ptr<double[]> i1data = in(1)->get_block(c2, c4, a3, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, c4, a3, a1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), c4.size(), a3.size(), a1.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*a1.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, x3, a1, x2);
}

void Task497::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index a1 = b(3);
  // tensor label: I724
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4, a3, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4, a3, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c2.size(), c4.size(), a3.size(), a1.size());
  }
  out()->put_block(odata, c2, c4, a3, a1);
}

void Task498::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x2 = b(3);
  // tensor label: I82
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, a1, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, a1, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, a1, x2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c4 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
      // tensor label: I784
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, x2)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), x2.size());
      dgemm_("T", "N", c2.size()*a1.size(), x3.size()*x2.size(), c4.size()*a3.size(),
             1.0, i0data_sorted, c4.size()*a3.size(), i1data_sorted, c4.size()*a3.size(),
             1.0, odata_sorted, c2.size()*a1.size());
    }
  }
  sort_indices<0,2,1,3,1,1,1,1>(odata_sorted, odata, c2.size(), a1.size(), x3.size(), x2.size());
  out()->put_block(odata, c2, x3, a1, x2);
}

void Task499::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I784
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, a3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), c4.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x3, c4, a3, x2);
}

#endif
