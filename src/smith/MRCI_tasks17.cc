//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks17.cc
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

#include <src/smith/MRCI_tasks17.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task800::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x5, x4)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x5.size(), x4.size());
    // tensor label: I149
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

void Task801::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I149
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x3.size(), c2.size());
  }
  out()->put_block(odata, x3, c2);
}

void Task802::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), x4.size());
      // tensor label: I1042
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, a3, c2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), a3.size(), c2.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->put_block(odata, x3, a1, x5, x4);
}

void Task803::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1042
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, a3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, a3, c2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), a1.size(), a3.size(), c2.size());
  }
  out()->put_block(odata, x3, a1, a3, c2);
}

void Task804::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x5, x4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x5.size(), x4.size());
      // tensor label: I1051
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), a1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->put_block(odata, x3, a1, x5, x4);
}

void Task805::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1051
  std::unique_ptr<double[]> odata = out()->move_block(x3, c3, a2, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c3, a2, a1);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x3.size(), c3.size(), a2.size(), a1.size());
  }
  out()->put_block(odata, x3, c3, a2, a1);
}

void Task806::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x5, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x5, x4)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x5.size(), x4.size());
      // tensor label: I1054
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a1, a3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a1, a3, c2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x3.size(), a1.size(), a3.size(), c2.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->put_block(odata, x3, a1, x5, x4);
}

void Task807::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1054
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, a3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, a3, c2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), a1.size(), a3.size(), c2.size());
  }
  out()->put_block(odata, x3, a1, a3, c2);
}

void Task808::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  // tensor label: I148
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I1081
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, c2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), c2.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<3,1,0,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, a1, x5, x4);
}

void Task809::Task_local::compute() {
  const Index a3 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1081
  std::unique_ptr<double[]> odata = out()->move_block(a3, x4, x3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x4, x3, c2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), x4.size(), x3.size(), c2.size());
  }
  out()->put_block(odata, a3, x4, x3, c2);
}

void Task810::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma50
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x4, x3, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I151
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

void Task811::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), x3.size());
    // tensor label: I152
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

void Task812::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  // tensor label: I152
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a2.size(), a1.size());
  }
  out()->put_block(odata, a2, a1);
}

void Task813::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, x4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, x4, a2)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), x4.size(), a2.size());
    // tensor label: I161
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

void Task814::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I161
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a2.size(), x3.size());
  }
  out()->put_block(odata, a2, x3);
}

void Task815::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I1075
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x4, x3)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x4.size(), x3.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, a1, x5, x4, x3);
}

void Task816::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1075
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), c2.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, a3, c2);
    sort_indices<2,3,0,1,1,1,-1,2>(i1data, odata, x4.size(), x3.size(), a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2, x4, x3);
}

void Task817::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I1078
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, x4, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, x4, x3)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), x4.size(), x3.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, a1, x5, x4, x3);
}

void Task818::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1078
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), c2.size(), x4.size(), x3.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x4, x3, a3, c2);
    sort_indices<2,3,0,1,1,1,1,1>(i1data, odata, x4.size(), x3.size(), a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2, x4, x3);
}

void Task819::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c3, a2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c3.size(), a2.size());
      // tensor label: I1090
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, a2, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, a2, x3)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x4.size(), c3.size(), a2.size(), x3.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<1,0,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, a1, x5, x4, x3);
}

void Task820::Task_local::compute() {
  const Index x4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x3 = b(3);
  // tensor label: I1090
  std::unique_ptr<double[]> odata = out()->move_block(x4, c3, a2, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c3, a2, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), c3.size(), a2.size(), x3.size());
  }
  out()->put_block(odata, x4, c3, a2, x3);
}

void Task821::Task_local::compute() {
  const Index a1 = b(0);
  const Index x5 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I151
  std::unique_ptr<double[]> odata = out()->move_block(a1, x5, x4, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x5, x4, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x5, x4, x3), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, a3)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), a3.size());
      // tensor label: I1111
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, a2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, a2, a1)]);
      sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size(), a2.size(), a1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), a3.size()*a2.size(),
             1.0, i0data_sorted, a3.size()*a2.size(), i1data_sorted, a3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->put_block(odata, a1, x5, x4, x3);
}

void Task822::Task_local::compute() {
  const Index a3 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1111
  std::unique_ptr<double[]> odata = out()->move_block(a3, x3, a2, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x3, a2, a1);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a3.size(), x3.size(), a2.size(), a1.size());
  }
  out()->put_block(odata, a3, x3, a2, a1);
}

void Task823::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: Gamma51
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    // tensor label: I154
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

void Task824::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
      // tensor label: I155
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

void Task825::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I155
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2);
}

void Task826::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I158
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

void Task827::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  // tensor label: I158
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, c2);
    sort_indices<0,1,1,1,2,1>(i0data, odata, a3.size(), c2.size());
  }
  out()->put_block(odata, a3, c2);
}

void Task828::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
        // tensor label: I1069
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, c2)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), c2.size());
        dgemm_("T", "N", a1.size(), x3.size(), c4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, c4.size()*a3.size()*c2.size(), i1data_sorted, c4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size());
  out()->put_block(odata, x3, a1);
}

void Task829::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1069
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, a3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, a3, c2);
    sort_indices<0,1,2,3,1,1,-4,1>(i0data, odata, x3.size(), c4.size(), a3.size(), c2.size());
  }
  out()->put_block(odata, x3, c4, a3, c2);
}

void Task830::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c4 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
        // tensor label: I1072
        std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c4, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c4, a3, c2)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x3.size(), c4.size(), a3.size(), c2.size());
        dgemm_("T", "N", a1.size(), x3.size(), c4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, c4.size()*a3.size()*c2.size(), i1data_sorted, c4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x3.size());
  out()->put_block(odata, x3, a1);
}

void Task831::Task_local::compute() {
  const Index x3 = b(0);
  const Index c4 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1072
  std::unique_ptr<double[]> odata = out()->move_block(x3, c4, a3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c4, a3, c2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x3.size(), c4.size(), a3.size(), c2.size());
  }
  out()->put_block(odata, x3, c4, a3, c2);
}

void Task832::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a4, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a4, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a4.size(), c2.size(), a3.size());
        // tensor label: I1099
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, a3, c2)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), c2.size());
        dgemm_("T", "N", x3.size(), a1.size(), a4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, a4.size()*a3.size()*c2.size(), i1data_sorted, a4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->put_block(odata, x3, a1);
}

void Task833::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1099
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, a3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1, a3, c2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, a4.size(), a1.size(), a3.size(), c2.size());
  }
  out()->put_block(odata, a4, a1, a3, c2);
}

void Task834::Task_local::compute() {
  const Index x3 = b(0);
  const Index a1 = b(1);
  // tensor label: I154
  std::unique_ptr<double[]> odata = out()->move_block(x3, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a4.size());
        // tensor label: I1102
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, a1, a3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, a1, a3, c2)]);
        sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), a1.size(), a3.size(), c2.size());
        dgemm_("T", "N", x3.size(), a1.size(), a4.size()*a3.size()*c2.size(),
               1.0, i0data_sorted, a4.size()*a3.size()*c2.size(), i1data_sorted, a4.size()*a3.size()*c2.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x3.size(), a1.size());
  out()->put_block(odata, x3, a1);
}

void Task835::Task_local::compute() {
  const Index a4 = b(0);
  const Index a1 = b(1);
  const Index a3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1102
  std::unique_ptr<double[]> odata = out()->move_block(a4, a1, a3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a1, a3, c2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, a4.size(), a1.size(), a3.size(), c2.size());
  }
  out()->put_block(odata, a4, a1, a3, c2);
}

void Task836::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& c2 : *range_[0]) {
        // tensor label: v2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, a1, x3, c2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, a1, x3, c2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x4.size(), a1.size(), x3.size(), c2.size());
        // tensor label: I1026
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x3, x4, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x3, x4, x0, x2, x1)]);
        sort_indices<2,1,0,3,4,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x3.size(), x4.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), c2.size()*x3.size()*x4.size(),
               1.0, i0data_sorted, c2.size()*x3.size()*x4.size(), i1data_sorted, c2.size()*x3.size()*x4.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task837::Task_local::compute() {
  const Index c2 = b(0);
  const Index x3 = b(1);
  const Index x4 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1026
  std::unique_ptr<double[]> odata = out()->move_block(c2, x3, x4, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x3, x4, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x3, x4, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: Gamma339
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x3, x5, x4, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x3, x5, x4, x0, x2, x1)]);
        sort_indices<0,1,3,2,4,5,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x3.size(), x5.size(), x4.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1027
        std::unique_ptr<double[]> i1data = in(1)->get_block(x7, x6, c2, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, x6, c2, x5)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), x6.size(), c2.size(), x5.size());
        dgemm_("T", "N", x3.size()*x4.size()*x0.size()*x2.size()*x1.size(), c2.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, x3.size()*x4.size()*x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x3.size(), x4.size(), x0.size(), x2.size(), x1.size(), c2.size());
  out()->put_block(odata, c2, x3, x4, x0, x2, x1);
}

void Task838::Task_local::compute() {
  const Index x7 = b(0);
  const Index x6 = b(1);
  const Index c2 = b(2);
  const Index x5 = b(3);
  // tensor label: I1027
  std::unique_ptr<double[]> odata = out()->move_block(x7, x6, c2, x5);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, c2, x5);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x7.size(), x6.size(), c2.size(), x5.size());
  }
  out()->put_block(odata, x7, x6, c2, x5);
}

void Task839::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma340
        std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x5, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x4, x5, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x4.size(), x5.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I1029
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a1, x5);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a1, x5)]);
        sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a1.size(), x5.size());
        dgemm_("T", "N", x0.size()*x2.size()*x1.size(), a1.size(), x4.size()*x3.size()*x5.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x5.size(), i1data_sorted, x4.size()*x3.size()*x5.size(),
               1.0, odata_sorted, x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x0.size(), x2.size(), x1.size(), a1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task840::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a1 = b(2);
  const Index x5 = b(3);
  // tensor label: I1029
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, a1, x5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x4, x3, a1, x5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x4, x3, a1, x5), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x5);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x5)]);
      sort_indices<0,2,1,3,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x5.size());
      // tensor label: I1030
      std::unique_ptr<double[]> i1data = in(1)->get_block(x4, c3, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, c3, x3, c2)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), c3.size(), x3.size(), c2.size());
      dgemm_("T", "N", a1.size()*x5.size(), x4.size()*x3.size(), c3.size()*c2.size(),
             1.0, i0data_sorted, c3.size()*c2.size(), i1data_sorted, c3.size()*c2.size(),
             1.0, odata_sorted, a1.size()*x5.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a1.size(), x5.size(), x4.size(), x3.size());
  out()->put_block(odata, x4, x3, a1, x5);
}

void Task841::Task_local::compute() {
  const Index x4 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1030
  std::unique_ptr<double[]> odata = out()->move_block(x4, c3, x3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, c3, x3, c2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x4.size(), c3.size(), x3.size(), c2.size());
  }
  out()->put_block(odata, x4, c3, x3, c2);
}

void Task842::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& x6 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a1, c2, x6);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a1, c2, x6)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a1.size(), c2.size(), x6.size());
        // tensor label: I1032
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x7, x0, x6, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x7, x0, x6, x2, x1)]);
        sort_indices<1,0,3,2,4,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x7.size(), x0.size(), x6.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), c2.size()*x7.size()*x6.size(),
               1.0, i0data_sorted, c2.size()*x7.size()*x6.size(), i1data_sorted, c2.size()*x7.size()*x6.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task843::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1032
  std::unique_ptr<double[]> odata = out()->move_block(c2, x7, x0, x6, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x0, x6, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x0, x6, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma341
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x5, x6, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x5, x6, x4, x3, x2, x1)]);
        sort_indices<2,4,5,0,1,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x5.size(), x6.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1033
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, c2, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, c2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), c2.size(), x4.size(), x3.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x2.size(), x1.size(), c2.size());
  out()->put_block(odata, c2, x7, x0, x6, x2, x1);
}

void Task844::Task_local::compute() {
  const Index x5 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1033
  std::unique_ptr<double[]> odata = out()->move_block(x5, c2, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), c2.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, c2, x4, x3);
}

void Task845::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x0 = b(2);
  const Index x6 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1032
  std::unique_ptr<double[]> odata = out()->move_block(c2, x7, x0, x6, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x0, x6, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x0, x6, x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma342
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x0, x3, x6, x5, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x0, x3, x6, x5, x4, x2, x1)]);
        sort_indices<2,4,5,0,1,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x0.size(), x3.size(), x6.size(), x5.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I1036
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, x4, x3, c2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, x4, x3, c2)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x5.size(), x4.size(), x3.size(), c2.size());
        dgemm_("T", "N", x7.size()*x0.size()*x6.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x0.size()*x6.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x0.size(), x6.size(), x2.size(), x1.size(), c2.size());
  out()->put_block(odata, c2, x7, x0, x6, x2, x1);
}

void Task846::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1036
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, x3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, c2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x5.size(), x4.size(), x3.size(), c2.size());
  }
  out()->put_block(odata, x5, x4, x3, c2);
}

void Task847::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& x7 : *range_[1]) {
      for (auto& x6 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, x7, x6);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, x7, x6)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), x7.size(), x6.size());
        // tensor label: I1044
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x7, x6, x0, x2, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x7, x6, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i1data, i1data_sorted, c2.size(), x7.size(), x6.size(), x0.size(), x2.size(), x1.size());
        dgemm_("T", "N", a1.size(), x0.size()*x2.size()*x1.size(), c2.size()*x7.size()*x6.size(),
               1.0, i0data_sorted, c2.size()*x7.size()*x6.size(), i1data_sorted, c2.size()*x7.size()*x6.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<0,1,2,3,1,1,1,1>(odata_sorted, odata, a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, a1, x0, x2, x1);
}

void Task848::Task_local::compute() {
  const Index c2 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x0 = b(3);
  const Index x2 = b(4);
  const Index x1 = b(5);
  // tensor label: I1044
  std::unique_ptr<double[]> odata = out()->move_block(c2, x7, x6, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x7, x6, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x7, x6, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma345
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x5, x0, x4, x3, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x5, x0, x4, x3, x2, x1)]);
        sort_indices<2,4,5,0,1,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x5.size(), x0.size(), x4.size(), x3.size(), x2.size(), x1.size());
        // tensor label: I1045
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, c2, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, c2, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), c2.size(), x4.size(), x3.size());
        dgemm_("T", "N", x7.size()*x6.size()*x0.size()*x2.size()*x1.size(), c2.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x0.size()*x2.size()*x1.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x0.size(), x2.size(), x1.size(), c2.size());
  out()->put_block(odata, c2, x7, x6, x0, x2, x1);
}

void Task849::Task_local::compute() {
  const Index x5 = b(0);
  const Index c2 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I1045
  std::unique_ptr<double[]> odata = out()->move_block(x5, c2, x4, x3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, c2, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x5.size(), c2.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, c2, x4, x3);
}

#endif
