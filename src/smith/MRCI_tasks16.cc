//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_tasks16.cc
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

#include <src/smith/MRCI_tasks16.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

void Task750::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a2 = b(2);
  const Index x5 = b(3);
  const Index c1 = b(4);
  const Index x4 = b(5);
  // tensor label: I894
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, a2, x5, c1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, a2, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, a2, x5, c1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c1, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c1, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c1.size(), x4.size());
    // tensor label: I895
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, x2, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, x2, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size(), x2.size(), a2.size());
    dgemm_("T", "N", x5.size()*c1.size()*x4.size(), x3.size()*x2.size()*a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*c1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c1.size(), x4.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, x3, x2, a2, x5, c1, x4);
}

void Task751::Task_local::compute() {
  const Index a3 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  // tensor label: I895
  std::unique_ptr<double[]> odata = out()->move_block(a3, x3, x2, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x3, x2, a2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), x3.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, a3, x3, x2, a2);
}

void Task752::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        for (auto& x4 : *range_[1]) {
          // tensor label: Gamma296
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x3, x4, x1, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x3, x4, x1, x0)]);
          sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x3.size(), x4.size(), x1.size(), x0.size());
          // tensor label: I900
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, x2, x5, c1, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, x2, x5, c1, x4)]);
          sort_indices<3,2,0,5,1,4,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), x2.size(), x5.size(), c1.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task753::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index c1 = b(4);
  const Index x4 = b(5);
  // tensor label: I900
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, x2, x5, c1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, a2, x2, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, a2, x2, x5, c1, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c1, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c1, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c1.size(), x4.size());
    // tensor label: I901
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, a3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), a3.size(), x2.size());
    dgemm_("T", "N", x5.size()*c1.size()*x4.size(), x3.size()*a2.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*c1.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c1.size(), x4.size(), x3.size(), a2.size(), x2.size());
  out()->put_block(odata, x3, a2, x2, x5, c1, x4);
}

void Task754::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I901
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, a3, x2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), a2.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, a3, x2);
}

void Task755::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& x5 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x5, x4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x5, x4)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x5.size(), x4.size());
        // tensor label: I915
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x5, x4, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x5, x4, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x5.size(), x4.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), c3.size()*x5.size()*x4.size(),
               1.0, i0data_sorted, c3.size()*x5.size()*x4.size(), i1data_sorted, c3.size()*x5.size()*x4.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task756::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I915
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x5, x4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x5, x4, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma240
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
      sort_indices<2,3,0,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
      // tensor label: I916
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, c3, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, c3, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c1.size(), c3.size(), x3.size(), x2.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, x5, x4, x1, x0);
}

void Task757::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I916
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c3, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, c1.size(), c3.size(), x3.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i1data = in(0)->get_block(x3, c3, c1, x2);
    sort_indices<2,1,0,3,1,1,1,2>(i1data, odata, x3.size(), c3.size(), c1.size(), x2.size());
  }
  {
    // tensor label: v2
    std::unique_ptr<double[]> i2data = in(0)->get_block(x3, x2, c1, c3);
    sort_indices<2,3,0,1,1,1,-1,1>(i2data, odata, x3.size(), x2.size(), c1.size(), c3.size());
  }
  out()->put_block(odata, c1, c3, x3, x2);
}

void Task758::Task_local::compute() {
  const Index c1 = b(0);
  const Index c3 = b(1);
  const Index x5 = b(2);
  const Index x4 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I915
  std::unique_ptr<double[]> odata = out()->move_block(c1, c3, x5, x4, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, c3, x5, x4, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, c3, x5, x4, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma4
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
      sort_indices<2,3,0,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I922
      std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, x2, c3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, x2, c3)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), c3.size());
      dgemm_("T", "N", x5.size()*x4.size()*x1.size()*x0.size(), c1.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x5.size()*x4.size()*x1.size()*x0.size());
    }
  }
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x1.size(), x0.size(), c1.size(), c3.size());
  out()->put_block(odata, c1, c3, x5, x4, x1, x0);
}

void Task759::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index c3 = b(3);
  // tensor label: I922
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, x2, c3);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, x2, c3);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c1.size(), x3.size(), x2.size(), c3.size());
  }
  out()->put_block(odata, c1, x3, x2, c3);
}

void Task760::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        for (auto& x3 : *range_[1]) {
          // tensor label: Gamma4
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x2, x3, x1, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x2, x3, x1, x0)]);
          sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x2.size(), x3.size(), x1.size(), x0.size());
          // tensor label: I924
          std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, a2, c1, x5, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, a2, c1, x5, x4)]);
          sort_indices<4,5,1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), a2.size(), c1.size(), x5.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task761::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  const Index x5 = b(4);
  const Index x4 = b(5);
  // tensor label: I924
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, a2, c1, x5, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x2, a2, c1, x5, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x2, a2, c1, x5, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a3, x5, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a3, x5, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a3.size(), x5.size(), x4.size());
    // tensor label: I925
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x3, x2, a2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x3, x2, a2)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x3.size(), x2.size(), a2.size());
    dgemm_("T", "N", c1.size()*x5.size()*x4.size(), x3.size()*x2.size()*a2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, c1.size()*x5.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c1.size(), x5.size(), x4.size(), x3.size(), x2.size(), a2.size());
  out()->put_block(odata, x3, x2, a2, c1, x5, x4);
}

void Task762::Task_local::compute() {
  const Index a3 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index a2 = b(3);
  // tensor label: I925
  std::unique_ptr<double[]> odata = out()->move_block(a3, x3, x2, a2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x3, x2, a2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a3.size(), x3.size(), x2.size(), a2.size());
  }
  out()->put_block(odata, a3, x3, x2, a2);
}

void Task763::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      for (auto& x5 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, x6, x5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, a2, x6, x5)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x7.size(), a2.size(), x6.size(), x5.size());
        // tensor label: I945
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x7, x6, x5, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x7, x6, x5, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i1data, i1data_sorted, c1.size(), x7.size(), x6.size(), x5.size(), x1.size(), x0.size());
        dgemm_("T", "N", a2.size(), c1.size()*x1.size()*x0.size(), x7.size()*x6.size()*x5.size(),
               1.0, i0data_sorted, x7.size()*x6.size()*x5.size(), i1data_sorted, x7.size()*x6.size()*x5.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task764::Task_local::compute() {
  const Index c1 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I945
  std::unique_ptr<double[]> odata = out()->move_block(c1, x7, x6, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x7, x6, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x7, x6, x5, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma312
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x4, x6, x5, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x4, x6, x5, x3, x2, x1, x0)]);
        sort_indices<1,4,5,0,2,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x4.size(), x6.size(), x5.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I946
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x7.size()*x6.size()*x5.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x7, x6, x5, x1, x0);
}

void Task765::Task_local::compute() {
  const Index c1 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I946
  std::unique_ptr<double[]> odata = out()->move_block(c1, x4, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x4, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, c1.size(), x4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c1, x4, x3, x2);
}

void Task766::Task_local::compute() {
  const Index c1 = b(0);
  const Index x7 = b(1);
  const Index x6 = b(2);
  const Index x5 = b(3);
  const Index x1 = b(4);
  const Index x0 = b(5);
  // tensor label: I945
  std::unique_ptr<double[]> odata = out()->move_block(c1, x7, x6, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x7, x6, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x7, x6, x5, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma313
        std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x2, x6, x5, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x2, x6, x5, x4, x3, x1, x0)]);
        sort_indices<1,4,5,0,2,3,6,7,0,1,1,1>(i0data, i0data_sorted, x7.size(), x2.size(), x6.size(), x5.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I949
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, c1, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, c1, x2)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), c1.size(), x2.size());
        dgemm_("T", "N", x7.size()*x6.size()*x5.size()*x1.size()*x0.size(), c1.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x7.size()*x6.size()*x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<5,0,1,2,3,4,1,1,1,1>(odata_sorted, odata, x7.size(), x6.size(), x5.size(), x1.size(), x0.size(), c1.size());
  out()->put_block(odata, c1, x7, x6, x5, x1, x0);
}

void Task767::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I949
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, c1, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), x3.size(), c1.size(), x2.size());
  }
  out()->put_block(odata, x4, x3, c1, x2);
}

void Task768::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x2, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a3, x2, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, a3.size(), x2.size(), c1.size(), a2.size());
      // tensor label: I951
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task769::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I951
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I952
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x2, x1, x0);
}

void Task770::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I952
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

void Task771::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x2, a3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x2, a3, a2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), x2.size(), a3.size(), a2.size());
      // tensor label: I954
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x2, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x2, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x2.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), a3.size()*x2.size(),
             1.0, i0data_sorted, a3.size()*x2.size(), i1data_sorted, a3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task772::Task_local::compute() {
  const Index a3 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I954
  std::unique_ptr<double[]> odata = out()->move_block(a3, x2, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x2, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x2, x1, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I955
        std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a3, x4, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a3, x4, x3)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, x5.size(), a3.size(), x4.size(), x3.size());
        dgemm_("T", "N", x2.size()*x1.size()*x0.size(), a3.size(), x5.size()*x4.size()*x3.size(),
               1.0, i0data_sorted, x5.size()*x4.size()*x3.size(), i1data_sorted, x5.size()*x4.size()*x3.size(),
               1.0, odata_sorted, x2.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x2, x1, x0);
}

void Task773::Task_local::compute() {
  const Index x5 = b(0);
  const Index a3 = b(1);
  const Index x4 = b(2);
  const Index x3 = b(3);
  // tensor label: I955
  std::unique_ptr<double[]> odata = out()->move_block(x5, a3, x4, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, x4, x3);
    sort_indices<0,1,2,3,1,1,-1,1>(i0data, odata, x5.size(), a3.size(), x4.size(), x3.size());
  }
  out()->put_block(odata, x5, a3, x4, x3);
}

void Task774::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a3, c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a3, c1, a2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a3.size(), c1.size(), a2.size());
      // tensor label: I993
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a2.size(), x1.size()*x0.size(), a3.size()*x5.size(),
             1.0, i0data_sorted, a3.size()*x5.size(), i1data_sorted, a3.size()*x5.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<0,2,3,1,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task775::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I993
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma240
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I994
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x5, x1, x0);
}

void Task776::Task_local::compute() {
  const Index a3 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I994
  std::unique_ptr<double[]> odata = out()->move_block(a3, x4, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x4, x3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, a3.size(), x4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, a3, x4, x3, x2);
}

void Task777::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I993
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I1000
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a3, x2)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x5, x1, x0);
}

void Task778::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1000
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, a3, x2);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x4.size(), x3.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x4, x3, a3, x2);
}

void Task779::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, a3)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), a3.size());
      // tensor label: I996
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x5, x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x5, x1, x0)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), x5.size(), x1.size(), x0.size());
      dgemm_("T", "N", a2.size()*c1.size(), x1.size()*x0.size(), a3.size()*x5.size(),
             1.0, i0data_sorted, a3.size()*x5.size(), i1data_sorted, a3.size()*x5.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x1.size(), x0.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task780::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I996
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  for (auto& x4 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: Gamma240
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x2, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x2, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x2.size(), x1.size(), x0.size());
        // tensor label: I997
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, x2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x5, x1, x0);
}

void Task781::Task_local::compute() {
  const Index a3 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I997
  std::unique_ptr<double[]> odata = out()->move_block(a3, x4, x3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x4, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, a3.size(), x4.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, a3, x4, x3, x2);
}

void Task782::Task_local::compute() {
  const Index a3 = b(0);
  const Index x5 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I996
  std::unique_ptr<double[]> odata = out()->move_block(a3, x5, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, x5, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, x5, x1, x0), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma252
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x2, x4, x3, x1, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x2, x4, x3, x1, x0)]);
        sort_indices<1,2,3,0,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x2.size(), x4.size(), x3.size(), x1.size(), x0.size());
        // tensor label: I1003
        std::unique_ptr<double[]> i1data = in(1)->get_block(x4, x3, a3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x4, x3, a3, x2)]);
        sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x4.size(), x3.size(), a3.size(), x2.size());
        dgemm_("T", "N", x5.size()*x1.size()*x0.size(), a3.size(), x4.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, x4.size()*x3.size()*x2.size(), i1data_sorted, x4.size()*x3.size()*x2.size(),
               1.0, odata_sorted, x5.size()*x1.size()*x0.size());
      }
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), x1.size(), x0.size(), a3.size());
  out()->put_block(odata, a3, x5, x1, x0);
}

void Task783::Task_local::compute() {
  const Index x4 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1003
  std::unique_ptr<double[]> odata = out()->move_block(x4, x3, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x4, x3, a3, x2);
    sort_indices<0,1,2,3,1,1,-1,2>(i0data, odata, x4.size(), x3.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, x4, x3, a3, x2);
}

void Task784::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        for (auto& x2 : *range_[1]) {
          // tensor label: Gamma338
          std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x3, x4, x2, x1, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x3, x4, x2, x1, x0)]);
          sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x3.size(), x4.size(), x2.size(), x1.size(), x0.size());
          // tensor label: I1023
          std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, x2, x5, a2, x4);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, x2, x5, a2, x4)]);
          sort_indices<3,1,5,2,0,4,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), x2.size(), x5.size(), a2.size(), x4.size());
          dgemm_("T", "N", x1.size()*x0.size(), c1.size()*a2.size(), x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, i0data_sorted, x3.size()*x2.size()*x5.size()*x4.size(), i1data_sorted, x3.size()*x2.size()*x5.size()*x4.size(),
                 1.0, odata_sorted, x1.size()*x0.size());
        }
      }
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c1.size(), a2.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task785::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x5 = b(3);
  const Index a2 = b(4);
  const Index x4 = b(5);
  // tensor label: I1023
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, x2, x5, a2, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x3, x2, x5, a2, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x3, x2, x5, a2, x4), 0.0);
  for (auto& a3 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, x4, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, x4, a3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), x4.size(), a3.size());
    // tensor label: I1024
    std::unique_ptr<double[]> i1data = in(1)->get_block(c1, x3, a3, x2);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, x3, a3, x2)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, c1.size(), x3.size(), a3.size(), x2.size());
    dgemm_("T", "N", x5.size()*a2.size()*x4.size(), c1.size()*x3.size()*x2.size(), a3.size(),
           1.0, i0data_sorted, a3.size(), i1data_sorted, a3.size(),
           1.0, odata_sorted, x5.size()*a2.size()*x4.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a2.size(), x4.size(), c1.size(), x3.size(), x2.size());
  out()->put_block(odata, c1, x3, x2, x5, a2, x4);
}

void Task786::Task_local::compute() {
  const Index c1 = b(0);
  const Index x3 = b(1);
  const Index a3 = b(2);
  const Index x2 = b(3);
  // tensor label: I1024
  std::unique_ptr<double[]> odata = out()->move_block(c1, x3, a3, x2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x3, a3, x2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), x3.size(), a3.size(), x2.size());
  }
  out()->put_block(odata, c1, x3, a3, x2);
}

void Task787::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      // tensor label: Gamma569
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x0.size());
      // tensor label: I1721
      std::unique_ptr<double[]> i1data = in(1)->get_block(x5, a2, c1, x4);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x5, a2, c1, x4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x5.size()*x4.size(),
             1.0, i0data_sorted, x5.size()*x4.size(), i1data_sorted, x5.size()*x4.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task788::Task_local::compute() {
  const Index x5 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x4 = b(3);
  // tensor label: I1721
  std::unique_ptr<double[]> odata = out()->move_block(x5, a2, c1, x4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x5.size(), a2.size(), c1.size(), x4.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x5, x4);
    sort_indices<2,1,0,3,1,1,-2,1>(i1data, odata, c1.size(), a2.size(), x5.size(), x4.size());
  }
  out()->put_block(odata, x5, a2, c1, x4);
}

void Task789::Task_local::compute() {
  const Index c1 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index a2 = b(3);
  // tensor label: I108
  std::unique_ptr<double[]> odata = out()->move_block(c1, x1, x0, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c1, x1, x0, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c1, x1, x0, a2), 0.0);
  for (auto& x7 : *range_[1]) {
    for (auto& x6 : *range_[1]) {
      // tensor label: Gamma573
      std::unique_ptr<double[]> i0data = in(0)->get_block(x7, x6, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x7, x6, x1, x0)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x7.size(), x6.size(), x1.size(), x0.size());
      // tensor label: I1729
      std::unique_ptr<double[]> i1data = in(1)->get_block(x7, a2, c1, x6);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x7, a2, c1, x6)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x7.size(), a2.size(), c1.size(), x6.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c1.size(), x7.size()*x6.size(),
             1.0, i0data_sorted, x7.size()*x6.size(), i1data_sorted, x7.size()*x6.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, c1, x1, x0, a2);
}

void Task790::Task_local::compute() {
  const Index x7 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x6 = b(3);
  // tensor label: I1729
  std::unique_ptr<double[]> odata = out()->move_block(x7, a2, c1, x6);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x7, a2, c1, x6);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x7.size(), a2.size(), c1.size(), x6.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x7, x6);
    sort_indices<2,1,0,3,1,1,-1,1>(i1data, odata, c1.size(), a2.size(), x7.size(), x6.size());
  }
  out()->put_block(odata, x7, a2, c1, x6);
}

void Task791::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index a1 = b(3);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2, x0, a1);
  {
    // tensor label: I144
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, x0, x2, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, a1.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x1, x2, x0, a1);
}

void Task792::Task_local::compute() {
  const Index a1 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I144
  std::unique_ptr<double[]> odata = out()->move_block(a1, x0, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, x0, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, x0, x2, x1), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x4 : *range_[1]) {
        // tensor label: Gamma48
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x0, x3, x4, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x0, x3, x4, x2, x1)]);
        sort_indices<0,2,3,1,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x0.size(), x3.size(), x4.size(), x2.size(), x1.size());
        // tensor label: I145
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

void Task793::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I145
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, x4)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), x4.size());
    // tensor label: I146
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

void Task794::Task_local::compute() {
  const Index x3 = b(0);
  const Index c2 = b(1);
  // tensor label: I146
  std::unique_ptr<double[]> odata = out()->move_block(x3, c2);
  {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c2);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x3.size(), c2.size());
  }
  out()->put_block(odata, x3, c2);
}

void Task795::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I145
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c3, x4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c3, x4)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c3.size(), x4.size());
      // tensor label: I1039
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, c3, a2, a1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, c3, a2, a1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), c3.size(), a2.size(), a1.size());
      dgemm_("T", "N", x5.size()*x4.size(), x3.size()*a1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x5.size()*x4.size());
    }
  }
  sort_indices<2,0,3,1,1,1,1,1>(odata_sorted, odata, x5.size(), x4.size(), x3.size(), a1.size());
  out()->put_block(odata, x3, x5, a1, x4);
}

void Task796::Task_local::compute() {
  const Index x3 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index a1 = b(3);
  // tensor label: I1039
  std::unique_ptr<double[]> odata = out()->move_block(x3, c3, a2, a1);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, c3, a2, a1);
    sort_indices<0,1,2,3,1,1,1,1>(i0data, odata, x3.size(), c3.size(), a2.size(), a1.size());
  }
  out()->put_block(odata, x3, c3, a2, a1);
}

void Task797::Task_local::compute() {
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index a1 = b(2);
  const Index x4 = b(3);
  // tensor label: I145
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, a1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, a1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, a1, x4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a1, c2, a3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a1, c2, a3)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x5.size(), a1.size(), c2.size(), a3.size());
      // tensor label: I1084
      std::unique_ptr<double[]> i1data = in(1)->get_block(a3, x4, x3, c2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, x4, x3, c2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), x4.size(), x3.size(), c2.size());
      dgemm_("T", "N", x5.size()*a1.size(), x4.size()*x3.size(), a3.size()*c2.size(),
             1.0, i0data_sorted, a3.size()*c2.size(), i1data_sorted, a3.size()*c2.size(),
             1.0, odata_sorted, x5.size()*a1.size());
    }
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), a1.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, x5, a1, x4);
}

void Task798::Task_local::compute() {
  const Index a3 = b(0);
  const Index x4 = b(1);
  const Index x3 = b(2);
  const Index c2 = b(3);
  // tensor label: I1084
  std::unique_ptr<double[]> odata = out()->move_block(a3, x4, x3, c2);
  {
    // tensor label: v2
    std::unique_ptr<double[]> i0data = in(0)->get_block(a3, x4, x3, c2);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, a3.size(), x4.size(), x3.size(), c2.size());
  }
  out()->put_block(odata, a3, x4, x3, c2);
}

void Task799::Task_local::compute() {
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
        // tensor label: Gamma49
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x3, x0, x2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x3, x0, x2, x1)]);
        sort_indices<0,1,2,3,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x3.size(), x0.size(), x2.size(), x1.size());
        // tensor label: I148
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

#endif
