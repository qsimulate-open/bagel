//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks18.cc
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


#include <src/smith/CASPT2_tasks18.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task850::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  // tensor label: I911
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
      // tensor label: I1068
      std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c1)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c1.size());
      dgemm_("T", "N", c3.size()*a2.size(), 1, a4.size()*c1.size(),
             1.0, i0data_sorted, a4.size()*c1.size(), i1data_sorted, a4.size()*c1.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size());
  out()->put_block(odata, c3, a2);
}

void Task851::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I1068
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, x0)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), x0.size());
      // tensor label: I1069
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", a4.size()*c1.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, a4.size()*c1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, a4.size(), c1.size());
  out()->put_block(odata, a4, c1);
}

void Task852::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1069
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task853::Task_local::compute() {
  const Index a4 = b(0);
  const Index c1 = b(1);
  // tensor label: I1068
  std::unique_ptr<double[]> odata = out()->move_block(a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, x1, x0)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), x1.size(), x0.size());
      // tensor label: I1075
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size());
      dgemm_("T", "N", c1.size()*a4.size(), 1, x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, c1.size()*a4.size());
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size());
  out()->put_block(odata, a4, c1);
}

void Task854::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I1075
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  {
    // tensor label: Gamma38
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<0,1,1,1,-1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x1, x0);
}

void Task855::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2);
  {
    // tensor label: I914
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), x1.size());
  }
  out()->put_block(odata, x1, x2);
}

void Task856::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  // tensor label: I914
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        for (auto& x0 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
          // tensor label: I915
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, c3, a2, c1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, c3, a2, c1)]);
          sort_indices<5,4,3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), c3.size(), a2.size(), c1.size());
          dgemm_("T", "N", 1, x2.size()*x1.size(), x0.size()*c3.size()*a2.size()*c1.size(),
                 1.0, i0data_sorted, x0.size()*c3.size()*a2.size()*c1.size(), i1data_sorted, x0.size()*c3.size()*a2.size()*c1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
  out()->put_block(odata, x2, x1);
}

void Task857::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  const Index c1 = b(5);
  // tensor label: I915
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, c3, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, c3, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, c3, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
    // tensor label: I916
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3, x2, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c1.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, c3, a2, c1);
}

void Task858::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I916
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  {
    // tensor label: Gamma282
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x0.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task859::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index c3 = b(3);
  const Index a2 = b(4);
  const Index c1 = b(5);
  // tensor label: I915
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, c3, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, c3, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, c3, a2, c1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
    // tensor label: I919
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3, x2, x1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c3.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c3.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c3.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, c3, a2, c1);
}

void Task860::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I919
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  {
    // tensor label: Gamma282
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x3, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x0.size(), x3.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task861::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  // tensor label: I914
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& c2 : *range_[0]) {
        for (auto& a3 : *range_[2]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, a3)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), a3.size());
          // tensor label: I1124
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x2, x1, a3, c2, a1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x2, x1, a3, c2, a1)]);
          sort_indices<0,5,4,3,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), x2.size(), x1.size(), a3.size(), c2.size(), a1.size());
          dgemm_("T", "N", 1, x2.size()*x1.size(), x0.size()*a3.size()*c2.size()*a1.size(),
                 1.0, i0data_sorted, x0.size()*a3.size()*c2.size()*a1.size(), i1data_sorted, x0.size()*a3.size()*c2.size()*a1.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
  out()->put_block(odata, x2, x1);
}

void Task862::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  const Index a1 = b(5);
  // tensor label: I1124
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
    // tensor label: I1125
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a3.size()*c2.size()*a1.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a3.size()*c2.size()*a1.size());
  }
  sort_indices<3,4,5,0,1,2,1,1,1,1>(odata_sorted, odata, a3.size(), c2.size(), a1.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a3, c2, a1);
}

void Task863::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1125
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task864::Task_local::compute() {
  const Index x0 = b(0);
  const Index x2 = b(1);
  const Index x1 = b(2);
  const Index a3 = b(3);
  const Index c2 = b(4);
  const Index a1 = b(5);
  // tensor label: I1124
  std::unique_ptr<double[]> odata = out()->move_block(x0, x2, x1, a3, c2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x2, x1, a3, c2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x2, x1, a3, c2, a1), 0.0);
  for (auto& x3 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
    sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
    // tensor label: I1128
    std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x0, x2, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x0, x2, x1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
    dgemm_("T", "N", a1.size()*c2.size()*a3.size(), x0.size()*x2.size()*x1.size(), x3.size(),
           1.0, i0data_sorted, x3.size(), i1data_sorted, x3.size(),
           1.0, odata_sorted, a1.size()*c2.size()*a3.size());
  }
  sort_indices<3,4,5,2,1,0,1,1,1,1>(odata_sorted, odata, a1.size(), c2.size(), a3.size(), x0.size(), x2.size(), x1.size());
  out()->put_block(odata, x0, x2, x1, a3, c2, a1);
}

void Task865::Task_local::compute() {
  const Index x3 = b(0);
  const Index x0 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I1128
  std::unique_ptr<double[]> odata = out()->move_block(x3, x0, x2, x1);
  {
    // tensor label: Gamma60
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x0.size(), x2.size(), x1.size());
  }
  out()->put_block(odata, x3, x0, x2, x1);
}

void Task866::Task_local::compute() {
  const Index c4 = b(0);
  const Index c1 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c4, c1);
  {
    // tensor label: I920
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c4.size(), c1.size());
  }
  out()->put_block(odata, c4, c1);
}

void Task867::Task_local::compute() {
  const Index c4 = b(0);
  const Index c1 = b(1);
  // tensor label: I920
  std::unique_ptr<double[]> odata = out()->move_block(c4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, c1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I921
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c4, a2, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c4, a2, c3)]);
        sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), c4.size(), a2.size(), c3.size());
        dgemm_("T", "N", c1.size(), c4.size(), x0.size()*a2.size()*c3.size(),
               1.0, i0data_sorted, x0.size()*a2.size()*c3.size(), i1data_sorted, x0.size()*a2.size()*c3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), c4.size());
  out()->put_block(odata, c4, c1);
}

void Task868::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I921
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4, a2, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
    // tensor label: I922
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c4.size()*a2.size()*c3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c4.size()*a2.size()*c3.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c4.size(), a2.size(), c3.size(), x0.size());
  out()->put_block(odata, x0, c4, a2, c3);
}

void Task869::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I922
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task870::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I921
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4, a2, c3), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
    // tensor label: I925
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a2.size()*c4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*a2.size()*c4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), c4.size(), x0.size());
  out()->put_block(odata, x0, c4, a2, c3);
}

void Task871::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I925
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task872::Task_local::compute() {
  const Index c4 = b(0);
  const Index c3 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(c4, c3);
  {
    // tensor label: I926
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c4.size(), c3.size());
  }
  out()->put_block(odata, c4, c3);
}

void Task873::Task_local::compute() {
  const Index c4 = b(0);
  const Index c3 = b(1);
  // tensor label: I926
  std::unique_ptr<double[]> odata = out()->move_block(c4, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I927
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c4, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c4, a2, c1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), c4.size(), a2.size(), c1.size());
        dgemm_("T", "N", c3.size(), c4.size(), x0.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, x0.size()*a2.size()*c1.size(), i1data_sorted, x0.size()*a2.size()*c1.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), c4.size());
  out()->put_block(odata, c4, c3);
}

void Task874::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I927
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4, a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
    // tensor label: I928
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c4.size()*a2.size()*c1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c4.size()*a2.size()*c1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c4.size(), a2.size(), c1.size(), x0.size());
  out()->put_block(odata, x0, c4, a2, c1);
}

void Task875::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I928
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,1,4>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task876::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I927
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4, a2, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
    // tensor label: I934
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a2.size()*c4.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a2.size()*c4.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), c4.size(), x0.size());
  out()->put_block(odata, x0, c4, a2, c1);
}

void Task877::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I934
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task878::Task_local::compute() {
  const Index a2 = b(0);
  const Index a4 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(a2, a4);
  {
    // tensor label: I929
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a4.size(), a2.size());
  }
  out()->put_block(odata, a2, a4);
}

void Task879::Task_local::compute() {
  const Index a4 = b(0);
  const Index a2 = b(1);
  // tensor label: I929
  std::unique_ptr<double[]> odata = out()->move_block(a4, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a2), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I930
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a4, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a4, c1)]);
        sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a4.size(), c1.size());
        dgemm_("T", "N", a2.size(), a4.size(), x0.size()*c3.size()*c1.size(),
               1.0, i0data_sorted, x0.size()*c3.size()*c1.size(), i1data_sorted, x0.size()*c3.size()*c1.size(),
               1.0, odata_sorted, a2.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a2.size(), a4.size());
  out()->put_block(odata, a4, a2);
}

void Task880::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  // tensor label: I930
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
    // tensor label: I931
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c3.size()*a4.size()*c1.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c3.size()*a4.size()*c1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, c3.size(), a4.size(), c1.size(), x0.size());
  out()->put_block(odata, x0, c3, a4, c1);
}

void Task881::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I931
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,4>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task882::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  // tensor label: I930
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a4, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a4, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a4, c1), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
    // tensor label: I937
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
    dgemm_("T", "N", c1.size()*a4.size()*c3.size(), x0.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, c1.size()*a4.size()*c3.size());
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), x0.size());
  out()->put_block(odata, x0, c3, a4, c1);
}

void Task883::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I937
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task884::Task_local::compute() {
  const Index x1 = b(0);
  const Index c1 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x1, c1);
  {
    // tensor label: I938
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c1.size());
  }
  out()->put_block(odata, x1, c1);
}

void Task885::Task_local::compute() {
  const Index x1 = b(0);
  const Index c1 = b(1);
  // tensor label: I938
  std::unique_ptr<double[]> odata = out()->move_block(x1, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I939
        std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c3)]);
        sort_indices<2,3,1,0,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c3.size());
        dgemm_("T", "N", c1.size(), x1.size(), x0.size()*a2.size()*c3.size(),
               1.0, i0data_sorted, x0.size()*a2.size()*c3.size(), i1data_sorted, x0.size()*a2.size()*c3.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size());
  out()->put_block(odata, x1, c1);
}

void Task886::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I939
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c3, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      // tensor label: I940
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x1, x0, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      dgemm_("T", "N", a2.size()*c3.size(), x1.size()*x0.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c3.size(), x1.size(), x0.size());
  out()->put_block(odata, x1, x0, a2, c3);
}

void Task887::Task_local::compute() {
  const Index x3 = b(0);
  const Index x1 = b(1);
  const Index x0 = b(2);
  const Index x2 = b(3);
  // tensor label: I940
  std::unique_ptr<double[]> odata = out()->move_block(x3, x1, x0, x2);
  {
    // tensor label: Gamma22
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
    sort_indices<0,1,2,3,1,1,1,4>(i0data, odata, x3.size(), x1.size(), x0.size(), x2.size());
  }
  out()->put_block(odata, x3, x1, x0, x2);
}

void Task888::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c3 = b(3);
  // tensor label: I939
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0, a2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0, a2, c3), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I946
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c3.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c3.size()*a2.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, c3.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x1, x0, a2, c3);
}

void Task889::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I946
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task890::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3);
  {
    // tensor label: I941
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c3.size());
  }
  out()->put_block(odata, x1, c3);
}

void Task891::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  // tensor label: I941
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I942
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, a2, c1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), a2.size(), c1.size());
        dgemm_("T", "N", c3.size(), x1.size(), x0.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, x0.size()*a2.size()*c1.size(), i1data_sorted, x0.size()*a2.size()*c1.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size());
  out()->put_block(odata, x1, c3);
}

void Task892::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I942
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a2, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<3,0,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I943
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", a2.size()*c1.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, a2.size()*c1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, a2.size(), c1.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a2, c1);
}

void Task893::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I943
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task894::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I942
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, a2, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, a2, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, a2, c1), 0.0);
  for (auto& x2 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I949
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, x2, x0, x1)]);
      sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      dgemm_("T", "N", c1.size()*a2.size(), x0.size()*x1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, c1.size()*a2.size());
    }
  }
  sort_indices<2,3,1,0,1,1,1,1>(odata_sorted, odata, c1.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, a2, c1);
}

void Task895::Task_local::compute() {
  const Index x3 = b(0);
  const Index x2 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I949
  std::unique_ptr<double[]> odata = out()->move_block(x3, x2, x0, x1);
  {
    // tensor label: Gamma12
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
    sort_indices<0,1,2,3,1,1,1,2>(i0data, odata, x3.size(), x2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x3, x2, x0, x1);
}

void Task896::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  // tensor label: den2
  std::unique_ptr<double[]> odata = out()->move_block(x1, a4);
  {
    // tensor label: I950
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), a4.size());
  }
  out()->put_block(odata, x1, a4);
}

void Task897::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  // tensor label: I950
  std::unique_ptr<double[]> odata = out()->move_block(x1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        for (auto& x0 : *range_[1]) {
          // tensor label: t2
          std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
          sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
          // tensor label: I951
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1, c1, a4, c3, a2);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1, c1, a4, c3, a2)]);
          sort_indices<2,5,4,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size(), c1.size(), a4.size(), c3.size(), a2.size());
          dgemm_("T", "N", 1, x1.size()*a4.size(), x0.size()*c1.size()*c3.size()*a2.size(),
                 1.0, i0data_sorted, x0.size()*c1.size()*c3.size()*a2.size(), i1data_sorted, x0.size()*c1.size()*c3.size()*a2.size(),
                 1.0, odata_sorted, 1);
        }
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), a4.size());
  out()->put_block(odata, x1, a4);
}

void Task898::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  const Index c3 = b(4);
  const Index a2 = b(5);
  // tensor label: I951
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1, c1, a4, c3, a2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1, c1, a4, c3, a2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1, c1, a4, c3, a2), 0.0);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I952
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
  sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
  dgemm_("T", "N", c1.size()*a4.size()*c3.size()*a2.size(), x0.size()*x1.size(), 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, c1.size()*a4.size()*c3.size()*a2.size());
  sort_indices<4,5,0,1,2,3,1,1,1,1>(odata_sorted, odata, c1.size(), a4.size(), c3.size(), a2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x1, c1, a4, c3, a2);
}

void Task899::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I952
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, x0.size(), x1.size());
  }
  out()->put_block(odata, x0, x1);
}

