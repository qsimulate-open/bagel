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
  const Index x3 = b(0);
  const Index x5 = b(1);
  const Index c1 = b(2);
  const Index x4 = b(3);
  // tensor label: I818
  std::unique_ptr<double[]> odata = out()->move_block(x3, x5, c1, x4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x3, x5, c1, x4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x3, x5, c1, x4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, a2, c1, x4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, a2, c1, x4)]);
    sort_indices<1,0,2,3,0,1,1,1>(i0data, i0data_sorted, x5.size(), a2.size(), c1.size(), x4.size());
    // tensor label: I819
    std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x3)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, a2.size(), x3.size());
    dgemm_("T", "N", x5.size()*c1.size()*x4.size(), x3.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, x5.size()*c1.size()*x4.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x5.size(), c1.size(), x4.size(), x3.size());
  out()->put_block(odata, x3, x5, c1, x4);
}

void Task851::Task_local::compute() {
  const Index a2 = b(0);
  const Index x3 = b(1);
  // tensor label: I819
  std::unique_ptr<double[]> odata = out()->move_block(a2, x3);
  {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, x3);
    sort_indices<0,1,1,1,-1,2>(i0data, odata, a2.size(), x3.size());
  }
  out()->put_block(odata, a2, x3);
}

void Task852::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I782
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      for (auto& x0 : *range_[1]) {
        for (auto& x1 : *range_[1]) {
          // tensor label: Gamma284
          std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x3, x2, x0, x1);
          std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x3, x2, x0, x1)]);
          sort_indices<1,2,3,4,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x3.size(), x2.size(), x0.size(), x1.size());
          // tensor label: I829
          std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3, x2, x1);
          std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3, x2, x1)]);
          sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size(), x2.size(), x1.size());
          dgemm_("T", "N", ci0.size(), 1, x0.size()*x3.size()*x2.size()*x1.size(),
                 1.0, i0data_sorted, x0.size()*x3.size()*x2.size()*x1.size(), i1data_sorted, x0.size()*x3.size()*x2.size()*x1.size(),
                 1.0, odata_sorted, ci0.size());
        }
      }
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task853::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I829
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c3, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c3.size(), x1.size());
    // tensor label: I830
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size());
    dgemm_("T", "N", x3.size()*x2.size()*x1.size(), x0.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x3.size()*x2.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task854::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  // tensor label: I830
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size());
      // tensor label: I831
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, x0.size()*c3.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->put_block(odata, x0, c3);
}

void Task855::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I831
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task856::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I829
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, c1, x1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), c1.size(), x1.size());
    // tensor label: I834
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c1.size());
    dgemm_("T", "N", x3.size()*x2.size()*x1.size(), x0.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x3.size()*x2.size()*x1.size());
  }
  sort_indices<3,0,1,2,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x1.size(), x0.size());
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task857::Task_local::compute() {
  const Index x0 = b(0);
  const Index c1 = b(1);
  // tensor label: I834
  std::unique_ptr<double[]> odata = out()->move_block(x0, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: f1
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2)]);
      sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size());
      // tensor label: I835
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
      dgemm_("T", "N", 1, x0.size()*c1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, 1);
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x0.size(), c1.size());
  out()->put_block(odata, x0, c1);
}

void Task858::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I835
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task859::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I829
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a2, c1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      // tensor label: I872
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task860::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I872
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2, c1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I873
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", x1.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x0, a2, c1, x1);
}

void Task861::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I873
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task862::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I829
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I876
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, x1)]);
      sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), c3.size()*a2.size(),
             1.0, i0data_sorted, c3.size()*a2.size(), i1data_sorted, c3.size()*a2.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task863::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index x1 = b(3);
  // tensor label: I876
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x1.size());
    // tensor label: I877
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", x1.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), c3.size(), a2.size());
  out()->put_block(odata, x0, c3, a2, x1);
}

void Task864::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I877
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task865::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I829
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: t2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, x3, x2)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), x3.size(), x2.size());
      // tensor label: I880
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, x1)]);
      sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), x1.size());
      dgemm_("T", "N", x3.size()*x2.size(), x0.size()*x1.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,0,1,3,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), x0.size(), x1.size());
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task866::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I880
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2, c1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), x1.size());
    // tensor label: I881
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", x1.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x0, a2, c1, x1);
}

void Task867::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I881
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task868::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I829
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      // tensor label: v2
      std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0, x1, a2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0, x1, a2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size(), x1.size(), a2.size());
      // tensor label: I1266
      std::unique_ptr<double[]> i1data = in(1)->get_block(x2, x3, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x2, x3, a2, c1)]);
      sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, x2.size(), x3.size(), a2.size(), c1.size());
      dgemm_("T", "N", x0.size()*x1.size(), x2.size()*x3.size(), a2.size()*c1.size(),
             1.0, i0data_sorted, a2.size()*c1.size(), i1data_sorted, a2.size()*c1.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task869::Task_local::compute() {
  const Index x2 = b(0);
  const Index x3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I1266
  std::unique_ptr<double[]> odata = out()->move_block(x2, x3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x2, x3, a2, c1);
}

void Task870::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  const Index x2 = b(2);
  const Index x1 = b(3);
  // tensor label: I829
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3, x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3, x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3, x2, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: h1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), x0.size());
    // tensor label: I1302
    std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c1, x2, x3);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c1, x2, x3)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), c1.size(), x2.size(), x3.size());
    dgemm_("T", "N", x0.size(), x1.size()*x2.size()*x3.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<0,3,2,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), x2.size(), x3.size());
  out()->put_block(odata, x0, x3, x2, x1);
}

void Task871::Task_local::compute() {
  const Index x1 = b(0);
  const Index c1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // tensor label: I1302
  std::unique_ptr<double[]> odata = out()->move_block(x1, c1, x2, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, c1, x1);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, x3.size(), x2.size(), c1.size(), x1.size());
  }
  out()->put_block(odata, x1, c1, x2, x3);
}

void Task872::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I782
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x3 : *range_[1]) {
      // tensor label: Gamma286
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x3);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x3)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x3.size());
      // tensor label: I837
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size());
      dgemm_("T", "N", ci0.size(), 1, x0.size()*x3.size(),
             1.0, i0data_sorted, x0.size()*x3.size(), i1data_sorted, x0.size()*x3.size(),
             1.0, odata_sorted, ci0.size());
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task873::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I837
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
        // tensor label: I838
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", x3.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size());
  out()->put_block(odata, x0, x3);
}

void Task874::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I838
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task875::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I837
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
        // tensor label: I841
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", x3.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size());
  out()->put_block(odata, x0, x3);
}

void Task876::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I841
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task877::Task_local::compute() {
  const Index ci0 = b(0);
  // tensor label: I782
  std::unique_ptr<double[]> odata = out()->move_block(ci0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(ci0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(ci0), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma288
      std::unique_ptr<double[]> i0data = in(0)->get_block(ci0, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(ci0, x0, x1)]);
      sort_indices<1,2,0,0,1,1,1>(i0data, i0data_sorted, ci0.size(), x0.size(), x1.size());
      // tensor label: I843
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
      sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
      dgemm_("T", "N", ci0.size(), 1, x0.size()*x1.size(),
             1.0, i0data_sorted, x0.size()*x1.size(), i1data_sorted, x0.size()*x1.size(),
             1.0, odata_sorted, ci0.size());
    }
  }
  sort_indices<0,1,1,1,1>(odata_sorted, odata, ci0.size());
  out()->put_block(odata, ci0);
}

void Task878::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I843
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c3.size(), x1.size());
        // tensor label: I844
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c4)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c4.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c4.size(), i1data_sorted, c3.size()*a2.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task879::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c4 = b(3);
  // tensor label: I844
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I845
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size(), c3.size(), a2.size());
  out()->put_block(odata, x0, c3, a2, c4);
}

void Task880::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I845
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task881::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I843
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c4, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c4, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c4.size(), x1.size());
        // tensor label: I848
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*a2.size()*c4.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c4.size(), i1data_sorted, c3.size()*a2.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task882::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c4 = b(3);
  // tensor label: I848
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, a2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, a2, c4), 0.0);
  for (auto& c1 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), c4.size());
    // tensor label: I849
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x0.size()*c3.size()*a2.size(), c1.size(),
           1.0, i0data_sorted, c1.size(), i1data_sorted, c1.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size(), c3.size(), a2.size());
  out()->put_block(odata, x0, c3, a2, c4);
}

void Task883::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I849
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task884::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I843
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c4 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c4, a2, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c4, a2, c1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c4.size(), a2.size(), c1.size(), x1.size());
        // tensor label: I852
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, c4)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), a2.size()*c1.size()*c4.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*c4.size(), i1data_sorted, a2.size()*c1.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task885::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I852
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2, c1, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I853
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x0, a2, c1, c4);
}

void Task886::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I853
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task887::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I843
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
        // tensor label: I856
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, c1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, c1, a4)]);
        sort_indices<1,3,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), c1.size(), a4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*c1.size()*a4.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*a4.size(), i1data_sorted, c3.size()*c1.size()*a4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task888::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  // tensor label: I856
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, c1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, c1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I857
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", a4.size(), x0.size()*c3.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size(), c3.size(), c1.size());
  out()->put_block(odata, x0, c3, c1, a4);
}

void Task889::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I857
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task890::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I843
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c4 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c4, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c4, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c4.size(), x1.size());
        // tensor label: I860
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a2, c1, c4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a2, c1, c4)]);
        sort_indices<2,1,3,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a2.size(), c1.size(), c4.size());
        dgemm_("T", "N", x1.size(), x0.size(), a2.size()*c1.size()*c4.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*c4.size(), i1data_sorted, a2.size()*c1.size()*c4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task891::Task_local::compute() {
  const Index x0 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index c4 = b(3);
  // tensor label: I860
  std::unique_ptr<double[]> odata = out()->move_block(x0, a2, c1, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a2, c1, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a2, c1, c4), 0.0);
  for (auto& c3 : *range_[0]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, c4);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, c4)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), c4.size());
    // tensor label: I861
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", c4.size(), x0.size()*a2.size()*c1.size(), c3.size(),
           1.0, i0data_sorted, c3.size(), i1data_sorted, c3.size(),
           1.0, odata_sorted, c4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size(), a2.size(), c1.size());
  out()->put_block(odata, x0, a2, c1, c4);
}

void Task892::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I861
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task893::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I843
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a4 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
        // tensor label: I864
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, c1, a4);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, c1, a4)]);
        sort_indices<2,3,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), c1.size(), a4.size());
        dgemm_("T", "N", x1.size(), x0.size(), c3.size()*c1.size()*a4.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*a4.size(), i1data_sorted, c3.size()*c1.size()*a4.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task894::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index c1 = b(2);
  const Index a4 = b(3);
  // tensor label: I864
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, c1, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c3, c1, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c3, c1, a4), 0.0);
  for (auto& a2 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a2);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, a2)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, a4.size(), a2.size());
    // tensor label: I865
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", a4.size(), x0.size()*c3.size()*c1.size(), a2.size(),
           1.0, i0data_sorted, a2.size(), i1data_sorted, a2.size(),
           1.0, odata_sorted, a4.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size(), c3.size(), c1.size());
  out()->put_block(odata, x0, c3, c1, a4);
}

void Task895::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I865
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task896::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I843
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    // tensor label: f1
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(a4, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, a4.size(), x1.size());
    // tensor label: I884
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, a4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, a4)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), a4.size());
    dgemm_("T", "N", x1.size(), x0.size(), a4.size(),
           1.0, i0data_sorted, a4.size(), i1data_sorted, a4.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task897::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I884
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
        // tensor label: I885
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", a4.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->put_block(odata, x0, a4);
}

void Task898::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I885
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task899::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I884
  std::unique_ptr<double[]> odata = out()->move_block(x0, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: I889
        std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", a4.size(), x0.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), x0.size());
  out()->put_block(odata, x0, a4);
}

