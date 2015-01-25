//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks12.cc
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

#include <src/smith/CASPT2_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task550::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I505
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, c1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I506
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size(), a2.size(), c1.size());
  out()->put_block(odata, c3, a2, c1, x1);
}

void Task551::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I506
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task552::Task_local::compute() {
  const Index a2 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(a2, a4);
  {
    // tensor label: I501
    std::unique_ptr<double[]> i0data = in(0)->get_block(a2, a4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a2.size(), a4.size());
  }
  out()->put_block(odata, a2, a4);
}

void Task553::Task_local::compute() {
  const Index a2 = b(0);
  const Index a4 = b(1);
  // tensor label: I501
  std::unique_ptr<double[]> odata = out()->move_block(a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a4), 0.0);
  for (auto& c3 : *range_[0]) {
    for (auto& c1 : *range_[0]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a4, c1, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a4, c1, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c3.size(), a4.size(), c1.size(), x1.size());
        // tensor label: I502
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        dgemm_("T", "N", a4.size(), a2.size(), c3.size()*c1.size()*x1.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*x1.size(), i1data_sorted, c3.size()*c1.size()*x1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size());
  out()->put_block(odata, a2, a4);
}

void Task554::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I502
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, c1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I503
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size(), a2.size(), c1.size());
  out()->put_block(odata, c3, a2, c1, x1);
}

void Task555::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I503
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task556::Task_local::compute() {
  const Index a2 = b(0);
  const Index a4 = b(1);
  // tensor label: I501
  std::unique_ptr<double[]> odata = out()->move_block(a2, a4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a4), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& x1 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, x1)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), x1.size());
        // tensor label: I508
        std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, c1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, c1, x1)]);
        sort_indices<2,0,3,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), c1.size(), x1.size());
        dgemm_("T", "N", a4.size(), a2.size(), c3.size()*c1.size()*x1.size(),
               1.0, i0data_sorted, c3.size()*c1.size()*x1.size(), i1data_sorted, c3.size()*c1.size()*x1.size(),
               1.0, odata_sorted, a4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a4.size(), a2.size());
  out()->put_block(odata, a2, a4);
}

void Task557::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x1 = b(3);
  // tensor label: I508
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, c1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, a2, c1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, a2, c1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I509
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c3, a2, c1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c3, a2, c1)]);
    sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x0.size(), c3.size(), a2.size(), c1.size());
    dgemm_("T", "N", x1.size(), c3.size()*a2.size()*c1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size(), a2.size(), c1.size());
  out()->put_block(odata, c3, a2, c1, x1);
}

void Task558::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I509
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task559::Task_local::compute() {
  const Index x1 = b(0);
  const Index c1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, c1);
  {
    // tensor label: I510
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c1.size());
  }
  out()->put_block(odata, x1, c1);
}

void Task560::Task_local::compute() {
  const Index x1 = b(0);
  const Index c1 = b(1);
  // tensor label: I510
  std::unique_ptr<double[]> odata = out()->move_block(x1, c1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c1), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I511
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c3, x1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c3, x1, x0)]);
        sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c3.size(), x1.size(), x0.size());
        dgemm_("T", "N", c1.size(), x1.size(), a2.size()*c3.size()*x0.size(),
               1.0, i0data_sorted, a2.size()*c3.size()*x0.size(), i1data_sorted, a2.size()*c3.size()*x0.size(),
               1.0, odata_sorted, c1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c1.size(), x1.size());
  out()->put_block(odata, x1, c1);
}

void Task561::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(a2, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma22
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x1, x0, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x1, x0, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x1.size(), x0.size(), x2.size());
      // tensor label: I512
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c3, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c3.size(), x2.size());
      dgemm_("T", "N", x1.size()*x0.size(), a2.size()*c3.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x1.size()*x0.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size(), a2.size(), c3.size());
  out()->put_block(odata, a2, c3, x1, x0);
}

void Task562::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index x2 = b(3);
  // tensor label: I512
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, c3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c3, x2);
    sort_indices<0,1,2,3,1,1,1,4>(i0data, odata, x3.size(), a2.size(), c3.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, c3, x2);
}

void Task563::Task_local::compute() {
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index x1 = b(2);
  const Index x0 = b(3);
  // tensor label: I511
  std::unique_ptr<double[]> odata = out()->move_block(a2, c3, x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c3, x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c3, x1, x0), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I518
      std::unique_ptr<double[]> i1data = in(1)->get_block(c3, a2, x3, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, a2, x3, x2)]);
      sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, c3.size(), a2.size(), x3.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), c3.size()*a2.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<3,2,1,0,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), c3.size(), a2.size());
  out()->put_block(odata, a2, c3, x1, x0);
}

void Task564::Task_local::compute() {
  const Index c3 = b(0);
  const Index a2 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I518
  std::unique_ptr<double[]> odata = out()->move_block(c3, a2, x3, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, x3, x2);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, c3.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, c3, a2, x3, x2);
}

void Task565::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3);
  {
    // tensor label: I513
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, x1.size(), c3.size());
  }
  out()->put_block(odata, x1, c3);
}

void Task566::Task_local::compute() {
  const Index x1 = b(0);
  const Index c3 = b(1);
  // tensor label: I513
  std::unique_ptr<double[]> odata = out()->move_block(x1, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, c3), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x0)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x0.size());
        // tensor label: I514
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x0, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x0, x1)]);
        sort_indices<2,0,1,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x0.size(), x1.size());
        dgemm_("T", "N", c3.size(), x1.size(), a2.size()*c1.size()*x0.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*x0.size(), i1data_sorted, a2.size()*c1.size()*x0.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x1.size());
  out()->put_block(odata, x1, c3);
}

void Task567::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x0 = b(2);
  const Index x1 = b(3);
  // tensor label: I514
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x0, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x2 : *range_[1]) {
      // tensor label: Gamma12
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x2, x0, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x2, x0, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x2.size(), x0.size(), x1.size());
      // tensor label: I515
      std::unique_ptr<double[]> i1data = in(1)->get_block(x3, a2, c1, x2);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x3, a2, c1, x2)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x3.size(), a2.size(), c1.size(), x2.size());
      dgemm_("T", "N", x0.size()*x1.size(), a2.size()*c1.size(), x3.size()*x2.size(),
             1.0, i0data_sorted, x3.size()*x2.size(), i1data_sorted, x3.size()*x2.size(),
             1.0, odata_sorted, x0.size()*x1.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x0.size(), x1.size(), a2.size(), c1.size());
  out()->put_block(odata, a2, c1, x0, x1);
}

void Task568::Task_local::compute() {
  const Index x3 = b(0);
  const Index a2 = b(1);
  const Index c1 = b(2);
  const Index x2 = b(3);
  // tensor label: I515
  std::unique_ptr<double[]> odata = out()->move_block(x3, a2, c1, x2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a2, c1, x2);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x3.size(), a2.size(), c1.size(), x2.size());
  }
  {
    // tensor label: t2
    std::unique_ptr<double[]> i1data = in(0)->get_block(c1, a2, x3, x2);
    sort_indices<2,1,0,3,1,1,1,2>(i1data, odata, c1.size(), a2.size(), x3.size(), x2.size());
  }
  out()->put_block(odata, x3, a2, c1, x2);
}

void Task569::Task_local::compute() {
  const Index x1 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, a4);
  {
    // tensor label: I522
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, x1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, a4.size(), x1.size());
  }
  out()->put_block(odata, x1, a4);
}

void Task570::Task_local::compute() {
  const Index a4 = b(0);
  const Index x1 = b(1);
  // tensor label: I522
  std::unique_ptr<double[]> odata = out()->move_block(a4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma16
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, x1);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, x1)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), x1.size());
    // tensor label: I523
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

void Task571::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I523
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
        // tensor label: I524
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

void Task572::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I524
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task573::Task_local::compute() {
  const Index x0 = b(0);
  const Index a4 = b(1);
  // tensor label: I523
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
        // tensor label: I527
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

void Task574::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I527
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x0);
    sort_indices<3,2,1,0,1,1,1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3, a2, c1);
}

void Task575::Task_local::compute() {
  const Index a1 = b(0);
  const Index x2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(a1, x2);
  {
    // tensor label: I528
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, a1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), a1.size());
  }
  out()->put_block(odata, a1, x2);
}

void Task576::Task_local::compute() {
  const Index x2 = b(0);
  const Index a1 = b(1);
  // tensor label: I528
  std::unique_ptr<double[]> odata = out()->move_block(x2, a1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, a1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, a1), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& x0 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x0, a1, c2, x1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i0data, i0data_sorted, x0.size(), a1.size(), c2.size(), x1.size());
        // tensor label: I529
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, x1, x2, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, x1, x2, x0)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), x1.size(), x2.size(), x0.size());
        dgemm_("T", "N", a1.size(), x2.size(), c2.size()*x1.size()*x0.size(),
               1.0, i0data_sorted, c2.size()*x1.size()*x0.size(), i1data_sorted, c2.size()*x1.size()*x0.size(),
               1.0, odata_sorted, a1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a1.size(), x2.size());
  out()->put_block(odata, x2, a1);
}

void Task577::Task_local::compute() {
  const Index c2 = b(0);
  const Index x1 = b(1);
  const Index x2 = b(2);
  const Index x0 = b(3);
  // tensor label: I529
  std::unique_ptr<double[]> odata = out()->move_block(c2, x1, x2, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, x1, x2, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, x1, x2, x0), 0.0);
  for (auto& x5 : *range_[1]) {
    for (auto& x4 : *range_[1]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: Gamma28
        std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, x1, x3, x2, x0);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x5, x4, x1, x3, x2, x0)]);
        sort_indices<0,1,3,2,4,5,0,1,1,1>(i0data, i0data_sorted, x5.size(), x4.size(), x1.size(), x3.size(), x2.size(), x0.size());
        // tensor label: I530
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

void Task578::Task_local::compute() {
  const Index x5 = b(0);
  const Index x4 = b(1);
  const Index c2 = b(2);
  const Index x3 = b(3);
  // tensor label: I530
  std::unique_ptr<double[]> odata = out()->move_block(x5, x4, c2, x3);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x5, x4, c2, x3);
    sort_indices<0,1,2,3,1,1,-1,4>(i0data, odata, x5.size(), x4.size(), c2.size(), x3.size());
  }
  out()->put_block(odata, x5, x4, c2, x3);
}

void Task579::Task_local::compute() {
  const Index c3 = b(0);
  const Index x2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c3, x2);
  {
    // tensor label: I531
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, c3);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), c3.size());
  }
  out()->put_block(odata, c3, x2);
}

void Task580::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I531
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c3), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a1, c2, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a1, c2, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), a1.size(), c2.size(), x3.size());
        // tensor label: I532
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
        sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
        dgemm_("T", "N", c3.size(), x2.size(), c2.size()*a1.size()*x3.size(),
               1.0, i0data_sorted, c2.size()*a1.size()*x3.size(), i1data_sorted, c2.size()*a1.size()*x3.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x2.size());
  out()->put_block(odata, x2, c3);
}

void Task581::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I532
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma29
      std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x3, x2, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x3, x2, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), x3.size(), x2.size(), x0.size());
      // tensor label: I533
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*a1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task582::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I533
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task583::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I531
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c3, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c3, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c3.size(), x3.size());
        // tensor label: I535
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x2, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x2.size(), x3.size());
        dgemm_("T", "N", c3.size(), x2.size(), c2.size()*a1.size()*x3.size(),
               1.0, i0data_sorted, c2.size()*a1.size()*x3.size(), i1data_sorted, c2.size()*a1.size()*x3.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x2.size());
  out()->put_block(odata, x2, c3);
}

void Task584::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // tensor label: I535
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x2, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma7
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I536
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<0,3,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x2.size()*x3.size(), c2.size()*a1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x2.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x3.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, a1, x2, x3);
}

void Task585::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I536
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task586::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I531
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c3), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c3, a2, c1, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c3, a2, c1, x3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, c3.size(), a2.size(), c1.size(), x3.size());
        // tensor label: I574
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x2, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x2, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x2.size(), x3.size());
        dgemm_("T", "N", c3.size(), x2.size(), a2.size()*c1.size()*x3.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*x3.size(), i1data_sorted, a2.size()*c1.size()*x3.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x2.size());
  out()->put_block(odata, x2, c3);
}

void Task587::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // tensor label: I574
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, x2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x2, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma7
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I575
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x2.size()*x3.size(), a2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x2.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x3.size(), a2.size(), c1.size());
  out()->put_block(odata, a2, c1, x2, x3);
}

void Task588::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I575
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task589::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I531
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& x3 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, x3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, x3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), x3.size());
        // tensor label: I577
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, c1, x2, x3);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, c1, x2, x3)]);
        sort_indices<1,0,3,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), c1.size(), x2.size(), x3.size());
        dgemm_("T", "N", c3.size(), x2.size(), a2.size()*c1.size()*x3.size(),
               1.0, i0data_sorted, a2.size()*c1.size()*x3.size(), i1data_sorted, a2.size()*c1.size()*x3.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x2.size());
  out()->put_block(odata, x2, c3);
}

void Task590::Task_local::compute() {
  const Index a2 = b(0);
  const Index c1 = b(1);
  const Index x2 = b(2);
  const Index x3 = b(3);
  // tensor label: I577
  std::unique_ptr<double[]> odata = out()->move_block(a2, c1, x2, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, c1, x2, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, c1, x2, x3), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma7
      std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x3, x1, x0);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x2, x3, x1, x0)]);
      sort_indices<2,3,0,1,0,1,1,1>(i0data, i0data_sorted, x2.size(), x3.size(), x1.size(), x0.size());
      // tensor label: I578
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, x0, a2, c1);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, x0, a2, c1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, x1.size(), x0.size(), a2.size(), c1.size());
      dgemm_("T", "N", x2.size()*x3.size(), a2.size()*c1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x2.size()*x3.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x3.size(), a2.size(), c1.size());
  out()->put_block(odata, a2, c1, x2, x3);
}

void Task591::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I578
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, x0, x1);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), x0.size(), x1.size());
  }
  out()->put_block(odata, x1, x0, a2, c1);
}

void Task592::Task_local::compute() {
  const Index x2 = b(0);
  const Index c3 = b(1);
  // tensor label: I531
  std::unique_ptr<double[]> odata = out()->move_block(x2, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, c3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c3, a2)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c3.size(), a2.size());
        // tensor label: I726
        std::unique_ptr<double[]> i1data = in(1)->get_block(a2, a1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, a1, x3, x2)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a2.size(), a1.size(), x3.size(), x2.size());
        dgemm_("T", "N", c3.size(), x2.size(), a2.size()*a1.size()*x3.size(),
               1.0, i0data_sorted, a2.size()*a1.size()*x3.size(), i1data_sorted, a2.size()*a1.size()*x3.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), x2.size());
  out()->put_block(odata, x2, c3);
}

void Task593::Task_local::compute() {
  const Index a2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I726
  std::unique_ptr<double[]> odata = out()->move_block(a2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a2, a1, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma60
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<1,3,0,2,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I727
      std::unique_ptr<double[]> i1data = in(1)->get_block(a2, x1, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a2, x1, a1, x0)]);
      sort_indices<3,1,0,2,0,1,1,1>(i1data, i1data_sorted, a2.size(), x1.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), a2.size()*a1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), a2.size(), a1.size());
  out()->put_block(odata, a2, a1, x3, x2);
}

void Task594::Task_local::compute() {
  const Index a2 = b(0);
  const Index x1 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I727
  std::unique_ptr<double[]> odata = out()->move_block(a2, x1, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, x1, a2);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), x1.size(), a2.size());
  }
  out()->put_block(odata, a2, x1, a1, x0);
}

void Task595::Task_local::compute() {
  const Index a1 = b(0);
  const Index a3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(a1, a3);
  {
    // tensor label: I543
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a1.size(), a3.size());
  }
  out()->put_block(odata, a1, a3);
}

void Task596::Task_local::compute() {
  const Index a1 = b(0);
  const Index a3 = b(1);
  // tensor label: I543
  std::unique_ptr<double[]> odata = out()->move_block(a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a3), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& c2 : *range_[0]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), x2.size());
        // tensor label: I544
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
        sort_indices<2,0,3,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
        dgemm_("T", "N", a3.size(), a1.size(), c2.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, c2.size()*x3.size()*x2.size(), i1data_sorted, c2.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size());
  out()->put_block(odata, a1, a3);
}

void Task597::Task_local::compute() {
  const Index c2 = b(0);
  const Index a1 = b(1);
  const Index x3 = b(2);
  const Index x2 = b(3);
  // tensor label: I544
  std::unique_ptr<double[]> odata = out()->move_block(c2, a1, x3, x2);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, a1, x3, x2)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, a1, x3, x2), 0.0);
  for (auto& x0 : *range_[1]) {
    for (auto& x1 : *range_[1]) {
      // tensor label: Gamma32
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x1, x2);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x1, x2)]);
      sort_indices<1,2,0,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x1.size(), x2.size());
      // tensor label: I545
      std::unique_ptr<double[]> i1data = in(1)->get_block(x1, c2, a1, x0);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x1, c2, a1, x0)]);
      sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, x1.size(), c2.size(), a1.size(), x0.size());
      dgemm_("T", "N", x3.size()*x2.size(), c2.size()*a1.size(), x1.size()*x0.size(),
             1.0, i0data_sorted, x1.size()*x0.size(), i1data_sorted, x1.size()*x0.size(),
             1.0, odata_sorted, x3.size()*x2.size());
    }
  }
  sort_indices<2,3,0,1,1,1,1,1>(odata_sorted, odata, x3.size(), x2.size(), c2.size(), a1.size());
  out()->put_block(odata, c2, a1, x3, x2);
}

void Task598::Task_local::compute() {
  const Index x1 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I545
  std::unique_ptr<double[]> odata = out()->move_block(x1, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, x1);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), x1.size());
  }
  out()->put_block(odata, x1, c2, a1, x0);
}

void Task599::Task_local::compute() {
  const Index a1 = b(0);
  const Index a3 = b(1);
  // tensor label: I543
  std::unique_ptr<double[]> odata = out()->move_block(a1, a3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a1, a3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a1, a3), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& x3 : *range_[1]) {
      for (auto& x2 : *range_[1]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, x3, x2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), x3.size(), x2.size());
        // tensor label: I553
        std::unique_ptr<double[]> i1data = in(1)->get_block(c2, a1, x3, x2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c2, a1, x3, x2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i1data, i1data_sorted, c2.size(), a1.size(), x3.size(), x2.size());
        dgemm_("T", "N", a3.size(), a1.size(), c2.size()*x3.size()*x2.size(),
               1.0, i0data_sorted, c2.size()*x3.size()*x2.size(), i1data_sorted, c2.size()*x3.size()*x2.size(),
               1.0, odata_sorted, a3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a3.size(), a1.size());
  out()->put_block(odata, a1, a3);
}

#endif
